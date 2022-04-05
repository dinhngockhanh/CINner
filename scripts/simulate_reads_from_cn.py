import numpy as np
import pandas as pd
from argparse import ArgumentParser


## Example usage:
## python3 scripts/simulate_reads_from_cn.py -i vignettes/EXPERIMENT-1A-SELECTIVE_cn_profiles_1.csv -g data/gc_map_500kb.csv 
## -m Min -o vignettes/EXPERIMENT-1A-SELECTIVE_read_profiles_1.csv
def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--cn_input', help='input: long-form csv containing chr, cell_id, state, etc')
    p.add_argument('-g', '--gc_input', nargs='?', help='long-form csv containing gc and mappability columns for each 500kb bin')
    p.add_argument('-s', '--sigma1', default=0.1, nargs='?', type=float, help='noise of read depth profiles')
    p.add_argument('-g1', '--gc_slope', default=1.2, nargs='?', type=float, help='slope of linear GC bias')
    p.add_argument('-g0', '--gc_int', default=0.0, nargs='?', type=float, help='intercept of linear GC bias')
    p.add_argument('-n', '--num_reads', default=1E6, nargs='?', type=int, help='number of reads per cell')
    p.add_argument('-m', '--minor_state_col', default=None, nargs='?', help='column representing true copy number state for the minor allele')
    p.add_argument('-o', '--output', help='same as input with read count and relevant simulation parameters added')

    return p.parse_args()


def model(true_CN, gc, sigma1, gc_slope, gc_int, num_reads, BAF=None):
    """Given the true copy number state, GC, and noise values for a cell, come up with that observed read count profile."""    
    # add gc bias to the true CN
    # Is a simple linear model sufficient here?
    observed_CN = true_CN * ((gc * gc_slope) + gc_int)
    
    # add some random noise to the observed copy number
    noisy_CN = np.random.gamma(observed_CN / sigma1, sigma1)
    
    # scale noisy_CN and then draw true read count from multinomial distribution
    noisy_CN_pval = noisy_CN / sum(noisy_CN)
    read_count = np.random.multinomial(num_reads, noisy_CN_pval)

    # draw read count for major and minor alleles if given a major allele fraction
    major_read_count = None
    minor_read_count = None
    if BAF is not None:
        minor_read_count = np.random.binomial(read_count, BAF)
        major_read_count = read_count - minor_read_count
    
    return read_count, major_read_count, minor_read_count


def simulate_cell(temp_df, cell_id, num_reads, gc_slope, gc_int, sigma1, cn_state_col='state', minor_state_col=None):    
    # two model modes: one with allele-specific read count and one without
    if minor_state_col is None:
        read_count, _, _ = model(
            temp_df[cn_state_col].values, temp_df['gc'].values,
            sigma1, gc_slope, gc_int, num_reads
        )
    else:
        BAF = temp_df[minor_state_col].values / temp_df[cn_state_col].values
        read_count, major_read_count, minor_read_count = model(
            temp_df[cn_state_col].values, temp_df['gc'].values,
            sigma1, gc_slope, gc_int, num_reads, BAF
        )
        temp_df['major_reads'] = major_read_count
        temp_df['minor_reads'] = minor_read_count
    
    # store simulated cell values in temp_df
    temp_df['reads'] = read_count
    temp_df['total_reads_per_cell'] = sum(read_count)
    temp_df['gc_slope'] = gc_slope
    temp_df['gc_int'] = gc_int
    temp_df['sigma1'] = sigma1
    temp_df['reads_per_million'] = (temp_df['reads'] / sum(temp_df['reads'].values)) * 1e6
    
    return temp_df


def simulate_reads_from_cn(cn, gc_slope, gc_int, sigma1, num_reads, cn_state_col='state', minor_state_col=None):
    """ Given a df with somatic CN states and true S-phase times for all cells, simulate read count with GC and replication bias. """
    sim_df = []
    
    for cell_id, chunk in cn.groupby('cell_id'):
        # TODO: draw gc_slope, gc_int, and sigma1 from distribution for each cell
        # instead of having them constant for all cells in a simulated dataset

        # simulate gc bias and to get read counts
        temp_df = simulate_cell(chunk, cell_id, num_reads, gc_slope, gc_int, sigma1, 
                                cn_state_col=cn_state_col, minor_state_col=minor_state_col)
        sim_df.append(temp_df)

    sim_df = pd.concat(sim_df, ignore_index=True)
    
    return sim_df


def main():
    argv = get_args()

    cn = pd.read_csv(argv.cn_input, dtype={'chr': str})

    # bring in GC values from reference genome
    if 'gc' not in cn.columns:
        gc = pd.read_csv(argv.gc_input, dtype={'chr': str})
        cn = pd.merge(cn, gc)

    print('simulating cells...')
    cn = simulate_reads_from_cn(cn, argv.gc_slope, argv.gc_int, argv.sigma1, argv.num_reads, cn_state_col='state', minor_state_col=argv.minor_state_col)
    print('cn.head()\n', cn.head())

    cn.to_csv(argv.output, index=False)



if __name__ == '__main__':
    main()
