#!/bin/bash

model_name=$1
n_simulations=$2
flag_parallel=$3
#-------------------------------Find ID's and cell counts of each sample
declare -a Sample_ID
declare -a Cell_count
INPUT="${model_name}-input-sampling.csv"
i=0
{
      read
      while IFS=, read -r next_Sample_ID next_Cell_count next_Age_sample next_T_sample
      do
           Sample_ID[$i]=$next_Sample_ID
           Cell_count[$i]=$next_Cell_count
           i=$((i+1))
      done
} < ${INPUT}
N_sample=${#Sample_ID[@]}
#------------------------------Function to run HMMcopy for a single cell
function HMMcopy_one_cell(){
      cell=$1
      sample_id=$2
      simulation=$3
      #-----Prepare variables for running HMMcopy
      var_cell_id="$sample_id-Library-${cell}-${cell}"
      var_input_filename="${model_name}_noisy_cn_profiles_long_${simulation}_$sample_id-Library-${cell}-${cell}.wig"
      #-----Correct readcounts for GC/mapp biases
      var_correction_filename="${model_name}_noisy_cn_profiles_long_${simulation}_$sample_id-Library-${cell}-${cell}_hmm_corrected.txt"
      docker run -v $PWD:$PWD -w $PWD quay.io/mondrianscwgs/hmmcopy:v0.0.45  \
            hmmcopy_utils correct_readcount  \
            --infile $var_input_filename  \
            --outfile $var_correction_filename  \
            --map_cutoff 0.9  \
            --gc_wig_file $var_gc_filename  \
            --map_wig_file $var_map_filename  \
            --cell_id $var_cell_id
      #-----Infer CN states with HMMcopy
      tempdir="output-$sample_id-Library-${cell}-${cell}"
      reads="reads-$sample_id-Library-${cell}-${cell}.csv.gz"
      metrics="metrics-$sample_id-Library-${cell}-${cell}.csv.gz"
      params="params-$sample_id-Library-${cell}-${cell}.csv.gz"
      segments="segments-$sample_id-Library-${cell}-${cell}.csv.gz"
      var_hmmcopy_tarball="$sample_id-Library-${cell}-${cell}.tar.gz"
      docker run -v $PWD:$PWD -w $PWD quay.io/mondrianscwgs/hmmcopy:v0.0.45  \
            hmmcopy_utils run_hmmcopy  \
            --corrected_reads $var_correction_filename  \
            --tempdir $tempdir  \
            --reads $reads \
            --metrics $metrics \
            --params $params \
            --segments $segments \
            --output_tarball $var_hmmcopy_tarball
            # --tempdir output  \
            # --reads reads.csv.gz \
            # --metrics metrics.csv.gz \
            # --params params.csv.gz \
            # --segments segments.csv.gz \
      #-----Move HMMcopy results outside
      cp "$tempdir/0/reads.csv" "${model_name}_noisy_cn_profiles_long_${simulation}_$sample_id-Library-${cell}-${cell}_hmm_reads.csv"
      cp "$tempdir/0/segs.csv" "${model_name}_noisy_cn_profiles_long_${simulation}_$sample_id-Library-${cell}-${cell}_hmm_segs.csv"
}
#----------------------------Run HMMcopy for each simulation/sample/cell
#-----Copy GC and Mappability WIG files to workspace
cp "${model_name}_gc.wig" "${model_name}/${model_name}_gc.wig"
cp "${model_name}_map.wig" "${model_name}/${model_name}_map.wig"
cd "${model_name}"
var_gc_filename="${model_name}_gc.wig"
var_map_filename="${model_name}_map.wig"
#-----Run HMMcopy for each simulation/sample/cell
if [[ $flag_parallel -eq 0 ]]
then
      #-----Run HMMcopy in sequential mode
      for (( simulation=1; simulation<=${n_simulations}; simulation++));
      do
            for (( sample=0; sample<${N_sample}; sample++ ));
            do
                  opt=${Sample_ID[$sample]}
                  temp="${opt%\"}"
                  temp="${temp#\"}"
                  sample_id=$temp
                  for (( cell=1; cell<=${Cell_count[$sample]}; cell++));
                  do
                        HMMcopy_one_cell $cell $sample_id $simulation
                  done
                  wait
            done
      done
else
      #-----Run HMMcopy in parallel mode
      for (( simulation=1; simulation<=${n_simulations}; simulation++));
      do
            for (( sample=0; sample<${N_sample}; sample++ ));
            do
                  opt=${Sample_ID[$sample]}
                  temp="${opt%\"}"
                  temp="${temp#\"}"
                  sample_id=$temp
                  for (( cell=1; cell<=${Cell_count[$sample]}; cell++));
                  do
                        HMMcopy_one_cell $cell $sample_id $simulation &
                  done
                  wait
            done
      done
fi
