#!/bin/sh

model_name=$1
simulation=$2
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
#-----------Move to workplace directory and prepare GC/MAPPABILITY files
cd "${model_name}"
var_gc_filename="${model_name}_gc_${simulation}.wig"
var_map_filename="${model_name}_map_${simulation}.wig"
#----------------------Run HMMcopy for each sample/cell based on true CN
for (( sample=0; sample<${N_sample}; sample++ ));
do
      opt=${Sample_ID[$sample]}
      temp="${opt%\"}"
      temp="${temp#\"}"
      sample_id=$temp
      for (( cell=1; cell<=${Cell_count[$sample]}; cell++));
      do
            #-----Prepare variables for running HMMcopy
            var_cell_id="$sample_id-Library-${cell}-${cell}"
            var_input_filename="${model_name}_noisy_cn_profiles_long_${simulation}_$var_cell_id.wig"
            if [[ -f "$var_input_filename" ]]; then
                  #-----Correct readcounts for GC/mapp biases
                  var_correction_filename="${model_name}_noisy_cn_profiles_long_${simulation}_${var_cell_id}_hmm_corrected.txt"
                  echo "PERFORMING CORRECT_READCOUNT FOR $var_cell_id IN SIMULATION $simulation"
                  docker run -v $PWD:$PWD -w $PWD quay.io/mondrianscwgs/hmmcopy:v0.0.45  \
                        hmmcopy_utils correct_readcount  \
                        --infile $var_input_filename  \
                        --outfile $var_correction_filename  \
                        --map_cutoff 0.9  \
                        --gc_wig_file $var_gc_filename  \
                        --map_wig_file $var_map_filename  \
                        --cell_id $var_cell_id
                  #-----Infer CN states with HMMcopy
                  tempdir="output-$simulation-$var_cell_id"
                  reads="reads-$simulation-$var_cell_id.csv.gz"
                  metrics="metrics-$simulation-$var_cell_id.csv.gz"
                  params="params-$simulation-$var_cell_id.csv.gz"
                  segments="segments-$simulation-$var_cell_id.csv.gz"
                  var_hmmcopy_tarball="$simulation-$var_cell_id.tar.gz"
                  echo "PERFORMING HMMCOPY FOR $var_cell_id IN SIMULATION $simulation"
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
                  cp "$tempdir/0/reads.csv" "${model_name}_noisy_cn_profiles_long_${simulation}_${var_cell_id}_hmm_reads.csv"
                  cp "$tempdir/0/segs.csv" "${model_name}_noisy_cn_profiles_long_${simulation}_${var_cell_id}_hmm_segs.csv"
            fi
      done
done
wait
#--------------------Run HMMcopy for each sample/cell based on neuvar CN
for (( sample=0; sample<${N_sample}; sample++ ));
do
      opt=${Sample_ID[$sample]}
      temp="${opt%\"}"
      temp="${temp#\"}"
      sample_id=$temp
      for (( cell=1; cell<=${Cell_count[$sample]}; cell++));
      do
            #-----Prepare variables for running HMMcopy
            var_cell_id="$sample_id-Library-${cell}-${cell}"
            var_input_filename="${model_name}_noisy_neuvar_cn_profiles_long_${simulation}_$var_cell_id.wig"
            if [[ -f "$var_input_filename" ]]; then
                  #-----Correct readcounts for GC/mapp biases
                  var_correction_filename="${model_name}_noisy_neuvar_cn_profiles_long_${simulation}_${var_cell_id}_hmm_corrected.txt"
                  echo "PERFORMING CORRECT_READCOUNT FOR $var_cell_id IN SIMULATION $simulation"
                  docker run -v $PWD:$PWD -w $PWD quay.io/mondrianscwgs/hmmcopy:v0.0.45  \
                        hmmcopy_utils correct_readcount  \
                        --infile $var_input_filename  \
                        --outfile $var_correction_filename  \
                        --map_cutoff 0.9  \
                        --gc_wig_file $var_gc_filename  \
                        --map_wig_file $var_map_filename  \
                        --cell_id $var_cell_id
                  #-----Infer CN states with HMMcopy
                  tempdir="output-$simulation-$var_cell_id"
                  reads="reads-$simulation-$var_cell_id.csv.gz"
                  metrics="metrics-$simulation-$var_cell_id.csv.gz"
                  params="params-$simulation-$var_cell_id.csv.gz"
                  segments="segments-$simulation-$var_cell_id.csv.gz"
                  var_hmmcopy_tarball="$simulation-$var_cell_id.tar.gz"
                  echo "PERFORMING HMMCOPY FOR $var_cell_id IN SIMULATION $simulation"
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
                  cp "$tempdir/0/reads.csv" "${model_name}_noisy_neuvar_cn_profiles_long_${simulation}_${var_cell_id}_hmm_reads.csv"
                  cp "$tempdir/0/segs.csv" "${model_name}_noisy_neuvar_cn_profiles_long_${simulation}_${var_cell_id}_hmm_segs.csv"
            fi
      done
done
wait
