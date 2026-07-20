#!/bin/bash

#SBATCH --nodes 1
#SBATCH --partition production
#SBATCH --account halla
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=4G
#SBATCH --time 60
#SBATCH --gres=disk:6G
#SBATCH --output=/work/halla/apex/disk1/sethhall/optics/apex_optics/slurm/%x.%a.%A.out
#SBATCH --error=/work/halla/apex/disk1/sethhall/optics/apex_optics/slurm/%x.%a.%A.err

source set_optics_replay_env.sh

cd /work/halla/apex/disk1/sethhall/optics/apex_optics

slurm_array_file="slurm/array/task_${SLURM_ARRAY_TASK_ID}.list"

PATH_SCRATCH="/scratch/slurm/${SLURM_JOB_ID}"

path_infile="${PATH_APEX_VOLATILE}/production/replay-${run0}-${run1}.root"
path_outfile="data/optics_replay/replay_${run0}_${run1}.root"

while read -r array_line
do

    # now, split the line into an array
    IFS='|' read -ra line_array <<< "${array_line}"

    input_file="${line_array[0]}"
    output_file="${line_array[1]}"

    temp_file="${PATH_SCRATCH}/temp.root"
    
    time root -l "scripts/invariant_mass_reco.C(\"${input_file}\", \"${output_file}\", false, \"${temp_file}\")"

done < <(cat "${slurm_array_file}")

