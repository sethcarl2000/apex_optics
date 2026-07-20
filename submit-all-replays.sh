#!/bin/bash

run_first=$1
run_last=$2

# list all replay files
run=${run_first}

task_id=-1

rm slurm/array/* 

while [[ ${run} -lt ${run_last} ]]
do
    
    task_id=$(( ${task_id} + 1 ))
    
    run1=$(( ${run} + 24 ))

    path_infile="/volatile/halla/apex/full_replay/production/replay-${run}-${run1}.root"
    path_outfile="data/optics_replay/replay-${run}-${run1}.root"

    echo "${path_infile}|${path_outfile}" > "slurm/array/task_${task_id}.list"

    run=$(( ${run} + 25 ))

done
    
cmd="sbatch --job-name=apex_optics_replay --partition=production --time=90 --array=0-${task_id} --mem=5G --cpus-per-task=10 submit-script.sh"

echo "${cmd}" 

eval $cmd
