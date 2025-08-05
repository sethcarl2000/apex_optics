#!/bin/bash

#SBATCH --nodes 1
#SBATCH --partition production
#SBATCH --account halla
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20 
#SBATCH --output=/volatile/halla/apex/mc/temp/slurm/%x.%j.out
#SBATCH --error=/volatile/halla/apex/mc/temp/slurm/%x.%j.err
#SBATCH --mem=3G  

########################################################################
#
# This job submits the script to 'train' a MLP for a given number of epochs.
#
#

root -l -b -q 'scripts/train_new_mlp.C(200, "data/mc/out_fp_L_production.root", "data/csv/mlp_sv_fp_554.dat", {5}, "plots/train_sv_fp_L_554.pdf")' 
