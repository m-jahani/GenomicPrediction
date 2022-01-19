#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for BayPass as job arrays
# ---------------------------------------------------------------------
#SBATCH --account=def-rieseber
#SBATCH --job-name=KfoldCrossValidation
#SBATCH --cpus-per-task=64
#SBATCH --mem=255000M
#SBATCH --time=02:59:00
#SBATCH --output=%j_%x.out
#SBATCH --array=1-46
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""

#Load the conda environment
source /home/mjahani/miniconda3/bin/activate rrblup

#input data
ROOT_DIR=/home/mjahani/scratch/GP
SCRIPT=${ROOT_DIR}/Kfold_cross_validation.R
GENOTYPE=${ROOT_DIR}/SAM.5maf.noheader.sed12.rrblup.csv
PHENOTYPE=${ROOT_DIR}/phenotype_common_georgia_corrected.csv

# Output directory
SAVE_DIR=${ROOT_DIR}/result

#run th escript
Rscript $SCRIPT $GENOTYPE $PHENOTYPE $SLURM_ARRAY_TASK_ID $SAVE_DIR 10


# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------
