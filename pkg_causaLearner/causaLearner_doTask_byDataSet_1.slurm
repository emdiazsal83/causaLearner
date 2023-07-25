#!/bin/bash
#SBATCH --array=53-100%20
#SBATCH --ntasks=1
# SBATCH --nodelist=nodo20
#SBATCH --cpus-per-task=56
#SBATCH --exclude=nodo17,nodo16
#SBATCH --job-name=R_caus_job
#SBATCH --output=log/R_caus_job_%A_%a.out

# Para usar el R del sistema
# module load R

# Run command
Rscript script_dag2-ME2-LINGAUSS_CMEMboot.R
