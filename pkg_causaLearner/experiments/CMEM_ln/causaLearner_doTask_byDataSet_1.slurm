#!/bin/bash
#SBATCH --array=67%16
#SBATCH --ntasks=1
# SBATCH --nodelist=nodo20
#SBATCH --cpus-per-task=56
#SBATCH --exclude=nodo16,nodo17
#SBATCH --job-name=R_caus_job
#SBATCH --output=log/R_caus_job_%A_%a.out

# Para usar el R del sistema
# module load R


# Run command
Rscript script_dag5-ME32-Cmplx_CMEMJ.R
