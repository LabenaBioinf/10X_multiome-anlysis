#!/bin/bash
#SBATCH --mem=800G
#SBATCH --job-name=FigR_Control
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --output=FigR_Control_gall.%j.out

srun singularity exec --bind /path/to/FigR/objects:/home/data SC_GEX-ATAC.sif Rscript FigR_networks_Control.R
