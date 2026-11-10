#!/bin/bash
#SBATCH --mem=1000G
#SBATCH --job-name=FigR_networks_MISC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --output=FigR_networks_MISC_all.%j.out

srun singularity exec --bind /path/to/FigR/objects:/home/data SC_GEX-ATAC.sif Rscript FigR_networks_MISC.R