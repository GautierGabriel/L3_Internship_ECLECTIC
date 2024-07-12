#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J full_quant
#SBATCH -o full_quant.out
#SBATCH -e full_quant.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
module load salmon/1.10.2
echo "Job ID : "$SLURM_JOB_ID
echo "START"
date '+%F %H:%M:%S'


salmon quantmerge --column numreads --quants /shared/ifbstor1/projects/malus_rnaseq/output/EXP3_full.quant_align/*  -o "/shared/ifbstor1/projects/malus_rnaseq/output/EXP3_full.quant_merge/EXP3_Salmon_quant.sf"

echo "DONE"
date '+%F %H:%M:%S'