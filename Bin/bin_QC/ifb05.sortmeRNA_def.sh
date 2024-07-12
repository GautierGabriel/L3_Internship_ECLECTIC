#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J full_sortmeRNA
#SBATCH -o full_sortmeRNA.out
#SBATCH -e full_sortmeRNA.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-119%30

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
module load sortmerna/4.3.6

# Getting the sample list and creating  bash array with it
mapfile list < /shared/projects/malus_rnaseq/bin/sample_list_2.txt
sample=${list[$SLURM_ARRAY_TASK_ID]}
sample=${sample::-1}
echo $sample

# Creating the output directories
mkdir /shared/projects/malus_rnaseq/output/EXP3_full.sortmerna/$sample

# Launching Sort Me RNA : ban RNA listed in ref
    sortmerna \
            -ref /shared/projects/malus_rnaseq/bin/sortmerna_database/smr_v4.3_default_db.fasta \
            -reads /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample/fpf.$sample"_1.fq.gz" \
            -reads /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample/fpf.$sample"_2.fq.gz" \
            -threads 8 \
            -workdir /shared/projects/malus_rnaseq/output/EXP3_full.sortmerna/$sample \
            -fastx \
            -aligned /shared/projects/malus_rnaseq/output/EXP3_full.sortmerna/$sample/rRNA_fpf.$sample \
            -other /shared/projects/malus_rnaseq/output/EXP3_full.sortmerna/$sample/str_fpf.$sample \
            -paired_in \
            -out2

# Delete unuse files
    rm -rf /shared/projects/malus_rnaseq/output/EXP3_full.sortmerna/$sample/idx
    rm -rf /shared/projects/malus_rnaseq/output/EXP3_full.sortmerna/$sample/kvdb
    rm -rf /shared/projects/malus_rnaseq/output/EXP3_full.sortmerna/$sample/readb
