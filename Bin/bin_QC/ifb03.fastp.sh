#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J EXP3_fastp
#SBATCH -o EXP3_fastp.out
#SBATCH -e EXP3_fastp.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-119%30

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
  
module load fastqc/0.12.1
module load fastp/0.23.1

# Getting the sample list and creating  bash array with it
mapfile list < /shared/projects/malus_rnaseq/bin/sample_list.txt
sample=${list[$SLURM_ARRAY_TASK_ID]}
sample=${sample::-1}
echo $sample

# Creating the output directories
mkdir /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample

# Launching fastp, trimming -f base, minimal lenght of -l
fastp -f 0 -l 100 -A -Q \
-i /shared/projects/malus_rnaseq/output/EXP3_full/$sample/fpf.$sample"_1.fq.gz" -I /shared/projects/malus_rnaseq/output/EXP3_full/$sample/fpf.$sample"_2.fq.gz" \
-o /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample/fpf.$sample"_1.fq.gz" -O /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample/fpf.$sample"_2.fq.gz" \
-h /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample/fastp.$sample.html -j /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample/fastp.$sample.json -R '$sample' &&\

# Launching QC on trimmed datas
fastqc -t 4 -f fastq -o /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample/fpf.$sample"_1.fq.gz" &&\
fastqc -t 4 -f fastq -o /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample /shared/projects/malus_rnaseq/output/EXP3_full_filter/$sample/fpf.$sample"_2.fq.gz" &&\

echo done