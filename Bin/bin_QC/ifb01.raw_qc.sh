#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J EXP3_raw_qc
#SBATCH -o EXP3_raw_qc.out
#SBATCH -e EXP3_raw_qc.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-120%16

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge  

module load fastqc/0.12.1

# Defining the working directory
mapfile list < /shared/projects/malus_rnaseq/bin/sample_list.txt

# Getting the sample list and creating  bash array with it
sample=${list[$SLURM_ARRAY_TASK_ID]}
sample=${sample::-1}
echo $sample

# Creating the output directories
#mkdir /shared/projects/malus_rnaseq/output/EXP3.fastp/$sample

# Launching QC
fastqc -t 4 -f fastq -o /shared/projects/malus_rnaseq/output/EXP3.fastp/$sample /shared/projects/malus_rnaseq/output/EXP3.fastp/$sample/fpf.$sample"_1.fq.gz" &&\ 
fastqc -t 4 -f fastq -o /shared/projects/malus_rnaseq/output/EXP3.fastp/$sample /shared/projects/malus_rnaseq/output/EXP3.fastp/$sample/fpf.$sample"_2.fq.gz" &&\ 
echo done