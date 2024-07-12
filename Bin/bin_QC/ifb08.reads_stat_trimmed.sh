#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J trimmed_reads_stat_fastp3
#SBATCH -o trimmed_reads_stat_fastp3.out
#SBATCH -e trimmed_reads_stat_fastp3.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

module purge
module load seqkit/2.1.0

mkdir -p /shared/projects/malus_rnaseq/output/EXP3_trimmed.reads_stat
cat /shared/projects/malus_rnaseq/bin/str_trimmed.list | \
xargs seqkit stats -j 10 -b -T > /shared/projects/malus_rnaseq/output/EXP3_trimmed.reads_stat/str_trimmed.stat

#cat /shared/projects/malus_rnaseq/bin/fastp_trimmed.list | \
#xargs seqkit stats -j 10 -b -T > /shared/projects/malus_rnaseq/output/EXP3_trimmed.reads_stat/fastp_trimmed.stat

