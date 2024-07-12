#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J trimmed_samtools_sort
#SBATCH -o trimmed_samtools_sort.out
#SBATCH -e trimmed_samtools_sort.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-119%40

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
module load samtools/1.18

mapfile list < sample_list_2.txt

sample=${list[$[$SLURM_ARRAY_TASK_ID]]}
sample=${sample::-1}
echo $sample
date '+%F %H:%M:%S'

mkdir /shared/projects/malus_rnaseq/output/EXP3_trimmed.samtools_sort/$sample

samtools sort -@ 8 /shared/projects/malus_rnaseq/output/EXP3_trimmed.STARmapping/$sample/$sample".Aligned.out.bam" -o /shared/projects/malus_rnaseq/output/EXP3_trimmed.samtools_sort/$sample/$sample".Aligned.out.sorted.bam"

echo "DONE"
date '+%F %H:%M:%S'