#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J full_qualimap
#SBATCH -o full_qualimap.out
#SBATCH -e full_qualimap.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-119%30

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
module load qualimap/2.2.2b 

mapfile list < sample_list_2.txt

sample=${list[$[$SLURM_ARRAY_TASK_ID]]}
sample=${sample::-1}
echo $sample

date '+%F %H:%M:%S'

mkdir /shared/projects/malus_rnaseq/output/EXP3_full.qualimap/$sample

unset DISPLAY

qualimap rnaseq -bam /shared/projects/malus_rnaseq/output/EXP3_full.samtools_sort/$sample/$sample".Aligned.out.sorted.bam" -gtf /shared/projects/malus_rnaseq/input/genome/STARindex/gene_models_20170612.gtf -outdir /shared/projects/malus_rnaseq/output/EXP3_full.qualimap/$sample/ --paired 

