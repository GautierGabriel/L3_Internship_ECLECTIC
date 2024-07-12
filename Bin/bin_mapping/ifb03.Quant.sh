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
#SBATCH --array=0-119%30
### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
module load salmon/1.10.2

mapfile list < sample_list_2.txt

sample=${list[$[$SLURM_ARRAY_TASK_ID]]}
sample=${sample::-1}
date '+%F %H:%M:%S'

mkdir /shared/projects/malus_rnaseq/output/EXP3_full.quant_align/$sample

echo $sample

salmon quant \
    -t /shared/projects/malus_rnaseq/input/genome/GDDH13_1-1_mrna.fasta \
    -l A -a /shared/projects/malus_rnaseq/output/EXP3_full.STARmapping/$sample/$sample".Aligned.toTranscriptome.out.bam" \
    -p 8 -o /shared/projects/malus_rnaseq/output/EXP3_full.quant_align/$sample

echo "DONE"
date '+%F %H:%M:%S'