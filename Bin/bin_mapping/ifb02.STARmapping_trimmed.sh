#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J trimmed_STARmapping
#SBATCH -o trimmed_STARmapping.out
#SBATCH -e trimmed_STARmapping.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --array=0-13%15
### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
module load star/2.7.11a

mapfile list < sample_list_3.txt

sample=${list[$[$SLURM_ARRAY_TASK_ID]]}
sample=${sample::-1}
echo $sample
date '+%F %H:%M:%S'

mkdir /shared/projects/malus_rnaseq/output/EXP3_trimmed.STARmapping/$sample

STAR --runThreadN 16 \
--genomeDir /shared/projects/malus_rnaseq/input/genome/STARindex/ \
--readFilesIn /shared/projects/malus_rnaseq/output/EXP3_trimmed.sortmerna/$sample/"str_fpf."$sample"_fwd.fq.gz" /shared/projects/malus_rnaseq/output/EXP3_trimmed.sortmerna/$sample/"str_fpf."$sample"_rev.fq.gz" \
--outFileNamePrefix /shared/projects/malus_rnaseq/output/EXP3_trimmed.STARmapping/$sample/$sample. \
--quantMode TranscriptomeSAM \
--twopassMode Basic \
--outSAMtype BAM Unsorted \
--readFilesCommand zcat \
--outSAMattributes NH HI AS NM MD \
--quantTranscriptomeBan Singleend 

echo "DONE"
date '+%F %H:%M:%S'