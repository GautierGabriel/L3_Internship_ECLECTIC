#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J full_sortmeRNA
#SBATCH -o full_sortmeRNA.out
#SBATCH -e full_sortmeRNA.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
module load star/2.7.11a



STAR --runMode genomeGenerate \
--genomeDir /shared/projects/malus_rnaseq/input/genome/STARindex \
--genomeSAindexNbases 13 \
--genomeFastaFiles /shared/projects/malus_rnaseq/input/genome/GDDH13_1-1_formatted.fasta \
--sjdbGTFfile /shared/projects/malus_rnaseq/input/genome/gene_models_20170612.gtf \
