#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J EXP3_full_multi_qc
#SBATCH -o EXP3_full_multi_qc.out
#SBATCH -e EXP3_full_multi_qc.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
module load multiqc/1.13

multiqc -n /shared/projects/malus_rnaseq/output/EXP3_full_filter.multi.qc/EXP3_Apple_fastp_report.html \
/shared/projects/malus_rnaseq/output/EXP3_full_filter/*/

