#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J qualimap_multi_qc
#SBATCH -o qualimap_multi_qc.out
#SBATCH -e qualimap_multi_qc.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

# modules loading
module purge
module load multiqc/1.13

multiqc -n /shared/projects/malus_rnaseq/output/EXP3_full/EXP3_full.multiqc_qualimap/EXP3_Apple_qualimap_report.html \
/shared/projects/malus_rnaseq/output/EXP3_full.qualimap/*/

