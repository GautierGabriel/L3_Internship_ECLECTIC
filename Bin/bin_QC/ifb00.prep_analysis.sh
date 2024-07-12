#!/bin/bash
#!/bin/bash   

### Job Requirments
#SBATCH --account=malus_rnaseq
#SBATCH -J F23A43000179
#SBATCH -o F23A43000179.out
#SBATCH -e F23A43000179.err

###  Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

### Email
#SBATCH --mail-user=gabriel.gautier@ens-lyon.fr
#SBATCH --mail-type=ALL

#Creating the sample list
ls /shared/projects/malus_rnaseq/input/EXP3 > sample_list.txt