#!/bin/bash   

### Job Requirements
#SBATCH --account=malus_rnaseq
#SBATCH -J get_mapping_stats
#SBATCH -o get_mapping_stats.out
#SBATCH -e get_mapping_stats.err

### Cluster Settings
#SBATCH --partition=fast 
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

output_file="mapping_trimmed.txt"

# Header for the output file
echo -e "Sample\tUniquely Mapped Reads\tAverage mapped length\tUnmapped Reads\tReads Mapped to Multiple Loci\tReads Mapped to Too Many Loci" > "$output_file"

for sample in $(find "/shared/ifbstor1/projects/malus_rnaseq/output/EXP3_trimmed.STARmapping/" -name "*.Log.final.out"); do
    sample_name=$(basename $(dirname "$sample"))
    
    # Extract information from log file 
    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" "$sample" | awk -F '|' '{print $2}' | xargs)
    average_mapped_length=$(grep "Average mapped length" "$sample" | awk -F '|' '{print $2}' | xargs)
    unmapped_reads=$(grep "Number of reads unmapped: too short" "$sample" | awk -F '|' '{print $2}' | xargs)
    reads_mapped_to_multiple_loci=$(grep "Number of reads mapped to multiple loci" "$sample" | awk -F '|' '{print $2}' | xargs)
    reads_mapped_to_too_many_loci=$(grep "Number of reads mapped to too many loci" "$sample" | awk -F '|' '{print $2}' | xargs)
    
    # Output the information to the file
    echo -e "$sample_name\t$uniquely_mapped_reads\t$average_mapped_length\t$unmapped_reads\t$reads_mapped_to_multiple_loci\t$reads_mapped_to_too_many_loci" >> "$output_file"
done

