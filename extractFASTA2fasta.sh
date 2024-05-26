#!/bin/bash
# This script extract individual fasta files from a combined one by a list as input.
# Define the input FASTA file and the file containing the list of desired IDs
/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/Prodigal
DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"

input_fasta="$DIR/data/Prodigal/DTU_2018_1050_1_MG_HO_5.scaf.min1000.faa"
ids_file="$DIR/data/hmmer_results/bla_class_A-NCBIFAM_IDs4alphafold2.txt"
output_folder="$DIR/data/fasta/bla_class_A-NCBIFAM2"

# Create the output folder if it doesn't exist
mkdir -p $output_folder

# Check if the input files exist
if [[ ! -f "$input_fasta" ]] || [[ ! -f "$ids_file" ]]; then
    echo "Error: Input file(s) not found."
    exit 1
fi
while read -r id; do
    awk -v id="$id" -v outdir="$output_folder" '
    /^>/ {
        printit = ($0 ~ ">"id" ");
    }
    printit {
        print > (outdir "/" id ".faa");
    }
    ' "$input_fasta"
done < "$ids_file"

echo "Selected FASTA entries have been extracted into $output_folder."