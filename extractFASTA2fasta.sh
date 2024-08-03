#!/bin/bash
# This script extract individual fasta files from a combined one by a list as input.
# Define the input FASTA file and the file containing the list of desired IDs
DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"

input_fasta="/home/projects/cge/data/projects/other/dist_decay_sewage/flankophile/input/DTU_2018_1050_1_MG_HO_5.scaf.min1000.fa"
ids_file="$DIR/temp/Prodigal/test_extracted_list.txt"
output_folder="$DIR/temp/Prodigal/test_assembly"

# Create the output folder if it doesn't exist
mkdir -p $output_folder

# Check if the input files exist
if [[ ! -f "$input_fasta" ]] || [[ ! -f "$ids_file" ]]; then
    echo "Error: Input file(s) not found."
    exit 1
fi
# Read each ID from the ids file
while read -r id; do
    echo "Processing ID: $id"
    # Use awk to extract the sequences and remove '*' characters
    awk -v id="$id" -v outdir="$output_folder" '
    BEGIN {
        printit = 0;
        outfname = "";
    }
    /^>/ {
        if (printit && outfname != "") close(outfname);
        printit = ($0 ~ ">"id"($|\\s|\\|)");  # Match the ID followed by end of line, space, or pipe
        if (printit) {
            outfname = outdir "/" id ".faa";
            print "Writing to " outfname > "/dev/stderr";
        }
    }
    printit {
        if ($0 !~ /^>/) gsub(/\*/, "");
        print >> outfname;
    }
    END {
        if (printit && outfname != "") close(outfname);
    }
    ' "$input_fasta"
done < "$ids_file"

echo "Selected FASTA entries have been extracted into $output_folder."