import re
import math
import sys

# Read the input file path and output file path from Snakemake
input_file = sys.argv[1]
output_file = sys.argv[2]

def read_hmm_profile(filepath):
    # Open the HMM profile file
    with open(filepath, 'r') as file:
        parsing = False
        results = []
        amino_acid_order = 'ACDEFGHIKLMNPQRSTVWY'  # Standard order of amino acids
        compo_values = []

        # Iterate over each line in the file
        for line in file:
            stripped_line = line.strip()
            # Start parsing when the HMM header is encountered
            if stripped_line.startswith('HMM'):
                parsing = True
                # Skip the header lines until the actual matrix starts
                next(file)  # Skip the line directly after 'HMM'
                continue

            # Stop parsing when the footer is reached
            if stripped_line.startswith('//'):
                break

            # Only parse lines when within the HMM matrix body
            if parsing:
                # Remove any trailing newlines and split by whitespace
                parts = stripped_line.split()
                
                if 'COMPO' in stripped_line:
                    # Extract the values after the keyword 'COMPO'
                    compo_parts = stripped_line.split()
                    compo_index = compo_parts.index('COMPO') + 1
                    compo_values = list(map(float, compo_parts[compo_index:compo_index + 20]))
                    print(f"COMPO values: {compo_values}")  # Debugging statement
                    continue
                
                # Check if the line contains emission probabilities
                # Ensure the first part is a digit and represents the position number
                if parts and parts[0].isdigit():
                    position = int(parts[0])  # Position number in the profile
                    emissions = list(map(float, parts[1:21]))  # Assuming 20 amino acids follow the position number
                    
                    # Check each emission rate against the threshold
                    for aa_index, emission in enumerate(emissions):
                        aa = amino_acid_order[aa_index]  # Map index to amino acid
                        compo = compo_values[aa_index] if compo_values else 3
                        x = round(math.exp(-emission) - math.exp(-compo), 3)
                        if x > 0.33:
                            results.append((position, aa, emission, compo, x))
                            
        return results

def write_results_to_tsv(results, output_filepath):
    # Write results to a TSV file with a header
    with open(output_filepath, 'w') as file:
        file.write("Position\tAmino Acid\tEmission Rate\tCOMPO\tProbability Above Background\n") #header
        for position, aa, emission, compo, x in results:
            file.write(f"{position}\t{aa}\t{emission}\t{compo}\t{x}\n")

# Read the HMM profile and filter based on emission rate
filtered_results = read_hmm_profile(input_file)

# Write the results to a TSV file
write_results_to_tsv(filtered_results, output_file)
