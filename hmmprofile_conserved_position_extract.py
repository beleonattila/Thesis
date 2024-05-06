import re

def read_hmm_profile(filepath):
    # Open the HMM profile file
    with open(filepath, 'r') as file:
        parsing = False
        results = []
        amino_acid_order = 'ACDEFGHIKLMNPQRSTVWY'  # Standard order of amino acids
        # AA_list  = {'C', 'D', 'E', 'H', 'K', 'N', 'R', 'S', 'Y'}  # Catalithic AA
        AA_list = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

        # Iterate over each line in the file
        for line in file:
            # Start parsing when the HMM header is encountered
            if line.startswith('HMM'):
                parsing = True
                # Skip the header lines until the actual matrix starts
                next(file)  # Skip the line directly after 'HMM'
                next(file)  # Skip the amino acid header line
                continue

            # Stop parsing when the footer is reached
            if line.startswith('//'):
                break

            # Only parse lines when within the HMM matrix body
            if parsing:
                # Remove any trailing newlines and split by whitespace
                parts = line.strip().split()
                
                # Check if the line contains emission probabilities
                # Ensure the first part is a digit and represents the position number
                if parts and parts[0].isdigit():
                    position = int(parts[0])  # Position number in the profile
                    emissions = list(map(float, parts[1:21]))  # Assuming 20 amino acids follow the position number
                    
                    # Check each emission rate against the threshold
                    for aa_index, emission in enumerate(emissions):
                        aa = amino_acid_order[aa_index]  # Map index to amino acid
                        if emission < 0.8 and aa in AA_list:
                            results.append((position, aa, emission))
                            
        return results

def write_results_to_tsv(results, output_filepath):
    # Write results to a TSV file with a header
    with open(output_filepath, 'w') as file:
        file.write("Position\tAmino Acid\tEmission Rate\n")  # Writing the header
        for position, aa, emission in results:
            file.write(f"{position}\t{aa}\t{emission}\n")  # Writing each result as a new line


def write_results_to_tsv(results, output_filepath):
    # Write results to a TSV file with a header
    with open(output_filepath, 'w') as file:
        file.write("Position\tAmino Acid\tEmission Rate\n")  # Writing the header
        for position, aa, emission in results:
            file.write(f"{position}\t{aa}\t{emission}\n")  # Writing each result as a new line

# Specify the path to your HMM profile
hmm_profile_path = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis//apps/NCBIfam-AMRFinder/HMM/bla_class_A-NCBIFAM.HMM'

output_tsv_path = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis//apps/NCBIfam-AMRFinder/bla_class_A-NCBIFAM_filtered_emissions.tsv'  # Define the output TSV file path

# Read the HMM profile and filter based on emission rate
filtered_results = read_hmm_profile(hmm_profile_path)

# Write the results to a TSV file
write_results_to_tsv(filtered_results, output_tsv_path)

print("Results have been written to:", output_tsv_path)
