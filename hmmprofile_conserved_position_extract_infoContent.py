import math

def read_hmm_profile(filepath, threshold):
    def calculate_information_content(frequencies, num_sequences):
        total_information = 0
        for freq in frequencies:
            if freq > 0:
                total_information += freq * math.log2(freq)
        R_i = math.log2(20) + total_information
        if num_sequences > 1:
            R_i_prime = R_i - ((20 - 1) / (2 * math.log(2) * num_sequences))
        else:
            R_i_prime = R_i
        return max(R_i_prime, 0)  # Ensure the information content is non-negative

    results = []
    amino_acid_order = 'ACDEFGHIKLMNPQRSTVWY'  # Standard order of amino acids

    with open(filepath, 'r') as file:
        parsing = False
        num_sequences = 0

        for line in file:
            if line.startswith('HMM'):
                parsing = True
                while True:
                    line = next(file).strip()
                    if line.startswith('NSEQ'):
                        num_sequences = int(line.split()[-1])
                    if line.startswith('HMM'):
                        next(file)  # Skip the amino acid header line
                        break
                continue

            if line.startswith('//'):
                break

            if parsing:
                parts = line.strip().split()
                if parts and parts[0].isdigit():
                    try:
                        position = int(parts[0])
                        emissions = list(map(float, parts[1:21]))

                        sum_emissions = sum(emissions)
                        if sum_emissions > 0:
                            frequencies = [emission / sum_emissions for emission in emissions]

                            information_content = calculate_information_content(frequencies, num_sequences)
                            
                            for aa_index, freq in enumerate(frequencies):
                                letter_height = freq * information_content * math.log2(20)  # Scale by max possible entropy
                                if letter_height >= threshold:
                                    aa = amino_acid_order[aa_index]
                                    results.append((position, aa, letter_height))
                    except ValueError:
                        continue  # Skip lines that cannot be processed correctly

    return results

def write_results_to_tsv(results, output_filepath):
    with open(output_filepath, 'w') as file:
        file.write("Position\tAmino Acid\tLetter Height\n")
        for position, aa, letter_height in results:
            file.write(f"{position}\t{aa}\t{letter_height:.4f}\n")

threshold = 0  # Set your threshold value here

# Specify the path to your HMM profile
hmm_profile_path = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/NCBIfam-AMRFinder/HMM/bla_class_A-NCBIFAM.HMM'

output_tsv_path = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/NCBIfam-AMRFinder/bla_class_A_main-NCBIFAM_filtered_emissions_new2.tsv'  # Define the output TSV file path

# Read the HMM profile and filter based on emission rate
filtered_results = read_hmm_profile(hmm_profile_path, threshold)

# Write the results to a TSV file
write_results_to_tsv(filtered_results, output_tsv_path)

print("Results have been written to:", output_tsv_path)
