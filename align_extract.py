def parse_stockholm_and_extract_snippets(input_file, output_file):
    with open(input_file, 'r') as file:
        sequences = {}
        rf_line = ""
        in_alignment_section = False

        # Read the file line by line
        for line in file:
            line = line.strip()
            if line.startswith("#=GC RF"):
                rf_line = line.split()[-1]  # Extract the RF annotation part
            elif line.startswith("//"):  # End of an alignment block
                in_alignment_section = False
            elif line.startswith("#"):  # Other headers not needed
                continue
            elif not line or line.startswith('//'):
                continue  # Skip empty lines or end of record
            else:
                parts = line.split()
                if len(parts) > 1:
                    seq_id = parts[0]
                    sequence = parts[1]
                    if seq_id in sequences:
                        sequences[seq_id] += sequence  # Append to existing sequence
                    else:
                        sequences[seq_id] = sequence  # Start a new sequence

    # Now extract the characters at positions marked by 'x' in the RF line
    x_positions = [i for i, char in enumerate(rf_line) if char == 'x']

    # Create snippets from sequences at positions marked by 'x'
    snippets = {seq_id: ''.join([seq[pos] for pos in x_positions if pos < len(seq)]) for seq_id, seq in sequences.items()}

    # Write the snippets to a file
    with open(output_file, 'w') as outfile:
        for seq_id, snippet in snippets.items():
            outfile.write(f">{seq_id}\n{snippet}\n")

# Example usage
input_filename = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/alig_bla_class_A-NCBIFAM.pfam'
output_filename = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_class_A-NCBIFAM_extracted_snippets_full.txt'
parse_stockholm_and_extract_snippets(input_filename, output_filename)
