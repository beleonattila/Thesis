import csv

def read_tsv(tsv_file):
    positions = {}
    with open(tsv_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            pos = int(row['Position']) - 1  # Convert to 0-based index
            aa = row['Amino Acid']
            positions[pos] = aa
    return positions

def read_fasta(fasta_file):
    fasta_sequences = {}
    with open(fasta_file, 'r') as file:
        identifier = None
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # New sequence
                if identifier:
                    fasta_sequences[identifier] = sequence
                identifier = line[1:]  # Remove '>'
                sequence = ""
            else:
                sequence += line
        if identifier:
            fasta_sequences[identifier] = sequence  # Add the last parsed sequence
    return fasta_sequences

def compare_amino_acids(fasta_sequences, positions):
    results = {}
    for identifier, sequence in fasta_sequences.items():
        total = len(positions)
        match_count = 0
        amino_acids = {}
        for pos, expected_aa in positions.items():
            if pos < len(sequence):
                fasta_aa = sequence[pos]
                amino_acids[pos] = fasta_aa
                if fasta_aa == expected_aa:
                    match_count += 1
            else:
                amino_acids[pos] = 'N/A'  # For positions not in the sequence
        ratio = match_count / total if total > 0 else 0
        results[identifier] = (ratio, amino_acids)
    return results

def write_output(results, positions, output_file):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        header = ['Fasta ID', 'Match Ratio'] + [f"AA Pos {pos+1}" for pos in positions.keys()]
        writer.writerow(header)
        for identifier, (ratio, aminos) in results.items():
            row = [identifier, f"{ratio:.2f}"] + [aminos[pos] if pos in aminos else '-' for pos in positions.keys()]
            writer.writerow(row)

# File paths
tsv_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis//apps/NCBIfam-AMRFinder/bla_class_A-NCBIFAM_filtered_emissions.tsv'
fasta_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_class_A-NCBIFAM_extracted_snippets_full.txt'
output_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_class_A-NCBIFAM_conservation_output_full.tsv'

# Process the files
positions = read_tsv(tsv_file)
fasta_sequences = read_fasta(fasta_file)
results = compare_amino_acids(fasta_sequences, positions)
write_output(results, positions, output_file)

