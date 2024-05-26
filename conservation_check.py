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

def read_domtblout(domtblout_file, hmm_profile):
    e_values = {}
    tlens = {}
    qlens = {}
    accs = {}
    align_lengths = {}
    
    with open(domtblout_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip comment lines
            fields = line.strip().split()
            query_id = fields[3]
            target_hmm = fields[0]
            if target_hmm == hmm_profile:
                e_value = float(fields[6])
                tlen = int(fields[2])
                qlen = int(fields[5])
                from_hmm = int(fields[15])
                to_hmm = int(fields[16])
                acc = float(fields[21])
                align_length = to_hmm - from_hmm + 1

                # Aggregate e-values, lengths, and accuracies
                if query_id not in e_values:
                    e_values[query_id] = e_value
                    tlens[query_id] = tlen
                    qlens[query_id] = qlen
                    align_lengths[query_id] = []
                    accs[query_id] = []
                
                align_lengths[query_id].append(align_length)
                accs[query_id].append(acc)

    # Calculate total alignment lengths and average accuracies
    total_align_lengths = {key: sum(value) for key, value in align_lengths.items()}
    average_accs = {key: sum(value) / len(value) if value else 0 for key, value in accs.items()}

    return e_values, tlens, qlens, total_align_lengths, average_accs



def compare_amino_acids(fasta_sequences, positions):
    results = {}
    c_AA_list = {'C', 'D', 'E', 'H', 'K', 'N', 'R', 'S', 'Y'}  # Catalytic AA
    for identifier, sequence in fasta_sequences.items():
        total = len(positions)
        match_count = 0
        catalytic_match_count = 0
        catalytic_total = 0
        amino_acids = {}
        
        for pos, expected_aas in positions.items():
            if pos < len(sequence):
                fasta_aa = sequence[pos]
                amino_acids[pos] = fasta_aa
                # Check if the amino acid at this position matches any of the expected amino acids
                if fasta_aa in expected_aas:
                    match_count += 1
                    # Check if the amino acid at this position is also a catalytic amino acid
                    if any(fasta_aa in c_AA_list for fasta_aa in expected_aas):
                        catalytic_match_count += 1
                # Count how many expected amino acids at this position are catalytic
                if any(e_aa in c_AA_list for e_aa in expected_aas):
                    catalytic_total += 1
            else:
                amino_acids[pos] = 'N/A'  # For positions not in the sequence

        # Calculate ratios
        overall_ratio = match_count / total if total > 0 else 0
        catalytic_ratio = catalytic_match_count / catalytic_total if catalytic_total > 0 else 0
        e_value = e_values.get(identifier, 'N/A')
        
        results[identifier] = (e_value, overall_ratio, catalytic_ratio, amino_acids)
    
    return results


def write_output(results, positions, tlens, qlens, total_align_lengths, average_accs, output_file):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        # Update the header with tlen and qlen
        header = ['Query_ID', 'E_value', 'tlen', 'qlen', 'Total_Align_Length', 'Align_Acc', 'Match_Ratio', 'Catalytic_Ratio'] + \
                 [f"{aa}{pos+1}" for pos, aa in positions.items()]
        writer.writerow(header)

        for identifier in results:
            e_value, overall_ratio, catalytic_ratio, aminos = results[identifier]
            tlen = tlens.get(identifier, 'N/A')  # Fetch tlen for the query
            qlen = qlens.get(identifier, 'N/A')  # Fetch qlen for the query
            total_align_length = total_align_lengths.get(identifier, 'N/A')
            average_acc = average_accs.get(identifier, 'N/A')

            # Format row and include all metrics
            row = [identifier, 
                   f"{e_value:.2e}" if isinstance(e_value, float) else e_value, 
                   tlen,
                   qlen,
                   total_align_length,
                   f"{average_acc:.2f}" if isinstance(average_acc, float) else average_acc,
                   f"{overall_ratio:.2f}",
                   f"{catalytic_ratio:.2f}"]  # Catalytic ratio added directly here
            row += [aminos.get(pos, '-') for pos in positions.keys()]  # Extend row with amino acids
            
            writer.writerow(row)



# File paths
hmm_profile = "bla_class_A-NCBIFAM"
domtblout_file_path = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_AMRFinder_scores.domtblout"
tsv_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis//apps/NCBIfam-AMRFinder/bla_class_A-NCBIFAM_filtered_emissions.tsv'
fasta_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_class_A-NCBIFAM_extracted_snippets_full_2.txt'
output_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_class_A-NCBIFAM_conservation_output_full_3.tsv'

# Process the files
e_values, tlens, qlens, total_align_lengths, average_accs = read_domtblout(domtblout_file_path, hmm_profile)
positions = read_tsv(tsv_file)
fasta_sequences = read_fasta(fasta_file)
results = compare_amino_acids(fasta_sequences, positions)
write_output(results, positions, tlens, qlens, total_align_lengths, average_accs, output_file)

