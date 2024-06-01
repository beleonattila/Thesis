import csv

# Groups of interchangeable amino acids
aa_groups = {
    'S': {'S', 'T', 'Y'},
    'T': {'S', 'T', 'Y'},
    'Y': {'S', 'T', 'Y'},
    'R': {'R', 'K'},
    'K': {'R', 'K'},
    'E': {'E', 'D'},
    'D': {'E', 'D'},
    'N': {'N', 'Q'},
    'Q': {'N', 'Q'}
}


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
    start_codon = {}  # Dictionary to store start info
    stop_codon = {}  # Dictionary to store stop info
    
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

                # Extract start and stop information from the query ID
                try:
                    partial = query_id.split('_partial=')[1]
                    start = partial[0]
                    stop = partial[1]
                except IndexError:
                    start = 'N/A'
                    stop = 'N/A'

                # Aggregate e-values, lengths, accuracies, start and stop info
                if query_id not in e_values:
                    e_values[query_id] = e_value
                    tlens[query_id] = tlen
                    qlens[query_id] = qlen
                    align_lengths[query_id] = []
                    accs[query_id] = []
                    start_codon[query_id] = start
                    stop_codon[query_id] = stop
                
                align_lengths[query_id].append(align_length)
                accs[query_id].append(acc)

    # Calculate total alignment lengths and average accuracies
    total_align_lengths = {key: sum(value) for key, value in align_lengths.items()}
    average_accs = {key: sum(value) / len(value) if value else 0 for key, value in accs.items()}

    return e_values, tlens, qlens, total_align_lengths, average_accs, start_codon, stop_codon


def read_active_site(filename):
    active_site = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                aa, pos = line[0], int(line[1:]) - 1  # Convert to 0-based index
                active_site[pos] = aa
    return active_site

def compare_amino_acids(fasta_sequences, positions, active_site, aa_groups, e_values):
    results = {}
    c_AA_list = {'C', 'D', 'E', 'H', 'K', 'N', 'R', 'S', 'Y'}  # Catalytic AA
    for identifier, sequence in fasta_sequences.items():
        total = len(positions)
        match_count = 0
        catalytic_match_count = 0
        catalytic_total = 0
        manual_match_count = 0
        flexible_match_count = 0
        manual_total = len(active_site)
        amino_acids = {}
        
        for pos, expected_aa in positions.items():
            if pos < len(sequence):
                fasta_aa = sequence[pos]
                amino_acids[pos] = fasta_aa
                
                # Check for exact match
                if fasta_aa == expected_aa:
                    match_count += 1
                    if fasta_aa in c_AA_list:
                        catalytic_match_count += 1

                # Check if the expected amino acid is catalytic
                if expected_aa in c_AA_list:
                    catalytic_total += 1

        # Check for matches in the manually annotated active site
        for pos, expected_aa in active_site.items():
            if pos < len(sequence):
                if sequence[pos] == expected_aa:
                    manual_match_count += 1
                # Check for flexible matching in active site
                if sequence[pos] in aa_groups.get(expected_aa, {expected_aa}):
                    flexible_match_count += 1

        overall_ratio = match_count / total if total > 0 else 0
        catalytic_ratio = catalytic_match_count / catalytic_total if catalytic_total > 0 else 0
        active_site_ratio = manual_match_count / manual_total if manual_total > 0 else 0
        interchangeable_active_site_ratio = flexible_match_count / manual_total if manual_total > 0 else 0
        e_value = e_values.get(identifier, 'N/A')
        
        results[identifier] = (e_value, overall_ratio, catalytic_ratio, active_site_ratio, interchangeable_active_site_ratio, amino_acids)
    
    return results


def write_output(results, positions, tlens, qlens, total_align_lengths, start_codon, stop_codon, average_accs, output_file):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        # Update the header with tlen and qlen
        header = ['Query_ID', 'E_value', 'tlen', 'qlen', 'Total_Align_Length', 'Start_codon', 'Stop_codon', 'Align_Acc', 'Match_Ratio', 'Catalytic_Ratio', 'Active_site_ratio', 'Interchangeable_active_site_ratio'] + \
                 [f"{aa}{pos+1}" for pos, aa in positions.items()]
        writer.writerow(header)

        for identifier in results:
            e_value, overall_ratio, catalytic_ratio, active_site_ratio, interchangeable_active_site_ratio, amino_acids = results[identifier]
            tlen = tlens.get(identifier, 'N/A')  # Fetch tlen for the query
            qlen = qlens.get(identifier, 'N/A')  # Fetch qlen for the query
            total_align_length = total_align_lengths.get(identifier, 'N/A')
            start = start_codon.get(identifier, 'N/A')
            stop = stop_codon.get(identifier, 'N/A')
            average_acc = average_accs.get(identifier, 'N/A')

            # Format row and include all metrics
            row = [identifier, 
                   f"{e_value:.2e}" if isinstance(e_value, float) else e_value,
                   tlen,
                   qlen,
                   total_align_length,
                   start,
                   stop,
                   f"{average_acc:.2f}" if isinstance(average_acc, float) else average_acc,
                   f"{overall_ratio:.2f}",
                   f"{catalytic_ratio:.2f}",  # Catalytic_aa ratio
                   f"{active_site_ratio:.2f}",  # Manually annotated active site ratio
                   f"{interchangeable_active_site_ratio:.2f}"]  # Manually annotated active ratio with interchangeable aa_groups
            row += [amino_acids.get(pos, '-') for pos in positions.keys()]  # Extend row with amino acids
            
            writer.writerow(row)



# File paths
active_site_path = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/NCBIfam-AMRFinder/bla_class_A-NCBIFAM_active_site.txt"
hmm_profile = "bla_class_A-NCBIFAM"
domtblout_file_path = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_AMRFinder_scores.domtblout"
tsv_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis//apps/NCBIfam-AMRFinder/bla_class_A-NCBIFAM_filtered_emissions.tsv'
fasta_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_class_A-NCBIFAM_extracted_snippets_full.txt'
output_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_class_A-NCBIFAM_conservation_output_full.tsv'

# Process the files
active_site = read_active_site(active_site_path)
e_values, tlens, qlens, total_align_lengths, average_accs, start_codon, stop_codon = read_domtblout(domtblout_file_path, hmm_profile)
positions = read_tsv(tsv_file)
fasta_sequences = read_fasta(fasta_file)
results = compare_amino_acids(fasta_sequences, positions, active_site, aa_groups, e_values)
write_output(results, positions, tlens, qlens, total_align_lengths, start_codon, stop_codon, average_accs, output_file)

