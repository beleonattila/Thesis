import csv
import sys
import os

# Groups of interchangeable amino acids
aa_groups = {
    'D': {'E', 'D'},
    'E': {'E', 'D'},
    'H': {'H'},
    'K': {'R', 'K'},
    'N': {'N', 'Q'},
    'Q': {'N', 'Q'},
    'R': {'R', 'K'},
    'S': {'S', 'T', 'Y'},
    'T': {'S', 'T', 'Y'},
    'Y': {'S', 'T', 'Y'}
}


def read_tsv(conserved_pos_tsv):
    positions = {}
    with open(conserved_pos_tsv, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            pos = int(row['Position']) - 1  # Convert to 0-based index
            aa = row['Amino Acid']
            positions[pos] = aa
    return positions

def read_alignment(alignment_file):
    alignment_sequences = {}
    with open(alignment_file, 'r') as file:
        identifier = None
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # New sequence
                if identifier:
                    alignment_sequences[identifier] = sequence
                identifier = line[1:]  # Remove '>'
                sequence = ""
            else:
                sequence += line
        if identifier:
            alignment_sequences[identifier] = sequence  # Add the last parsed sequence
    return alignment_sequences

def read_domtblout(domtblout_file, hmm_profile):
    e_values = {}
    tlens = {}
    qlens = {}
    domains_info = {}
    
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
                domain_num = fields[10]
                domain_id = fields[9]
                align_length = to_hmm - from_hmm + 1

                if query_id not in e_values:
                    e_values[query_id] = e_value
                    tlens[query_id] = tlen
                    qlens[query_id] = qlen
                    domains_info[query_id] = []
                
                domains_info[query_id].append({
                    "acc": acc,
                    "align_length": align_length,
                    "domain_num": domain_num,
                    "domain_id": domain_id
                })

    return e_values, tlens, qlens, domains_info


def read_active_site(filename):
    active_site = {}
    if filename and os.path.exists(filename):  # Check if the file exists
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if line:
                    aa, pos = line[0], int(line[1:]) - 1  # Convert to 0-based index
                    active_site[pos] = aa
    return active_site

def compare_amino_acids(alignment_sequences, positions, active_site, aa_groups, e_values):
    results = {}
    c_AA_list = {'D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y'}  # Catalytic AA
    for identifier, sequence in alignment_sequences.items():
        total = len(positions)
        match_count = 0
        catalytic_match_count = 0
        catalytic_total = 0
        manual_match_count = 0
        flexible_match_count = 0
        flexible_match_count_active_site = 0
        manual_total = len(active_site)
        amino_acids = {}
        
        for pos, expected_aa in positions.items():
            if pos < len(sequence):
                alignment_aa = sequence[pos]
                amino_acids[pos] = alignment_aa
                
                # Check for exact match
                if alignment_aa == expected_aa:
                    match_count += 1
                    if alignment_aa in c_AA_list:
                        catalytic_match_count += 1
                        
                # Check if the expected amino acid is catalytic
                if expected_aa in c_AA_list:
                    catalytic_total += 1
                    
                    # Check for flexible matching in active site
                    if alignment_aa in aa_groups.get(expected_aa, {expected_aa}):
                        flexible_match_count += 1

        # Check for matches in the manually annotated active site
        for pos, expected_aa in active_site.items():
            if pos < len(sequence):
                if sequence[pos] == expected_aa:
                    manual_match_count += 1
                # Check for flexible matching in active site
                if sequence[pos] in aa_groups.get(expected_aa, {expected_aa}):
                    flexible_match_count_active_site += 1

        overall_ratio = match_count / total if total > 0 else 0
        catalytic_ratio = catalytic_match_count / catalytic_total if catalytic_total > 0 else 0
        interchangeable_catalytic_ratio = flexible_match_count / catalytic_total if catalytic_total > 0 else 'N/A'
        active_site_ratio = manual_match_count / manual_total if manual_total > 0 else 'N/A'
        interchangeable_active_site_ratio = flexible_match_count_active_site / manual_total if manual_total > 0 else 'N/A'
        
        results[identifier] = (e_values.get(identifier, 'N/A'), overall_ratio, catalytic_ratio, interchangeable_catalytic_ratio, active_site_ratio, interchangeable_active_site_ratio, amino_acids)
    
    return results


def write_output(results, positions, tlens, qlens, domains_info, output_file):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        # Update the header with tlen and qlen
        header = ['Query_ID', 'E_value', 'tlen', 'qlen', 'Align_Acc', 'Align_Length', 'Domain_Num', 'Domain_ID', 'Match_Ratio', 'Catalytic_Ratio','Interchangeable_Catalytic_Ratio' , 'Active_site_ratio', 'Interchangeable_active_site_ratio'] + \
                 [f"{aa}{pos+1}" for pos, aa in positions.items()]
        writer.writerow(header)

        for identifier in results:
            e_value, overall_ratio, catalytic_ratio, interchangeable_catalytic_ratio, active_site_ratio, interchangeable_active_site_ratio, amino_acids = results[identifier]
            tlen = tlens.get(identifier, 'N/A')  # Fetch tlen for the query
            qlen = qlens.get(identifier, 'N/A')  # Fetch qlen for the query

            # Write individual entries for each domain
            for domain_info in domains_info.get(identifier, []):
                row = [
                    identifier,
                    f"{e_value:.2e}" if isinstance(e_value, float) else e_value,
                    tlen,
                    qlen,
                    f"{domain_info['acc']:.2f}" if isinstance(domain_info['acc'], float) else domain_info['acc'],
                    domain_info['align_length'],
                    domain_info['domain_num'],
                    domain_info['domain_id'],
                    f"{overall_ratio:.2f}",
                    f"{catalytic_ratio:.2f}",  # Catalytic_aa ratio
                    f"{interchangeable_catalytic_ratio:.2f}" if isinstance(interchangeable_catalytic_ratio, float) else interchangeable_catalytic_ratio,
                    f"{active_site_ratio:.2f}" if isinstance(active_site_ratio, float) else active_site_ratio,  # Manually annotated active site ratio
                    f"{interchangeable_active_site_ratio:.2f}" if isinstance(interchangeable_active_site_ratio, float) else interchangeable_active_site_ratio  # Manually annotated active ratio with interchangeable aa_groups
                ]
                row += [amino_acids.get(pos, '-') for pos in positions.keys()]  # Extend row with amino acids
                
                writer.writerow(row)





# File paths
# active_site_file = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/NCBIfam-AMRFinder/bla_subclass_B1-NCBIFAM_active_site.txt"
# hmm_profile = "bla_subclass_B1-NCBIFAM"
# domtblout_file = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_AMRFinder_scores.domtblout"
# conserved_pos_tsv = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/NCBIfam-AMRFinder/bla_subclass_B1-NCBIFAM_filtered_emissions.tsv'
# alignment_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_subclass_B1-NCBIFAM_extracted_snippets_full.txt'
# output_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results/bla_subclass_B1-NCBIFAM_conservation_output_full_v4.tsv'


# Process the files
def main():
    # This part contains the code that should only run when the script is executed directly
    if len(sys.argv) < 5:
        print("Usage: python conservation_check.py <domtblout_file> <conserved_pos_tsv> <alignment_file> <output_file> [<active_site_file>]")
        sys.exit(1)
    
    domtblout_file = sys.argv[1]
    conserved_pos_tsv = sys.argv[2]
    alignment_file = sys.argv[3]
    output_file = sys.argv[4]
    manual_annotation_folder = sys.argv[5] if len(sys.argv) > 5 else None
    
    hmm_profile  = os.path.basename(conserved_pos_tsv).split('.')[0]
    if manual_annotation_folder:
        active_site_file = os.path.join(manual_annotation_folder, hmm_profile + "_active_site.txt")
        if not os.path.exists(active_site_file):
            active_site_file = None
    else:
        active_site_file = None

    positions = read_tsv(conserved_pos_tsv)
    alignment_sequences = read_alignment(alignment_file)
    e_values, tlens, qlens, domains_info = read_domtblout(domtblout_file, hmm_profile)
    active_site = read_active_site(active_site_file)
    results = compare_amino_acids(alignment_sequences, positions, active_site, aa_groups, e_values)
    write_output(results, positions, tlens, qlens, domains_info, output_file)

if __name__ == "__main__":
    main()