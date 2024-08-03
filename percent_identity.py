def parse_clustal(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Remove header and empty lines
    lines = [line.strip() for line in lines if line.strip() and not line.startswith('CLUSTAL')]
    
    sequences = {}
    
    for line in lines:
        if not line[0].isspace():
            parts = line.split(maxsplit=1)
            key = parts[0]
            seq = parts[1] if len(parts) > 1 else ''
            
            # Skip lines of alignment annotation
            if set(seq).issubset(set("* . :")):
                continue
            
            if key not in sequences:
                sequences[key] = seq
            else:
                sequences[key] += seq
    
    return sequences

def calculate_percent_identity(sequences):
    if len(sequences) != 2:
        print(f"Error: Expected exactly 2 sequences, but found {len(sequences)}")
        for key, sequence in sequences.items():
            print(f"Sequence {key}: {sequence[:50]}... (length {len(sequence)})")
        raise ValueError("Exactly two sequences are required for percent identity calculation.")
    
    seq1, seq2 = sequences.values()
    if len(seq1) != len(seq2):
        print("Error: Sequences must be of the same length.")
        print(f"Length of sequence 1: {len(seq1)}")
        print(f"Length of sequence 2: {len(seq2)}")
        raise ValueError("Sequences must be of the same length.")
    
    total_length = len(seq1)
    match_count = 0
    aligned_residues = 0
    
    # Calculate percent identity for the full length and aligned residues
    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            aligned_residues += 1
            if a == b:
                match_count += 1
    
    if aligned_residues == 0:
        print("Error: No aligned residues found.")
        raise ValueError("Aligned length is zero.")
    
    percent_identity_noGap = (match_count / aligned_residues) * 100
    percent_identity = (match_count / total_length) * 100
    
    # Trim gaps from both ends for the trimmed alignment percent identity
    start_index = 0
    end_index = total_length
    
    # Find the first non-gap character
    while start_index < total_length and (seq1[start_index] == '-' or seq2[start_index] == '-'):
        start_index += 1
    
    # Find the last non-gap character
    while end_index > start_index and (seq1[end_index - 1] == '-' or seq2[end_index - 1] == '-'):
        end_index -= 1
    
    trimmed_seq1 = seq1[start_index:end_index]
    trimmed_seq2 = seq2[start_index:end_index]
    
    trimmed_match_count = sum(1 for a, b in zip(trimmed_seq1, trimmed_seq2) if a == b)
    trimmed_length = len(trimmed_seq1)
    
    if trimmed_length == 0:
        print("Error: Trimmed alignment length is zero.")
        raise ValueError("Trimmed alignment length is zero.")
    
    percent_identity_trimmed = (trimmed_match_count / trimmed_length) * 100
    
    return percent_identity, percent_identity_noGap, percent_identity_trimmed, aligned_residues


# Example usage
file_path = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/results/AlphaFold/AMR_best/DTU_2018_1050_1_MG_HO_5_NODE_38153_length_3794_cov_1.031906_3/msas/2zo4.aln'
sequences = parse_clustal(file_path)
percent_identity, percent_identity_noGap, percent_identity_trimmed, aligned_residues = calculate_percent_identity(sequences)
print(f"Percent Identity: {percent_identity:.2f}%")
print(f"Percent Identity without gaps: {percent_identity_noGap:.2f}%")
print(f"Percent Identity alignment length: {percent_identity_trimmed:.2f}%")
print(f"Number of Aligned Residues: {aligned_residues}")