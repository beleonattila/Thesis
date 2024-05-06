def check_positions(tsv_file, alignment_positions):
    matched_positions = 0
    with open(tsv_file, 'r') as file:
        next(file)  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                pos = int(parts[0])  # Assuming the position is the first column
                if pos in alignment_positions:
                    matched_positions += 1
    return matched_positions

# Example call
tsv_file = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis//apps/NCBIfam-AMRFinder/bla_class_B_core_filtered_emissions.tsv'
alignment_positions = 'alignment_positions.txt'
matched_count = check_positions(tsv_file, alignment_positions)
print(f"Number of matching positions: {matched_count}")
