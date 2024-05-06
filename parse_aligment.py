import re

def parse_hmmscan_output(filename, output_filename):
    with open(filename, 'r') as file, open(output_filename, 'w') as outfile:
        for line in file:
            # Match lines that contain sequence alignments with positions
            if line.strip() and '==' not in line and any(char.isdigit() for char in line):
                parts = line.strip().split()
                # Ensure it's the sequence line with positions at the end
                if len(parts) > 2 and parts[0].isdigit() and parts[-1].isdigit():
                    query_start = int(parts[1])
                    query_seq = parts[2]
                    query_end = int(parts[-1])

                    # Write alignment info and iterate over each character in the sequence part
                    pos = query_start
                    for char in query_seq:
                        if char != '-':
                            outfile.write(f"Position {pos}: {char}\n")
                            pos += 1
                        else:
                            pos += 1

# Example usage
input_hmmscan_file = 'ali_test2_bla_B_core.txt'
output_positions_file = 'alignment_positions.txt'
parse_hmmscan_output(input_hmmscan_file, output_positions_file)
