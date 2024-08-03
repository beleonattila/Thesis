import sys

def parse_domtblout(domtblout_file, output_file):
    with open(domtblout_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_query_id = None
        current_target_id = None

        for line in infile:
            if line.startswith('#'):
                continue

            parts = line.strip().split()
            if len(parts) < 18:
                continue

            query_id = parts[3]
            target_id = parts[0]

            if query_id != current_query_id:
                # New query ID encountered
                current_query_id = query_id
                current_target_id = target_id
                outfile.write(line)
            elif query_id == current_query_id and target_id == current_target_id:
                # Same query ID and target ID as previous line
                outfile.write(line)
            else:
                # Different target ID, skip this line
                continue

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_domains.py <input_domtblout_file> <output_file>")
        sys.exit(1)
        
    domtblout_file = sys.argv[1]
    output_file = sys.argv[2]
    parse_domtblout(domtblout_file, output_file)
