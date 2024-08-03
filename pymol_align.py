import pymol
from pymol import cmd

# Initialize PyMOL without GUI
pymol.finish_launching(['pymol', '-cq'])

# File paths
BASE_DIR = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/alphafold/"
prot_1 = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/alphafold/class_B1/DTU_2018_1050_1_MG_HO_5_NODE_89_length_76240_cov_2.248170_46_partial=00/ranked_0.pdb"
prot_2 = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/alphafold/class_B1/DTU_2018_1050_1_MG_HO_5_NODE_89_length_76240_cov_2.248170_46_partial=00/4ysk_B.pdb"

# Load protein structure files
cmd.load(prot_1, "obj1")
cmd.load(prot_2, "obj2")

# Perform sequence alignment
alignment_result = cmd.align("obj1", "obj2")

# Define dictionaries for residues
model1_residues = {}
model2_residues = {}

# Update global dictionaries
globals()['model1_residues'] = model1_residues
globals()['model2_residues'] = model2_residues

# Extract sequences and alignment information
def three_to_one(resn):
    code = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    return code.get(resn, 'X')

# Use PyMOL's iterate function to populate the dictionaries
cmd.iterate("obj1 and name CA", "model1_residues[resi] = resn", space=globals())
cmd.iterate("obj2 and name CA", "model2_residues[resi] = resn", space=globals())

# Initialize lists for aligned sequences
model1_sequence = []
model2_sequence = []

# Get alignment information
alignment = cmd.get_raw_alignment("obj1", "obj2")

# Build aligned sequences based on the alignment
for pair in alignment:
    resi1, resi2 = pair
    if resi1 in model1_residues:
        model1_sequence.append(three_to_one(model1_residues[resi1]))
    else:
        model1_sequence.append('-')

    if resi2 in model2_residues:
        model2_sequence.append(three_to_one(model2_residues[resi2]))
    else:
        model2_sequence.append('-')

# Join the sequences into strings
model1_sequence = ''.join(model1_sequence)
model2_sequence = ''.join(model2_sequence)

# Save the alignment in CLUSTAL format
clustal_file = "alignment.clustal"
with open(clustal_file, "w") as f:
    f.write("CLUSTAL\n\n")
    line_length = 60
    for i in range(0, len(model1_sequence), line_length):
        f.write("obj1     {}\n".format(model1_sequence[i:i+line_length]))
        f.write("obj2     {}\n".format(model2_sequence[i:i+line_length]))
        f.write("\n")

# Save RMSD result
with open("rmsd_output.txt", "w") as f:
    f.write("RMSD: {}\n".format(alignment_result[0]))

# Quit PyMOL
cmd.quit()