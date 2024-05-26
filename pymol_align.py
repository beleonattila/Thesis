import pymol

# File paths
prot_1 = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/alphafold/results/DTU_2018_1050_1_MG_HO_5_NODE_213565_length_1407_cov_1.423437_3/ranked_0.pdb"
prot_2 = '/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/PDB/class_A_hit_test/1n9b.pdb'

# Load protein structure files
pymol.cmd.load(prot_1, "obj1")
pymol.cmd.load(prot_2, "obj2")

# Perform sequence alignment
alignment_result = pymol.cmd.align("obj1", "obj2")

# Save aligned structures to a new PDB file
pymol.cmd.save("aligned_structures.pdb", "obj1 or obj2")

with open("rmsd_output.txt", "w") as f:
    f.write("RMSD: {}\n".format(alignment_result[0]))


# Quit PyMOL
pymol.cmd.quit()
