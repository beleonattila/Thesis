#!/bin/bash

PROJECT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"

hmmscan --notextw $PROJECT_DIR/apps/NCBIfam-AMRFinder/bla_class_B_core.hmm/bla_class_B_core-NCBIFAM.HMM $PROJECT_DIR/data/hmmer_results/test2_bla_B_core.fasta | ali_test2_bla_B_core.txt

# /home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/NCBIfam-AMRFinder