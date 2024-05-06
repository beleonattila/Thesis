#!/bin/bash

module purge
module load tools
module load hmmer/3.4

# Path to your combined HMM file
COMBINED_HMM="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm"

# Directory where you want to save extracted HMM profiles
OUTPUT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/NCBIfam-AMRFinder/AMRFinder-bla.hmm"
mkdir -p ${OUTPUT_DIR}

# Path to your file with unique target names
TARGET_NAMES="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/NCBIfam-AMRFinder/unique__AMRFinder_targets.txt"

# Loop through each target name and extract the corresponding HMM profile
while read TARGET; do
    hmmfetch ${COMBINED_HMM} "${TARGET}" > "${OUTPUT_DIR}/${TARGET}.hmm"
done < ${TARGET_NAMES}
