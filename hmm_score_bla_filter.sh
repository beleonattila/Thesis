#!/bin/bash

PROJECT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"
RESFAMS_SCORES="hmm_ResFams_scores.domtblout"
AMRFINDER_SCORES="hmm_AMRFinder_scores.domtblout"

grep -Fwf $PROJECT_DIR/apps/ResFams/only_bla_ResFams_targets.txt $PROJECT_DIR/data/hmmer_results/$RESFAMS_SCORES > $PROJECT_DIR/data/hmmer_results/filtered_ResFams_scores.domtblout
grep -Fwf $PROJECT_DIR/apps/NCBIfam-AMRFinder/only_bla_AMRFinder_targets.txt $PROJECT_DIR/data/hmmer_results/$AMRFINDER_SCORES > $PROJECT_DIR/data/hmmer_results/filtered_AMRFinder_scores.domtblout
