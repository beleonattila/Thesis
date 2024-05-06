#!/bin/bash

PROJECT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"
RESFAMS_SCORES="hmm_ResFams_scores.domtblout"
AMRFINDER_SCORES="hmm_AMRFinder_scores.domtblout"

awk '!seen[$4]++' $PROJECT_DIR/data/hmmer_results/$RESFAMS_SCORES > $PROJECT_DIR/data/hmmer_results/top_score_$RESFAMS_SCORES
awk '!seen[$4]++' $PROJECT_DIR/data/hmmer_results/$AMRFINDER_SCORES > $PROJECT_DIR/data/hmmer_results/top_score_$AMRFINDER_SCORES
