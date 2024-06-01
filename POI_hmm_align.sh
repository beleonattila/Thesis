#!/bin/bash
PROJECT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"
HMM_profile="bla_class_A-NCBIFAM"

awk -v profile="$HMM_profile" '$1 == profile && !/^#/ {print $4}' $PROJECT_DIR/data/hmmer_results/bla_AMRFinder_scores.domtblout | awk '!seen[$0]++' > $PROJECT_DIR/data/hmmer_results/$HMM_profile"_seq_IDs.txt"
awk '/^>/{print s;if(s!="")print "";s="";print;next}{s=s""$0}END{print s}' $PROJECT_DIR/data/Prodigal/DTU_2018_1050_1_MG_HO_5.scaf.min1000.faa | grep -A 1 -E -f <(sed 's/$/[^[:alnum:]_]/' $PROJECT_DIR/data/hmmer_results/$HMM_profile"_seq_IDs.txt") | grep -v "^--$" > $PROJECT_DIR/data/hmmer_results/$HMM_profile"_seq_IDs.faa"


module purge
module load tools
module load hmmer/3.4

hmmalign --outformat pfam $PROJECT_DIR/apps/NCBIfam-AMRFinder/HMM/$HMM_profile".HMM" $PROJECT_DIR/data/hmmer_results/$HMM_profile"_seq_IDs.faa" > $PROJECT_DIR/data/hmmer_results/"alig_"$HMM_profile".pfam"

