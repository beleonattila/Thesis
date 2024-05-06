#!/bin/bash
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=cge -A cge
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N hmmer_ResFams_PanRes
### Output files (comment out the next 2 lines to get the job name used instead)
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 1 hour)
#PBS -l walltime=3:00:00
##

module purge
module load tools
module load hmmer/3.4


cd $PBS_O_WORKDIR

PROJECT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"
OUT_DIR="PanRes_scan"

### mkdir $PROJECT_DIR/data/$OUT_DIR

hmmscan --cpu 40 --domtblout $PROJECT_DIR/data/$OUT_DIR/PanRes_ResFams_scores.domtblout $PROJECT_DIR/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm  $PROJECT_DIR/apps/PanRes_v1_0_1/panres_aa.fa
### hmmscan --cpu 40 --domtblout $PROJECT_DIR/data/$OUT_DIR/sample1__AMRFinder_scores.domtblout $PROJECT_DIR/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm $PROJECT_DIR/data/prokka_results/PROKKA.faa
### hmmscan --cpu 40 --domtblout $PROJECT_DIR/data/$OUT_DIR/sample1_ResFams_scores.domtblout $PROJECT_DIR/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm $PROJECT_DIR/data/prokka_results/PROKKA.faa
