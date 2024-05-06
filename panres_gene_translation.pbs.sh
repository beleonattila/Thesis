#!/bin/bash
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=cge -A cge
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N prodigal_PanRes
### Output files (comment out the next 2 lines to get the job name used instead)
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=5
### Memory
#PBS -l mem=16gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 1 hour)
#PBS -l walltime=1:00:00
##

module purge
module load tools
module load prodigal/2.6.3


cd $PBS_O_WORKDIR

INPUT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/apps/PanRes_v1_0_1"
OUT_DIR=$INPUT_DIR/PanRes_prodigal

mkdir -p $OUT_DIR

prodigal -i $INPUT_DIR/panres_genes.fa -o $OUT_DIR/gene_coords.gbk -a $OUT_DIR/protein_translations.faa