#!/bin/bash
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=cge -A cge
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N prokka_test
### Output files (comment out the next 2 lines to get the job name used instead)
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 1 hour)
#PBS -l walltime=24:00:00
##

module purge
module load tools
module load perl/5.36.1
module load ncbi-blast/2.15.0+
module load prokka/1.14.5

cd $PBS_O_WORKDIR

raw_dir="/home/projects/cge/data/projects/other/dist_decay_sewage/flankophile/input"
PROJECT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"
sample=DTU_2018_1050_1_MG_HO_5.scaf.min1000.fa

prokka --outdir $PROJECT_DIR/data/prokka_result --force --centre cge --compliant --cpus 0 --metagenome parallel --gnu --prefix PROKKA $raw_dir/$sample
