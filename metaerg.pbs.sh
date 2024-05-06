#!/bin/bash
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=cge -A cge
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N MetaErg_test
### Output files (comment out the next 2 lines to get the job name used instead)
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 1 hour)
#PBS -l walltime=2:00:00
##

module purge
module load tools
module load anaconda2/4.4.0
module load minpath/1.4
module load jdk/22
module load minced/0.4.2
module load hmmer/3.4
module load diamond/2.1.8
module load aragorn/1.2.36
module load perl
module load ncbi-blast/2.15.0+
module load prodigal/2.6.3
module load signalp/6.0h
module load tmhmm/2.0c
module load metaerg/1.2.3


cd $PBS_O_WORKDIR

### raw_dir="/home/projects/cge/data/projects/other/dist_decay_sewage/flankophile/input"
PROJECT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"
### sample=DTU_2018_1050_1_MG_HO_5.scaf.min1000.fa

raw_dir="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data"
sample=testAssembly2.fsa

metaerg.pl --outdir $PROJECT_DIR/data/metaerg_result --force --cpus 40 --prefix MetaErg $raw_dir/$sample

