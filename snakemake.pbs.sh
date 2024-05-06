#!/bin/bash
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=cge -A cge
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N snakemake_job
### Output files (comment out the next 2 lines to get the job name used instead)
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 1 hour)
#PBS -l walltime=48:00:00
##

module purge
module load tools
module load snakemake/8.4.2

cd $PBS_O_WORKDIR

snakemake --unlock
snakemake --cores 40 --max-status-checks-per-second 0.001 --rerun-incomplete
snakemake --dag | dot -Tpng > dag.png
