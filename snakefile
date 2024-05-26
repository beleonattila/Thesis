# Define the path to the project folder
PROJECT_DIR = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"

# Define the path to the input metagenomic assembly file
ASSEMBLY_PATH = "/home/projects/cge/data/projects/other/dist_decay_sewage/flankophile/input"
ASSEMBLY_NAME = "DTU_2018_1050_1_MG_HO_5.scaf.min1000.fa"
# ASSEMBLY = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/testAssembly.fsa"

# Define the output directory for hmmer
HMM_OUTDIR = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results"

from datetime import datetime
import os

ASSEMBLY_BASE_NAME = os.path.splitext(ASSEMBLY_NAME)[0]

# Get current time with precision up to minutes and format it as a string
current_time = datetime.now().strftime("%Y-%m-%d_%H-%M")

rule all:
    input:
        f"{HMM_OUTDIR}/bla_ResFams_scores.domtblout",
        f"{HMM_OUTDIR}/bla_AMRFinder_scores.domtblout"

rule panres_bla_ID_extract:
    input:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_annotations.tsv"
    output:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_bla_ids.txt"
    shell:
        """
        awk -F'\t' '$2 == "class" && $3 == "beta_lactam" {{print $1}}' {input} | sort -u > {output}
        """

rule PanRes_fasta_bla_filter:
    input:
        fasta=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_genes.fa",
        bla_ids=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_bla_ids.txt"
    output:
        fasta=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_bla_genes.fa"
    shell:
        """
		seqtk subseq {input.fasta} {input.bla_ids} > {output.fasta}
        """
		
rule panres_protein_translation:
    input:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_bla_genes.fa"
    output:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_bla_proteins.faa"
    shell:
        """
        module load tools
        module load prodigal/2.6.3
        mkdir -p f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal"
        prodigal -i {input} -p meta -a {output}
        """

rule hmmpress_ResFams:
    input:
        hmm=f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm"
    output:
        f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm.h3p",
        f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm.h3m",
        f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm.h3i",
        f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm.h3f"
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmpress {input.hmm}
        """
        
rule hmmpress_AMRFinder:
    input:
        hmm=f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm"
    output:
        f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm.h3p",
        f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm.h3m",
        f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm.h3i",
        f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm.h3f"
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmpress {input.hmm}
        """

rule hmmscan_PanRes_ResFams:
    input:
        fasta_file=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_bla_proteins.faa",
        isPressed=f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm.h3m"
    output:
        score=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_hmm_ResFams_scores.domtblout"
    resources:
        cores=18
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmscan --cpu {resources.cores} --domtblout {output.score} {PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm {input.fasta_file}
        """

rule hmmscan_PanRes_AMRFinder:
    input:
        fasta_file=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_bla_proteins.faa",
        isPressed=f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm.h3m"
    output:
        score=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_hmm_AMRFinder_scores.domtblout"
    resources:
        cores=18
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmscan --cpu {resources.cores} --domtblout {output.score} {PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm {input.fasta_file}
        """

rule PanRes_hmm_id_extract:
    input:
        ResFams=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_hmm_ResFams_scores.domtblout",
        AMR=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_hmm_AMRFinder_scores.domtblout"
    output:
        ResFams=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_ResFams_bla_HMM_ids.txt",
        AMR=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_AMR_bla_HMM_ids.txt"
    shell:
        """
        awk '$1 !~ /^#/ && $7 < 1e-10 {{print $1}}' {input.ResFams} | sort | uniq > {output.ResFams}
        awk '$1 !~ /^#/ && $7 < 1e-10 {{print $1}}' {input.AMR} | sort | uniq > {output.AMR}
        """

rule metagenom_prodigal:
    input:
        f"{ASSEMBLY_PATH}/{ASSEMBLY_NAME}"
    output:
        f"{PROJECT_DIR}/data/Prodigal/{ASSEMBLY_BASE_NAME}.faa"
    shell:
        """
        module load tools
        module load prodigal/2.6.3
        mkdir -p f"{PROJECT_DIR}/data/Prodigal"
        prodigal -i {input} -p meta -a {output}
        """

rule hmmscan_ResFams:
    input:
        fasta_file=f"{PROJECT_DIR}/data/Prodigal/{ASSEMBLY_BASE_NAME}.faa",
        isPressed=f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm.h3m"
    output:
        score=f"{HMM_OUTDIR}/hmm_ResFams_scores.domtblout"
    resources:
        cores=15
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmscan --cpu {resources.cores} --domtblout {output.score} {PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm {input.fasta_file}
        """

rule hmmscan_AMRFinder:
    input:
        fasta_file=f"{PROJECT_DIR}/data/Prodigal/{ASSEMBLY_BASE_NAME}.faa",
        isPressed=f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm.h3m"
    output:
        score=f"{HMM_OUTDIR}/hmm_AMRFinder_scores.domtblout"
    resources:
        cores=15
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmscan --cpu {resources.cores} --domtblout {output.score} {PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm {input.fasta_file}
        """


rule bla_filter_sample_hmm_scores:
    input:
        ResFams_score=f"{HMM_OUTDIR}/hmm_ResFams_scores.domtblout",
        AMR_score=f"{HMM_OUTDIR}/hmm_AMRFinder_scores.domtblout",
        ResFams_ids=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_ResFams_bla_HMM_ids.txt",
        AMR_ids=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_AMR_bla_HMM_ids.txt"
    output:
        ResFams_score=f"{HMM_OUTDIR}/bla_ResFams_scores.domtblout",
        AMR_score=f"{HMM_OUTDIR}/bla_AMRFinder_scores.domtblout"
    shell:
        """
        grep -Fwf {input.ResFams_ids} {input.ResFams_score} > {output.ResFams_score}
        grep -Fwf {input.AMR_ids} {input.AMR_score} > {output.AMR_score}
        """