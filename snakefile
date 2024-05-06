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
        f"{HMM_OUTDIR}/POI_ResFams_scores.domtblout",
        f"{HMM_OUTDIR}/POI_AMRFinder_scores.domtblout"

rule panres_bla_filter:
    input:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_annotations.tsv"
    output:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_bla_ids.txt"
    shell:
        """
        awk -F'\t' '$2 == "class" && $3 == "beta_lactam" {{print $1}}' {input} | sort -u > {output}
        """

rule panres_protein_translation:
    input:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_genes.fa"
    output:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/protein_translations.faa"
    shell:
        """
        module load tools
        module load prodigal/2.6.3
        mkdir -p f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal
        prodigal -i {input} -a {output}
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
        fasta_file=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/protein_translations.faa",
        isPressed=f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm.h3m"
    output:
        score=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/hmm_ResFams_scores.domtblout"
    resources:
        cores=10
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmscan --cpu {resources.cores} --domtblout {output.score} {PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm {input.fasta_file}
        """

rule hmmscan_PanRes_AMRFinder:
    input:
        fasta_file=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/protein_translations.faa",
        isPressed=f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm.h3m"
    output:
        score=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/hmm_AMRFinder_scores.domtblout"
    resources:
        cores=10
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmscan --cpu {resources.cores} --domtblout {output.score} {PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm {input.fasta_file}
        """

rule PanRes_id_restore:
    input:
        ResFams_score=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/hmm_ResFams_scores.domtblout",
        AMR_score=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/hmm_AMRFinder_scores.domtblout"
    output:
        orig_id_ResFams_score=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/hmm_ResFams_orig_id_scores.domtblout",
        orig_id_AMR_score=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/hmm_AMRFinder_orig_id_scores.domtblout"
    shell:
        """
        awk '{{sub(/..$/, "", $4)}} 1' {input.ResFams_score} > {output.orig_id_ResFams_score}
        awk '{{sub(/..$/, "", $4)}} 1' {input.AMR_score} > {output.orig_id_AMR_score}
        """

rule PanRes_hmmScore_bla_filter:
    input:
        ResFams=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/hmm_ResFams_orig_id_scores.domtblout",
        AMR=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/hmm_AMRFinder_orig_id_scores.domtblout",
        bla_id=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_bla_ids.txt"
    output:
        ResFams=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/only_bla_ResFams_scores.domtblout",
        AMR=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/only_bla_AMR_scores.domtblout"
    shell:
        """
        grep -Fwf {input.bla_id} {input.ResFams} > {output.ResFams}
        grep -Fwf {input.bla_id} {input.AMR} > {output.AMR}
        """

rule PanRes_hmm_id_extract:
    input:
        ResFams=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/only_bla_ResFams_scores.domtblout",
        AMR=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/only_bla_AMR_scores.domtblout"
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

rule confident_hmm_hit_list:
    input:
        ResFams_score=f"{HMM_OUTDIR}/hmm_ResFams_scores.domtblout",
        AMR_score=f"{HMM_OUTDIR}/hmm_AMRFinder_scores.domtblout"
    output:
        ResFams=f"{HMM_OUTDIR}/ResFams_confident_hit_list.txt",
        AMR=f"{HMM_OUTDIR}/AMR_confident_hit_list.txt"
    shell:
        """
        awk '!seen[$4]++' {input.ResFams_score} | awk '!/^#/ && $7 < 1e-5 {{print $4}}' > {output.ResFams}
        awk '!seen[$4]++' {input.AMR_score} | awk '!/^#/ && $7 < 1e-5 {{print $4}}' > {output.AMR}
        """
        
rule remove_confident_hits_from_hmm_scores:
    input:
        ResFams_score=f"{HMM_OUTDIR}/hmm_ResFams_scores.domtblout",
        AMR_score=f"{HMM_OUTDIR}/hmm_AMRFinder_scores.domtblout",
        ResFams_conf_list=f"{HMM_OUTDIR}/ResFams_confident_hit_list.txt",
        AMR_conf_list=f"{HMM_OUTDIR}/AMR_confident_hit_list.txt"
    output:
        ResFams_score=f"{HMM_OUTDIR}/hmm_ResFams_low_confidence_scores.domtblout",
        AMR_score=f"{HMM_OUTDIR}/hmm_AMRFinder_low_confidence_scores.domtblout"
    shell:
        """
        awk 'FNR==NR {{ ids[$1]; next }} !($4 in ids)' {input.ResFams_conf_list} {input.ResFams_score} > {output.ResFams_score}
        awk 'FNR==NR {{ ids[$1]; next }} !($4 in ids)' {input.AMR_conf_list} {input.AMR_score} > {output.AMR_score}
        """

rule bla_filter_sample_hmm_scores:
    input:
        ResFams_score=f"{HMM_OUTDIR}/hmm_ResFams_low_confidence_scores.domtblout",
        AMR_score=f"{HMM_OUTDIR}/hmm_AMRFinder_low_confidence_scores.domtblout",
        ResFams_ids=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_ResFams_bla_HMM_ids.txt",
        AMR_ids=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_AMR_bla_HMM_ids.txt"
    output:
        ResFams_score=f"{HMM_OUTDIR}/POI_ResFams_scores.domtblout",
        AMR_score=f"{HMM_OUTDIR}/POI_AMRFinder_scores.domtblout"
    shell:
        """
        grep -Fwf {input.ResFams_ids} {input.ResFams_score} > {output.ResFams_score}
        grep -Fwf {input.AMR_ids} {input.AMR_score} > {output.AMR_score}
        """

rule hmm_frequency:
    input:
        ResFams_score=f"{HMM_OUTDIR}/POI_ResFams_scores.domtblout",
        AMR_score=f"{HMM_OUTDIR}/POI_AMRFinder_scores.domtblout"
    output:
        ResFams_freq=f"{HMM_OUTDIR}/POI_ResFams_freq.txt",
        AMR_freq=f"{HMM_OUTDIR}/POI_AMRFinder_freq.txt"
    shell:
        """
        awk '!/^#/ {print $1}' {input.AMR_score} | sort | uniq -c | sort > {output.AMR_freq}
        awk '!/^#/ {print $1}' {input.ResFams_score} | sort | uniq -c | sort > {output.ResFams_freq}
        """
