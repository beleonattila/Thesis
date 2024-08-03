# Define the path to the project folder
PROJECT_DIR = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis"

# Define the path to the input metagenomic assembly file
ASSEMBLY_PATH = "/home/projects/cge/data/projects/other/dist_decay_sewage/flankophile/input"
ASSEMBLY_NAME = "DTU_2018_1050_1_MG_HO_5.scaf.min1000.fa"
# ASSEMBLY_PATH = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/temp/Prodigal/test_assembly"
# ASSEMBLY_NAME = "test_assembly.fa"

# Define the output directory for hmmer
# HMM_OUTDIR = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results_test"

# Define the results directory
RESULT_DIR = "/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/results"


from datetime import datetime
import os
	
ASSEMBLY_BASE_NAME = os.path.splitext(ASSEMBLY_NAME)[0]

# Get current time with precision up to minutes and format it as a string
current_time = datetime.now().strftime("%Y-%m-%d_%H-%M")

def get_AMR_profiles(wildcards):
    checkpoint_output = checkpoints.extract_AMR_profiles.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        profiles = [line.strip() for line in f]
    return profiles

def get_ResFams_profiles(wildcards):
    checkpoint_output = checkpoints.extract_ResFams_profiles.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        profiles = [line.strip() for line in f]
    return profiles


rule all:
    input:
        AMR=f"{RESULT_DIR}/HMM/POI/AMR_combined_scores.tsv",
        ResFams=f"{RESULT_DIR}/HMM/POI/ResFams_conbined_scores.tsv"

rule panres_bla_ID_extract:
    input:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_annotations.tsv"
    output:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/panres_bla_ids.txt"
    shell:
        """
        python3 panres_annotation_beta_lactam_filter.py {input} {output}
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
        mkdir -p {PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal
        prodigal -i {input} -q -a {output}
        """

rule panres_multi_ORF_protein_id_list:
    input:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_bla_proteins.faa"
    output:
        f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/multi_domain_ids.txt"
    shell:
        r"""
        awk '/^>/ {{id = substr($1, 2); if (id ~ /_[^1]$/) {{print id}}}}' {input} | sed 's/_[0-9]*$//' | sort | uniq > {output}
        """

rule panres_multi_ORF_protein_filter:
    input:
        fasta=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_bla_proteins.faa",
        multi_domain_ids=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/multi_domain_ids.txt"
    output:
        filtered_fasta=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_bla_proteins_filtered.faa",
        temp_file=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/no_endings.faa"
    shell:
        r"""
        # Remove endings from the original FASTA file
        sed -E 's/(>[^ ]+)_([0-9]+)/\1/' {input.fasta} > {output.temp_file}
        
        # Remove sequences with the base IDs from multi_domain_ids.txt
        awk 'BEGIN {{while ((getline < "{input.multi_domain_ids}") > 0) ids[$1]=1}} /^>/ {{flag=0; id=substr($1, 2); if (ids[id]) flag=1}} !flag' {output.temp_file} > {output.filtered_fasta}
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
        fasta_file=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_bla_proteins_filtered.faa",
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
        fasta_file=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_bla_proteins_filtered.faa",
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

rule PanRes_best_hit_filter:
    input:
        ResFams=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_hmm_ResFams_scores.domtblout",
        AMR=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/PanRes_hmm_AMRFinder_scores.domtblout"
    output:
        ResFams=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_ResFams_best_hit.domtblout",
        AMR=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_AMR_bla_best_hit.domtblout"
    shell:
        """
		sort -k4,4 -k7,7g {input.ResFams} | awk '!seen[$4]++' > {output.ResFams}
        sort -k4,4 -k7,7g {input.AMR} | awk '!seen[$4]++' > {output.AMR}
        """

rule PanRes_hmm_id_extract:
    input:
        ResFams=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_ResFams_best_hit.domtblout",
        AMR=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_AMR_bla_best_hit.domtblout"
    output:
        ResFams=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_ResFams_bla_HMM_ids.txt",
        AMR=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_AMR_bla_HMM_ids.txt"
    shell:
        """
        awk '/^#/ {{next}} {{print $1}}' {input.ResFams} | sort | uniq > {output.ResFams}
		awk '/^#/ {{next}} {{print $1}}' {input.AMR} | sort | uniq > {output.AMR}
        """

rule sample_prodigal:
    input:
        f"{ASSEMBLY_PATH}/{ASSEMBLY_NAME}"
    output:
        f"{RESULT_DIR}/Prodigal/{ASSEMBLY_BASE_NAME}.faa"
    shell:
        """
        module load tools
        module load prodigal/2.6.3
        mkdir -p {RESULT_DIR}/Prodigal
        prodigal -i {input} -p meta -a {output}
        """

rule sample_partial_AA_sequence_filter:
    input:
        f"{RESULT_DIR}/Prodigal/{ASSEMBLY_BASE_NAME}.faa"
    output:
        f"{RESULT_DIR}/Prodigal/{ASSEMBLY_BASE_NAME}_complete_only.faa"
    shell:
        """
        # Remove partials and '*' (stop codon) simbol at the end of sequences
        awk '/^>/ {{ out=(match($0, /partial=00/) ? "{output}" : "/dev/null") }} {{ print >> out }}' {input}
        """

rule sample_hmmscan_ResFams:
    input:
        fasta_file=f"{RESULT_DIR}/Prodigal/{ASSEMBLY_BASE_NAME}_complete_only.faa",
        profile=f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm",
        isPressed=f"{PROJECT_DIR}/apps/ResFams/Resfams-full.hmm/Resfams-full.hmm.h3m"
    output:
        score=f"{RESULT_DIR}/HMM/hmm_ResFams_scores.domtblout"
    resources:
        cores=18
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmscan --cpu {resources.cores} --domtblout {output.score} {input.profile} {input.fasta_file}
        """

rule sample_hmmscan_AMRFinder:
    input:
        fasta_file=f"{RESULT_DIR}/Prodigal/{ASSEMBLY_BASE_NAME}_complete_only.faa",
        profile=f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm",
        isPressed=f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/AMRFinder-full.hmm/AMRFinder-full.hmm.h3m"
    output:
        score=f"{RESULT_DIR}/HMM/hmm_AMRFinder_scores.domtblout"
    resources:
        cores=18
    shell:
        """
        module load tools
        module load hmmer/3.4
        hmmscan --cpu {resources.cores} --domtblout {output.score} {input.profile} {input.fasta_file}
        """

rule best_hit_filter_sample_AMR_scores:
    input:
        f"{RESULT_DIR}/HMM/hmm_AMRFinder_scores.domtblout"
    output:
        f"{RESULT_DIR}/HMM/best_hit_AMRFinder_scores.domtblout"
    shell:
        """
        python3 filter_best_hits.py {input} {output}
        """
		
rule best_hit_filter_sample_ResFams_scores:
    input:
        f"{RESULT_DIR}/HMM/hmm_ResFams_scores.domtblout"
    output:
        f"{RESULT_DIR}/HMM/best_hit_ResFams_scores.domtblout"
    shell:
        """
        python3 filter_best_hits.py {input} {output}
        """

rule bla_filter_sample_AMR_score:
    input:
        AMR_score=f"{RESULT_DIR}/HMM/best_hit_AMRFinder_scores.domtblout",
        AMR_ids=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_AMR_bla_HMM_ids.txt"
    output:
        f"{RESULT_DIR}/HMM/bla_AMRFinder_scores.domtblout"
    shell:
        """
        grep -Fwf {input.AMR_ids} {input.AMR_score} > {output}
        """
		
rule bla_filter_sample_ResFams_score:
    input:
        ResFams_score=f"{RESULT_DIR}/HMM/best_hit_ResFams_scores.domtblout",
        ResFams_ids=f"{PROJECT_DIR}/apps/PanRes_v1_0_1/PanRes_prodigal/Panres_ResFams_bla_HMM_ids.txt"
    output:
        f"{RESULT_DIR}/HMM/bla_ResFams_scores.domtblout"
    shell:
        """
        grep -Fwf {input.ResFams_ids} {input.ResFams_score} > {output}
        """

rule bla_eValue_threshold_sample_AMR_score:
    input:
        f"{RESULT_DIR}/HMM/bla_AMRFinder_scores.domtblout"
    output:
        f"{RESULT_DIR}/HMM/POI_AMRFinder_scores.domtblout"
    shell:
        """
        awk 'BEGIN {{FS=" "}} !/^#/ && $7 > 1e-7 {{print $0}}' {input} > {output}
        """

rule bla_eValue_threshold_sample_ResFams_score:
    input:
        ResFams_score=f"{RESULT_DIR}/HMM/bla_ResFams_scores.domtblout"
    output:
        ResFams_score=f"{RESULT_DIR}/HMM/POI_ResFams_scores.domtblout"
    shell:
        """
        awk 'BEGIN {{FS=" "}} !/^#/ && $7 > 1e-7 {{print $0}}' {input} > {output}
        """

checkpoint extract_AMR_profiles:
    input:
        f"{RESULT_DIR}/HMM/POI_AMRFinder_scores.domtblout"
    output:
        "AMR_profiles.txt"
    run:
        profiles = set()
        with open(input[0]) as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.split()
                    profile = fields[0]
                    profiles.add(profile)
        with open(output[0], 'w') as out:
            for profile in profiles:
                out.write(profile + '\n')

checkpoint extract_ResFams_profiles:
    input:
        f"{RESULT_DIR}/HMM/POI_ResFams_scores.domtblout"
    output:
        "ResFams_profiles.txt"
    run:
        profiles = set()
        with open(input[0]) as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.split()
                    profile = fields[0]
                    profiles.add(profile)
        with open(output[0], 'w') as out:
            for profile in profiles:
                out.write(profile + '\n')

rule extract_seq_ids_amr:
    input:
        domtblout = f"{RESULT_DIR}/HMM/POI_AMRFinder_scores.domtblout"
    output:
        seq_ids = f"{RESULT_DIR}/HMM/POI/seq_IDs/AMR/{{profile}}_seq_IDs.txt"
    shell:
        """
        awk -v profile="{wildcards.profile}" '$1 == profile && !/^#/ {{print $4}}' {input.domtblout} | awk '!seen[$0]++' > {output.seq_ids}
        """

rule extract_seq_ids_resfams:
    input:
        domtblout = f"{RESULT_DIR}/HMM/POI_ResFams_scores.domtblout"
    output:
        seq_ids = f"{RESULT_DIR}/HMM/POI/seq_IDs/ResFams/{{profile}}_seq_IDs.txt"
    shell:
        """
        awk -v profile={wildcards.profile} '$1 == profile && !/^#/ {{print $4}}' {input.domtblout} | awk '!seen[$0]++' > {output.seq_ids}
        """
		
rule extract_sequences_amr:
    input:
        seq_ids = f"{RESULT_DIR}/HMM/POI/seq_IDs/AMR/{{profile}}_seq_IDs.txt",
        fasta = f"{RESULT_DIR}/Prodigal/{ASSEMBLY_BASE_NAME}_complete_only.faa"
    output:
        sequences = f"{RESULT_DIR}/HMM/POI/fasta/AMR/{{profile}}_seq_IDs.faa"
    shell:
        """
        awk '/^>/{{print s;if(s!="")print "";s="";print;next}}{{s=s""$0}}END{{print s}}' {input.fasta} | grep -A 1 -E -f <(sed 's/$/[^[:alnum:]_]/' {input.seq_ids}) | grep -v "^--$" > {output.sequences}
        """

rule extract_sequences_resfams:
    input:
        seq_ids = f"{RESULT_DIR}/HMM/POI/seq_IDs/ResFams/{{profile}}_seq_IDs.txt",
        fasta = f"{RESULT_DIR}/Prodigal/{ASSEMBLY_BASE_NAME}_complete_only.faa"
    output:
        sequences = f"{RESULT_DIR}/HMM/POI/fasta/ResFams/{{profile}}_seq_IDs.faa"
    shell:
        """
        awk '/^>/{{print s;if(s!="")print "";s="";print;next}}{{s=s""$0}}END{{print s}}' {input.fasta} | grep -A 1 -E -f <(sed 's/$/[^[:alnum:]_]/' {input.seq_ids}) | grep -v "^--$" > {output.sequences}
        """

rule align_sequences_amr:
    input:
        hmm = f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/HMM/{{profile}}.HMM",
        sequences = f"{RESULT_DIR}/HMM/POI/fasta/AMR/{{profile}}_seq_IDs.faa"
    output:
        alignment = f"{RESULT_DIR}/HMM/POI/alignment/AMR/{{profile}}.pfam"
    shell:
        """
        module purge
        module load tools
        module load hmmer/3.4
        hmmalign --outformat pfam {input.hmm} {input.sequences} > {output.alignment}
        """

rule align_sequences_resfams:
    input:
        hmm = f"{PROJECT_DIR}/apps/ResFams/HMM/{{profile}}.hmm",
        sequences = f"{RESULT_DIR}/HMM/POI/fasta/ResFams/{{profile}}_seq_IDs.faa"
    output:
        alignment = f"{RESULT_DIR}/HMM/POI/alignment/ResFams/{{profile}}.pfam"
    shell:
        """
        module purge
        module load tools
        module load hmmer/3.4
        hmmalign --outformat pfam {input.hmm} {input.sequences} > {output.alignment}
        """

rule alignment_profile_trim_AMR:
    input:
        f"{RESULT_DIR}/HMM/POI/alignment/AMR/{{profile}}.pfam"
    output:
        f"{RESULT_DIR}/HMM/POI/alignment/AMR/{{profile}}_trim.pfam"
    shell:
        """
        python3 align_extract.py {input} {output}
        """

rule alignment_profile_trim_ResFams:
    input:
        f"{RESULT_DIR}/HMM/POI/alignment/ResFams/{{profile}}.pfam"
    output:
        f"{RESULT_DIR}/HMM/POI/alignment/ResFams/{{profile}}_trim.pfam"
    shell:
        """
        python3 align_extract.py {input} {output}
        """

rule HMM_profile_conserved_site_extraction_AMR:
    input:
        f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/HMM/{{profile}}.HMM"
    output:
        f"{RESULT_DIR}/HMM/POI/conserved_sites/AMR/{{profile}}.tsv"
    shell:
        """
        python3 hmmprofile_conserved_position_extract.py {input} {output}
        """

rule HMM_profile_conserved_site_extraction_ResFams:
    input:
        f"{PROJECT_DIR}/apps/ResFams/HMM/{{profile}}.hmm"
    output:
        f"{RESULT_DIR}/HMM/POI/conserved_sites/ResFams/{{profile}}.tsv"
    shell:
        """
        python3 hmmprofile_conserved_position_extract.py {input} {output}
        """

rule create_conservation_score_table_AMR:
    input:
        domtblout_file = f"{RESULT_DIR}/HMM/POI_AMRFinder_scores.domtblout",
        conserved_pos_tsv = f"{RESULT_DIR}/HMM/POI/conserved_sites/AMR/{{profile}}.tsv",
        alignment_file = f"{RESULT_DIR}/HMM/POI/alignment/AMR/{{profile}}_trim.pfam"
    output:
        f"{RESULT_DIR}/HMM/POI/conservation_score/AMR/{{profile}}.tsv"
    params:
        manual_annotation_folder = f"{PROJECT_DIR}/apps/NCBIfam-AMRFinder/active_sites"
    shell:
        """
        python3 conservation_check.py {input.domtblout_file} {input.conserved_pos_tsv} {input.alignment_file} {output} {params.manual_annotation_folder}
        """

rule create_conservation_score_table_ResFams:
    input:
        domtblout_file = f"{RESULT_DIR}/HMM/POI_ResFams_scores.domtblout",
        conserved_pos_tsv = f"{RESULT_DIR}/HMM/POI/conserved_sites/ResFams/{{profile}}.tsv",
        alignment_file = f"{RESULT_DIR}/HMM/POI/alignment/ResFams/{{profile}}_trim.pfam"
    output:
        f"{RESULT_DIR}/HMM/POI/conservation_score/ResFams/{{profile}}.tsv"
    params:
        manual_annotation_folder = f"{PROJECT_DIR}/apps/ResFams/active_sites"
    shell:
        """
        python3 conservation_check.py {input.domtblout_file} {input.conserved_pos_tsv} {input.alignment_file} {output} {params.manual_annotation_folder}
        """

rule combine_conservation_score_table_AMR:
    input:
        expand(f"{RESULT_DIR}/HMM/POI/conservation_score/AMR/{{profile}}.tsv", profile=get_AMR_profiles)
    output:
        f"{RESULT_DIR}/HMM/POI/AMR_combined_scores.tsv"
    shell:
        """
        # Extract and concatenate the first 13 columns of all input TSV files
        head -n 1 {input[0]} | awk '{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}}' > {output}
        for file in {input}; do
            tail -n +2 "$file" | awk '{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}}' >> {output}
        done
        """

rule combine_conservation_score_table_ResFams:
    input:
        expand(f"{RESULT_DIR}/HMM/POI/conservation_score/ResFams/{{profile}}.tsv", profile=get_ResFams_profiles)
    output:
        f"{RESULT_DIR}/HMM/POI/ResFams_conbined_scores.tsv"
    shell:
        """
        # Extract and concatenate the first 13 columns of all input TSV files
        head -n 1 {input[0]} | awk '{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}}' > {output}
        for file in {input}; do
            tail -n +2 "$file" | awk '{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}}' >> {output}
        done
        """
