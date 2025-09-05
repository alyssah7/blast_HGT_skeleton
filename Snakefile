
## USER defined variables (in theory do not need to be touched)
### CSV GOES IN MY PROJECT DIRECTORY
# each reference genome must exist in the top level file named "{referenceID}.gbff" # GENOME.FASTA?

''' GLOBAL '''
# Global variables: In theory do not need to be changed
import sys
SCRIPTS_DIRECTORY = "/nfs/tamilab001/c3ddb-scratch-mit_lieberman/scripts" # for all the other functions like the few matlab ones from Aro
PROJECT_SCRIPTS_DIRECTORY = "./scripts" # this is bc so many scripts have the same/similar names in matlab and python versions -- hopefully this will keep them properly referenced
REF_GENOME_DIRECTORY = "./reference" # this is different for Vedanta? maybe not
CURRENT_DIRECTORY = os.getcwd()
sys.path.insert(0, SCRIPTS_DIRECTORY)
sys.path.insert(0, PROJECT_SCRIPTS_DIRECTORY)
from gus_helper_functions import *
from itertools import compress
from pathlib import Path

spls = "all_samples_rm_low_cov.csv" 

''' VARIABLES '''
# User defined variables: Make sure these are right before run!
flag="mapping" #options are 'all' (mapping+case), 'mapping', 'case', 'assembly', 'bracken'
maxFQ = -30 # purity threshold (from mapping quality) for including position-v3/genome_bowtie2.1.bt2", # already exists bc she already did bowtie with this ref and same ref for all samples

''' PRE-SNAKEMAKE '''
# Extract info from samples.csv
# Format: Path,Sample,FileName,Reference,Group,Outgroup
# Required fields for each mode:
    # all: Path,Sample,FileName,Reference,Group,Outgroup
    # mapping: Path,Sample,FileName,Reference,Outgroup
    # case: Path,Sample,Reference,Group,Outgroup
    # assembly: Path,Sample,FileName,Reference
    # bracken: Path,Sample,FileName,Reference
[PATH_ls, SAMPLE_ls, FILENAME_ls, REF_Genome_ls, GROUP_ls, OUTGROUP_ls] = read_samples_CSV(spls)
# Write sample_info.csv for each sample
split_samplesCSV(PATH_ls, SAMPLE_ls, FILENAME_ls, REF_Genome_ls, GROUP_ls, OUTGROUP_ls)
strain_ls = list(range(1,9))
UNIQ_GROUP_ls = set(GROUP_ls)

''' FUNCTIONS '''
# generate all unique pairs (q,s) with q < s
pairs = [(q, s) for i, q in enumerate(SAMPLE_ls) for s in SAMPLE_ls[i+1:]]
# Collect unique q values from pairs
Q_LIST = sorted({q for q, s in pairs})

input_all = []

# input_all.append("process_results/longest_contigs_all.csv")
input_all.append("process_results/filtered_blast_all.tsv")
input_all.append("process_results/final_filters/RE_filtered_blast_all.tsv")
input_all.append("genomad/all_done.txt")

rule all:
    input:
        # expand(
        #     "results/blast_{q}__{s}_out.tsv",
        #     zip,
        #     q=[q for q, s in pairs],
        #     s=[s for q, s in pairs]
        # )
        input_all

rule blastn:
    group: "group0",
    input:
        query_fasta   = "/orcd/data/tami/003/projects/abh_vedanta/mobile_phage_gainloss_followups/assemble_classify_nonVE303_receivers/Assembly/genomes/{q}.fasta",
        subject_fasta = "/orcd/data/tami/003/projects/abh_vedanta/mobile_phage_gainloss_followups/assemble_classify_nonVE303_receivers/Assembly/genomes/{s}.fasta",
    output:
        tsv = "results/blast_{q}__{s}_out.tsv",
    params:
        pident=99.9,
    conda:
        "envs/ncbi_blast.yaml",
    shell:
        """
        blastn -query {input.query_fasta} -subject {input.subject_fasta} \
               -perc_identity {params.pident} -outfmt 6 -out {output.tsv}
        """

rule combine_blast:
    input:
        lambda wc: expand(
            "results/blast_{q}__{s}_out.tsv",
            q=[wc.q],
            s=[s for q, s in pairs if q == wc.q]
        )
    output:
        "process_results/blast_{q}.tsv"
    shell:
        r"""
        mkdir -p process_results
        > {output}
        for f in {input}; do
            q=$(basename "$f" | sed -E 's/^blast_([^_]+)__([^_]+)_out\.tsv$/\1/')
            s=$(basename "$f" | sed -E 's/^blast_([^_]+)__([^_]+)_out\.tsv$/\2/')
            awk -v q="$q" -v s="$s" '{{print $0 "\t" q "\t" s}}' "$f" >> {output}
        done
        """


rule longest_contig_summary_all:
    input:
        fasta_files = expand(
            "/orcd/data/tami/003/projects/abh_vedanta/mobile_phage_gainloss_followups/assemble_classify_nonVE303_receivers/Assembly/genomes/{q}.fasta",
            q=SAMPLE_ls
        )
    output:
        csv = "process_results/longest_contigs_all.csv"
    run:
        import re
        from pathlib import Path

        out_path = Path(output.csv)
        out_path.parent.mkdir(exist_ok=True, parents=True)

        with open(out_path, "w") as out:
            out.write("sample,longest_length,longest_cov,total_length\n")
            for f in input.fasta_files:
                sample = Path(f).stem
                longest_len = 0
                coverage = 0.0
                total_length = 0
                with open(f) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            m_len = re.search(r"length_(\d+)", line)
                            m_cov = re.search(r"cov_([\d.]+)", line)
                            if m_len and m_cov:
                                l = int(m_len.group(1))
                                c = float(m_cov.group(1))
                                total_length += l
                                if l > longest_len:
                                    longest_len = l
                                    coverage = c
                out.write(f"{sample},{longest_len},{coverage},{total_length}\n")

rule filter_blast:
    input:
        "process_results/blast_{q}.tsv",
    params:
        length=10000,
        cov=3,
        contig_info="process_results/longest_contigs_all.csv",
        lenthresh=20000,
        kmercov=10,
        fraccov=0.2,
    output:
        "process_results/filtered_blast_{q}.tsv",
    conda:
        "envs/updated_py_for_snakemake.yaml",        
    shell:
        '''
        python scripts/filter_blast_hits.py -s {input} -l {params.length} -c {params.cov} -o {output} -g {params.contig_info} -e {params.lenthresh} -k {params.kmercov} -f {params.fraccov}

        '''

rule cat_filtered_blast:
    input:
        expand("process_results/filtered_blast_{q}.tsv",q=Q_LIST),
    output:
        "process_results/filtered_blast_all.tsv",   
    shell:
        '''
        cat {input} > {output}
        '''

rule refilter_blast_hits:
    input:
        contig_info="process_results/longest_contigs_all.csv",
        blast_results = "process_results/filtered_blast_all.tsv",
    params:
        outdir="process_results/final_filters",
        meta="2025-04_combined_isolate_metadata.csv",
        num_events_thresh=10, # isolate pairs with more than this many events will be deemed same-strain
    output:
        refiltered_tsv="process_results/final_filters/RE_filtered_blast_all.tsv",
        isolates_of_interest="process_results/final_filters/isolates_of_interest.txt",
    conda:
        "envs/py_for_GL.yaml",
    shell:
        '''
        mkdir -p {params.outdir}
        python scripts/second_filter_blast_hits.py -i {input.blast_results} -o {params.outdir} -m {params.meta} -l {params.num_events_thresh}
        '''

# Use a checkpoint so Snakemake can read the subsample file later
checkpoint parse_subsamples:
    input:
        "process_results/final_filters/isolates_of_interest.txt",
    output:
        directory("checkpoints/subsamples")
    run:
        subs = [l.strip() for l in open(input[0])]
        outdir = Path(output[0])
        outdir.mkdir(parents=True, exist_ok=True)
        with open(outdir / "samples.txt", "w") as out:
            out.write("\n".join(subs))

def get_subsamples(wildcards):
    outdir = Path(checkpoints.parse_subsamples.get().output[0])
    subs_file = outdir / "samples.txt"
    with open(subs_file) as f:
        return [l.strip() for l in f]

rule genomad:
    input:
        fasta = "/orcd/data/tami/003/projects/abh_vedanta/mobile_phage_gainloss_followups/assemble_classify_nonVE303_receivers/Assembly/genomes/{sample}.fasta",
    params:
        outdir="genomad/{sample}",
        db="/orcd/data/tami/003/databases/genomad_db_v1.6/genomad_db",
    conda:
        "genomad"
    output:
        "genomad/{sample}/genome_summary/{sample}_virus_genes.tsv"
    shell:
        """
        genomad end-to-end --cleanup {input.fasta} {params.outdir} {params.db}
        """

rule genomad_collect:
    input:
        lambda wildcards: expand(
            "genomad/{sample}/genome_summary/{sample}_virus_genes.tsv",
            sample=get_subsamples(wildcards)
        )
    output:
        "genomad/all_done.txt"
    shell:
        "touch {output}"

