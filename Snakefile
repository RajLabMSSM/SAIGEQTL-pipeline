# SAIGEQTL mapping pipeline
# Beomjin Jang
# Corrected for Snakemake execution and TypeError

import glob
import pandas as pd
import os

# --- Configuration ---
# These would typically be in a separate config.yaml file
# config = {
#     "mode": "cis",
#     "dataCode": "example_data",
#     "R_version": "R/4.2.0",
#     "plink2_version": "plink2/2.3",
#     "MAF": 0.01,
#     "sampleKey": "path/to/sample_key.txt",
#     "phenoMeta": "path/to/pheno_meta.txt",
#     "bfile": "path/to/genotypes",
#     "phenotypes": "path/to/phenotypes.txt",
#     "covarColList": "batch,sex,age,PC1,PC2,PC3,PEER1,PEER2,PEER3",
#     "sampleCovarColList": "batch,sex,age,PC1,PC2,PC3",
#     "sampleIDColinphenoFile": "participant_id",
#     "threads": 8
# }

mode = config["mode"]
dataCode = config["dataCode"]

print(" * SAIGEQTL-pipeline *")
print(" Beomjin Jang 2025 ")
print(" * Data code is : %s " % dataCode)
print(" * Mode selected is: %s" % mode)

# --- Modules and Software ---
R_VERSION = config.get("R_version", "R/4.2.0")
PLINK_VERSION = config.get("PLINK_version", "plink2/2.3")
TABIX_VERSION = config.get("TABIX_version", "tabix/0.2.6")
BCFTOOLS_VERSION = config.get("BCFTOOLS_version", "bcftools/1.9")
SAIGEQTL_bin = config.get("SAIGEQTL_bin", "/sc/arion/projects/ad-omics/data/software/SAIGEQTL/docker/saigeqtl.sif")

## MAF - minor allele frequency - default = 0.01
MAF = config.get("MAF", 0.01)

# --- Chunking Logic ---
chunk_factor = 15 # 

# Load chromosome weights to calculate the number of chunks per chromosome
chunk_df = pd.read_csv("scripts/gencode_v38_chr_chunk_weights.tsv", sep="\t") # 
chunk_df = chunk_df.assign(chunk=chunk_df['ceil'] * chunk_factor) # 
chunk_dict = chunk_df.set_index("chr").T.to_dict()

# Define the chromosomes to be processed in the workflow
# The original file was set to run only chr22 for testing
chromosomes = ["chr" + str(i) for i in range(1, 23)]

# --- Input/Output Configuration ---
sample_key = config["sampleKey"]
phenoMeta = config["phenoMeta"]
genofile = config["bfile"]
phenotypes = config["phenotypes"]
covarColList = config["covarColList"]
sampleCovarColList = config["sampleCovarColList"]
sampleIDColinphenoFile = config["sampleIDColinphenoFile"]
threads = config["threads"]

SAIGEQTL_folder = "results/" + dataCode + "/"
SAIGEQTL_tmp_folder = SAIGEQTL_folder + "SAIGEQTL_tmp/"
prefix = SAIGEQTL_folder + dataCode

# Create output directories if they don't exist
if not os.path.exists(SAIGEQTL_folder):
    os.mkdir(SAIGEQTL_folder)
if not os.path.exists(SAIGEQTL_tmp_folder):
    os.mkdir(SAIGEQTL_tmp_folder)

# --- Target File Generation ---
# Generate the full list of target files that Snakemake needs to create.
all_chunk_files = []
for chrom in chromosomes:
    if chrom in chunk_dict:
        num_chunks = chunk_dict[chrom]['chunk']
        # [FIX] Cast num_chunks to int() to prevent TypeError in range()
        for i in range(1, int(num_chunks) + 1):
            all_chunk_files.append(
                f"{SAIGEQTL_tmp_folder}SAIGEQTL_{chrom}_chunk_{i}"
            )

# --- Snakemake Rules ---
rule all:
    input:
        SAIGEQTL_folder + dataCode + "_full_assoc.tsv.gz"

#1. Make sample key 
rule getParticipants:
    output:
        txt = prefix + "_participants.txt"
    run:
        #sample_key = sample_key
        sk = pd.read_csv(sample_key, sep = "\t")
        participants = sk[["participant_id"]]
        participants.to_csv(output.txt, index = False, header = False, sep = "\t")

# 2 Convert VCF to plink, remove multi-allelic SNPs, blacklisted regions of genome, individuals not in sample key, variants with nan allele frequency, and maintain allele-order to prevent allele-flipping
rule prepare_PLINK:
    input:
        participants = prefix + "_participants.txt"
    output:
        bed = prefix + "_genotypes.bed",
        bim = prefix + "_genotypes.bim",
        fam = prefix + "_genotypes.fam",
        log = prefix + "_genotypes.log",
    params:
        stem = prefix + "_genotypes",
        blacklist = "scripts/Lifted_HighLDregion_hg38_RK_12_12_19.bed"
    run:
        bfile = genofile
        shell("""
            ml {PLINK_VERSION}
            plink2 --make-bed \
                   --output-chr M \
                   --max-alleles 2 \
                   --geno 0.1 \
                   --mac 20 \
                   --maf 0.01 \
                   --keep-allele-order \
                   --keep {input.participants} \
                   --exclude range {params.blacklist} \
                   --allow-extra-chr \
                   --bfile {bfile} \
                   --out {params.stem}_step1

            plink2 --make-bed \
                   --output-chr chrM \
                   --max-alleles 2 \
                   --geno 0.1 \
                   --maf 0.01 \
                   --max-maf 0.9975 \
                   --keep-allele-order \
                   --keep {input.participants} \
                   --exclude range {params.blacklist} \
                   --allow-extra-chr \
                   --bfile {bfile} \
                   --out {params.stem}

        #rm {params.stem}_step1.log
        """)

# 3. QTL mapping with SAIGEQTL
rule run_SAIGEQTL:
    input:
        bed = prefix + "_genotypes.bed",
        bim = prefix + "_genotypes.bim",
        fam = prefix + "_genotypes.fam",
    params:
        script="scripts/run_SAIGEQTL.R",
        phenoMeta=phenoMeta,
        geno_prefix = prefix + "_genotypes",
        phenotypes=phenotypes,
        covarColList=covarColList,
        sampleCovarColList=sampleCovarColList,
        sampleIDColinphenoFile=sampleIDColinphenoFile,
        SAIGEQTL_tmp_folder=SAIGEQTL_tmp_folder,
        threads=threads
    output:
        SAIGEQTL_tmp_folder + "SAIGEQTL_{chrom}_chunk_{CHUNK}"
    run:
        max_chunk = chunk_dict[wildcards.chrom]['chunk']
        shell(
            "ml {R_VERSION}; "
            "Rscript {params.script} "
            " --chrom {wildcards.chrom} "
            " --pheno_meta {params.phenoMeta} "
            " --geno_file {params.geno_prefix} " # 
            " --pheno_file {params.phenotypes} "
            " --cov_cols {params.covarColList} "
            " --sample_cov_cols {params.sampleCovarColList} "
            " --sample_cols {params.sampleIDColinphenoFile} "
            " --prefix {params.SAIGEQTL_tmp_folder} "
            " --threads {params.threads} "
            " --SAIGEQTL {SAIGEQTL_bin} "
            " -i {wildcards.CHUNK} " # 
            " -n {max_chunk} "
        )

# 4. Collate top QTL mapping with SAIGEQTL
rule collate_top_chrom:
    input:
        all_chunk_files
    params:
        script="scripts/collate_top_chrom.R",
        prefix=SAIGEQTL_tmp_folder
    output:
        SAIGEQTL_folder + dataCode + "_top_assoc.tsv.gz"
    run:
        shell(
            "ml {R_VERSION}; "
            "Rscript {params.script} "
            " --output_file {output} "
            " --prefix {params.prefix} "
        )


#5. Collate SAIGEQTL results per chromosome
rule SAIGEQTL_collate:
    input:
        SAIGEQTL_folder + dataCode + "_top_assoc.tsv.gz"
    output:
        SAIGEQTL_folder + "{chrom}_full_assoc.tsv.gz"
    params:
        script="scripts/collate_SAIGEQTL.R",
        prefix=SAIGEQTL_tmp_folder
    run:
        shell(
            "ml {R_VERSION}; "
            "Rscript {params.script} "
            " --prefix {params.prefix} "
            " --chrom {wildcards.chrom} "
            " --output_file {output}"
        )

#6. FullCollate SAIGEQTL results
rule fullCollate:
    input:
        expand(SAIGEQTL_folder + "{chrom}_full_assoc.tsv.gz", chrom=chromosomes)
    output:
        gz=SAIGEQTL_folder + dataCode + "_full_assoc.tsv.gz",
        tbi=SAIGEQTL_folder + dataCode + "_full_assoc.tsv.gz.tbi"
    params:
        tsv=SAIGEQTL_folder + dataCode + "_full_assoc.tsv",
        header=SAIGEQTL_folder + "chr1_full_assoc.tsv_header",
        prefix=SAIGEQTL_folder,
    shell:
        """
        set +o pipefail
        ml {BCFTOOLS_VERSION}
        ml {TABIX_VERSION}
        echo "Sorting and collating all chromosome results"

        # Write the header to the final output TSV
        cat {params.header} > {params.tsv}

        # All chromosome files append to the Final file
        zcat {input} >> {params.tsv}

        # Sort the concatenated file and append to the final output with the header
        sort -k2,2V -k3,3n {params.tsv} > {params.tsv}.sorted

        # Clean up the temporary concatenated file
        #rm {params.tsv}.sorted {params.header}
        #
        #echo "Compressing and indexing full associations"
        bgzip -f -c {params.tsv}.sorted > {output.gz}
        tabix -S 1 -s 2 -b 3 -e 3 {output.gz}
        """
