# SAIGEQTL mapping pipeline
# Beomjin Jang
# Corrected for Snakemake execution and TypeError

import glob
import pandas as pd
import os

# --- Configuration ---
# These would typically be in a separate config.yaml file
# config = {
#     "outFolder": "results/",
#     "dataCode": "Mic",
#     "R_version": "R/4.2.0",
#     "PLINK_version": "plink2/2.3",
#     "TABIX_version": "tabix/0.2.6",
#     "BCFTOOLS_version": "bcftools/1.9",
#     "SAIGEQTL_bin": "/path/to/SAIGEQTL.sif",
#     "sampleKey": "path/to/sample_key.txt",
#     "phenoMeta": "path/to/phenotype_metadata.tsv",
#     "bfile": "path/to/genotype_prefix",
#     "phenotypes": "path/to/pheno_feature.tsv",
#     "covarColList": "sex,age,PC1,PC2,PC3,PC4,pct_mt,InferredCov1",
#     "sampleCovarColList": "sex,age,PC1,PC2,PC3,PC4",
#     "sampleIDColinphenoFile": "participant_id",
#     "offsetCol": "log_total_n",
#     "covariateFile": "path/to/pheno_cov_PEER.txt",
#     "chr_type": "chr1",
#     "threads": 4
# }

dataCode = config["dataCode"]

print(" * SAIGEQTL-pipeline *")
print(" Beomjin Jang 2026 ")
print(" * Data code is : %s " % dataCode)

# --- Modules and Software ---
R_VERSION = config.get("R_version", "R/4.2.0")
PLINK_VERSION = config.get("PLINK_version", "plink2/2.3")
TABIX_VERSION = config.get("TABIX_version", "tabix/0.2.6")
BCFTOOLS_VERSION = config.get("BCFTOOLS_version", "bcftools/1.9")
SAIGEQTL_bin = config.get("SAIGEQTL_bin", "/path/to/saigeqtl.sif")

## MAF - minor allele frequency - default = 0.01
MAF = config.get("MAF", 0.01)

# --- Chunking Logic ---
chunk_factor = 15 # 

# Load chromosome weights to calculate the number of chunks per chromosome
chunk_df = pd.read_csv("scripts/gencode_v38_chr_chunk_weights.tsv", sep="\t") # 
chunk_df = chunk_df.assign(chunk=chunk_df['ceil'] * chunk_factor) # 
chunk_dict = chunk_df.set_index("chr").T.to_dict()

# Define the chromosomes to be processed in the workflow
chromosomes = [str(chrom) for chrom in config.get("chromosomes", ["chr" + str(i) for i in range(1, 23)])]

# Limit each chromosome's chunks to the number of features in phenoMeta.
meta = pd.read_csv(config["phenoMeta"], sep="\t", header=None,
                   names=["chr", "start", "end", "feature"], dtype={"chr": str})
feature_counts = meta.groupby("chr").size()
for chrom in chromosomes:
    weight_chrom = chrom if chrom.startswith("chr") else "chr" + chrom
    num_chunks = int(chunk_dict[weight_chrom]['chunk'])
    chunk_dict[chrom] = dict(chunk_dict[weight_chrom])
    chunk_dict[chrom]['chunk'] = min(num_chunks, int(feature_counts.get(chrom, 0)))

# --- Input/Output Configuration ---
sample_key = config["sampleKey"]
phenoMeta = config["phenoMeta"]
genofile = config["bfile"]
phenotypes = config["phenotypes"]
covarColList = config["covarColList"]
sampleCovarColList = config["sampleCovarColList"]
sampleIDColinphenoFile = config["sampleIDColinphenoFile"]
offsetCol = config["offsetCol"]

threads = config["threads"]

SAIGEQTL_folder = "results/" + dataCode + "/"
SAIGEQTL_tmp_folder = SAIGEQTL_folder + "SAIGEQTL_tmp/"
genotype_folder = SAIGEQTL_folder + "genotype/"
prefix = genotype_folder + dataCode

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

pval_columns = ["pval_main"]
top_assoc_files = [
    SAIGEQTL_folder + dataCode + "_top_assoc_" + col + ".tsv.gz"
    for col in pval_columns
]

step4_prefix = "step4_" + dataCode + "_"
step4_outputs = [
    SAIGEQTL_folder + step4_prefix + col + suffix
    for col in pval_columns
    for suffix in ("_egene.tsv", "_egene_genes.txt")
]

# --- Snakemake Rules ---
rule all:
    input:
        step4_outputs,
        top_assoc_files,
        SAIGEQTL_folder + dataCode + "_full_assoc.tsv.gz",
        SAIGEQTL_folder + dataCode + "_full_assoc.tsv.gz.tbi",

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
        participants = prefix + "_participants.txt",
        source_bed = genofile + ".bed",
        source_bim = genofile + ".bim",
        source_fam = genofile + ".fam"
    output:
        bed = prefix + "_genotypes.bed",
        bim = prefix + "_genotypes.bim",
        fam = prefix + "_genotypes.fam",
        step1bed = prefix + "_genotypes_step1.bed",
        step1bim = prefix + "_genotypes_step1.bim",
        step1fam = prefix + "_genotypes_step1.fam",
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
                   --maf 0.05 \
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
                   --maf 0.05 \
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
        step1bed = prefix + "_genotypes_step1.bed",
        step1bim = prefix + "_genotypes_step1.bim",
        step1fam = prefix + "_genotypes_step1.fam",
        metadata = phenoMeta,
        phenotype = phenotypes,
        runner = "scripts/run_SAIGEQTL.R",
    params:
        script="scripts/run_SAIGEQTL.R",
        phenoMeta=phenoMeta,
        geno_prefix = prefix + "_genotypes",
        phenotypes=phenotypes,
        covarColList=covarColList,
        sampleCovarColList=sampleCovarColList,
        sampleIDColinphenoFile=sampleIDColinphenoFile,
        offset_col = offsetCol,
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
            " --cov_cols '{params.covarColList}' "
            " --sample_cov_cols '{params.sampleCovarColList}' "
            " --sample_cols '{params.sampleIDColinphenoFile}' "
            " --offset_col '{params.offset_col}' "
            " --prefix {params.SAIGEQTL_tmp_folder} "
            " --threads {params.threads} "
            " --SAIGEQTL {SAIGEQTL_bin} "
            " -i {wildcards.CHUNK} " # 
            " -n {max_chunk} "
        )

# 4. Combine Step 2 results using the local script, after all chunks finish.
rule concat_step2_results:
    input:
        chunks = all_chunk_files,
        metadata = phenoMeta,
        script = "scripts/concat_step2_results.py"
    output:
        step2_concat = SAIGEQTL_folder + "step3_input.txt"
    params:
        folder = SAIGEQTL_tmp_folder + "step_02/",
        maf_min = config.get("acat_maf_min", 0.05),
        maf_max = config.get("acat_maf_max", 0.95),
    shell:
        "ml python;  python3 {input.script:q} "
        " --input-dir {params.folder:q} --output {output.step2_concat:q} "
        " --file-pattern 'SAIGEQTL_*.txt' "
        " --gene-regex '^SAIGEQTL_[^_]+_(?P<gene>.+)[.]txt$' "
        " --maf-min {params.maf_min} --maf-max {params.maf_max} "

# 5. Calculate gene-level p-values using the copied SAIGEQTL Step 3 script.
rule step3_gene_pvalue:
    input:
        txt = SAIGEQTL_folder + "step3_input.txt",
        script = "scripts/step3_gene_pvalue.R"
    output:
        SAIGEQTL_folder + "step3_longformat.txt"
    params:
        folder = SAIGEQTL_folder,
    shell:
        "ml {R_VERSION};  Rscript {input.script:q} "
        " --input {input.txt:q} --outdir {params.folder:q}"

# 6. Apply FDR using local Step 4; no dynamic context is required.
rule step4_get_egenes:
    input:
        txt = SAIGEQTL_folder + "step3_longformat.txt",
        script = "scripts/step4_get_egenes.R"
    output:
        step4_outputs
    params:
        folder = SAIGEQTL_folder,
        fdr = config.get("fdr", 0.05),
        prefix = step4_prefix,
    shell:
        "ml {R_VERSION};  Rscript {input.script:q} "
        " --input {input.txt:q} --outdir {params.folder:q} "
        " --prefix {params.prefix:q} --fdr {params.fdr}"

# 7. Collate top associations after Step 4.
rule collate_top_chrom:
    input:
        step4 = step4_outputs,
        step3 = SAIGEQTL_folder + "step3_longformat.txt",
        script = "scripts/collate_top_chrom.R"
    params:
        prefix = SAIGEQTL_tmp_folder,
        output_prefix = SAIGEQTL_folder + dataCode + "_top_assoc",
        step4_prefix = SAIGEQTL_folder + step4_prefix
    output:
        top_assoc_files
    run:
        shell(
            "ml {R_VERSION}; "
            "ml {TABIX_VERSION}; "
            "Rscript {input.script:q} "
            " --output_file {params.output_prefix:q} "
            " --prefix {params.prefix:q} "
            " --step3 {input.step3:q} "
            " --step4_prefix {params.step4_prefix:q} "
        )

# 8. Collate full SAIGEQTL results per chromosome.
rule SAIGEQTL_collate:
    input:
        top = top_assoc_files,
        metadata = phenoMeta,
        script = "scripts/collate_SAIGEQTL.R"
    output:
        gz = SAIGEQTL_folder + "{chrom}_full_assoc.tsv.gz",
        header = temp(SAIGEQTL_folder + "{chrom}_full_assoc.tsv_header")
    params:
        prefix = SAIGEQTL_tmp_folder
    run:
        shell(
            "ml {R_VERSION}; "
            "ml {TABIX_VERSION}; "
            "Rscript {input.script:q} "
            " --prefix {params.prefix:q} "
            " --chrom {wildcards.chrom} "
            " --pheno_meta {input.metadata:q} "
            " --output_file {output.gz:q} "
        )

# 9. Merge, sort, bgzip and index all chromosome results.
rule fullCollate:
    input:
        assoc = expand(SAIGEQTL_folder + "{chrom}_full_assoc.tsv.gz",
                       chrom=[chrom for chrom in chromosomes if chunk_dict[chrom]['chunk'] > 0]),
        headers = expand(SAIGEQTL_folder + "{chrom}_full_assoc.tsv_header",
                         chrom=[chrom for chrom in chromosomes if chunk_dict[chrom]['chunk'] > 0])
    output:
        gz = SAIGEQTL_folder + dataCode + "_full_assoc.tsv.gz",
        tbi = SAIGEQTL_folder + dataCode + "_full_assoc.tsv.gz.tbi"
    params:
        tsv = SAIGEQTL_folder + dataCode + "_full_assoc.tsv",
        sorted = SAIGEQTL_folder + dataCode + "_full_assoc.tsv.sorted",
        header = lambda wildcards, input: input.headers[0]
    shell:
        """
        ml {TABIX_VERSION}
        echo "Sorting and collating all chromosome results"

        # Combine chromosome bodies (each input has no header).
        gzip -cd {input.assoc:q} > {params.tsv:q}

        # Sort by chromosome and position, then add one header.
        LC_ALL=C sort -k2,2V -k3,3n {params.tsv:q} > {params.sorted:q}
        cat {params.header:q} {params.sorted:q} | bgzip -c > {output.gz:q}

        # Index the final file.
        tabix -f -S 1 -s 2 -b 3 -e 3 {output.gz:q}

        # Remove uncompressed intermediate files.
        rm {params.tsv:q} {params.sorted:q}
        """
