# SAIGEQTL-pipeline

Beomjin Jang

Towfique Raj Lab

Mount Sinai, New York

2025

Snakemake pipeline for running [SAIGE-QTL](https://pubmed.ncbi.nlm.nih.gov/38798318/).

### Input data :

The inputs required by the user are a *sampleKey* matching RNA samples to DNA samples, a *bfile* containing the genotype information.

The *sample key* must have two columns named **participant_id** listing the names of the genotype IDs and **sample_id** listing the name of the phenotype sample IDs.

The *genotypes* must be in bfile format with all chromsomes in a single file with prefix 'chr'.

The *phenoMeta* must be a tab-separated table with feature information including chr, start and end.(chr, start, end, feature)

The *Phenotype* can be either space or tab-delimited with a header. It is required that the file contains one column for sample IDs and one more columns for the phenotype. It may contain columns for covariates.

The *covarColList* should be covariate columns in *Phenotype*

The *sampleCovarColList* should be covariate columns for samples in *Phenotype*, which have the same values for all cells from the same individual.

The *sampleIDColinphenoFile* should be genotype IDs column in *Phenotype* 

Please take a look at the format details on [SAIGE-QTL-doc](https://weizhou0.github.io/SAIGE-QTL-doc/).



## Config options

* threads 

Number of threads to use for runSAIGEQTL rule mclapply parallelization


