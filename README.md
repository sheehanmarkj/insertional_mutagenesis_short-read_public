# insertional\_mutagenesis\_short-read

## Short read viral integration processing pipeline
##### This pipeline is intended to take fastq files from Illumina paired-end target-enrichment sequencing (TES) and test for integration between a host and viral vector.
Job scripts are written for LSF job scheduler and will need adjusting for other systems

## Required input

- Sample fastq files, with files named \{sample\}\_1.fastq.gz and \{sample\}\_2.fastq.gz
- Host fasta, bwa indexed
- ./viral\_genome/viral\_genome.fa: Vector genome fasta, bwa indexed
- Directory structure as contained in this repo
- Optionally, ./data/secondary/meta/Meta.csv metadata file

## Running analysis
```{bash}

parent=/path/to/parent
project=project

cd ${parent}/${project}

my_step0=${parent}/${project}/functions/step0.viralign_mask.sh
$my_step0 \
-i ${parent}/${project}/data/raw_fastq \
-e fastq.gz \
-r paired \
-v ${parent}/${project}/viral_genome/viral_genome.fa

my_step1=${parent}/${project}/functions/step1.trimalign.sh
$my_step1 \
-g cyno \
-i ${parent}/${project}/data/raw_fastq \
-e fastq.gz \
-r paired \
-o ${parent}/${project}/data/trimalign

my_step2=${parent}/${project}/functions/step2.ms.extract_update.sh
${my_step2} \
-i ${parent}/${project}/data/trimalign \
-g cyno \
-r paired

my_step3=${parent}/${project}/functions/step3.merge_update.sh
${my_step3} \
-g cyno \
-i ${parent}/${project}/data/trimalign \
-f ${parent}/${project}/data/raw_fastq \
-r paired


```

``` {r}

# Folders and files
parent <- "/path-to-parent/"
project <- "project"
source(file.path(function_dir,"functions/cigarParse.R"))
source(file.path(function_dir,"functions/loadIS.R"))

meta.data <-read_csv(file.path(parent,project,"/secondary/meta/Meta.csv"))

# Generating secondary files
loadIS(meta = meta.data,
       base.dir = paste0(parent,"/",project,"/secondary/"), # directory of secondary analysis
       pattern.is = "sorted.bam$", # suffix pattern of is bam files
       set = project, # any string name you want to give to this project
       do.is = FALSE,
       do.vir.origin = FALSE,
       do.read.info = TRUE)

```


## Results

* project >
  * data >
	* secondary >
  	* is/
	* reads\_info/
	* vir\_orig/
	* project\_read.rds: R object with candidate IS reads and info 


### Proceed to figure generation and other secondary analyses

