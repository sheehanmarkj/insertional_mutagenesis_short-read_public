#!/bin/bash
module load ib R/4.1.2
module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

#set -e
#set -o pipefail

# Creating a usage sentence
usage="To use this viral alignment for initial sample reads:
$(basename $0) [-s sample-name] [-i path-to-orig-fastq] [-e fastq|fq|fastq.gz|fq.gz] [-r single|paired] [-v path-to-viral-genome] [-w this is genewerk data] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi

while getopts 's:i:e:r:m:v:w' OPTION; do
        case "$OPTION" in
                s) sample="$OPTARG"
                printf -- "\033[32m The selected rerun sample is $OPTARG. \033[0m\n"
                ;;
                i) indir="$OPTARG"
                printf -- "\033[32m The original fastq files are at $OPTARG. \033[0m\n"
                ;;
                e) suffix="$OPTARG"
                printf -- "\033[32m The suffix of your files is $OPTARG. \033[0m\n"
                ;;
                r) reads="$OPTARG"
                printf -- "\033[32m Your readtype is $OPTARG end. \033[0m\n"
                ;;
                v) vir_genome="$OPTARG"
                printf -- "\033[32m The viral genome is $OPTARG. \033[0m\n"
                ;;
                w) gw=1
                printf -- "\033[32m You've indicated that this is genewerk data. \033[0m\n"
                ;;
        esac
done

# Checking for missings

if [ "x" == "x$sample" ]; then
  printf -- "\033[31m [-s sample] option is required. Please indicate the sample name to rerun.\033[0m\n"
  exit
fi

if [ "x" == "x$indir" ]; then
  printf -- "\033[31m [-i path-to-data-dir] option is required. Please select the directory that contains your fastq files.\033[0m\n"
  exit
fi

if [ "x" == "x$suffix" ]; then
  printf -- "\033[31m [-e fastq|fq|fastq.gz|fq.gz] option is required. Please indicate the suffix of your original fastq files.\033[0m\n"
  exit
fi

if [ "x" == "x$reads" ]; then
  printf -- "\033[31m [-r single|paired] option is required. Please select the read type of your sequencng.\033[0m\n"
  exit
fi

if [ "x" == "x$vir_genome" ]; then
  printf -- "\033[31m [-v path-to-viral-genome] option is required. Please direct this option to a bwa indexed viral vector genome.\033[0m\n"
  exit
fi


function_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Defining tools
my_bwa=/path-to-program/bwa-0.7.17/bwa
my_stools=/path-to-program/samtools-1.9/samtools
my_sblaster=/path-to-program/samblaster/samblaster
## Here we're using my_clipm to match at least 30 bp of viral sequence
my_clipm=${function_dir}/samclip2_Pfmod2
my_bedtools=/path-to-program/bedtools2/bin/bedtools
my_makebed=${function_dir}/makeYourBed.R
my_seqtk=/path-to-program/seqtk/seqtk


mkdir -p ${indir}/full_viral/

cd $indir
proj_dir=$(pwd | sed 's/\/data\/raw_fastq//')

f=${sample}

## First, remove any sample files from previous run - don't delete the raw fastq

rm ${indir}/full_viral/${f}.*
rm ${indir}/full_viral/${f}_*
rm ${indir}/${f}*.bam
rm ${indir}/${f}*.fasta; rm ${indir}/${f}*.qual
rm ${indir}/${f}_*trimmed*


if [ "$reads" == "single" ]; then

  # Defining fastq to use
  fq1=${f}.${suffix}

  if [ "$gw" == 1 ]; then
    # Defining RG code
    rg=$(echo \"\@RG\\tID:${f}\\tPL:Illumina\\tPU:x\\tLB:single\\tSM:${f}\")
    
    # Running pipeline
    bsub -q medium -app large -J "v.int.viralign-${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/full_viral/${f}.v.int.viralign.log" -e "${indir}/full_viral/${f}.v.int.viralign.err" \
    "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -m 40 -q 30,30 ${indir}/$fq1 | gzip > ${indir}/${f}_trimmed.fastq.gz
    ${my_bwa} mem ${vir_genome} -t 4 -R ${rg} ${indir}/${f}_trimmed.fastq.gz | \
    ${my_stools} view -h -F 4 - > ${indir}/full_viral/${f}_unfiltered.bam
    ${my_stools} view -h -q 30 ${indir}/full_viral/${f}_unfiltered.bam | \
    ${my_clipm} --max 350 --match 30 --ref ${vir_genome} | \
    ${my_stools} sort -nT ${indir}/full_viral/${f}.prefix.bam - | \
    ${my_stools} fixmate -m - - | \
    ${my_stools} sort -T ${indir}/full_viral/${f}.postfix.bam - | \
    ${my_stools} markdup -r - -O BAM ${indir}/full_viral/${f}.bam"

  else

    # Defining RG code
    rg=$(echo \"\@RG\\tID:${f}\\tPL:Illumina\\tPU:x\\tLB:single\\tSM:${f}\")

    # Running pipeline
    bsub -q medium -app large -J "v.int.viralign-${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/full_viral/${f}.v.int.viralign.log" -e "${indir}/full_viral/${f}.v.int.viralign.err" \
    "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -m 40 -q 30,30 ${indir}/$fq1 | gzip > ${indir}/${f}_trimmed.fastq.gz
    ${my_bwa} mem ${vir_genome} -t 4 -R ${rg} ${indir}/${f}_trimmed.fastq.gz | \
    ${my_stools} view -h -F 4 - > ${indir}/full_viral/${f}_unfiltered.bam
    ${my_stools} view -h -q 30 ${indir}/full_viral/${f}_unfiltered.bam | \
    ${my_clipm} --max 120 --match 30 --ref ${vir_genome} | \
    ${my_stools} sort -nT ${indir}/full_viral/${f}.prefix.bam - | \
    ${my_stools} fixmate -m - - | \
    ${my_stools} sort -T ${indir}/full_viral/${f}.postfix.bam - | \
    ${my_stools} markdup -r - -O BAM ${indir}/full_viral/${f}.bam"
  fi

	# generating masked fastq for alignment to host
	bsub -q medium -app large -w 'ended("v.int.viralign-'${f}'")' -J "mask_vir_aligned_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
    -o "${indir}/full_viral/${f}.mask_vir_aligned.log" -e "${indir}/full_viral/${f}.mask_vir_aligned.err" \
	"$my_stools view ${indir}/full_viral/${f}_unfiltered.bam | awk '{print \$1"\"\\t\""\$6"\"\\t\""\$2}' | Rscript ${my_makebed} | awk '\$2 != \$3 {print \$0}' | sed 's/^ //' > ${indir}/full_viral/${f}_vir_align.bed
	#zcat -f ${indir}/${f}_trimmed.fastq.gz | paste - - - - | tee >(cut -f 1-2 | tr "\"\\t\"" "\"\\n\"" | sed 's/@/>/' > ${indir}/${f}.fasta) | cut -f 3-4 | tr "\"\\t\"" "\"\\n\"" > ${indir}/${f}.qual
	#recommend switching to tee split implemented above, have to get LSF to agree to use bash instead of bourne shell, SO
	zcat -f ${indir}/${f}_trimmed.fastq.gz | paste - - - - | cut -f 1-2 | tr "\"\\t\"" "\"\\n\"" | sed 's/@/>/' > ${indir}/${f}.fasta
	zcat -f ${indir}/${f}_trimmed.fastq.gz | paste - - - - | cut -f 3-4 | tr "\"\\t\"" "\"\\n\"" > ${indir}/${f}.qual
	$my_bedtools maskfasta -fi ${indir}/${f}.fasta -bed ${indir}/full_viral/${f}_vir_align.bed -fo ${indir}/${f}_masked.fasta
	$my_seqtk seq -l0 ${indir}/${f}_masked.fasta > ${indir}/${f}.temp; mv ${indir}/${f}.temp ${indir}/${f}_masked.fasta
	paste <(cat ${indir}/${f}_masked.fasta | paste - - ) <(cat ${indir}/${f}.qual | paste - - ) | tr "\"\\t"" "\"\\n"" | sed 's/>/@/' | gzip > ${indir}/${f}_trimmed_masked.fastq.gz
	rm ${indir}/full_viral/${f}_unfiltered.bam ${indir}/full_viral/${f}_vir_align.bed ${indir}/${f}.fasta ${indir}/${f}.qual ${indir}/${f}_masked.fasta"


else

  # Submitting paired end script
  
  bsub -q long -app large -J "v.int.viralign-${f}" -n 8,8 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
  -o "${indir}/full_viral/${f}.v.int.viralign.log" -e "${indir}/full_viral/${f}.v.int.viralign.err" \
  "${function_dir}/step0.2.paired.sh \
  -g ${vir_genome} \
  -i ${indir} \
  -e ${suffix} \
  -o ${indir}/full_viral \
  -f ${f} >& ${indir}/full_viral/${f}.0.2.txt"
fi


## Submit output check, currently only implemented for paired data
if [ "$reads" == "paired" ]; then
    bsub -q express -app small -J "step0.check" \
    -w 'ended("v.int.viralign*")' \
    "${function_dir}/check_step0.sh \
    -d ${proj_dir}"
fi


# Exit
printf -- '\n';
exit 0;


