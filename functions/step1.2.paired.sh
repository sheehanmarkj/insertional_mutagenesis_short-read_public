#!/bin/bash

module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

set -e
set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this stringent alignment pipeline please specify the following:
$(basename $0) [-g cyno|mouse] [-i path-to-data-dir] [-e fastq|fq|fastq.gz|fq.gz] [-o path-do-output-dir] [-f sample name] [-r rg tag] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi


while getopts 'g:i:o:f:r:' OPTION; do
        case "$OPTION" in
                g) genome="$OPTARG"
                printf -- "\033[32m The genome you've chosen is $OPTARG. \033[0m\n"
                ;;
                i) indir="$OPTARG"
                printf -- "\033[32m Your input directory is: $OPTARG. \033[0m\n"
                ;;
                o) outdir="$OPTARG"
                printf -- "\033[32m Your output directory is $OPTARG. \033[0m\n"
                ;;
                f) f="$OPTARG"
                printf -- "\033[32m The sample name is $OPTARG. \033[0m\n"
                ;;
        esac
done

# Checking for missings

if [ "x" == "x$genome" ]; then
  printf -- "\033[31m [-g cyno|mouse] option is required. Please select the genome to use.\033[0m\n"
  exit
fi

if [ "x" == "x$indir" ]; then
  printf -- "\033[31m [-i path-to-data-dir] option is required. Please select the directory that contains your fastq files.\033[0m\n"
  exit
fi

if [ "x" == "x$outdir" ]; then
  printf -- "\033[31m [-o path-do-output-dir] option is required. Please select the output directory.\033[0m\n"
  exit
fi


if [ "x" == "x$f" ]; then
  printf -- "\033[31m [-f sample name] is required. Please put in a sample name. \033[0m\n"
  exit
fi

function_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Defining tools
my_bwa=/path-to-program/bwa-0.7.17/bwa
my_stools=/path-to-program/samtools-1.9/samtools
my_extract=/path-to-program/extract_reads.py


# Defining RG code
#rg=$(echo \@RG\\tID:${f}\\tPL:Illumina\\tPU:x\\tLB:$reads\\tSM:${f})

# Defining fqs

fq1=${f}_1_trimmed_masked.fastq.gz
fq2=${f}_2_trimmed_masked.fastq.gz
fq_unmasked=${f}_trimmed.fastq.gz


# Define temp sort directory
proj_dir=$(echo $indir | sed 's/\/data\/raw_fastq//')
sort_dir=${proj_dir}/data/temp_sort_dirs/${f}/
mkdir -p ${sort_dir}



# Code to run
rg=$(echo \"\@RG\\tID:${f}\\tPL:Illumina\\tPU:x\\tLB:$reads\\tSM:${f}\")

bsub -q long -app large -J "v.int.ta-${f}_unmasked" -n 16,16 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
-o "${outdir}/${f}.v.int.1_unmasked.log" -e "${outdir}/${f}.v.int.1_unmasked.err" \
"cutadapt --interleaved -m 40 ${indir}/${fq_unmasked} | \
${my_bwa} mem ${genome} -t 16 -R ${rg} -p - | \
${my_stools} view -h -u - | \
${my_stools} sort -nT ${outdir}/${f}_unmasked.prefix.bam - | \
${my_stools} fixmate -m - - | \
${my_stools} sort -T ${outdir}/${f}_unmasked.postfix.bam - | \
${my_stools} markdup -S - -O BAM ${outdir}/${f}_unmasked_temp.bam

${my_stools} view -F 1024 ${outdir}/${f}_unmasked_temp.bam | cut -f 1 | sort -T ${sort_dir} -u > ${outdir}/${f}_nondupReads.txt
${my_stools} view -b -F 1024 ${outdir}/${f}_unmasked_temp.bam > ${outdir}/${f}_unmasked.bam
rm ${outdir}/${f}_unmasked_temp.bam"

bsub -q medium -app large -w 'ended("v.int.ta-'${f}'*")' -J "v.int.ta-${f}_remDup" -n 16,16 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
-o "${outdir}/${f}.v.int.1_remDup.log" -e "${outdir}/${f}.v.int.1_remDup.err" \
"python ${my_extract} -b ${outdir}/${f}_temp.bam -n ${outdir}/${f}_nondupReads.txt -o ${outdir}/${f}.bam
rm ${outdir}/${f}_temp.bam"


rg=$(echo \@RG\\tID:${f}\\tPL:Illumina\\tPU:x\\tLB:$reads\\tSM:${f})
paste <(zcat -f ${indir}/$fq1 | paste - - - -) \
      <(zcat -f ${indir}/$fq2 | paste - - - -) | \
    tr '\t' '\n' | \
    cutadapt --interleaved -m 40 - | \
    ${my_bwa} mem ${genome} -t 16 -R ${rg} -p - | \
    ${my_stools} view -h -u - | \
    ${my_stools} sort -nT ${outdir}/${f}.prefix.bam - | \
    ${my_stools} fixmate -m - - | \
    ${my_stools} sort -T ${outdir}/${f}.postfix.bam - | \
    ${my_stools} markdup -r - -O BAM ${outdir}/${f}_temp.bam


# Exit
printf -- '\n';
exit 0;
