#!/bin/bash
module load ib R/4.1.2
module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

set -e
set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this stringent alignment pipeline please specify the following:
$(basename $0) [-g viral-genome] [-i path-to-data-dir] [-e fastq|fq|fastq.gz|fq.gz] [-o path-to-output-dir] [-f sample name] [-r rg tag] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi


while getopts 'g:i:e:o:f:r:' OPTION; do
        case "$OPTION" in
                g) genome="$OPTARG"
                printf -- "\033[32m The genome you've chosen is $OPTARG. \033[0m\n"
                ;;
                i) indir="$OPTARG"
                printf -- "\033[32m Your input directory is: $OPTARG. \033[0m\n"
                ;;
                e) suffix="$OPTARG"
                printf -- "\033[32m The suffix of your files is $OPTARG. \033[0m\n"
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
  printf -- "\033[31m [-g viral-genome] option is required. Please select the genome to use.\033[0m\n"
  exit
fi

if [ "x" == "x$indir" ]; then
  printf -- "\033[31m [-i path-to-data-dir] option is required. Please select the directory that contains your fastq files.\033[0m\n"
  exit
fi

if [ "x" == "x$suffix" ]; then
  printf -- "\033[31m [-e fastq|fq|fastq.gz|fq.gz] option is required. Please select the suffix of your fastq files.\033[0m\n"
  exit
fi

if [ "x" == "x$outdir" ]; then
  printf -- "\033[31m [-o path-to-output-dir] option is required. Please select the output directory.\033[0m\n"
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
## Using my_clipm for viral alignments
my_clipm=${function_dir}/samclip2_Pfmod2
my_bedtools=/path-to-program/bedtools2/bin/bedtools
my_makebed=${function_dir}/makeYourBed.R
my_seqtk=/path-to-program/seqtk/seqtk


# Defining RG code
rg=$(echo \@RG\\tID:${f}\\tPL:Illumina\\tPU:x\\tLB:paired\\tSM:${f})

# Defining fqs

fq1=${f}_1.${suffix}
fq2=${f}_2.${suffix}

# Code to run

paste <(zcat -f ${indir}/$fq1 | paste - - - -) \
      <(zcat -f ${indir}/$fq2 | paste - - - -) | \
    tr '\t' '\n' | \
    cutadapt --interleaved -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -m 40 -q 30,30 - | gzip > ${indir}/${f}_trimmed.fastq.gz

${my_bwa} mem ${genome} -t 4 -R ${rg} -p ${indir}/${f}_trimmed.fastq.gz | \
    ${my_stools} view -h -F 4 - > ${outdir}/${f}_unfiltered.bam
${my_stools} view -h -q 30 ${outdir}/${f}_unfiltered.bam | \
    ${my_clipm} --max 120 --match 30 --ref ${genome} | \
    ${my_stools} sort -nT ${outdir}/${f}.prefix.bam - | \
    ${my_stools} fixmate -m - - | \
    ${my_stools} sort -T ${outdir}/${f}.postfix.bam - | \
    ${my_stools} markdup -r - -O BAM ${outdir}/${f}.bam

# generate beds for first and second in pair
$my_stools view -f 64 ${outdir}/${f}_unfiltered.bam | awk '{print $1"\t"$6"\t"$2}' | Rscript ${my_makebed} | awk '$2 != $3 {print $0}' | sed 's/^ //' > ${outdir}/${f}_vir_align_1.bed
$my_stools view -f 128 ${outdir}/${f}_unfiltered.bam | awk '{print $1"\t"$6"\t"$2}' | Rscript ${my_makebed} | awk '$2 != $3 {print $0}' | sed 's/^ //' > ${outdir}/${f}_vir_align_2.bed

## split fastqs into seq and qual
zcat -f ${indir}/${f}_trimmed.fastq.gz | paste - - - - - - - - | \
tee >(cut -f 1-4 | tee >(cut -f 1-2 | tr "\t" "\n" | sed 's/@/>/' > ${indir}/${f}_1.fasta) | cut -f 3-4 | tr "\t" "\n" > ${indir}/${f}_1.qual) | \
cut -f 5-8 | tee >(cut -f 1-2 | tr "\t" "\n" | sed 's/@/>/' > ${indir}/${f}_2.fasta) | cut -f 3-4 | tr "\t" "\n" > ${indir}/${f}_2.qual


## mask fastas
$my_seqtk seq -l0 -M ${outdir}/${f}_vir_align_1.bed -n N ${indir}/${f}_1.fasta > ${indir}/${f}_1_masked.fasta
$my_seqtk seq -l0 -M ${outdir}/${f}_vir_align_2.bed -n N ${indir}/${f}_2.fasta > ${indir}/${f}_2_masked.fasta
#$my_bedtools maskfasta -fi ${indir}/${f}_1.fasta -bed ${outdir}/${f}_vir_align_1.bed -fo ${indir}/${f}_1_masked.fasta
#$my_bedtools maskfasta -fi ${indir}/${f}_2.fasta -bed ${outdir}/${f}_vir_align_2.bed -fo ${indir}/${f}_2_masked.fasta
#$my_seqtk seq -l0 ${indir}/${f}_1_masked.fasta > ${indir}/${f}.temp; mv ${indir}/${f}.temp ${indir}/${f}_1_masked.fasta #sed -i ':a;N;$!ba;s/\([ACGTN]\)\n\([ACGTN]\)/\1\2/g' ${indir}/${f}_1_masked.fasta
#$my_seqtk seq -l0 ${indir}/${f}_2_masked.fasta > ${indir}/${f}.temp; mv ${indir}/${f}.temp ${indir}/${f}_2_masked.fasta #sed -i ':a;N;$!ba;s/\([ACGTN]\)\n\([ACGTN]\)/\1\2/g' ${indir}/${f}_2_masked.fasta

## restore quality scores, gzip
paste <(cat ${indir}/${f}_1_masked.fasta | paste - - ) <(cat ${indir}/${f}_1.qual | paste - - ) | tr '\t' '\n' | sed 's/^>/@/' | gzip > ${indir}/${f}_1_trimmed_masked.fastq.gz
paste <(cat ${indir}/${f}_2_masked.fasta | paste - - ) <(cat ${indir}/${f}_2.qual | paste - - ) | tr '\t' '\n' | sed 's/^>/@/' | gzip > ${indir}/${f}_2_trimmed_masked.fastq.gz

## cleanup?
#rm ${indir}/${f}_*.fasta ${indir}/${f}_*.qual ${outdir{/${f}_vir_align_*.bed
rm ${outdir}/${f}_unfiltered.bam ${outdir}/${f}_vir_align_1.bed ${outdir}/${f}_vir_align_2.bed ${indir}/${f}_1.fasta ${indir}/${f}_1.qual ${indir}/${f}_2.fasta ${indir}/${f}_2.qual ${indir}/${f}_1_masked.fasta ${indir}/${f}_2_masked.fasta

# Exit
printf -- '\n';
exit 0;
