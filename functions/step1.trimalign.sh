#!/bin/bash

module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

set -e
set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this stringent alignment pipeline please specify the following:
$(basename $0) [-g cyno|mouse] [-i path-to-data-dir] [-e fastq|fq|fastq.gz|fq.gz] [-r single|paired] [-o path-do-output-dir] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi


while getopts 'g:i:e:r:o:' OPTION; do
        case "$OPTION" in
                g) species="$OPTARG"
                printf -- "\033[32m The genome you've chosen is $OPTARG. \033[0m\n"
                ;;
                i) indir="$OPTARG"
                printf -- "\033[32m Your input directory is: $OPTARG. \033[0m\n"
                ;;
                e) suffix="$OPTARG"
                printf -- "\033[32m The suffix of your files is $OPTARG. \033[0m\n"
                ;;
                r) reads="$OPTARG"
                printf -- "\033[32m Your readtype is $OPTARG end. \033[0m\n"
                ;;
                o) outdir="$OPTARG"
                printf -- "\033[32m Your output directory is $OPTARG. \033[0m\n"
                ;;
        esac
done

# Checking for missings

if [ "x" == "x$species" ]; then
  printf -- "\033[31m [-g cyno|mouse] option is required. Please select the species from which this data was sequenced.\033[0m\n"
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

if [ "x" == "x$reads" ]; then
  printf -- "\033[31m [-r single|paired] option is required. Please select the read type of your sequencng.\033[0m\n"
  exit
fi


if [ "x" == "x$outdir" ]; then
  printf -- "\033[31m [-o path-do-output-dir] option is required. Please select the output directory.\033[0m\n"
  exit
fi

# Check error status of previous script, and prepare sort directory
proj_dir=$(echo $indir | sed 's/\/data\/raw_fastq//')

sort_dir=${proj_dir}/data/temp_sort_dirs
mkdir -p ${sort_dir}

if [ "$reads" == "paired" ]; then
  if [ -s "${proj_dir}/step0_check.txt" ]; then
    printf -- "Possible issue, check ${proj_dir}/step0_check.txt for more information. Remove this file before continuing this step. Exiting."
    exit
  fi
fi


# Starting the process

mkdir -p ${outdir}

## Defining genome
if [ "$species" == "cyno" ]; then
  genome=/path-to-cyno-genome/Macaca_fascicularis.Macaca_fascicularis_5.0.dna.toplevel.fa
  elif [ "$species" == "human" ]; then
  genome=/path-to-human-genome/Homo_sapiens.GRCh38.dna.toplevel.fa
  elif [ "$species" == "dog" ]; then
  genome=/path-to-dog-genome/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa
  elif [ "$species" == "rat" ]; then
  genome=/path-to-rat-genome/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
  else
  genome=/path-to-mouse-genome/Mus_musculus.GRCm38.dna.toplevel.fa
fi

function_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Defining tools
my_bwa=/path-to-program/bwa-0.7.17/bwa
my_stools=/path-to-program/samtools-1.9/samtools
my_sblaster=/path-to-program/samblaster/samblaster
my_bedtools=/path-to-program/bedtools2/bin/bedtools
my_extract=/path-to-program/extract_reads.py



## Creating a sample list

cd $indir

#if [ "$reads" == "paired" ]; then
#  ls -1 *${suffix} | grep -v trimmed.${suffix} | sed -e "s/\_1.*//g" | sed -e "s/\_2.*//g" |  uniq > $indir/sample_list.txt
#  else
#  ls -1 *${suffix} | grep -v trimmed.${suffix} | sed "s/\.$suffix//g"| uniq > $indir/sample_list.txt
#fi

list=$indir/sample_list.txt

# Checking the total number of samples

tot=$(cat $indir/sample_list.txt | wc -l)

if [ "$tot" -eq "0" ]; then
        printf -- "\033[31m There were no samples of your specified suffix detected in your data directory. Please check if you have your fastq.gz or fq.gz files in the sample directory. If so, please make sure that you have specified fastq or fq in the -e option above.\033]0m\n"
fi

# Trimming and aligning samples to correct host

for f in $(cat $list); do
    # Defining RG code
    rg=$(echo \"\@RG\\tID:${species}_${f}\\tPL:Illumina\\tPU:x\\tLB:$reads\\tSM:${species}_${f}\")

  if [ "$reads" == "single" ]; then

    # Defining fastq to use
    fq1=${f}_trimmed_masked.fastq.gz
	fq1_unmasked=${f}_trimmed.fastq.gz

    # Running pipeline
    bsub -q long -app large -J "v.int.ta-${f}" -n 16,16 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
    -o "${outdir}/${f}.v.int.1.log" -e "${outdir}/${f}.v.int.1.err" \
    "cutadapt -m 40 ${indir}/$fq1 | \
    ${my_bwa} mem ${genome} -t 16 -R ${rg} - | \
    ${my_stools} view -h -u - | \
    ${my_stools} sort -nT ${outdir}/${f}.prefix.bam - | \
    ${my_stools} fixmate -m - - | \
    ${my_stools} sort -T ${outdir}/${f}.postfix.bam - | \
    ${my_stools} markdup -r - -O BAM ${outdir}/${f}_temp.bam"

    bsub -q long -app large -J "v.int.ta-${f}_unmasked" -n 16,16 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
    -o "${outdir}/${f}.v.int.1_unmasked.log" -e "${outdir}/${f}.v.int.1_unmasked.err" \
    "cutadapt -m 40 ${indir}/$fq1_unmasked | \
    ${my_bwa} mem ${genome} -t 16 -R ${rg} - | \
    ${my_stools} view -h -u - | \
    ${my_stools} sort -nT ${outdir}/${f}_unmasked.prefix.bam - | \
    ${my_stools} fixmate -m - - | \
    ${my_stools} sort -T ${outdir}/${f}_unmasked.postfix.bam - | \
    ${my_stools} markdup -S - -O BAM ${outdir}/${f}_unmasked_temp.bam

    ${my_stools} view -F 1024 ${outdir}/${f}_unmasked_temp.bam | cut -f 1 | sort -u > ${outdir}/${f}_nondupReads.txt
    ${my_stools} view -b -F 1024 ${outdir}/${f}_unmasked_temp.bam > ${outdir}/${f}_unmasked.bam
    rm ${outdir}/${f}_unmasked_temp.bam"

    bsub -q medium -app large -w 'ended("v.int.ta-'${f}'*")' -J "v.int.${f}_remDup" -n 16,16 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
    -o "${outdir}/${f}.v.int.1_remDup.log" -e "${outdir}/${f}.v.int.1_remDup.err" \
    "python ${my_extract} -b ${outdir}/${f}_temp.bam -n ${outdir}/${f}_nondupReads.txt -o ${outdir}/${f}.bam
    rm ${outdir}/${f}_temp.bam"


    else

    # Submitting paired end script

    bsub -q long -app large -J "v.int.ta-${f}" -n 16,16 -M 64GB -R "span[hosts=1] rusage[mem=64GB]" \
    -o "${outdir}/${f}.v.int.1.log" -e "${outdir}/${f}.v.int.1.err" \
    "${function_dir}/step1.2.paired.sh \
    -g ${genome} \
    -i ${indir} \
    -o ${outdir} \
    -f ${f} >& $outdir/err1.2.txt"
  fi
done

## Submit output check, currently only implemented for paired data
#### Runs too soon, not holding on the v.int.ta* from jobs run in paired script. Moved to start of step2.
#if [ "$reads" == "paired" ]; then
#    bsub -q express -app small -J "step1.check" \
#    -w 'ended("v.int.ta*")' \
#    "${function_dir}/check_step1.sh \
#    -d ${proj_dir}"
#fi


# Exit
printf -- '\n';
exit 0;
