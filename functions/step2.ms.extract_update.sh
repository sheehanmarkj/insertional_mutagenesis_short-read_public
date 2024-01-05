#!/bin/bash

module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

#set -e
#set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this post-alignment mate and softclip extraction pipeline please specify the following:
$(basename $0) [-i path-to-data-dir] [-g cyno|mouse] [-r single|paired] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi


while getopts 'i:g:r:' OPTION; do
        case "$OPTION" in
                i) indir="$OPTARG"
                printf -- "\033[32m Your input directory is: $OPTARG. \033[0m\n"
                ;;
                g) species="$OPTARG"
                printf -- "\033[32m The genome you've chosen is $OPTARG. \033[0m\n"
                ;;                
		r) reads="$OPTARG"
                printf -- "\033[32m Your readtype is $OPTARG end. \033[0m\n"
                ;;
        esac
done

# Checking for missings

if [ "x" == "x$indir" ]; then
  printf -- "\033[31m [-i path-to-data-dir] option is required. Please select the directory that contains your fastq files.\033[0m\n"
  exit
fi

if [ "x" == "x$species" ]; then
  printf -- "\033[31m [-g cyno|mouse] option is required. Please select the species from which this data was sequenced.\033[0m\n"
  exit
fi

if [ "x" == "x$reads" ]; then
  printf -- "\033[31m [-r single|paired] option is required. Please select the read type of your sequencng.\033[0m\n"
  exit
fi

# Check error status of previous script
proj_dir=$(echo $indir | sed 's/\/data\/trimalign//')


if [ "$reads" == "paired" ]; then
  grep "xit code" ${proj_dir}/data/trimalign/*.v.int.*.log > ${proj_dir}/step1_check.txt
  list=${proj_dir}/data/raw_fastq/sample_list.txt
  for f in $(cat $list); do
    if [ ! -f ${proj_dir}/data/trimalign/${f}_unmasked.bam ]; then printf "${f}_unmasked.bam not found\n" >> ${proj_dir}/step1_check.txt; fi
    if [ ! -f ${proj_dir}/data/trimalign/${f}.bam ]; then printf "${f}.bam not found\n" >> ${proj_dir}/step1_check.txt; fi
  done
  if [ -s "${proj_dir}/step1_check.txt" ]; then
    printf -- "Possible issue, check ${proj_dir}/step1_check.txt for more information. Remove this file before continuing this step. Exiting."
    exit
  fi
fi

# Starting the process
cd $indir

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
## Defining tools
my_stools=/path-to-program/samtools-1.9/samtools
my_bedtools=/path-to-program/bedtools2/bin/bedtools
# Here my-clip is used to extract reads reads that have clips, so appropriate to use the softclip filter, rather than the match filter
my_clipm=${function_dir}/samclip2_Pfmod2

# Making directories
mkdir -p ${indir}/mate
mkdir -p ${indir}/softclip

# Defining directories
mate_out=${indir}/mate
soft_out=${indir}/softclip

list=${proj_dir}/data/raw_fastq/sample_list.txt

for sample in $(cat $list); do
   if [ "$reads" == "paired" ]; then

     # Define temp sort directory
     sort_dir=${proj_dir}/data/temp_sort_dirs/${sample}/

     # Singly unmapped reads (remember to never include q30 flag here)
     bsub -q medium -app large -J "v.int.mates_${sample}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
     -o "${mate_out}/${sample}.v.int.2.log" -e "${mate_out}/${sample}.v.int.2.err" \
     "${my_stools} view -f 68 -F 8 ${indir}/${sample}.bam | cut -f 1 | sort -T ${sort_dir} | uniq > ${mate_out}/${sample}_unaligned_reads_1_masked.txt
	${my_stools} view -f 132 -F 8 ${indir}/${sample}.bam | cut -f 1 | sort -T ${sort_dir} | uniq > ${mate_out}/${sample}_unaligned_reads_2_masked.txt
	${my_stools} view -f 68 -F 8 ${indir}/${sample}_unmasked.bam | cut -f 1 | sort -T ${sort_dir} | uniq > ${mate_out}/${sample}_unaligned_reads_1_unmasked.txt
	${my_stools} view -f 132 -F 8 ${indir}/${sample}_unmasked.bam | cut -f 1 | sort -T ${sort_dir} | uniq > ${mate_out}/${sample}_unaligned_reads_2_unmasked.txt
	cat ${mate_out}/${sample}_unaligned_reads_1_masked.txt ${mate_out}/${sample}_unaligned_reads_1_unmasked.txt | sort -T ${sort_dir} | uniq -d > ${mate_out}/${sample}_unaligned_reads_1.txt
	cat ${mate_out}/${sample}_unaligned_reads_2_masked.txt ${mate_out}/${sample}_unaligned_reads_2_unmasked.txt | sort -T ${sort_dir} | uniq -d > ${mate_out}/${sample}_unaligned_reads_2.txt
	${my_stools} view -h -f 8 -F 4 -q 30 ${indir}/${sample}.bam | ${my_clipm} --max 120 --match 30 --ref ${genome} | \
        $my_stools sort > ${mate_out}/${sample}_mapped_mates.bam"


     # Softclipped reads
     bsub -q medium -app large -J "v.int.soft_${sample}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
     -o "${soft_out}/${sample}.v.int.2.log" -e "${soft_out}/${sample}.v.int.2.err" \
     "$my_stools view -h ${indir}/${sample}.bam |\
     $my_clipm \
     --match 30 \
     --invert \
     --max 30 \
     --ref $genome |\
     $my_stools sort > ${soft_out}/${sample}_softclipped.bam #### 

     # Getting read names from bam
     ${my_stools} view -f 64 ${soft_out}/${sample}_softclipped.bam | cut -f 1 | sort -T ${sort_dir} | uniq > ${soft_out}/${sample}_softclipped_reads_1_masked.txt
     ${my_stools} view -f 128 ${soft_out}/${sample}_softclipped.bam | cut -f 1 | sort -T ${sort_dir} | uniq > ${soft_out}/${sample}_softclipped_reads_2_masked.txt
     
     ## repeat, unmasked
     $my_stools view -h ${indir}/${sample}_unmasked.bam |\
     $my_clipm \
     --match 30 \
     --max 30 \
     --ref $genome |\
     $my_stools sort > ${soft_out}/${sample}_softclipped_unmasked.bam

     # Getting read names from bam
     ${my_stools} view -f 64 ${soft_out}/${sample}_softclipped_unmasked.bam | cut -f 1 | sort -T ${sort_dir} | uniq > ${soft_out}/${sample}_softclipped_reads_1_unmasked.txt
     ${my_stools} view -f 128 ${soft_out}/${sample}_softclipped_unmasked.bam | cut -f 1 | sort -T ${sort_dir} | uniq > ${soft_out}/${sample}_softclipped_reads_2_unmasked.txt
     
     # Intersection 
     cat ${soft_out}/${sample}_softclipped_reads_1_masked.txt ${soft_out}/${sample}_softclipped_reads_1_unmasked.txt ${soft_out}/${sample}_softclipped_reads_1_unmasked.txt | sort -T ${sort_dir} | uniq -u > ${soft_out}/${sample}_softclipped_reads_1.txt

     cat ${soft_out}/${sample}_softclipped_reads_2_masked.txt ${soft_out}/${sample}_softclipped_reads_2_unmasked.txt ${soft_out}/${sample}_softclipped_reads_2_unmasked.txt | sort -T ${sort_dir} | uniq -u > ${soft_out}/${sample}_softclipped_reads_2.txt"


   elif [ "$reads" == "single" ]; then

     # Softclipped reads
     bsub -q medium -app large -J "v.int.soft_${sample}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
     -o "${soft_out}/${sample}.v.int.2.log" -e "${soft_out}/${sample}.v.int.2.err" \
     "$my_stools view -h ${indir}/${sample}.bam |\
     $my_clipm \
     --match 30 \
     --invert \
     --max 30 \
     --ref $genome |\
     $my_stools sort > ${soft_out}/${sample}_softclipped.bam

     # Getting read names from bam
     ${my_stools} view ${soft_out}/${sample}_softclipped.bam | cut -f 1 | sort | uniq > ${soft_out}/${sample}_softclipped_reads_masked.txt"

     # Softclipped reads, unmasked alignment - because of alignments to regions with seq overlap to virus, reversing logic and EXCLUDING found reads (i.e. alignments to host with fewer than 30 clipped bases)
     bsub -q medium -app large -J "v.int.soft_unmasked${sample}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
     -o "${soft_out}/${sample}.v.int.2_unmasked.log" -e "${soft_out}/${sample}.v.int.2_unmasked.err" \
     "$my_stools view -h ${indir}/${sample}_unmasked.bam |\
     $my_clipm \
     --match 30 \
     --max 30 \
     --ref $genome |\
     $my_stools sort > ${soft_out}/${sample}_softclipped_unmasked.bam

     # Getting read names from bam
     ${my_stools} view ${soft_out}/${sample}_softclipped_unmasked.bam | cut -f 1 | sort | uniq > ${soft_out}/${sample}_softclipped_reads_unmasked.txt"


	# Intersection of softclipped reads
     bsub -q medium -app large -w 'ended("v.int.soft*'${sample}'")' -J "v.int.soft_intersect${sample}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
     -o "${soft_out}/${sample}.v.int.2_intersect.log" -e "${soft_out}/${sample}.v.int.2_intersect.err" \
     "cat ${soft_out}/${sample}_softclipped_reads_masked.txt ${soft_out}/${sample}_softclipped_reads_unmasked.txt ${soft_out}/${sample}_softclipped_reads_unmasked.txt | sort | uniq -u > ${soft_out}/${sample}_softclipped_reads.txt"


   fi

done

## Submit output check, currently only implemented for paired data
if [ "$reads" == "paired" ]; then
    bsub -q express -app small -J "step2.check" \
    -w 'ended("v.int*")' \
    "${function_dir}/check_step2.sh \
    -d ${proj_dir}"
fi


# Exit
printf -- '\n';
exit 0;
