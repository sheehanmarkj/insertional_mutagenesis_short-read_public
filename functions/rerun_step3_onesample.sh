#!/bin/bash

module load eb/2017  GCC/5.4.0-2.27  OpenMPI/2.0.0 cutadapt/1.9.1-Python-2.7.12

#set -e
#set -o pipefail

#############################################################################
# Creating a usage sentence
usage="To use this viral alignment for initial sample, mate and softclip reads:
$(basename $0) [-s sample-name] [-g cyno|mouse|human] [-i path-to-folder-containing-mate-and-softclip-folders] [-f path-to-folder-containing-initial-fastqs] [-r single|paired] [-w this is genewerk data] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi

while getopts 's:g:i:f:r:w' OPTION; do
        case "$OPTION" in
                s) sample="$OPTARG"
                printf -- "\033[32m The sample you've chosen is $OPTARG. \033[0m\n"
                ;;
                g) species="$OPTARG"
                printf -- "\033[32m The genome you've chosen is $OPTARG. \033[0m\n"
                ;;
                i) indir="$OPTARG"
                printf -- "\033[32m The sample directory is $OPTARG. \033[0m\n"
                ;;
                f) fqdir="$OPTARG"
                printf -- "\033[32m The fastq directory is $OPTARG. \033[0m\n"
                ;;
                r) reads="$OPTARG"
                printf -- "\033[32m Your readtype is $OPTARG end. \033[0m\n"
                ;;
                w) gw=1
                printf -- "\033[32m You've indicated that this is genewerk data. \033[0m\n"
                ;;
        esac
done

# Checking for missings
if [ "x" == "x$sample" ]; then
  printf -- "\033[31m [-s sample-name] option is required. Please select the sample to rerun.\033[0m\n"
  exit
fi

if [ "x" == "x$species" ]; then
  printf -- "\033[31m [-g cyno|mouse|human] option is required. Please select the species from which this data was sequenced.\033[0m\n"
  exit
fi

if [ "x" == "x$indir" ]; then
  printf -- "\033[31m [-i path-to-alignment-dir] option is required. Please select the directory that contains your host-alignment files.\033[0m\n"
  exit
fi

if [ "x" == "x$fqdir" ]; then
  printf -- "\033[31m [-f path-to-fastq-dir] option is required. Please select the directory that contains your intiial fastq files.\033[0m\n"
  exit
fi

if [ "x" == "x$reads" ]; then
  printf -- "\033[31m [-r single|paired] option is required. Please select the read type of your sequencng.\033[0m\n"
  exit
fi

# Check error status of previous script
proj_dir=$(echo $indir | sed 's/\/data\/trimalign//')

if [ "$reads" == "paired" ]; then
  if [ -s "${proj_dir}/step2_check.txt" ]; then
    printf -- "Possible issue, check ${proj_dir}/step2_check.txt for more information. Remove this file before continuing this step. Exiting."
    exit
  fi
fi


## Defining genome
if [ "$species" == "cyno" ]; then
  genome=/path-to-cyno-genome/Macaca_fascicularis.Macaca_fascicularis_5.0.dna.toplevel.fa
  elif [ "$species" == "human" ]; then
  genome=/path-to-human-genome/Homo_sapiens.GRCh38.dna.toplevel.fa
  elif [ "$species" == "rat" ]; then
  genome=/path-to-rat-genome/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
  else
  genome=/path-to-mouse-genome/Mus_musculus.GRCm38.dna.toplevel.fa
fi

function_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Defining tools
my_stools=/path-to-program/samtools-1.9/samtools
my_sblaster=/path-to-program/samblaster/samblaster
my_bedtools=/path-to-program/bedtools2/bin/bedtools
my_clipm=${function_dir}/samclip2_Pfmod2
my_extract=/path-to-program/extract_reads.py



# First remove any residual sample files from previous run
rm ${proj_dir}/data/temp_sort_dirs/${sample}/*
rm ${indir}/mate/${sample}.mate_vir_bam_extract.log
rm ${indir}/mate/${sample}.mate_vir_bam_extract.err
rm ${fqdir}/full_viral/${sample}_1.bam
rm ${indir}/mate/${sample}_unmapped_mate_1_temp.bam
rm ${indir}/mate/${sample}_unmapped_mate_1.bam
rm ${fqdir}/full_viral/${sample}_2.bam
rm ${indir}/mate/${sample}_unmapped_mate_2_temp.bam
rm ${indir}/mate/${sample}_unmapped_mate_2.bam
rm ${indir}/mate/${sample}_unaligned_mapped_filtered_reads_info.txt
rm ${indir}/mate/${sample}_unaligned_mapped_filtered_reads.txt
rm ${indir}/mate/${sample}_unmapped_mate.bam
rm ${indir}/mate/${sample}_viral_mate.bam
rm ${indir}/mate/samtools*; rm ${indir}/softclip/samtools*; rm ${indir}/samtools*
rm ${indir}/softclip/${sample}_softclipped_vir_1_temp.bam
rm ${indir}/softclip/${sample}_softclipped_vir_2_temp.bam
rm ${indir}/softclip/${sample}_softclipped_vir_1.bam
rm ${indir}/softclip/${sample}_softclipped_vir_2.bam
rm ${indir}/softclip/${sample}_softclipped_mapped_filtered_reads_info.txt
rm ${indir}/softclip/${sample}_unmapped_softclip.bam
rm ${indir}/softclip/${sample}_viral_softclip.bam
rm ${indir}/${sample}_merged.bam
rm ${indir}/${sample}_sorted.bam
rm ${indir}/${sample}_reads_info.txt


# next part is extracting the reads from viral bam
# Starting with mates
if [ "$reads" == "paired" ]; then
	nam=${sample}
	sort_dir=${proj_dir}/data/temp_sort_dirs/${nam}/

        cd ${indir}/mate

        ## extracting virally mapped reads and filtering
	
        bsub -q short -app large -J "v.int.3.mate_vir_bam_extract_${nam}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
                -o "${indir}/mate/${nam}.mate_vir_bam_extract.log" -e "${indir}/mate/${nam}.mate_vir_bam_extract.err" \
                "
                #nam=${sample}
                #sort_dir=${proj_dir}/data/temp_sort_dirs/${nam}/
                ${my_stools} view -b -f 64 ${fqdir}/full_viral/${nam}.bam | ${my_stools} sort > ${fqdir}/full_viral/${nam}_1.bam
                python ${my_extract} -b ${fqdir}/full_viral/${nam}_1.bam -n ${indir}/mate/${nam}_unaligned_reads_1.txt -o ${indir}/mate/${nam}_unmapped_mate_1_temp.bam
                ${my_stools} sort ${indir}/mate/${nam}_unmapped_mate_1_temp.bam > ${indir}/mate/${nam}_unmapped_mate_1.bam; rm ${indir}/mate/${nam}_unmapped_mate_1_temp.bam
                ${my_stools} view -b -f 128 ${fqdir}/full_viral/${nam}.bam | ${my_stools} sort > ${fqdir}/full_viral/${nam}_2.bam
                python ${my_extract} -b ${fqdir}/full_viral/${nam}_2.bam -n ${indir}/mate/${nam}_unaligned_reads_2.txt -o ${indir}/mate/${nam}_unmapped_mate_2_temp.bam
                ${my_stools} sort ${indir}/mate/${nam}_unmapped_mate_2_temp.bam > ${indir}/mate/${nam}_unmapped_mate_2.bam; rm ${indir}/mate/${nam}_unmapped_mate_2_temp.bam
                ${my_stools} view ${indir}/mate/${nam}_unmapped_mate_1.bam | cut -f 1,3,4,6,9 | sort -T ${sort_dir} | uniq > ${indir}/mate/${nam}_unaligned_mapped_filtered_reads_info.txt
                ${my_stools} view ${indir}/mate/${nam}_unmapped_mate_2.bam | cut -f 1,3,4,6,9 | sort -T ${sort_dir} | uniq >> ${indir}/mate/${nam}_unaligned_mapped_filtered_reads_info.txt
                cut -f 1 ${indir}/mate/${nam}_unaligned_mapped_filtered_reads_info.txt | sort -T ${sort_dir} | uniq > ${indir}/mate/${nam}_unaligned_mapped_filtered_reads.txt
                ${my_stools} merge ${indir}/mate/${nam}_unmapped_mate.bam ${indir}/mate/${nam}_unmapped_mate_1.bam ${indir}/mate/${nam}_unmapped_mate_2.bam
                python ${my_extract} -b ${indir}/mate/${nam}_mapped_mates.bam -n ${indir}/mate/${nam}_unaligned_mapped_filtered_reads.txt -o ${indir}/mate/${nam}_viral_mate.bam
                "


        # Extracting softclips
        cd ${indir}/softclip/

        bsub -q short -app large -w 'ended("v.int.3.mate_vir_bam_extract_*")' -J "v.int.3.softclip_vir_bam_extract_${nam}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
                -o "${indir}/softclip/${nam}.softclip_vir_bam_extract.log" -e "${indir}/softclip/${nam}.softclip_vir_bam_extract.err" \
                "
                #nam=${sample}
                #sort_dir=${proj_dir}/data/temp_sort_dirs/${nam}/
                python ${my_extract} -b ${fqdir}/full_viral/${nam}_1.bam -n ${indir}/softclip/${nam}_softclipped_reads_1.txt -o ${indir}/softclip/${nam}_softclipped_vir_1_temp.bam
                python ${my_extract} -b ${fqdir}/full_viral/${nam}_2.bam -n ${indir}/softclip/${nam}_softclipped_reads_2.txt -o ${indir}/softclip/${nam}_softclipped_vir_2_temp.bam
                ${my_stools} sort ${indir}/softclip/${nam}_softclipped_vir_1_temp.bam > ${indir}/softclip/${nam}_softclipped_vir_1.bam; rm ${indir}/softclip/${nam}_softclipped_vir_1_temp.bam
                ${my_stools} sort ${indir}/softclip/${nam}_softclipped_vir_2_temp.bam > ${indir}/softclip/${nam}_softclipped_vir_2.bam; rm ${indir}/softclip/${nam}_softclipped_vir_2_temp.bam
                ${my_stools} view ${indir}/softclip/${nam}_softclipped_vir_1.bam | cut -f 1,3,4,6,9 | sort -T ${sort_dir} | uniq > ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads_info.txt
                ${my_stools} view ${indir}/softclip/${nam}_softclipped_vir_2.bam | cut -f 1,3,4,6,9 | sort -T ${sort_dir} | uniq >> ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads_info.txt
                cut -f 1 ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads_info.txt | sort -T ${sort_dir} | uniq > ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads.txt
                ${my_stools} merge ${indir}/softclip/${nam}_unmapped_softclip.bam ${indir}/softclip/${nam}_softclipped_vir_1.bam ${indir}/softclip/${nam}_softclipped_vir_2.bam
                python ${my_extract} -b ${indir}/softclip/${nam}_softclipped.bam -n ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads.txt -o ${indir}/softclip/${nam}_viral_softclip.bam
                "


elif [ "$reads" == "single" ]; then

        # Extracting softclips
        cd ${indir}/softclip/

        nam=${sample}

        bsub -q medium -app large -J "v.int.3.softclip_vir_bam_extract_${nam}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
                -o "${indir}/softclip/${nam}.softclip_vir_bam_extract.log" -e "${indir}/softclip/${nam}.softclip_vir_bam_extract.err" \
                "python ${my_extract} -b ${fqdir}/full_viral/${nam}.bam -n ${indir}/softclip/${nam}_softclipped_reads.txt -o ${indir}/softclip/${nam}_softclipped_vir_temp.bam
                ${my_stools} sort ${indir}/softclip/${nam}_softclipped_vir_temp.bam > ${indir}/softclip/${nam}_softclipped_vir.bam; rm ${indir}/softclip/${nam}_softclipped_vir_temp.bam"


        bsub -q medium -app large -w 'ended("v.int.3.softclip_vir_bam_extract*")' -J "v.int.3.softclip_read_info_extract_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
                -o "${indir}/softclip/${nam}.softclip_read_info_extract.log" -e "${indir}/softclip/${nam}.softclip_read_info_extract.err" \
                "${my_stools} view ${indir}/softclip/${nam}_softclipped_vir.bam | cut -f 1,3,4,6,9 | sort | uniq > ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads_info.txt
                cut -f 1 ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads_info.txt | sort | uniq > ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads.txt
                mv ${indir}/softclip/${nam}_softclipped_vir.bam ${indir}/softclip/${nam}_unmapped_softclip.bam"

        bsub -q medium -app large -w 'ended("v.int.3.softclip_read_info_extract*")' -J "v.int.3.softclip_bam_extract_${nam}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
                -o "${indir}/softclip/${nam}.softclip_bam_extract.log" -e "${indir}/softclip/${nam}.softclip_bam_extract.err" \
                "python ${my_extract} -b ${indir}/softclip/${nam}_softclipped.bam -n ${indir}/softclip/${nam}_softclipped_mapped_filtered_reads.txt -o ${indir}/softclip/${nam}_viral_softclip.bam"


fi



# Merging the two while still in the softclip directory - single end stream has no mate bam
cd ${indir}/softclip/

base=${sample}
# Define temp sort directory
sort_dir=${proj_dir}/data/temp_sort_dirs/${base}/

if [ "$reads" == "paired" ]; then

        bsub -q short -app large -w 'ended("*_extract_*")' -J "v.int.3.merge_ms_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
            -o "${indir}/${base}.merge.ms.log" -e "${indir}/${base}.merge.ms.err" \
            "$my_stools merge -f ${indir}/${base}_merged.bam ${indir}/mate/${base}_viral_mate.bam ${indir}/softclip/${base}_viral_softclip.bam
            $my_stools sort ${indir}/${base}_merged.bam -o ${indir}/${base}_sorted.bam
	    ${my_stools} view ${indir}/${base}_sorted.bam | cut -f 1,3,4,6,9 | sort -T ${sort_dir} | uniq > ${indir}/${base}_reads_info.txt"

elif [ "$reads" == "single" ]; then

        bsub -q medium -app large -w 'ended("*_extract_*")' -J "v.int.3.merge_ms_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
            -o "${indir}/${base}.merge.ms.log" -e "${indir}/${base}.merge.ms.err" \
            "cp ${indir}/softclip/${base}_viral_softclip.bam ${indir}/${base}_merged.bam
            $my_stools sort ${indir}/${base}_merged.bam -o ${indir}/${base}_sorted.bam
	    ${my_stools} view ${indir}/${base}_sorted.bam | cut -f 1,3,4,6,9 | sort -T ${sort_dir} | uniq > ${indir}/${base}_reads_info.txt"

fi

#bsub -q express -app large -w 'ended("v.int.3.merge_*")' -J "v.int.3.host_read_info_extract_${f}" -n 4,4 -M 16GB -R "span[hosts=1] rusage[mem=16GB]" \
#    -o "${indir}/${nam}.host_read_info_extract.log" -e "${indir}/${nam}.host_read_info_extract.err" \
#    "${my_stools} view ${indir}/${base}_sorted.bam | cut -f 1,3,4,6,9 | sort -T ${sort_dir} | uniq > ${indir}/${base}_reads_info.txt"



## Submit output check, currently only implemented for paired data
if [ "$reads" == "paired" ]; then
    bsub -q express -app small -J "step3.check" \
    -w 'ended("v.int.3*")' \
    "${function_dir}/check_step3.sh \
    -d ${proj_dir}"
fi



# Exit
printf -- '\n';
exit 0;

