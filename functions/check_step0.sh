#!/bin/bash

#############################################################################
# Creating a usage sentence
usage="To use this step0 error check please specify the following:
$(basename $0) [-d project-directory] [-h] usage info"

## Printing it when there are no arguments
if [ "$#" == "0" ] || [ "$1" == "--help" ]; then
        printf -- "\033[33m$usage. \033[0m\n" >&2
        exit
fi


while getopts 'd:' OPTION; do
        case "$OPTION" in
                d) proj_dir="$OPTARG"
                printf -- "\033[32m Your project directory is: $OPTARG. \033[0m\n"
                ;;
        esac
done



grep "xit code" ${proj_dir}/data/raw_fastq/full_viral/*.log > ${proj_dir}/step0_check.txt

list=${proj_dir}/data/raw_fastq/sample_list.txt
for f in $(cat $list); do
	if [ ! -f ${proj_dir}/data/raw_fastq/${f}_trimmed.fastq.gz ]; then printf "${f}_trimmed.fastq.gz not found\n" >> ${proj_dir}/step0_check.txt; fi
	if [ ! -f ${proj_dir}/data/raw_fastq/${f}_1_trimmed_masked.fastq.gz ]; then printf "${f}_1_trimmed_masked.fastq.gz not found\n" >> ${proj_dir}/step0_check.txt; fi
	if [ ! -f ${proj_dir}/data/raw_fastq/${f}_2_trimmed_masked.fastq.gz ]; then printf "${f}_2_trimmed_masked.fastq.gz not found\n" >> ${proj_dir}/step0_check.txt; fi
	if [ ! -f ${proj_dir}/data/raw_fastq/full_viral/${f}.bam ]; then printf "${f}.bam not found\n" >> ${proj_dir}/step0_check.txt; fi
done


# Exit
printf -- '\n';
exit 0;

