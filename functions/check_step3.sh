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

grep "xit code" ${proj_dir}/data/trimalign/*merge.ms.log > ${proj_dir}/step3_check.txt
#grep "xit code" ${proj_dir}/data/trimalign/*host_read_info_extract.log >> ${proj_dir}/step3_check.txt
grep "xit code" ${proj_dir}/data/trimalign/mate/*.mate*.log >> ${proj_dir}/step3_check.txt
grep "xit code" ${proj_dir}/data/trimalign/softclip/*.softclip*.log >> ${proj_dir}/step3_check.txt


list=${proj_dir}/data/raw_fastq/sample_list.txt

for f in $(cat $list); do
	if [ ! -f ${proj_dir}/data/trimalign/${f}_sorted.bam ]; then printf "${f}_sorted.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/${f}_merged.bam ]; then printf "${f}_merged.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/${f}_reads_info.txt ]; then printf "${f}_reads_info.txt not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/mate/${f}_unmapped_mate_1.bam ]; then printf "${f}_unmapped_mate_1.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/mate/${f}_unmapped_mate_2.bam ]; then printf "${f}_unmapped_mate_2.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/mate/${f}_unmapped_mate.bam ]; then printf "${f}_unmapped_mate.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/mate/${f}_viral_mate.bam ]; then printf "${f}_viral_mate.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/mate/${f}_unaligned_mapped_filtered_reads.txt ]; then printf "${f}_unaligned_mapped_filtered_reads.txt not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/mate/${f}_unaligned_mapped_filtered_reads_info.txt ]; then printf "${f}_unaligned_mapped_filtered_reads_info.txt not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/softclip/${f}_softclipped_vir_1.bam ]; then printf "${f}_softclipped_vir_1.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/softclip/${f}_softclipped_vir_2.bam ]; then printf "${f}_softclipped_vir_2.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/softclip/${f}_unmapped_softclip.bam ]; then printf "${f}_unmapped_softclip.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/softclip/${f}_viral_softclip.bam ]; then printf "${f}_viral_softclip.bam not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/softclip/${f}_softclipped_mapped_filtered_reads.txt ]; then printf "${f}_softclipped_mapped_filtered_reads.txt not found\n" >> ${proj_dir}/step3_check.txt; fi
	if [ ! -f ${proj_dir}/data/trimalign/softclip/${f}_softclipped_mapped_filtered_reads_info.txt ]; then printf "${f}_softclipped_mapped_filtered_reads_info.txt not found\n" >> ${proj_dir}/step3_check.txt; fi
done


if [ ! -s "${proj_dir}/step3_check.txt" ]; then
	# Moving to sub-directories
	indir=${proj_dir}/data/trimalign
	outdir=${proj_dir}/data/secondary
	
	## Sorted IS files to secondary analysis
	mkdir -p ${outdir}/is
	cp ${indir}/*sorted.bam ${outdir}/is/
	
	## Moving viral info reads
	mkdir -p ${outdir}/reads_info
	cp ${indir}/*reads_info.txt ${outdir}/reads_info
	cp ${indir}/mate/*reads_info.txt ${outdir}/reads_info
	cp ${indir}/softclip/*reads_info.txt ${outdir}/reads_info
	
	## Vir origins
	mkdir -p ${outdir}/vir_orig/
	cp ${indir}/mate/*unmapped_mate.bam ${outdir}/vir_orig/
	cp ${indir}/softclip/*unmapped_softclip.bam ${outdir}/vir_orig/
fi


# Exit
printf -- '\n';
exit 0;

