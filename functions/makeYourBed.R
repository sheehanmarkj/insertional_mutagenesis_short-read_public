#!/usr/bin/env Rscript
## take read\tcigar\tFLAG as input
input<-file('stdin', 'r')
while( length( row <- readLines(input, n=1) ) > 0) {
	read = strsplit(row,'\t')[[1]][1]
	cigar = strsplit(row,'\t')[[1]][2]
	flags = strsplit(row,'\t')[[1]][3]
	unrolled_cigar = strsplit(cigar, "(?<=[A-Z])", perl=TRUE)[[1]]
	if ( intToBits(flags)[5] == 01 ) { unrolled_cigar = rev(unrolled_cigar)  }
	loc_init = 0; loc_end = loc_init
	output=character()
	for ( i in 1:length(unrolled_cigar) ) {
		chunk = unrolled_cigar[i]
		letter = gsub("[0-9]","",chunk); nbase = as.numeric(gsub("[A-Z]","",chunk))
		if ( letter == "I" ) {
			output = c(output, paste0(read,"\t",loc_init,"\t",loc_end,"\n")) ## output accumulated matching bases, resume next line from end
			loc_init = loc_end + nbase # skip n bases 
			loc_end = loc_init
		} else if ( letter == "D" ) {
			next # do nothing, consumes ref
		} else if ( letter == "S" ) {
			output = c(output, paste0(read,"\t",loc_init,"\t",loc_end,"\n")) ## output accumulated matching bases, resume next line from end
			loc_init = loc_end + nbase # skip n bases 
			loc_end = loc_init
		} else if ( letter == "H" ) {
			loc_end = loc_end + nbase # maps elsewhere in same genome, so mask these too (e.g. vector rearrangement)
		} else if ( letter == "M" ) {
			loc_end = loc_end + nbase
		} else {
			stop(paste0("Unrecognized letter ",letter," in CIGAR string, read ",read))
		}
	}
	output = c(output, paste0(read,"\t",loc_init,"\t",loc_end,"\n"))
	cat(output) #write(output, stdout())
}

