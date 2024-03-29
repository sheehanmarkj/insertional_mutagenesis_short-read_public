library(GenomicRanges)
library(chromstaR)
library(GenomicFeatures)
library(pasillaBamSubset)
library(GenomicAlignments)
library(ggbio)
library(ORFik)
library(tidyverse)
library(exomeCopy)
library(RMariaDB)
library(pool)
library(data.table)
library(DOSE)
library(knitr)

loadIS = function(meta,
                   base.dir,
                   pattern.is = "sorted.bam$",
                   set,
                   do.is = TRUE,
                   do.vir.origin = FALSE,
                   do.read.info = FALSE,
                   sample.type.in.name = FALSE){
  
  # Section1: IS read and agglomeration
  if(do.is == TRUE){
    files<-list.files(path = paste0(base.dir,"/is"), pattern = pattern.is, full.names = TRUE)
    the.set <- set
    
    # Creating an empty list to add for each file
    is.list<-GRangesList()
    vir.list<-GRangesList()
    
    # Loading in the bam files that have IS reads for each sample
    for(i in files){
      # Extracting sample id from file name
      nam<-gsub(paste0('\\_',pattern.is),'',i)
      nam2<-gsub(".*is\\/","",nam)
      print(nam2)
      
      # accounting for spike samples
      #if(nam2 %like% "spike"){
       # nam2<-gsub("spike","",nam2)
      #}
      
      if(sample.type.in.name == TRUE){
        stype<-str_to_lower(unlist(str_split(nam2,pattern="_"))[2])
      }
      
      # Loading GRanges object for the current IS file
      is.ranges <- tryCatch(
        readBamFileAsGRanges(i, min.mapq = 0),
        error = function(cond){
          message(paste0("The file has no reads: ", nam2))
          return(c())
        } 
      )
      
      print(paste0("Done with is.ranges for ",nam2))
      
      if(class(is.ranges)=="GRanges"){
        # Subsetting metadata to current sample
        meta2<-meta %>% filter(sample_id == nam2)
        
        # Adding those metadata to GRanges object
        is.ranges$sample_id<-meta2$sample_id
        is.ranges$sex<-meta2$sex
        is.ranges$dose_group<-meta2$dose_group
        is.ranges$tissue<-meta2$tissue
        is.ranges$set<-set
        is.ranges$spike<-ifelse(nam %like% "spike",1,0)
        if(sample.type.in.name == TRUE){
          is.ranges$stype<-stype
        }
        
        
        # Adding to the list
        is.list<-c(is.list, GRangesList(is.ranges))
        names(is.list)[length(is.list)]<-nam2
        
      }
    }
    
    # Save as RDS to read later
    base.nam <- gsub(".*\\/(.*).*",'\\1',
                     gsub("\\/$","",base.dir))
    saveRDS(is.list, paste0(base.dir,base.nam,"_is.rds"))
  }
  
  # Section 2: Loading in viral origins
  if(do.vir.origin == TRUE){
    print("###############################")
    print("Starting viral origins section")
    
    ## Loading in files and 
    files<-list.files(path = paste0(base.dir,"/vir_orig"), pattern = ".bam$", full.names = TRUE)
    
    vir.list<-GRangesList()
    
    for(i in files){
      # Extracting sample id from file name
      nam<-gsub(paste0(base.dir,"/vir_orig/"),"",i)
      nam1<-gsub('\\_softclipped.*','',nam)
      nam2<-gsub('\\_unmapped.*','',nam1)
      
      print(paste0("Working on viral origins of", nam2))
      
      # Loading in the GRanges object for the current file
      vir.ranges <- tryCatch(
        readBamFileAsGRanges(i, min.mapq = 30),
        error = function(cond){
          message(paste0("The file has no reads: ", i))
          return(c())
        } 
      )
      
      if(class(vir.ranges)=="GRanges"){
        # Subsetting metadata to current sample
        meta2<-meta %>% filter(sample_id == nam2)
        
        # Adding those metadata to GRanges object
        vir.ranges$sample_id<-meta2$sample_id
        vir.ranges$sex<-meta2$sex
        vir.ranges$dose_group<-meta2$dose_group
        vir.ranges$tissue<-meta2$tissue
        vir.ranges$set<-set
        
        # Differentiating between mates vs softclip
        test.t<-gsub(".*softclipped.*","softclipped",i)
        test<-gsub(".*unmapped.*","mate",test.t)
        
        if(test == "softclipped"){
          vir.ranges$source<-"soft"
        }else{
          vir.ranges$source<-"mate"
        }
        
        # Adding to the list
        vir.list<-c(vir.list, GRangesList(vir.ranges))
        names(vir.list)[length(vir.list)]<-nam2
      }
    }
    
    # Save as RDS to read later
    base.nam <- gsub(".*\\/(.*).*",'\\1',
                     gsub("\\/$","",base.dir))
    saveRDS(vir.list, paste0(base.dir,base.nam,"_vir.rds"))
  }
  
  # Section 3: read info
  if(do.read.info == TRUE){
    # Section 3: Read info loading
    files<-list.files(path = paste0(base.dir,"/reads_info"), pattern = "\\_reads\\_info.txt$", full.names = TRUE)
    files<-files[grep("softclipped|unaligned",files,invert=TRUE)]
	
    read.list<-data.frame(read = c(),
                          chromosome = c(),
                          start = c(),
                          cigar = c(),
                          direction = c(),
                          source = c(),
                          type = c(),
                          width = c(),
                          sample_id = c())
    
    for(i in files){
      # Extracting sample id from file name
      nam<-gsub(paste0(base.dir,"/reads_info/"),"",i)
      nam1<-gsub('\\_reads\\_info.*','',nam)
      
      # Message
      print(paste0("Reading read info for: ", nam1))
      
      # Reading in files
      host.all <- if(file.exists(paste0(base.dir,"reads_info/",nam1,"_reads_info.txt"))){
        fread(paste0(base.dir,"reads_info/",nam1,"_reads_info.txt"),header=FALSE)
      }else{data.table(source = c(), type = c())}
      virus.soft <- if(file.exists(paste0(base.dir,"reads_info/",nam1,"_softclipped_mapped_filtered_reads_info.txt"))){
        fread(paste0(base.dir,"reads_info/",nam1,"_softclipped_mapped_filtered_reads_info.txt"),header=FALSE)
      }else{data.table(source = c(), type = c())}
      virus.mate <- if(file.exists(paste0(base.dir,"reads_info/",nam1,"_unaligned_mapped_filtered_reads_info.txt"))){
        fread(paste0(base.dir,"reads_info/",nam1,"_unaligned_mapped_filtered_reads_info.txt"),header=FALSE)
      }else{data.table(source = c(), type = c())}
      
      # Adding source annotation
      host.all$source <- "host"
      virus.soft$source <- "virus"
      virus.mate$source <- "virus"
      
      # adding type
      host.all$type <- NA
      virus.soft$type <- "soft"
      virus.mate$type <- "mate"
      
      # Removing missing files
      if(dim(host.all)[2]<6){host.all <- NULL}
      if(dim(virus.soft)[2]<6){virus.soft <- NULL}
      if(dim(virus.mate)[2]<6){virus.mate <- NULL}
      
      # Restricting viral reads to ones found in host
      if(!is.null(virus.soft)){
        virus.soft <- virus.soft %>% filter(V1 %in% host.all$V1)
      }
      if(!is.null(virus.mate)){
        virus.mate <- virus.mate %>% filter(V1 %in% host.all$V1)
      }
      
      # Binding and naming
      read.info <- rbind(host.all, virus.soft, virus.mate)
      if(!is.null(read.info) & dim(as.data.frame(read.info))[1]>0){
        colnames(read.info) <- c("read","chromosome","start","cigar","direction","source","type")
        
        # Adding cigar info and naming
        print("Parsing cigar")
        read.info<-cbind(read.info, width = cigar.parse(read.info$cigar)$tot.map)
        read.info$sample_id<-nam1
        
        # Removing all reads if host reads don't pass QC
        if(is.null(host.all)){read.info <- NULL}
        
        # Adding to list
        read.list<-rbind(read.list, read.info)
      }
    }
    
    # Save as RDS to read later
    base.nam <- gsub(".*\\/(.*)\\/.*\\/.*",'\\1',
                     gsub("\\/$","",base.dir))
    saveRDS(read.list, paste0(base.dir,base.nam,"_read.rds"))
    
  }
  print(paste0("##########################################"))
  print(paste0("Reading of ", set," completed successfully"))
}
