# This script analyses next-generation sequencing data.  The FASTQ files from various NGS samples were previously aligned 
# to a set of "reference" sequences using Magic BLAST.  These alignments are used as the input to this script, which will 
# determine the number of reads for each of the reference sequences.  Across biological replicates of similar samples, 
# the average reads and error are determined.

closeAllConnections()
rm(list=ls())

library(Rsamtools)
library(stringr)

#Function list:
#merge_dfs: Utility function to take a list of dataframes and merge them by a column.
#preprocess: Take a BAM file input, pull out the desired data into a dataframe, and calculate other data.
#removedata: Take the dataframe and filter out undesired entries.
#getfreqs: Convert a BAM-derived dataframe into a dataframe of sequence counts.

#procsample: Wrapper that reads in BAM files (ext, nc1, nc2) and then runs (1) preprocess and (2) removedata.
#            Returns the a list of the three count tables (one per BAM file).
#add_frequencies: Take the list of 3 count dataframes and calculates the total counts and frequencies across all three.
#calculate_freqs: Wrapper that runs procsample, getfreqs, and add_frequencies for a sample.
#analyze_sampleset: Function to load the data for a specific sample and return the total frequencies

setwd("~/path/to/datafiles")

#This function takes a list of dataframes and a column name common to all of these dataframes.
#The dfs are iteratively merged using the designated column to do so.
#The merged df is returned.
merge_dfs <- function(dflist, by) {
  #Set a first_df flag
  first_df <- TRUE
  #Loop over elements in the df list.
  for (curdf in names(dflist)) {
    #Is this the first df?
    if (first_df) {
      #Set mgdf as the first df, reset the first_df flag, and continue loop.
      mgdf <- dflist[[curdf]]
      first_df <- FALSE
      next
    }
    
    #Otherwise, mgdf should become the merge of mgdf and dflist[[curdf]]
    mgdf <- merge(mgdf, dflist[[curdf]], by=by)
  }
  
  return(mgdf)
}

#Function that takes a raw bam file and converts it to a dataframe containing the following values:
#rname, cigar, expcigar (ie. cigar numbers in a list), Scount (number of 'S' in the cigar string), qname
preprocess <- function(datalist) {
  #Create a dataframe with the desired columns.
  datadf <- data.frame(datalist$rname, datalist$cigar, datalist$flag, as.character(datalist$cigar), str_count(datalist$cigar, pattern = 'S'), datalist$qname)
  #Give the df columns names.
  colnames(datadf) <- c("rname", "cigar", "flag", "expcigar", "Scount", "qname")
  
  #Convert the expanded cigar column into a list split by characters, and then turn the numbers numeric.
  datadf$expcigar <- strsplit(as.character(datadf$expcigar), '[A-Z]')
  datadf$expcigar <- sapply(datadf$expcigar, as.numeric)
  
  return(datadf)
}

#Function to subset a dataframe dataset.  All NA records are removed, as are alignments
#that do not meet the quality standard.
removedata <- function(datadf) {
  #Remove the NA records from datadf
  subdata <- subset(datadf, !is.na(rname))
  
  #Remove data with bad cigar string characters in them (I,D,N,H,P,X).
  subdata <- subdata[grep("[IDNHPX]", subdata$cigar, invert = TRUE),]
  
  #Remove data with flag values of 4 (unmapped), 256 (secondary alignment), or 272 (secondary alignment + reverse complement).
  subdata <- subset(subdata, flag != 4 & flag != 256 & flag != 272)
  
  #Select rows that meet our condition re: leading and trailing 'S'
  keeprows <- apply(subdata, 1, function(ngsdata) {
    #If Scount is <2, keep the data.
    if (ngsdata$Scount < 2) {return(TRUE)}
    #If the first and last skips are less than 70, then keep.  Otherwise, drop.
    if (ngsdata$expcigar[1] < 70 & ngsdata$expcigar[length(ngsdata$expcigar)] < 70) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  
  #Subset the data using the keeprows values
  subdata <- subdata[keeprows,]
  
  return(subdata)
}

#Function to convert a dataframe with rname in it to a frequency table.
getfreqs <- function(dataset) {
  freqtable <- as.data.frame(table(dataset$rname))
  colnames(freqtable) <- c('design', 'Freq')
  
  return(freqtable)
}

#Function to open the three files for each sample, process them,
procsample <- function(sample) {
  cursample <- sprintf("%02d", sample)
  samplenames <- c("extendedFrags", "notCombined_1", "notCombined_2")
  
  #Create data, which contains a list of the three fastq alignments for each sample.
  data <- sapply(samplenames, function(sample) {
    filepath <- paste0("aligned/BY", cursample, ".", sample, ".bam")
    message(paste("Reading", filepath))
    dataset <- scanBam(filepath)
    message(paste("Done reading", filepath))
    return(dataset)
  })
  
  message(paste("Preprocessing Sample", cursample))
  fulldata <- lapply(data, preprocess)
  message(paste("Filtering Sample", cursample))
  subdata <- lapply(fulldata, removedata)
  return(subdata)
}

#Function to take the three frequency files, merge them to one dataframe, and then add the frequencies.
add_frequencies <- function(infreqs, sampleid='') {
  #Merge the first two frequency tables.
  merged_freqs <- merge(infreqs$extendedFrags, infreqs$notCombined_1, by='design')
  #Merge that with the last frequency table.
  merged_freqs <- merge(merged_freqs, infreqs$notCombined_2, by='design')
  
  merged_freqs$SeqCount <- merged_freqs$Freq.x + merged_freqs$Freq.y + merged_freqs$Freq
  merged_freqs <- merged_freqs[,c(1,5)]

  #Calculate the total number of sequences.
  total_seqs <- sum(merged_freqs$SeqCount)
  
  #If the specific percentage of an individual sequence is greater than 1%, make a new column where it is re-set to 0.
  #This is to be able to consider the frequency data without the overwhelmingly dominant sequences skewing the data.
  SeqCount_small <- sapply(merged_freqs$SeqCount, function(X, coltotal) {
    if (X/coltotal > 0.01) {
      return(0)
    } else {
      return(X)
    }
  }, coltotal = total_seqs)
  
  #Calculate the percentage of total sequences present for each sequence.
  merged_freqs$SeqPct <- merged_freqs$SeqCount / total_seqs
  merged_freqs$small_SeqPct <- SeqCount_small / sum(SeqCount_small, na.rm = TRUE)

  colnames(merged_freqs) <- c('design', paste0('SeqCount', sampleid), paste0('SeqPct', sampleid), paste0('small_SeqPct', sampleid))
  
  return(merged_freqs)
}

#Function to load the data for a specific sample and return the total frequencies
calculate_freqs <- function(sample) {
  #Return a list of the data for the sample.  The list contains the three datasets (extended,
  #nc1, nc2).  Each dataset is a dataframe containing the extracted data from preprocess.
  data <- procsample(sample)
  
  #Generate the frequency tables for the individual datasets (extended, nc1, nc2).
  freqs <- lapply(data, getfreqs)
  #Double the counts for extendedFrags to reflect that it is two combined sequences.
  freqs$extendedFrags$Freq <- freqs$extendedFrags$Freq * 2
  
  #Sum the frequencies to get a master frequency table.
  freqsum <- add_frequencies(freqs, sample)
  
  return(freqsum)
}

#Function to load a set of samples, generate frequencies, and calculate average/stdev.
analyze_sampleset <- function(sample_list, sample_name = "") {
  freq_list <- lapply(sample_list, calculate_freqs)
  names(freq_list) <- paste0('S', sample_list)

  #Merge the dataframes in freq_list into a single dataframe.
  mgfreqs <- merge_dfs(freq_list, by = 'design')
  
  #Figure out which columns are either design or a SeqPct column.
  pctcols <- grep("^SeqPct", colnames(mgfreqs))
  descols <- grep("design", colnames(mgfreqs))
  targetcols <- c(descols, pctcols)
  
  #Copy over targetcols to a new variable, replacing 0 with NA.
  pctna <- replace(mgfreqs[targetcols], mgfreqs[targetcols] == 0, NA)
  
  #Remove pctcols from mgfreqs.
  mgfreqs <- mgfreqs[,-pctcols]
  
  #Merge pctna with mgfreqs by design.
  mgfreqs <- merge(mgfreqs, pctna, by = 'design')
  
  #Calculate the average frequency (percentage) and stdev for each design.
  #Generate the new column names
  colname.avr <- paste0("avr_", sample_name)
  colname.sd <- paste0("stdev_", sample_name)
  colname.pcterr <- paste0("pcterror_", sample_name)
  colname.smallavr <- paste0("small_avr_", sample_name)
  colname.smallsd <- paste0("small_stdev_", sample_name)
  colname.smallpcterr <- paste0("small_pcterror_", sample_name)
  
  #For robustness, re-determine pctcols.
  pctcols <- grep("^SeqPct", colnames(mgfreqs))
  #Calculate the mean on a per-row basis.
  mgfreqs[[colname.avr]] <- apply(mgfreqs, 1, function(design, targetcols = pctcols) {
    mean(as.numeric(design[targetcols]), na.rm = TRUE)
    })
  smallpctcols <- grep("^small_SeqPct", colnames(mgfreqs))
  mgfreqs[[colname.smallavr]] <- apply(mgfreqs, 1, function(design, targetcols = smallpctcols) {
    mean(as.numeric(design[targetcols]), na.rm = TRUE)
  })
  
  #For robustness, re-determine pctcols.
  pctcols <- grep("^SeqPct", colnames(mgfreqs))
  mgfreqs[[colname.sd]] <- apply(mgfreqs, 1, function(design, targetcols = pctcols) {
    sd(as.numeric(design[targetcols]), na.rm = TRUE)
  })
  smallpctcols <- grep("^small_SeqPct", colnames(mgfreqs))
  mgfreqs[[colname.smallsd]] <- apply(mgfreqs, 1, function(design, targetcols = smallpctcols) {
    sd(as.numeric(design[targetcols]), na.rm = TRUE)
  })
  
  #Calculate the percent error based on stdev
  mgfreqs[[colname.pcterr]] <- mgfreqs[[colname.sd]] / mgfreqs[[colname.avr]]
  mgfreqs[[colname.smallpcterr]] <- mgfreqs[[colname.smallsd]] / mgfreqs[[colname.smallavr]]
  
  return(mgfreqs)
}

sample_list <- list(
  c(1:4), #"MTX0",
  c(5,7), #"MTX30",
  c(9:11), #"MTX50",
  c(12:14), #"MTX100",
  c(15:16), #"FullLibrary",
  c(17:18)  #"PCRProduct"
)

sample_names <- c(
  "MTX0",
  "MTX30",
  "MTX50",
  "MTX100",
  "FullLibrary",
  "PCRProduct"
)
 
names(sample_list) <- sample_names

all_freqs <- mapply(analyze_sampleset, sample_list = sample_list, sample_name = sample_names, SIMPLIFY = FALSE)

freq_table <- merge_dfs(all_freqs, by = 'design')

outfile <- "full_frequency_table_small.csv"
write.csv(freq_table, file = outfile, quote = FALSE)
