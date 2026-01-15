# Load libraries
library("HaplotypR")
library("ShortRead")
library("pwalign")

# Define output directory
outputDir <- "/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/CSP_TRAP_AMA1_ROUND3_No_CV109AorB"
# Create output directory
if(!dir.exists(outputDir))
  dir.create(outputDir, recursive=T)

# Specify number of pools to analyse
Poolnum <- 20

# Specify the pools that have CSP/TRAP amplicons in them
select <- "CVPool"

# Specify mismatch options
minMMrate <- 0.5
minOccGen <- 2

# Specify haplotype options
minCov <- 3
detectionLimit <- 1/200
minOccHap <- 1
minCovSample <- 25

# source edited funs
source("/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/deplexBySample_ext.R")
source("/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/checkBarcode_ext.R")
source("/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/createHaplotypeTable_ext.R")


##############Run demultiplexing by sample 
# Define a function to perform demultiplexing
demultiplex <- function(pattern, outputDir, sampleNum) {
  # Set input file path
  primerFile <- paste0("/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/markerFile_CSP_TRAP_AMA1.txt")
  sampleFile <- paste0("/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/samplefiles/", "sampleFile_", pattern, ".txt")
  fnBarcodeF <- paste0("/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/barcode_Fwd.fasta")
  fnBarcodeR <- paste0("/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/barcode_Rev.fasta")
  
  # Set reads
  all_reads <- list.files(
    "/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/CV_pf_DIV_amplicons_raw_fastq/",
    pattern = paste0("^", pattern, "_R[12]_001\\.fastq\\.gz$"),
    full.names = TRUE
  )

  read1 <- all_reads[grepl("_R1_", all_reads)]
  read2 <- all_reads[grepl("_R2_", all_reads)]

  if (length(read1) != 1 || length(read2) != 1) {
    stop(paste0(
      "Could not find exactly one R1 and one R2 for ", pattern,
      ". Found R1: ", paste(read1, collapse = ","),
      " | Found R2: ", paste(read2, collapse = ",")
    ))
  }

  # Create output subdirectory
  outDeplexSample <- file.path(outputDir, paste0("dePlexSample", sampleNum))
  dir.create(outDeplexSample, showWarnings = FALSE)
  
  # Demultiplex by samples - no SNPs or deletions allowed
  dePlexSample <- deplexSampleExtended(read1, read2, fnBarcodeF, fnBarcodeR, outDeplexSample, max.mismatch = 0, with.indels = FALSE)
  
  sampleTab <- read.delim(sampleFile, stringsAsFactors = FALSE)
  dePlexSample <- renameDemultiplexedFiles(sampleTab, dePlexSample)
  
  write.table(dePlexSample, file.path(outputDir, paste0("demultiplex", sampleNum, "SampleSummary.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
  dePlexSample <- na.omit(dePlexSample)
  # Save the dePlexSample dataframe to the R environment with a variable name based on sampleNum
  assign(paste0("dePlexSample", sampleNum), dePlexSample, envir = .GlobalEnv)
  # Save the primerFile to the R environment
  assign("primerFile", primerFile, envir = .GlobalEnv)
}

# Set the folder path to fastq files
folder_path <- "/mnt/storage11/sophie/plasmodium_cape_verde/DIV_amplicons/CV_pf_DIV_amplicons_raw_fastq"

# List all the files in the folder
file_list <- list.files(folder_path, full.names = TRUE)
file_list <- file_list[grep(select, file_list)]

# Define a function to extract the sample name from the file name
extract_sample_name <- function(file_name) {
  # Extract the pool ID from the filename, e.g.
  # "CVPool7_R1_001.fastq.gz" -> "CVPool7"
  file_name <- basename(file_name)
  sample_name <- sub("_R[12]_001\\.fastq\\.gz$", "", file_name)
  return(sample_name)
}
# Define function to run demultiplex for each pool
generate_and_execute_demultiplex_lines <- function(file_list, outputDir) {
  unique_names <- character(0)
  
  for (i in seq_along(file_list)) {
    file_name <- basename(file_list[i])
    sample_name <- extract_sample_name(file_name)
    unique_name <- sample_name  # e.g. "CVPool7"
    
    # Skip duplicate pools
    if (unique_name %in% unique_names) {
      next
    }
    
    # Add this pool to the processed set
    unique_names <- union(unique_names, unique_name)
    
    # Build command
    demultiplex_line <- paste0("demultiplex(\"", unique_name, "\", outputDir, ", i, ")")
    cat(demultiplex_line, "\n")
    
    # Execute
    eval(parse(text = demultiplex_line))
  }
}

# Call the function with file_list and outputDir
generate_and_execute_demultiplex_lines(file_list, outputDir)

## Merge dePlexSample data frames from all pools

# Create an empty list to store the dataframes
dePlexSampleList <- mget(ls(pattern = "^dePlexSample[0-9]+$"))
dePlexSample <- do.call(rbind, dePlexSampleList)

write.table(
  dePlexSample,
  file.path(outputDir, "demultiplexSampleSummary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

## Demultiplex by marker
# create output subdirectory
outDeplexMarker <- file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)

# process each marker
markerTab <- read.delim(primerFile, stringsAsFactors=F)
dePlexMarker <- demultiplexByMarker(dePlexSample, markerTab, outDeplexMarker)

# save summary table
write.table(
  dePlexMarker,
  file.path(outputDir, "demultiplexMarkerSummary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#remove lines without sampleID
dePlexMarker <- na.omit(dePlexMarker)


## Merge paired reads (R1/R2)

# create output subdirectory
outProcFiles <- file.path(outputDir, "processedReads")
dir.create(outProcFiles)

postfix <- "_merge"
refSeq <- DNAStringSet(markerTab$ReferenceSequence)
names(refSeq) <- markerTab$MarkerID


# Write a small FASTA per marker for reference

lapply(seq_along(refSeq), function(i) {
  writeFasta(
    refSeq[i],
    file.path(outputDir, paste0(names(refSeq)[i], postfix, ".fasta"))
  )
})

procReadsMerge <- mergeAmpliconReads(
  fastqFileR1 = as.character(dePlexMarker$FileR1),
  fastqFileR2 = as.character(dePlexMarker$FileR2),
  outputDir   = outProcFiles,
  method      = "vsearch"
)

###
procReads <- cbind(
  dePlexMarker[, c("SampleID", "SampleName", "BarcodePair", "MarkerID")],
  procReadsMerge
)

write.table(
  procReads,
  file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ======================================================================================================================================================
# DOWN-SAMPLING
# ======================================================================================================================================================

library(ShortRead)

# Count reads helper (unchanged)
count_fastq_reads <- function(fq_gz) {
  con <- gzfile(fq_gz, open = "r")
  on.exit(close(con))
  n_lines <- 0L
  repeat {
    chunk <- readLines(con, n = 100000)
    if (length(chunk) == 0) break
    n_lines <- n_lines + length(chunk)
  }
  as.integer(n_lines / 4L)
}

## Create downsampling function

downsample_fastq_gz <- function(in_fq, out_fq,
                                target_reads,
                                seed = 1L,
                                verbose = TRUE) {

  if (!file.exists(in_fq)) stop(sprintf("Input FASTQ not found: %s", in_fq))
  if (is.null(target_reads) || target_reads < 1) stop("target_reads must be >= 1")

  set.seed(seed)

  # Count reads
  total <- count_fastq_reads(in_fq)
  if (verbose) message(sprintf("  %s has %d reads", basename(in_fq), total))

  # Zero reads → do nothing
  if (total == 0L) {
    if (verbose) message("    (0 reads) — skipping, no output file created")
    return(invisible(0L))
  }

  # Below target → copy file
  if (total < target_reads) {
    if (verbose) message(sprintf("    (<%d) — copying instead of downsampling", target_reads))
    if (!file.copy(in_fq, out_fq, overwrite = TRUE))
      stop(sprintf("Failed to copy %s -> %s", in_fq, out_fq))
    return(invisible(total))
  }

  # Load all reads at once (fast for small FASTQs)
  fq <- ShortRead::readFastq(in_fq)

  # Randomly sample target reads
  idx <- sample(seq_len(length(fq)), target_reads)

  ShortRead::writeFastq(
    fq[idx],
    file = out_fq,
    compress = TRUE
  )

  invisible(target_reads)
}

## ---- Downsample merged FASTQs----

downsampled_dir <- file.path(outputDir, "downsampled_processed_reads")
dir.create(downsampled_dir, showWarnings = FALSE, recursive = TRUE)

downsample_target_reads <- 1000L   # change if needed
downsample_seed         <- 42L

kept_vec  <- integer(nrow(procReads))
out_files <- character(nrow(procReads))

message("Downsampling FASTQs into: ", downsampled_dir)

for (i in seq_len(nrow(procReads))) {

  in_fq  <- procReads$ReadFile[i]
  base   <- sub("\\.fastq\\.gz$", "", basename(in_fq))
  out_fq <- file.path(downsampled_dir, paste0(base, ".downsampled.fastq.gz"))

  kept_vec[i] <- downsample_fastq_gz(
    in_fq,
    out_fq,
    target_reads = downsample_target_reads,
    seed = downsample_seed,
    verbose = TRUE
  )

  out_files[i] <- if (file.exists(out_fq)) out_fq else NA

  message(sprintf("[%d/%d] %s -> %s (kept %d reads)",
                  i, nrow(procReads),
                  basename(in_fq), basename(out_fq), kept_vec[i]))
}

# Update procReads
procReads$ReadFile <- out_files
procReads$numRead  <- vapply(procReads$ReadFile,
                             function(f) if (!is.na(f)) count_fastq_reads(f) else 0L,
                             integer(1))
# Include low-depth filter
procReads <- procReads[procReads$numRead > 10, ]

# =============================================
# Added downsampling chunk has ended

## =============================================
## Calculate mismatch rate and call SNPs

# process each marker
snpLst <- lapply(markerTab$MarkerID, function(marker){
  # Calculate mismatch rate
  seqErrLst <- calculateMismatchFrequencies(as.character(procReads[procReads$MarkerID == marker, "ReadFile"]),
                                            refSeq[marker],
                                            method ="pairwiseAlignment", # c("pairwiseAlignment","compareDNAString"),
                                            minCoverage=100L)
  names(seqErrLst) <- procReads[procReads$MarkerID == marker, "SampleID"]
  seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, file.path(outputDir, sprintf("mismatchRate_rate_%s%s.txt", marker, postfix)), sep="\t", row.names=F)

  # Call SNPs
  potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start=snp, width=1))
  }))
  snps <- data.frame(Chr=marker, Pos=potSNP, Ref=snpRef, Alt="N", stringsAsFactors=F)
  write.table(snps, file=file.path(outputDir, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt",
                                                      minMMrate*100, minOccGen, marker, postfix)),
              row.names=F, col.names=T, sep="\t", quote=F)

  # Plot mismatch rate and SNP calls
  png(file.path(outputDir, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png",
                                   minMMrate*100, minOccGen, marker, postfix)),
      width=1500 , height=600)
  matplot(seqErr, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
          ylab="Mismatch Rate", xlab="Base Position", main=marker, cex.axis=2, cex.lab=2)
  abline(v=snps[,"Pos"], lty=2, col="grey")
  abline(h=minMMrate, lty=1, col="red")
  dev.off()

  return(snps)
})

names(snpLst) <- markerTab$MarkerID

## Call haplotypes

# call final haplotypes
finalTab <- createFinalHaplotypTableExtended(
  outputDir = outputDir,
  sampleTable = procReads,
  markerTable = markerTab,
  referenceSeq = refSeq,
  snpList = snpLst,
  postfix = postfix,
  minHaplotypCoverage = minCov, # ≥X reads per haplotype
  minReplicate = minOccHap, # halotype must recur in ≥X replicates/samples
  detectability = detectionLimit, # ≥X% within-sample frequency
  minSampleCoverage = minCovSample # ≥X reads per sample to be kept
)

CSP_finalTab<-finalTab[["csp"]]
TRAP_finalTab<-finalTab[["trap"]]
AMA1_finalTab <- finalTab[["ama1"]]

# Write raw per-marker tables
write.table(CSP_finalTab,  file.path(outputDir, "CSP_finalTab.txt"),  sep="\t", row.names=FALSE, quote=FALSE)
write.table(TRAP_finalTab, file.path(outputDir, "TRAP_finalTab.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(AMA1_finalTab, file.path(outputDir, "AMA1_finalTab.txt"), sep="\t", row.names=FALSE, quote=FALSE)


# Filter: keep haplotypes present in BOTH technical replicates of a specimen
library(dplyr)
artefacts <- c("Noise","Indels","Singelton","Chimera")  # fixed typo; handle NA separately

filter_replicates <- function(df) {
  df %>%
    filter(!is.na(Haplotype),
           !Haplotype %in% artefacts,
           Reads > 0) %>%
    mutate(Specimen = sub("_Rep[0-9]+.*$", "", SampleID)) %>%
    group_by(Specimen, Haplotype) %>%
    filter(n_distinct(SampleID) >= 2) %>%   # seen in ≥2 replicates for that Specimen
    ungroup()
}

# Read back from outputDir (matches where you wrote them), filter, write results
CSP_filtered  <- read.delim(file.path(outputDir, "CSP_finalTab.txt"),  stringsAsFactors = FALSE) |> filter_replicates()
TRAP_filtered <- read.delim(file.path(outputDir, "TRAP_finalTab.txt"), stringsAsFactors = FALSE) |> filter_replicates()
AMA1_filtered <- read.delim(file.path(outputDir, "AMA1_finalTab.txt"), stringsAsFactors = FALSE) |> filter_replicates()

write.table(CSP_filtered,  file.path(outputDir, "CSP_finalTab_passed_technical_replicates.txt"),  sep="\t", row.names=FALSE, quote=FALSE)
write.table(TRAP_filtered, file.path(outputDir, "TRAP_finalTab_passed_technical_replicates.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(AMA1_filtered, file.path(outputDir, "AMA1_finalTab_passed_technical_replicates.txt"), sep="\t", row.names=FALSE, quote=FALSE)