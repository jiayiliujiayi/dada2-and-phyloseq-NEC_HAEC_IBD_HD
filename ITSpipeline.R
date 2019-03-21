# 需要显示版本信息
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")


# 更换目录
# CHANGE ME to the directory containing the fastq files after unzipping.
getwd()
path <- "/Users/JiayiLiu/Documents/compare/raw_data_processing/eachsample"
# list.files返回指定目录中的文件名
list.files(path)
# pattern指定返回指定类型的文件
# full.names = TRUE返回带有路径的完整文件名
# 返回测序正向文件完整文件名
fnFs <- sort(list.files(path, pattern = "1.fq.gz", full.names = TRUE))
# 返回测序反向文件完整文件名
fnRs <- sort(list.files(path, pattern = "2.fq.gz", full.names = TRUE))


# Identify primers --------------------------------------------------------
# 识别primer
#define your forward primer sequence
#定义上游下游引物的序列
FWD <- "ACTCCTACGGGAGGCAGCAG" ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVHHHTWTCTAAT"

#verify the existence in your data
#确认上下游引物是否存在于你的测序数据中
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#pre-过滤掉模糊碱基(N)
#pre-filter ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#count the number of times the primers appear in the forward and reverse read, 
#while considering all possible primer orientations.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


# Remove Primers ----------------------------------------------------------
#use cutadapt
cutadapt <- "/Users/JiayiLiu/miniconda2/bin/cutadapt" #去cutadapt的网站上有安装教程
system2(cutadapt, args = "--version")
###If the above command succesfully executed, 
#R has found cutadapt and you are ready to continue following along.

#create output filenames for the cutadapt-ed files & define the parameters 
# critical parameters = primers, and they need to be in the right orientation, i.e. the FWD primer should have been matching the forward-reads in its forward orientation, and the REV primer should have been matching the reverse-reads in its forward orientation. 
# Warning: A lot of output will be written to the screen by cutadapt!
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

##sanity check: count the number of primers in the data
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# DADA2 pipeline ----------------------------------------------------------
## read the names of cutadapt-ed files
## & string manipulation to get forward and reverse fastq files
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "1.fq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "2.fq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format:
sample.names <- gsub(pattern = "raw.split.",
                    replacement = "", 
                    basename(fnFs))
sample.names <- gsub(".1.fq.gz", "", sample.names)
head(sample.names)

##Inspect read quality profiles
###visualize the quality of forward reads
plotQualityProfile(cutFs[1:2])
###visualize the quality of reverse reads
plotQualityProfile(cutRs[1:2])

##Filter and trim
###Assigning the filenames for the output of the filtered reads
### to be stored as .fq.gz files
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
# note begin ===
### note: use standard filtering paraments: 
### maxN=0 (DADA2 requires sequences contain no Ns), truncQ = 2, rm.phix = TRUE and maxEE=2
### maxEE = maximum number of “expected errors” allowed in a read
# note end ===
###start filter & trim
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)


# Learn the Error Rates ---------------------------------------------------
# note begin ===
#Please ignore all the “Not all sequences were the same length.” messages 
#in the next couple sections. 
#We know they aren’t, and it’s OK!
# note end ===

## learn the error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

## plot error rates as a "sanity check"
plotErrors(errF, nominalQ = TRUE)


# Dereplicate identical reads ---------------------------------------------

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
## Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# Sample Inference --------------------------------------------------------

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


# Merge paired reads ------------------------------------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


# Construct Sequence Table ------------------------------------------------
# notes begin ===
# consruct a sequence variant table (ASV) table
# notes end ===
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths: no need for ITS pipeline 
table(nchar(getSequences(seqtab)))


# Remove chimeras ---------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))


# Track reads through the pipeline ----------------------------------------
# notes begin ===
#inspect the the number of reads that made it through each step in the pipeline 
#to verify everything worked as expected.
# notes end ===

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)


# Assign taxonomy ---------------------------------------------------------

# notes begin ===
## utilize the silva reference 
## website: https://zenodo.org/record/1172783#.XJDfEC2B0Wo
# notes end ===

taxa <- assignTaxonomy(seqtab.nochim, "silva/silva_nr_v132_train_set.fa.gz", multithread=TRUE)


#add species
time1 <- Sys.time()
taxa <- addSpecies(taxa, "silva/silva_species_assignment_v132.fa.gz")
time2 <- Sys.time()
time = time2 - time1 

taxa.print <- taxa 
# 另存物种注释变量，去除序列名，只显示物种信息
# Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# save files --------------------------------------------------------------

setwd(path)
dir.create("output")
setwd("output")
seqtable.taxa.plus <- cbind('#seq'=rownames(taxa), t(seqtab.nochim), taxa)
# ASV表格导出
write.table(seqtab.nochim, "dada2_counts.txt", sep="\t", quote=F, row.names = T)
# 带注释文件的ASV表格导出
write.table(seqtable.taxa.plus , "dada2_counts.taxon.species.txt", sep="\t", quote=F, row.names = F)
# track文件保存
write.table(track , "dada2_track.txt", sep="\t", quote=F, row.names = F)
# save RDdata
save(list = ls(all=TRUE), file = "dada2.RData")
