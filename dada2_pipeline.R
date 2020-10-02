###DADA-2 16s-rRna analysis pipeline
#!/usr/bin/env Rscript

# https://benjjneb.github.io/dada2/tutorial.html

# Functions
printTocmd <- function(s){
	cat(paste(format(Sys.time(), "%Y-%m-%d__%X --->"),s,"\n"))
}
usage <- function() {
	cat("\nUSAGE: Rscriptdada2_pipeline.R <path (path contains file names 'fastq' with input files)> <refDB (silva,refseq,greengenes)> <ncores (2-16, empty = single-thread)>\n\n")
}
exitWithUsage <- function(err) {
	stop(paste(format(Sys.time(), "%Y.%m.%d-%X --->"),err,usage()), call.=FALSE)
}

printTocmd("Loading required packages...")

library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ape)
library(biomformat)

# Set script args
args = commandArgs(trailingOnly=TRUE)
path   <- args[1]
DB     <- args[2]

# Set config vars
fastq_path      <- paste0(path,"/fastq")
out_path        <- paste0(path,"/out/")
fastq_files     <- sort(list.files(fastq_path, pattern=".fastq", full.names=TRUE))
sample.names    <- sapply(strsplit(basename(fastq_files),".fastq"),`[`,1)
plots_path 	    <- paste0(out_path,"/plots")
fastqFilt_path  <- paste0(fastq_path,"/filtered")
fastqFilt_files <- file.path(fastqFilt_path, paste0(sample.names, "_filt.fastq.gz"))
availableDBs    <- c("silva","refseq","greengenes")


# Arguement consistency checks
if (!file.exists(path)){
	exitWithUsage("Input directory doesn't exist")
}
if (!file.exists(fastq_path)){
	exitWithUsage("Directory named 'fastq' containing input fastq files doesn't exist.")
}
if(!exists("DB") || !(DB %in% availableDBs)) {
	exitWithUsage("You must select a reference Database for taxonomy profiling among 'silva','refseq' or 'greengenes'")
}

printTocmd("Analysis has started.")

#QC
unlink(plots_path, force=TRUE)
dir.create(plots_path, showWarnings=FALSE)
png(filename=paste0(plots_path,"/QC.png"),  width=1200, height=1200)
plotQualityProfile(fastq_path, aggregate=TRUE)
dev.off()

printTocmd("Quality filtering fastQ files...")
unlink(list.files(fastqFilt_path), force=TRUE, recursive=TRUE)

names(fastqFilt_files) <- sample.names
out <- filterAndTrim(
	fastq_files, 
	fastqFilt_files, 
	trimLeft=15, # IonTorrent-specific suggested by authors
	# trimRight=40, # Trim low quality bases @read ends
	# truncLen=50, # Too many short reads are discarded with this one
	minLen=50, # Remove smaller than min allowed reads
	maxN=0, 
	maxEE=4,  # Lenient Expected Error - Default: 2
	truncQ=2, 
	rm.phix=TRUE,
	compress=TRUE, 
	multithread=16
	)

# Seq error rates
printTocmd("Denoising fastQ files...")
errF <- learnErrors(fastqFilt_files, multithread=16)
png(filename=paste0(plots_path,"/noise.png"),  width=1200, height=1200)
plotErrors(errF, nominalQ=TRUE)
dev.off()

dadaFs <- dada(fastqFilt_files, err=errF, multithread=16)

# Seq abundance table
printTocmd("Creating OTU abundance table...")
seqtab <- makeSequenceTable(as.list(dadaFs))

# Remove chimeras
printTocmd("Removing chimeric reads...")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=16, verbose=TRUE)
saveRDS(seqtab.nochim, file=paste0(out_path,"seqtab.rds"))

# Calculate remaining sequencies
getN <- function(x) sum(getUniques(x))
printTocmd("Writing reads filtering report...")
track <- cbind(out, sapply(dadaFs, function(x) sum(getUniques(x))), rowSums(seqtab.nochim))
colnames(track) <- c("input", "QC-filtered", "denoised", "chimeras-filtered")
rownames(track) <- sample.names
write.table(track, file = paste0(out_path,"data_trunc.tsv"))

if (DB == "silva") {
	printTocmd("Assigning taxonomy using Silva v132...") # ~6.5h runtime
	taxa <- assignTaxonomy(
		seqs=seqtab.nochim,
		refFasta=paste0("/home/makis/tools/metagen/DBs/silva_nr_v132_train_set.fa.gz"), 
		tryRC=TRUE,
		multithread=FALSE
	)
	# assign species (extra in Silva)
	taxa <- addSpecies(
		taxtab=taxa, 
		refFasta=paste0("/home/makis/tools/metagen/DBs/silva_species_assignment_v132.fa.gz")
	)
}

if (DB == "greengenes") {
	printTocmd("Assigning taxonomy using GreenGenes v13.8...") # Needs >50bp reads
	taxa <- assignTaxonomy(
		seqs=seqtab.nochim,
		refFasta=paste0("/home/makis/tools/metagen/DBs/gg_13_8_train_set_97.fa.gz"), 
		tryRC=TRUE,
		multithread=FALSE
	)
	# special formatting in GreenGenes database
	taxa <- gsub("[a-z]__", "", taxa)
}

if (DB == "refseq") {
	printTocmd("Assigning taxonomy using RefSeq...")
	taxa <- assignTaxonomy(
		seqs=seqtab.nochim,
		refFasta=paste0("/home/makis/tools/metagen/DBs/RefSeq-RDP16S_v2_May2018.fa.gz"), 
		tryRC=TRUE,
		multithread=FALSE
	)
}

samples.out <- sample.names
organ <- sapply(strsplit(samples.out, "_"), `[`, 1)
samdf <- data.frame(Subject=samples.out, Organ=organ)
rownames(samdf) <- samples.out

# PhyloSeq object
printTocmd("Creating the phyloseq object...")
ps <- phyloseq(
		otu_table(
			seqtab.nochim, 
			taxa_are_rows=FALSE),
		sample_data(samdf),
		tax_table(taxa)
		)

random_tree = rtree(
	ntaxa(ps), 
	rooted=TRUE, 
	tip.label=taxa_names(ps)
	)

ps <- merge_phyloseq(ps,random_tree)

# TESTING
# ps <- readRDS("~/silva_ps.rds")

printTocmd("Exporting phyloseq object .rds file...")
saveRDS(ps, file=paste0(out_path,DB,"_ps.rds"))

# Convert ASV to FASTA format with sequences as IDs
printTocmd("Generating functional annotation fasta file...")
uniquesToFasta(
	unqs=seqtab.nochim,
	fout=paste0(out_path,"asv.fasta"), 
	ids=colnames(seqtab.nochim)
)
# Convert ASV to biom format
printTocmd("Generating functional annotation biom file ...")
st.biom <- make_biom(t(seqtab.nochim))
write_biom(st.biom, paste0(out_path,"my.biom"))
