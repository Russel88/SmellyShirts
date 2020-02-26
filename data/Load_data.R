library(dada2)
library(foreach)
library(phyloseq)
library(ape)
library(ips)

path <- paste0(getwd(),"/Fastq.trim")

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))

# Extract sample names
sample.names <- gsub("Order150-","",sapply(strsplit(fnFs, "_"), `[`, 1))

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Filter and learn error 
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, compress=TRUE, trimLeft = 6)
errF <- learnErrors(filtFs)
errR <- learnErrors(filtRs)

# Derep
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
  
names(derepFs) <- sample.names
names(derepRs) <- sample.names
  
# Inference
dadaFs <- dada(derepFs, err=errF, multithread = 7L)
dadaRs <- dada(derepRs, err=errR, multithread = 7L)

# Merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 5)

# SeqTable
seqtab <- makeSequenceTable(mergers)

# Chimeras
seqtab.nc <- removeBimeraDenovo(seqtab, method="consensus", multithread = 7L, verbose=TRUE)

# Trim by length
newtab <- as.data.frame(t(seqtab.nc[,nchar(colnames(seqtab.nc)) > 380]))

# Remove singletons
newtab2 <- newtab[rowSums(newtab) > 1, ]

# Taxonomy
set.seed(100) # Initialize random number generator for reproducibility
taxa <- assignTaxonomy(rownames(newtab2), "rdp_train_set_16_trim.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "rdp_species_assignment_16_trim.fa.gz", allowMultiple = TRUE, verbose=TRUE)

# Remove unclassified
unc <- rownames(taxa[is.na(taxa[,3]),])
newtab3 <- newtab2[!rownames(newtab2) %in% unc,]

taxa <- taxa[rownames(taxa) %in% rownames(newtab3),]

# Fix taxonomy
for(i in 1:nrow(taxa)){
  if(is.na(taxa[i,2])){
    taxa[i,2:7] <- paste0("Kingdom.",taxa[i,1])
  } else {
    if(is.na(taxa[i,3])){
      taxa[i,3:7] <- paste0("Phylum.",taxa[i,2])
    } else {
      if(is.na(taxa[i,4])){
        taxa[i,4:7] <- paste0("Class.",taxa[i,3])
      } else {
        if(is.na(taxa[i,5])){
          taxa[i,5:7] <- paste0("Order.",taxa[i,4])
        } else {
          if(is.na(taxa[i,6])){
            taxa[i,6:7] <- paste0("Family.",taxa[i,5])
          } else {
            if(is.na(taxa[i,7])){
              taxa[i,7] <- paste0("Genus.",taxa[i,6])
            } else {
              taxa[i,7] <- paste0(taxa[i,6],".",taxa[i,7])
            }
          }
        }
      }
    }
  }
}

# Fix naming
seqs <- rownames(newtab3)
rownames(newtab3) <- paste0("ASV_",1:nrow(newtab3))
rownames(taxa) <- paste0("ASV_",1:nrow(taxa))
colnames(newtab3) <- sample.names

# Sequence objects
ref <- Biostrings::DNAStringSet(seqs)
names(ref) <- rownames(taxa)
seq.dna <- as.DNAbin(ref)
write.fas(seq.dna, file = "seqs.fasta")

### Run the Bowtie2.txt script in UNIX terminal externally to map ASVs to human genome
human <- read.table("mapped.txt")$V1

# Remove potential human seqs
ref <- ref[!names(ref) %in% human]
seq.dna <- seq.dna[!names(seq.dna) %in% human]
newtab4 <- newtab3[!rownames(newtab3) %in% human,]
taxa <- taxa[!rownames(taxa) %in% human,]

# Reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nc), colSums(newtab), colSums(newtab2), colSums(newtab3), colSums(newtab4))
colnames(track) <- c("denoised", "merged", "tabled", "nonchim", "noshort", "nosingleton", "nounclassified", "nohuman")
colSums(track)/colSums(track)[1]

# Make tree
align <- mafft(seq.dna, path = "C:/mafft/mafft-win/mafft.bat", method = "globalpair", maxiterate = 1000)
write.fas(align, file = "alignment.fasta")
system2("C:/Program Files/FastTree/FastTree.exe", args = c("-nt", "-out tree","alignment.fasta"))
tree <- read.tree("tree")

# Sample data
samp <- read.table("sample_data.csv",sep = ";", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

## Fix fabric
samp[samp$Fabric == "cot","Fabric"] <- "Cotton"
samp[samp$Fabric == "pol","Fabric"] <- "Polyester"
samp[is.na(samp$First_Tshirt),"First_Tshirt"] <- "none"
samp[samp$First_Tshirt == "cot","First_Tshirt"] <- "Cotton"
samp[samp$First_Tshirt == "pol","First_Tshirt"] <- "Polyester"

## Fix treatment
samp[samp$Treatment == "dirt","Treatment"] <- "Worn"
samp[samp$Treatment == "wash","Treatment"] <- "Worn + Washed"
samp[samp$Treatment == "prewash","Treatment"] <- "Unworn"
samp[samp$Treatment == "prewash.wash","Treatment"] <- "Unworn + Washed"

samp$Worn <- NA
samp[samp$Treatment %in% c("Worn","Worn + Washed"),"Worn"] <- "Yes"
samp[samp$Treatment %in% c("Unworn","Unworn + Washed"),"Worn"] <- "No"

samp$Washed <- NA
samp[samp$Treatment %in% c("Unworn + Washed","Worn + Washed"),"Washed"] <- "Yes"
samp[samp$Treatment %in% c("Unworn","Worn"),"Washed"] <- "No"

## Fridge incubation 
samp$Fridge <- NA
samp[samp$Fabric == samp$First_Tshirt,"Fridge"] <- "Fridge Incubation"
samp[samp$Fabric != samp$First_Tshirt,"Fridge"] <- "Direct Process"
samp[samp$Fabric == "none","Fridge"] <- NA


##### qPCR ######
qpcr1 <- read.table("qpcr_123.csv",sep = ";", header = TRUE, stringsAsFactors = FALSE)
qpcr2 <- read.table("qpcr_456.csv",sep = ";", header = TRUE, stringsAsFactors = FALSE)
qpcr3 <- read.table("qpcr_789.csv",sep = ";", header = TRUE, stringsAsFactors = FALSE)
qpcr4 <- read.table("qpcr_101112.csv",sep = ";", header = TRUE, stringsAsFactors = FALSE)

qpcr <- list(qpcr1, qpcr2, qpcr3, qpcr4)

qpcr.df <- foreach(i = 1:length(qpcr), .combine = rbind) %do% {
  
  # Subset samples and columns
  qpcr.sub <- qpcr[[i]]
  ## TOM/tom: Empty, VAND/vand: Water, JR1/2: Other samples not part of this study
  qpcr.x <- qpcr.sub[!qpcr.sub$Sample.Name %in% c("standard   ","TOM","tom","VAND","vand","JR1","JR2"),c("Sample.Name","Concentration")]
  qpcr.x <- qpcr.x[qpcr.x$Concentration != "-            ",]

  # Take the median
  qpcr.x$Concentration <- as.numeric(qpcr.x$Concentration)
  qpcr.agg <- aggregate(Concentration ~ Sample.Name, data = qpcr.x, FUN = median)
  
  # Mark those below water controls
  qpcr.ctrl <- qpcr.sub[qpcr.sub$Sample.Name %in% c("VAND","vand"),c("Sample.Name","Concentration")]
  max.ctrl <- max(as.numeric(qpcr.ctrl$Concentration))
  qpcr.agg$Low <- "No"
  qpcr.agg[qpcr.agg$Concentration <= max.ctrl,"Low"] <- "Yes"
  
  return(qpcr.agg)
}

# Put on sample names
plateToList <- read.table("PlateToList.csv",sep = ";", header = TRUE, stringsAsFactors = FALSE)
qpcr.final <- merge(qpcr.df, plateToList, by.x = "Sample.Name", by.y = "Pos")[,2:4]

# Merge with sample data
samp.x <- merge(samp, qpcr.final, by = "Sample", all.x = TRUE, sort = FALSE)
samp.o <- samp[order(samp$Sample),]
samp.xo <- samp.x[order(samp.x$Sample),]
all(samp.o$Sample == samp.xo$Sample)
rownames(samp.xo) <- rownames(samp.o)

# Make phyloseq
otu <- otu_table(newtab4, taxa_are_rows = TRUE)
tax <- tax_table(taxa)
phy <- phy_tree(tree)
samdat <- sample_data(samp.xo)

SS <- phyloseq(otu, tax, samdat, phy, ref)

# No controls
SS.nc <- subset_samples(SS, Type == "sample")
SS.nc <- prune_taxa(taxa_sums(SS.nc) > 0, SS.nc)

# 16S rRNA gene Copy number correction
rrn <- read.table("rrnDB-5.4_pantaxa_stats_RDP.txt", header = TRUE, sep = "\t")

tax <- as.data.frame(unclass(tax_table(SS.nc)))
tax$ASV <- rownames(tax)

## Iteratively search after copy number
tax.1 <- merge(tax, rrn[rrn$rank == "genus" & rrn$name %in% tax$Genus,c("name","mean")], 
               by.x = "Genus", by.y = "name", all.x = TRUE, sort = FALSE)

tax.2 <- merge(tax.1, rrn[rrn$rank == "family" & rrn$name %in% tax.1$Family,c("name","mean")], 
               by.x = "Family", by.y = "name", all.x = TRUE, sort = FALSE)

tax.3 <- merge(tax.2, rrn[rrn$rank == "order" & rrn$name %in% tax.2$Order,c("name","mean")], 
               by.x = "Order", by.y = "name", all.x = TRUE, sort = FALSE)

tax.4 <- merge(tax.3, rrn[rrn$rank == "class" & rrn$name %in% tax.3$Class,c("name","mean")], 
               by.x = "Class", by.y = "name", all.x = TRUE, sort = FALSE)

tax.5 <- merge(tax.4, rrn[rrn$rank == "phylum" & rrn$name %in% tax.4$Phylum,c("name","mean")], 
               by.x = "Phylum", by.y = "name", all.x = TRUE, sort = FALSE)

tax.6 <- merge(tax.5, rrn[rrn$rank == "domain" & rrn$name %in% tax.5$Kingdom,c("name","mean")], 
               by.x = "Kingdom", by.y = "name", all.x = TRUE, sort = FALSE)

rownames(tax.6) <- tax.6$ASV

sum(!is.na((tax.1[,ncol(tax.1)])))/nrow(tax)
sum(!is.na((tax.2[,ncol(tax.2)])))/nrow(tax)
sum(!is.na((tax.3[,ncol(tax.3)])))/nrow(tax)
sum(!is.na((tax.4[,ncol(tax.4)])))/nrow(tax)
sum(!is.na((tax.5[,ncol(tax.5)])))/nrow(tax)
sum(!is.na((tax.6[,ncol(tax.6)])))/nrow(tax)

## Pick the copy number
copy <- apply(tax.6[,9:14], 1, function(x) x[!is.na(x)][1])
tax.final <- cbind(tax.6[,1:7],copy)
tax.final <- tax.final[order(rownames(tax.final)), ]

# Absolute abundances
## Order data
otu.new <- otu_table(SS.nc)
otu.new <- otu.new[order(rownames(otu.new)), ]

samp.xo.nc <- samp.xo[samp.xo$Type == "sample", ]
samp.xo.nc <- samp.xo.nc[order(rownames(samp.xo.nc)), ]

all(rownames(otu.new) == rownames(tax.final))
all(colnames(otu.new) == rownames(samp.xo.nc))

## Absolute abundance calculation
otu.new <- apply(otu.new, 2, function(x) x/sum(x)) # Relative abundance
otu.new2 <- t(t(otu.new) * samp.xo.nc$Concentration) # Multiply by qPCR concentration
otu.new3 <- otu.new2 / tax.final$copy # Divide by copy numbers
otu.new4 <- ceiling(otu.new3) # Round up

SSa <- phyloseq(otu_table(otu.new4, taxa_are_rows = TRUE), 
                   tax_table(as.matrix(tax.final)), 
                   sample_data(samp.xo.nc), 
                   phy, ref)

# No Low biomass or read coverage
low <- gsub("-", ".", names(sample_sums(SS)[sample_sums(SS) < 3716]))

SSx.nc <- subset_samples(SS.nc, Low == "No" & !Sample %in% low)
SSax <- subset_samples(SSa, Low == "No" & !Sample %in% low)
SSx <- subset_samples(SS, (Low == "No" | is.na(Low)) & !Sample %in% low)

SSx.nc <- prune_taxa(taxa_sums(SSx.nc) > 0, SSx.nc)
SSax <- prune_taxa(taxa_sums(SSax) > 0, SSax)
SSx <- prune_taxa(taxa_sums(SSx) > 0, SSx)

# Absolute no controls
SSa.nc <- subset_samples(SSa, Guy != "CON")
SSax.nc <- subset_samples(SSax, Guy != "CON")

SSa.nc <- prune_taxa(taxa_sums(SSa.nc) > 0, SSa.nc)
SSax.nc <- prune_taxa(taxa_sums(SSax.nc) > 0, SSax.nc)

# x means no low reads/biomass samples
# a means absolute abundances
# .nc means no controls
# Save
save(SS, SSx, SSa, SSax, SS.nc, SSx.nc, SSa.nc, SSax.nc, file = "physeq.RData")

rm(derepFs, derepRs)
save.image("dada2.RData")
