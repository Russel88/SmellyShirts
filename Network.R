# Packages
library(phyloseq)
library(DAtest)

# Load data
load("data/physeq.RData")

# Subset only dirty
SSax.nc.d <- subset_samples(SSax.nc, Washed == "No")

# Remove low abundant
SSax.nc.d <- preDA(SSax.nc.d, min.abundance = 0.001)

# Normalize to Fabric and CLR
subC <- subset_samples(SSax.nc.d, Fabric == "Cotton")
asvC <- as.data.frame(otu_table(subC))
asvC2 <- t(apply(asvC, 1, function(x) x-mean(x)))
otu_table(subC) <- otu_table(asvC2, taxa_are_rows = TRUE)

subP <- subset_samples(SSax.nc.d, Fabric == "Polyester")
asvP <- as.data.frame(otu_table(subP))
asvP2 <- t(apply(asvP, 1, function(x) x-mean(x)))
otu_table(subP) <- otu_table(asvP2, taxa_are_rows = TRUE)

SSax.nc.dN <- merge_phyloseq(subP, subC)
SSax.nc.dN <- prune_taxa(taxa_names(SSax.nc.dN)[taxa_names(SSax.nc.dN) != "Others"],SSax.nc.dN)

asv <- as.data.frame(unclass(otu_table(SSax.nc.dN)))
asv <- as.matrix(asv)

cors <- cor(t(asv), method = "spearman")
diag(cors) <- NA

corsX <- cors
corsX[abs(corsX) < 0.6] <- NA
sum(!is.na(corsX))

hist(corsX)

# Export
rowCol <- expand.grid(rownames(corsX), colnames(corsX))
labs <- rowCol[as.vector(upper.tri(corsX,diag = FALSE)),]
edges <- cbind(labs, corsX[upper.tri(corsX,diag = FALSE)])
colnames(edges) <- c("Source","Target","Weight")
edges$Sign <- "Pos"
edges <- edges[!is.na(edges$Weight),]
edges[edges$Weight < 0,"Sign"] <- "Neg"
edges$Weight <- abs(edges$Weight)
edges[edges$Sign == "Neg", "Weight"] <- (1 - edges[edges$Sign == "Neg", "Weight"])/10

write.table(edges, file = "edge.csv", quote = FALSE, row.names = FALSE, sep = ";")

tax <- as.data.frame(unclass(tax_table(SSax.nc)))
tax <- tax[rownames(tax) %in% unique(c(as.character(edges$Source), as.character(edges$Target))),]
out.n <- data.frame(Id = rownames(tax),
                    Label = tax$Species)
out.n$Label <- as.character(out.n$Label)
out.n[str_count(out.n$Label,"/") > 2,"Label"] <- gsub("\\..*","\\ sp\\.",out.n[str_count(out.n$Label,"/") > 2,"Label"])
out.n[grep("Genus",out.n$Label),"Label"] <- paste(gsub("Genus.","",out.n[grep("Genus",out.n$Label),"Label"]),"sp.")
out.n$Label <- gsub("sp","sp.",gsub("\\."," ",out.n$Label))
out.n$Label <- paste(out.n$Id, out.n$Label)
outn2 <- merge(out.n, tax[, "Genus", drop=FALSE], by = "row.names")[, -1]

write.table(outn2, file = "nodes.csv", quote = FALSE, row.names = FALSE, sep = ";")

# Save
save.image("RData/Network.RData")

# 10 most abundant
test <- subset_samples(SSax.nc, Washed == "No")
testr <- transform_sample_counts(test, function(x) x/sum(x))
test2 <- as.data.frame(unclass(otu_table(testr)))
test3 <- apply(test2, 1, mean)
test4 <- test3[order(test3, decreasing = TRUE)]

sum(test4[1:10])
sum(test4)

tax_table(SSax)[names(test4[1:10])]
