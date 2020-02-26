# Packages
library(ggplot2)
library(DAtest)
library(phyloseq)
library(foreach)

# Load data
load("data/physeq.RData")

## Treatment
# Make pairing variable
sample_data(SSax.nc)$Pair <- paste(sample_data(SSax.nc)$Guy,sample_data(SSax.nc)$Fabric,sample_data(SSax.nc)$Orientation, sep = "_")

# Group low occuring 
SSax.nc.x <- preDA(SSax.nc, min.samples = 10)

# Test
set.seed(42)
test <- testDA(SSax.nc.x, predictor = "Washed", paired = "Pair", relative = FALSE, effectSize = 10)
summary(test)

# Run test
res.lli <- DA.lli(SSax.nc.x, predictor = "Washed", paired = "Pair", relative = FALSE)

# Plot
# Subset to only balanced
SSax.nc.b <- subset_samples(SSax.nc, Pair %in% names(table(Pair)[table(Pair) == 2]))
SSax.nc.bx <- preDA(SSax.nc.b, min.samples = 10)
phydf <- psmelt(SSax.nc.bx)
phydfX <- phydf[order(phydf$Washed, phydf$Pair), ]

feature_sub <- res.lli[res.lli$pval.adj < 0.01 & abs(res.lli$logFC) > 1.5, "Feature"]

changedf <- foreach(i = feature_sub, .combine = rbind) %do% {
    data.frame(ASV = i,
               FC = log2((phydfX[phydfX$OTU == i, "Abundance"][36:70]+1)/(phydfX[phydfX$OTU == i, "Abundance"][1:35]+1)))
}

namedf <- foreach(i = feature_sub, .combine = rbind) %do% {
    data.frame(ASV = i,
               Name = unique(phydfX[phydfX$OTU == i, "Species"]))
}

namedf$Name <- gsub("\\.", " ", namedf$Name)
namedf[namedf$ASV == "ASV_12", "Name"] <- "Streptococcus sp."
namedf[namedf$ASV == "ASV_30", "Name"] <- "Anaerococcus sp."
namedf[namedf$ASV == "ASV_103", "Name"] <- "Staphylococcus sp."

changedf$ASV <- factor(changedf$ASV, levels = c("ASV_103", "ASV_73", "ASV_30", "ASV_19", "ASV_12", "ASV_8"))

p <- ggplot(changedf, aes(ASV, FC)) +
    theme_bw() +
    geom_hline(yintercept = 0, size = 1) +
    geom_boxplot() +
    coord_flip() +
    geom_text(aes(x = ASV, y = 0.5, label = Name), data = namedf, hjust = 0, vjust = -0.8, size = 3.5) +
    xlab(NULL) +
    ylab("Log2 fold change (Washed/Unwashed)")
p
ggsave(filename = "Figures/DA_Wash.pdf",plot = p,width=14,height=10,units="cm")

# Save data
save.image(file="RData/DA.wash.RData")
