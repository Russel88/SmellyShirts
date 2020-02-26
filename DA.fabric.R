# Packages
library(ggplot2)
library(DAtest)
library(phyloseq)

# Load data
load("data/physeq.RData")

# Genus level
SSax.nc.g <- tax_glom(SSax.nc, "Genus")

# Remove low abundant
SSax.nc.gx <- preDA(SSax.nc.g, min.samples = 10)

# Test
set.seed(42)
test <- testDA(SSax.nc.gx, predictor = "Fabric", relative = FALSE)
summary(test)

# Run test
res.per <- DA.per(SSax.nc.gx, predictor = "Fabric", relative = FALSE)

# Plot
phydf <- psmelt(SSax.nc.gx)

feature_sub <- res.per[res.per$pval.adj < 0.05 & abs(res.per$log2FC) > 3, "Feature"]
phydfX <- phydf[phydf$OTU %in% feature_sub, ]
abund <- rowSums(otu_table(SSax.nc.gx))[order(rowSums(otu_table(SSax.nc.gx)), decreasing = TRUE)]
phydfX <- phydfX[order(match(phydfX$OTU, res.per[order(res.per$log2FC), "Feature"])), ]
phydfX$OTU <- factor(phydfX$OTU, levels = unique(phydfX$OTU), labels = unique(phydfX$Genus))

p <- ggplot(phydfX, aes(OTU, Abundance+1, fill = Fabric)) +
    theme_bw() +
    scale_y_log10() +
    geom_boxplot() +
    coord_flip() +
    xlab(NULL) +
    ylab("Abundance + 1")
p
ggsave(filename = "Figures/DA_Fabric.pdf",plot = p,width=14,height=10,units="cm")

# Save data
save.image(file="RData/DA.fabric.RData")
