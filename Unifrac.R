# Packages
library(ggplot2)
library(vegan)
library(phyloseq)
library(ape)

# Load data
load("data/physeq.RData")

# Root
tree_old <- phy_tree(SSax.nc)
tree_new <- root(tree_old, outgroup = "ASV_3720", resolve.root = TRUE)
phy_tree(SSax.nc) <- phy_tree(tree_new)

# UniFrac and average
UF <- UniFrac(SSax.nc, weighted = TRUE)

# Extract sample_data
Samp <- sample_data(SSax.nc)
Samp$Block <- paste(Samp$Guy, Samp$Orientation)

# Subject specific effects
set.seed(1)
adonis(UF ~ Guy + Fabric + Washed + Fridge + Guy:Fabric + Guy:Washed, 
        data = data.frame(Samp), permutations = 999)

# Plot PCoA
PCOA <- cmdscale(UF, eig = TRUE)

# Collect in data.frame
df <- data.frame(MDS1 = c(PCOA$points[,1]),
                 MDS2 = c(PCOA$points[,2]),
                 Sample = Samp$Sample)

# Merge with all sample data
df <- merge(df, unclass(Samp), by = "Sample")

# Plot
p <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Washed)) + 
  geom_point() + 
  theme_bw() + 
  facet_grid(.~Fabric) +
  stat_ellipse(level = 0.9, geom = "polygon",alpha=0.1) +
  xlab(paste0("PCoA 1 (",round((PCOA$eig/sum(PCOA$eig))[1]*100,1),"%)")) +
  ylab(paste0("PCoA 2 (",round((PCOA$eig/sum(PCOA$eig))[2]*100,1),"%)"))
p
ggsave(filename = "Figures/PCOA.UFa.Fabric.Wash.pdf",plot = p,width=16,height=8,units="cm")

px <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Guy)) + 
  geom_point() + 
  theme_bw() + 
  stat_ellipse(level = 0.9, geom = "polygon",alpha=0.1) +
  xlab(paste0("PCoA 1 (",round((PCOA$eig/sum(PCOA$eig))[1]*100,1),"%)")) +
  ylab(paste0("PCoA 2 (",round((PCOA$eig/sum(PCOA$eig))[2]*100,1),"%)")) +
  scale_color_discrete(labels = LETTERS[1:10], name = "Subject")
px
ggsave(filename = "Figures/PCOA.UFa.Guy.pdf",plot = px,width=14,height=10,units="cm")

# pRDA
test <- dbrda(UF ~ Samp$Fabric + Samp$Washed + Condition(Samp$Guy))
summary(test)
set.seed(1)
anova(test, by = "terms", model = "reduced")
RsquareAdj(test)

rdadf <- test$CCA$wa
rdadf <- merge(rdadf, Samp, by = "row.names")
rdape <- summary(test)$cont$importance[, c(1, 2)]
rat <- summary(test)$cont$importance[1, 2] / summary(test)$cont$importance[1, 1]

pp <- ggplot(data = rdadf, aes(dbRDA1, dbRDA2)) +
    theme_bw() +
    geom_point(aes(colour = Fabric, shape = Washed),
               size = 2) +
    stat_ellipse(aes(colour = Fabric, linetype = Washed),
                  level = 0.9, geom = "polygon", alpha=0.1) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    coord_equal(ratio = rat) +
    xlab(paste0("dbRDA1 (",round(rdape[2, 1]*100,1),"%)")) +
    ylab(paste0("dbRDA2 (",round(rdape[2, 2]*100,1),"%)"))
pp
ggsave(filename = "Figures/dbRDA.pdf",plot = pp,width=14,height=9,units="cm")

# Save data
save.image(file="RData/Unifrac.RData")

