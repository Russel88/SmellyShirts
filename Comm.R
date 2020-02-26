# Packages
library(phyloseq)
library(ggplot2)
library(DAtest)
library(eulerr)

# Load data
load("data/physeq.RData")

# Genus
SSax.g <- tax_glom(SSax, "Genus")

# Test Acinetobacter
test <- subset_samples(SSax.g, Washed == "No" & Fabric == "Cotton")
test <- subset_taxa(test, Genus == "Acinetobacter")
dftest <- psmelt(test)
t.test(log10(Abundance) ~ Treatment, data = dftest)

# Group less abundant
SSax.gp <- preDA(SSax.g, min.abundance = 0.03)

df <- psmelt(SSax.gp)
df$Genus <- as.character(df$Genus)
df[is.na(df$Genus),"Genus"] <- "Others"

df$Washed <- factor(df$Washed, labels = c("Unwashed","Washed"))

df$Genus <- factor(df$Genus, levels = c("Corynebacterium","Staphylococcus","Acinetobacter","Propionibacterium","Moraxella","Others"))

df$Guy <- as.character(df$Guy)
df[df$Guy == "CON", "Guy"] <- "Unworn"

with(df[df$Guy == "SI",], table(Fabric, Washed))

p <- ggplot(df, aes(Guy, Abundance, fill = Genus)) +
  theme_bw() +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(Fabric ~ Washed) +
  scale_x_discrete(labels = c(LETTERS[1:10],"Unworn")) +
  scale_fill_manual(values = c("red3","blue3","gold2","green4","black","darkgrey")) +
  geom_vline(xintercept = 10.5, size = 1) +
  xlab("Subject") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(p, file = "Figures/Rel.pdf", width = 16, height = 12, units = "cm")

# Venn diagram
# Worn
TS.uw <- taxa_sums(subset_samples(SSax, Worn == "No"))
TS.uw <- names(TS.uw[TS.uw > 0])

TS.w <- taxa_sums(subset_samples(SSax, Worn == "Yes"))
TS.w <- names(TS.w[TS.w > 0])

# Wash
TS.uwa <- taxa_sums(subset_samples(SSax, Washed == "No"))
TS.uwa <- names(TS.uwa[TS.uwa > 0])

TS.wa <- taxa_sums(subset_samples(SSax, Washed == "Yes"))
TS.wa <- names(TS.wa[TS.wa > 0])

# Fabric
TS.c <- taxa_sums(subset_samples(SSax, Fabric == "Cotton"))
TS.c <- names(TS.c[TS.c > 0])

TS.p <- taxa_sums(subset_samples(SSax, Fabric == "Polyester"))
TS.p <- names(TS.p[TS.p > 0])

fit <- euler(list(Cotton_Unworn = intersect(TS.c, TS.uw),
                  Polyester_Worn = intersect(TS.p, TS.w),
                  Cotton_Worn = intersect(TS.c, TS.w)),
             shape = "ellipse")
plot(fit, quantities = TRUE)

fitw <- euler(list(Cotton_Unwashed = intersect(TS.w, intersect(TS.c, TS.uwa)),
                  Polyester_Unwashed = intersect(TS.w, intersect(TS.p, TS.uwa)),
                  Cotton_Washed = intersect(TS.w, intersect(TS.c, TS.wa)),
                  Polyester_Washed = intersect(TS.w, intersect(TS.p, TS.wa))),
             shape = "ellipse")
plot(fitw, quantities = TRUE)

# Abundance dependent on core versus fabric-specific 
df <- psmelt(SSax)
df$Washed <- factor(df$Washed, labels = c("Unwashed","Washed"))

df$Guy <- as.character(df$Guy)
df[df$Guy == "CON", "Guy"] <- "Unworn"

df$Prevalence <- "Shared"
df[df$OTU %in% setdiff(TS.uw, TS.w), "Prevalence"] <- "Unworn cotton only"
df[df$OTU %in% setdiff(setdiff(intersect(TS.c, TS.w), TS.uw), TS.p), "Prevalence"] <- "Worn cotton only"
df[df$OTU %in% setdiff(setdiff(intersect(TS.p, TS.w), TS.uw), TS.c), "Prevalence"] <- "Worn polyester only"

p <- ggplot(df, aes(Guy, Abundance, fill = Prevalence)) +
    theme_bw() +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(Fabric ~ Washed) +
    scale_x_discrete(labels = c(LETTERS[1:10],"Unworn")) +
    geom_vline(xintercept = 10.5, size = 1) +
    xlab("Subject") +
    ylab("Relative Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(p, file = "Figures/Rel_Shared.pdf", width = 16, height = 10, units = "cm")
ggsave(p, file = "Figures/Rel_Shared.png", width = 16, height = 10, units = "cm")

# Save
save.image("RData/Rel.RData")