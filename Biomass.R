# Packages
library(ggplot2)
library(COEF) # github.com/Russel88/COEF
library(nlme)
library(phyloseq)

# Load data
load("data/physeq.RData")

# Treatment variable and medians
samp <- as.data.frame(unclass(sample_data(SS)))
samp$Treatmentx <- gsub(" \\+ ","\n\\+\n",samp$Treatment)
sampm <- aggregate(Concentration ~ Treatmentx + Fabric, data = samp[samp$Fabric != "none",], FUN = median)

# Plot
p <- ggplot(samp[samp$Fabric != "none",], aes(Treatmentx, Concentration)) +
  theme_bw() +
  geom_hline(yintercept = max(samp[samp$Low == "Yes","Concentration"], na.rm = TRUE)) +
  geom_jitter(width = 0.2) +
  geom_crossbar(data = sampm, aes(ymin = Concentration, ymax = Concentration), 
                width = 0.5) +
  scale_y_log10(labels = fancy_scientific, breaks = c(1e4,1e5,1e6,1e7)) +
  ylab("16S rRNA gene copies") +
  xlab(NULL) +
  facet_grid(. ~ Fabric)
p

ggsave(filename = "Figures/qPCR.pdf",plot = p,width=12,height=8,units="cm")

# Statistics
samp$Worn <- "No"
samp[samp$Treatment %in% c("Worn","Worn + Washed"),"Worn"] <- "Yes"
samp$Block <- paste(samp$Guy, samp$Orientation)

# Tests for cotton
t.test(log10(Concentration) ~ Worn, data = samp[samp$Fabric == "Cotton",])
t.test(log10(Concentration) ~ Treatment, data = samp[samp$Fabric == "Cotton" & samp$Worn == "No",])

# Mixed-effect model for Worn samples only
fit <- lme(log10(Concentration) ~ Fabric + Treatment, random = ~ 1|Block, data = samp[samp$Worn == "Yes",], method = "ML")
plot(fit)
summary(fit)
intervals(fit)

# Test random effect
fit_noR <- lm(log10(Concentration) ~ Fabric + Treatment, data = samp[samp$Worn == "Yes",])

anova(fit, fit_noR)
drop1(fit)

# Save data
save.image(file="RData/Biomass.RData")
