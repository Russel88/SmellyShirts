# Packages
library(ggplot2)
library(vegan)
library(nlme)
library(multcomp)
library(foreach)
library(phyloseq)

# Load data
load("data/physeq.RData")

# Rarefaction curve
tab <- otu_table(SSx)
subs <- c(1e2*1:300)
rars <- c()
for(i in subs){
    rars <- c(rars, rarefy(t(tab), sample = i))
}
df <- as.data.frame(cbind(rars,unlist(lapply(subs, function(x) rep(x, ncol(tab))))))
df$Sample <- gsub("\\.", "-", rep(rownames(df)[1:88], 300))

# Subset
colS <- colSums(tab)
df2 <- foreach(i = names(colS), .combine = rbind) %do% {
    df[df$Sample == i & df$V2 < colS[names(colS) == i],]
}

# Add sample data
samp <- sample_data(SSx)
samp$Sample <- rownames(samp)
df2 <- merge(df2, unclass(samp), by = "Sample")
df2$Fabric <- as.character(df2$Fabric)
df2[df2$Fabric == "none", "Fabric"] <- "Control"

# 20 first library size from low to high
lowss <- sample_sums(SSx)[order(sample_sums(SSx))][1:20]
lowss

## Plot it
p.r <- ggplot(df2, aes(x = V2,y = rars, group = Sample)) + 
  geom_line() + 
  theme_bw() +
  facet_grid(. ~ Fabric, scales = 'free_y') +
  geom_vline(xintercept = max(lowss[lowss <= 4000])) +
  ylab("Observed richness") +
  xlab("Read depth") +
  theme(legend.title = element_blank(),
        legend.background = element_blank())
p.r

ggsave("Figures/Rarefaction.pdf", p.r, units = "cm", width = 20, height = 8)
ggsave("Figures/Rarefaction.png", p.r, units = "cm", width = 20, height = 8)

# Richness
rich <- estimate_richness(SSx, measures = c("Observed"))
rownames(rich) <- gsub("\\.","-",rownames(rich))
rich <- merge(rich, as.data.frame(sample_data(SSx)), by = "row.names")
rich <- rich[rich$Fabric != "none",]
rich$Treatment <- gsub(" ","\\\n",rich$Treatment)
richo <- aggregate(Observed ~ Fabric + Treatment, data = rich, FUN = median)

p1 <- ggplot(rich, aes(Treatment, Observed)) + 
  theme_bw() +
  facet_grid(. ~ Fabric, scales = 'free_x', space = "free_x") +
  geom_jitter(width = 0.2) +
  geom_crossbar(data = richo, aes(ymin = Observed, ymax = Observed), width = 0.5) +
  ylab("Observed richness") +
  xlab(NULL)
p1

ggsave("Figures/Richness.pdf",p1,units="cm",width=12,height=8)

# Stats
rich$Block <- paste(rich$Guy, rich$Orientation)
rich$Comb <- interaction(rich$Washed,rich$Fabric)

fit.worn <- lme(log2(Observed) ~ 0 + Comb, random = ~ 1|Block, 
                 data = rich[rich$Fabric != "Control" & 
                               rich$Treatment %in% c("Worn","Worn\n+\nWashed"),])
summary(fit.worn)
plot(fit.worn)
qqnorm(resid(fit.worn))
qqline(resid(fit.worn))

M <- rbind(c(1,1,-1,-1), # Cotton vs Polyester
           c(-1,1,0,0), # Wash effect in cotton
           c(0,0,-1,1)) # Wash effect in polyester
MultComp <- glht(fit.worn,linfct=mcp(Comb=M))
confint(MultComp)

# Cotton only
t.test(log2(Observed) ~ Worn, data = rich[rich$Fabric == "Cotton",])

# Unworn only
t.test(log2(Observed) ~ Washed, data = rich[rich$Worn == "No",])

# Save data
save.image(file="RData/Alpha.RData")
