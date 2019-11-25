###diversity analysis
library(vegan)
library(phyloseq)

library(shiny)
library(tidyverse)
library(phyloseq)
library(plotly)
#install.packages('matrixStats')
library(matrixStats)
library(viridis)
phyloseq_import <- import_mothur(mothur_shared_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared", mothur_constaxonomy_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")
ETN_KP_yc_import <- import_mothur(mothur_shared_file = "data/ETN_KP_yc.trim.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared", mothur_constaxonomy_file = "data/ETN_KP_yc.trim.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy")
ETN_KP_wc_import <- import_mothur(mothur_shared_file = "data/ETN_KP_wc.trim.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared", mothur_constaxonomy_file = "data/ETN_KP_wc.trim.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy")
cave2_import <- import_mothur(mothur_shared_file = "data/cave2.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared", mothur_constaxonomy_file = "data/cave2.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy")
surface3_import <- import_mothur(mothur_shared_file = "data/surface3.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared", mothur_constaxonomy_file = "data/surface3.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")
map_lava <- read.csv("data/lava_tubes_meta.csv")

rownames(map_lava) <- map_lava$sample.id #assign rownames to be sample.id names
map_phy <- sample_data(map_lava)

map_KP <- read.csv("data/ETN_KP_yc_meta.csv")
rownames(map_KP) <- map_KP$sample.id #assign rownames to be sample.id names
map_phy_KP <- sample_data(map_KP)

map_KPwc <- read.csv("data/ETN_KP_wc_meta.csv")
rownames(map_KPwc) <- map_KPwc$sample.id #assign rownames to be sample.id names
map_phy_KPwc <- sample_data(map_KPwc)

map_cave2 <- read.csv("data/cave2_meta.csv")
rownames(map_cave2) <- map_cave2$sample.id #assign rownames to be sample.id names
map_phy_cave2 <- sample_data(map_cave2)

map_surface3 <- read.csv("data/surface3_meta.csv")
rownames(map_surface3) <- map_surface3$sample.id #assign rownames to be sample.id names
map_phy_surface3 <- sample_data(map_surface3)
#merged metadata into phyloseq object
sample_names(cave2_import) <- "cave2"
sample_names(surface3_import) <- "surface3"



moth_merge <- merge_phyloseq(phyloseq_import,map_phy)

moth_merge_KP <- merge_phyloseq(ETN_KP_yc_import, map_phy_KP)
moth_merge_KPwc <- merge_phyloseq(ETN_KP_wc_import, map_phy_KPwc)
moth_merge_cave2 <- merge_phyloseq(cave2_import, map_phy_cave2)

moth_merge_surface3 <- merge_phyloseq(surface3_import, map_phy_surface3)
sample_names(moth_merge_surface3)
moth_merge_KPLT <- merge_phyloseq(moth_merge_KP, moth_merge, moth_merge_KPwc, moth_merge_cave2, moth_merge_surface3)

head(sample_data(moth_merge_KPLT))

colnames(tax_table(moth_merge_KPLT)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

ps.rarefied = rarefy_even_depth(moth_merge_KPLT, rngseed=1, sample.size=0.9*min(sample_sums(moth_merge_KPLT)), replace=F)

moth_phylum <- ps.rarefied %>%
  tax_glom(taxrank = "Phylum", NArm =FALSE) %>% #agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)}) %>% #transform to rel abund
  psmelt() %>%
filter(Abundance > 0.01) %>%
arrange(Phylum)

#use moth_merge_KPLT from the project_graphingmothur.R
plot_richness(moth_merge_KPLT, x = "new.name", color = NULL, shape = NULL,
              title = NULL, scales = "free_y", nrow = 1, shsi = NULL,
              measures = NULL, sortby = NULL)

ggsave(filename = "plots/alpha_diversityALL.png", height = 4, width = 10, units = "in", dpi = 500)


#############creating rarefaction curve
####make OTU table into dataframe
OTU1 = as(otu_table(moth_merge_KPLT), "matrix")
# # Coerce to data.frame
OTUdf = as.data.frame(OTU1)

data(moth_merge_KPLT, package = "vegan")
mmKPLT <- OTUdf[1:7, ]

col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)
rarefaction <- with(pars[1:7,],
                    rarecurve(t(otu_table(moth_merge_KPLT)), col = col, lty = lty, step=50, cex=0.5))
png("plots/rarefaction_curve.png", width = 480, height = 480)

# rarefy without replacement
ps.rarefied = rarefy_even_depth(moth_merge_KPLT, rngseed=1, sample.size=0.9*min(sample_sums(moth_merge_KPLT)), replace=F)
# `set.seed(1)` was used to initialize repeatable random subsampling.
# Please record this for your records so others can reproduce.
# Try `set.seed(1); .Random.seed` for the full vector
# ...
# 79004OTUs were removed because they are no longer 
# present in any sample after random subsampling
# 
# ...
plot_bar(ps.rarefied, fill="Rank2")
plot <- plot_bar(ps.rarefied, fill="Phylum")
print(plot)



ggplot(moth_phylum, aes(x = new.name, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=TRUE, option="B") 

#########################################################################################
#**nMDS**
otus <- as.data.frame(otu_table(moth_merge_KPLT))
tax <- as.data.frame(tax_table(moth_merge_KPLT))

tax_grouped <- tax %>%
  group_by()


otus = read.csv(file = "EPIPHYTE_OTU_TABLE.csv", header = TRUE, check.names = FALSE, row.names = 1)
metadata = read.csv(file = "EPIPHYTE_METADATA.csv", header = TRUE, check.names = FALSE, row.names = 1)
good_samples <- colnames(otus[(colSums(decostand(otus,"pa")) >= 1)])     # decostand(x,"pa") counts presence/absence
otus = otus[,good_samples] #creates object from good_samples of just the otus?
metadata = metadata[good_samples,]
t_otus <- as.data.frame(t(otus))
min_depth = min(colSums(otus))
t_otus_rarefied <- as.data.frame(round(rrarefy(t_otus, min_depth)))

#transform data by sqrt: de-emphasizes large abundances relative to very rare taxa
sqrt_t_otus_rarefied = sqrt(t_otus_rarefied)

#creates a matrix that compares distances between abundances of each OTU between sites
#many ways to do this: must use a method... 
#options = euclidian distance, manhattan, Bray-Curtis = indices = c(...)
rank.totus <- rankindex(as.matrix(sqrt_t_otus_rarefied), t_otus_rarefied, indices = c("bray", "euclid", "manhattan", "horn"), method = "spearman")

print(paste("The highest rank was given by the", names(sort(rank.totus, decreasing = TRUE)[1]), "method."))
#[1] "The highest rank was given by the bray method."

#this tells us bray gives the best data, use "bray" method to determine distances: otus_dist
#create distance matrix
scale(tt_otus_rarefied,scale = TRUE,center = TRUE)
otus_dist = as.matrix((vegdist(t_otus_rarefied, "bray")))

#perform NMDS
NMDS = metaMDS(otus_dist)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Host = metadata$Host_Taxon, Location = metadata$Location)

head(NMDS) #just axes and spread of data

#adds ellipses 1 stdev from center of data? 
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Location)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "NMDS Plot")
