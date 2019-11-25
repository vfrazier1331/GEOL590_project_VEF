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

#makes an otu_table file and a tax_table file
###################################################

# tax1 = as(tax_table(phyloseq_import), "matrix")
# 
# taxdf = as.data.frame(tax1)
# OTU1 = as(otu_table(phyloseq_import), "matrix")
# # Coerce to data.frame
# 
# OTUdf = as.data.frame(OTU1)
# merged_df <- cbind(OTUdf,taxdf)
# print(merged_df)
#  


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


# #make a dataframe with a column for the read counts of each sample
# sample_sum_df <- data.frame(sum = sample_sums(moth_merge))
# 
# #plots the distribution of sample read count in histogram
# ggplot(sample_sum_df, aes(x = sum)) +
#   geom_histogram(color = "black", fill = "green", binwidth = 2500) +
#   ggtitle("Distribution of sample sequencing depth") +
#   xlab("Read counts") +
#   theme(axis.title.y = element_blank())

#sub-sample to rarefy samples with different sequencing depths
ps.rarefied_phylum = rarefy_even_depth(moth_merge_KPLT, rngseed=1, sample.size=0.9*min(sample_sums(moth_merge_KPLT)), replace=F)
# Try `set.seed(1); .Random.seed` for the full vector
# ...
# 79004OTUs were removed because they are no longer 
# present in any sample after random subsampling
# 
# ...
#make stacked barplot
#prune out low abundance taxa (below 1%)
moth_phylum <- moth_merge_KPLT %>%
  tax_glom(taxrank = "Phylum") %>% #agglomerate at phylum level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>% #transform to rel abund
  psmelt() %>%
  filter(Abundance > 1.0) %>%
  arrange(desc(Phylum))

#########RAREFIED VERSION
moth_phylum_rare <- ps.rarefied_phylum %>%
  tax_glom(taxrank = "Phylum") %>% #agglomerate at phylum level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>% #transform to rel abund
  psmelt() %>%
  filter(Abundance > 1.0) %>%
  arrange(desc(Phylum))
  
moth_genus <- moth_merge_KPLT %>%
  tax_glom(taxrank = "Genus") %>% #agglomerate at phylum level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>% #transform to rel abund
  psmelt() %>%
  filter(Abundance > 1.0) %>%
  arrange(desc(Genus))

# Sum remaining taxa with a relative abundance < 1% and make a new dataframe
Remainders_phy <- (moth_phylum) %>%
  group_by(sample.id,new.name) %>% 
  summarise(Abundance = (100-sum(Abundance))) %>% 
  as.data.frame()
Remainders_phy$Order<-"_Phylums < 1%"
Remainders_phy$Phylum<-"Phylum < 1%"
Remainders_phy$Class<-"_Phylums < 1%"
Remainders_phy$Kingdom<-"Bacteria"

#####REMAINDERS RAREFIED
Remainders_phy_rare <- (moth_phylum_rare) %>%
  group_by(sample.id,new.name) %>% 
  summarise(Abundance = (100-sum(Abundance))) %>% 
  as.data.frame()
Remainders_phy_rare$Order<-"_Phylums < 1%"
Remainders_phy_rare$Phylum<-"Phylum < 1%"
Remainders_phy_rare$Class<-"_Phylums < 1%"
Remainders_phy_rare$Kingdom<-"Bacteria"
# Join dataframes
moth_phy_rarebarchart <- full_join(moth_phylum_rare,Remainders_phy_rare)
moth_phy_barchart <- full_join(moth_phylum,Remainders_phy)

# reorder based on phylogeny
moth_phy_barchart <- moth_phy_barchart[with(moth_phy_barchart, order(Phylum,Class,Order)),]
moth_phy_rarebarchart <- moth_phy_rarebarchart[with(moth_phy_rarebarchart, order(Phylum,Class,Order)),]

# make plot with labels for each Phylum (easier to read)
ggplot(moth_phy_barchart, aes(x = new.name, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option="D") +
  theme_minimal() +
  geom_text(aes(label = Phylum), color = "white", size = 3, hjust = 0.5, vjust = 3, position =     "stack") +
  ggtitle(label = "Sub-aerial Biofilms: Relative Abundance of Taxa by Phylum") +
  theme(axis.title.x = element_blank(), legend.position = "none")

ggsave(filename = "plots/phylum_multiple_env.png", height = 8, width = 9, units = "in", dpi = 500)

#######BARCHART WITH RAREFIED DATA
ggplot(moth_phy_rarebarchart, aes(x = new.name, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option="D") +
  theme_minimal() +
  geom_text(aes(label = Phylum), color = "white", size = 3, hjust = 0.5, vjust = 3, position =     "stack") +
  ggtitle(label = "Sub-aerial Biofilms: Relative Abundance of Taxa by Phylum") +
  theme(axis.title.x = element_blank(), legend.position = "none")

ggsave(filename = "plots/phylum_multiple_env_rarefied.png", height = 8, width = 9, units = "in", dpi = 500)


# phylum_colors <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
#                    "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
#                    "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
# ggplot(moth_phylum, aes(x = new.name, y = Abundance, fill = Phylum)) +
#   geom_bar(stat = "identity") +
#   scale_fill_viridis(discrete = TRUE, option="B") 

