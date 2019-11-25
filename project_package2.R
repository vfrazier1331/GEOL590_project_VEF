library(shiny)
library(tidyverse)
library(phyloseq)
library(plotly)
install.packages('matrixStats')
library(matrixStats)
phyloseq_import <- import_mothur(mothur_shared_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared", mothur_constaxonomy_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")
#makes an otu_table file and a tax_table file
###################################################
##MAKE data frame
##Convert otu_table into a dataframe + DON'T need to transpose

OTU1 = as(otu_table(phyloseq_import), "matrix")
# Coerce to data.frame

OTUdf = as.data.frame(OTU1)
#contains the number of each OTU in each sample (A-C)
class(OTUdf)
##convert tax_table to dataframe

tax1 = as(tax_table(phyloseq_import), "matrix")

taxdf = as.data.frame(tax1)



#########################################################
#MERGE DF: Now we have a tax and otu dataframes... need to merge them into a single dataframe to graph

merged_df <- cbind(OTUdf,taxdf)
print(merged_df)

########################################################
#FORMAT DF: re-order and re-name columns

merged_df$otu.number <- seq.int(nrow(merged_df))

colnames(merged_df)[colnames(merged_df)=="Rank1"] <- "Kingdom"
colnames(merged_df)[colnames(merged_df)=="Rank2"] <- "Phylum"
colnames(merged_df)[colnames(merged_df)=="Rank3"] <- "Class"
colnames(merged_df)[colnames(merged_df)=="Rank4"] <- "Order"
colnames(merged_df)[colnames(merged_df)=="Rank5"] <- "Family"
colnames(merged_df)[colnames(merged_df)=="Rank6"] <- "Genus"

merged_df_avg <- cbind(merged_df,IVMean=rowMeans(merged_df[1:3], na.rm=TRUE))

merged_df_rmsing <- merged_df_avg %>% filter(IVMean >= 2)

count_samples <- ncol(OTUdf)
count_samples

#make three new columns calculating relative abundance/percentages for each sampled 1:count_samples
merged_df_percent <- merged_df_rmsing %>% 
  mutate_at(.funs = funs(percent_ab = (./sum(merged_df_rmsing[1:count_samples]))*100), .vars = c(1:count_samples))

# count_dataframe <- ncol(merged_df_percent)
# last_col = count_dataframe-count_samples+1
# print(last_col)


#finding mean relative abundances
percent_avg <- cbind(merged_df_percent,mean_ab=rowMeans(merged_df_percent[last_col:count_dataframe], na.rm=TRUE))
print(percent_avg)

grouped_phylum <- percent_avg %>%
  group_by(percent_avg$Phylum) %>%
  summarize(n())
grouped_phylum



#grouping those with abundances less than 0.01%
###compile taxa by order (filtering out low abundance taxa)
# transformed <- t(percent_avg)
# print(transformed)
# count_rows <- nrow(transformed)
# count_rows
# sample_rownum <- count_rows-count_samples
# as.data.frame(sample_rownum)
# #####Now that the df is all fucked up, need to make new df using subset, name columns using data from rows...
# phylum = transformed['Phylum',]
# order = transformed['Order',]
# rel_ab <- transformed[sample_rownum:count_rows,]
# 
# df <- rbind(phylum, rel_ab)
# df_order <- rbind(order, rel_ab)
# colnames(df) <- as.character(unlist(df[1,]))
# colnames(df_order) <- as.character(unlist(df_order[1,]))
# df_phylum <- df[-1,]
# df_ord <- df_order[-1,]


p <- ggplot(data=df_phylum, aes(x=Sample, y=Abundance, fill=phylum))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))



compile_order <- percent_avg %>%
  tax_glom(taxrank = "Order") %>%            #agglomerate taxa at order level
  psmelt() %>%                              #melt phyloseq object to long format 
  filter(mean_ab > 0.05) %>%                 #filter our orders below 0.05%
  arrange(desc(Phylum), desc(Class), desc(Order))


write.csv("data/")

proteo_order <- plot_ly(merged_df_rmsing,
                        labels = ~Order,
                        parents = ~Phylum,
                        values = ~A,
                        type = 'sunburst',
                        branchvalues = 'total'
)
print(proteo_order)


# Make the plot

pie <- plot_ly(merged_df_rmsing, labels = ~Phylum, values = ~A, type = 'pie',
               textposition = 'inside', 
               showlegend = FALSE,
               textinfo = 'label+percent') %>%
  layout(title = 'Community Composition',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

print(pie)

bar <-  ggplot(merged_df_rmsing, aes_string(x = input$xvar, y = input$yvar, fill = input$color)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45, hjust=1))

print(bar)