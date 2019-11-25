library(shiny)
library(tidyverse)
library(phyloseq)
library(plotly)


devtools::create("taxplots")
##Already have a function for importing mothur files (phyloseq)

mothur_files_import <- function(shared_file, taxonomy_file) {
  #load packages... 
  library(tidyverse)
  library(phyloseq)
  
  phyloseq_obj <- import_mothur(mothur_shared_file = shared_file, mothur_constaxonomy_file = taxonomy_file)
  phyloseq_obj
}

phyloseq_import <- mothur_files_import(shared_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared", taxonomy_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")

#             mothur_constaxonomy_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")
#makes an otu_table file and a tax_table file
###################################################mothur_shared_file = input$shared_file, mothur_constaxonomy_file = input$tax_file
#MAKE DF

###Create function for concatenating and organizing dataframe

#input into this function needs to be a phyloseq class object....
make_otudf <- function(phyloseq_obj) {
  library(tidyverse)
  library(phyloseq)
  OTU1 = as(otu_table(phyloseq_obj), "matrix")
  # Coerce to data.frame
  OTUdf = as.data.frame(OTU1)
  OTUdf
}
#make a dataframe containing ONLY the OTU abundance data from the dataset using the make_otudf
lavatubes_otu_table <- make_otudf(phyloseq_import)

#make a dataframe concatenating abundance and taxonomic identification for one dataset containing one or more samples
make_taxabundance_df <- function(phyloseq_obj) {
library(tidyverse)
library(phyloseq)
#x = phyloseq_obj made with mothur_files_import (phyloseq object class)
#y = csv file with sample names, number of samples, metadata

##Convert otu_table into a dataframe + DON'T need to transpose
OTU1 = as(otu_table(phyloseq_obj), "matrix")
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
#contains the number of each OTU in each sample 
##convert tax_table to dataframe
tax1 = as(tax_table(phyloseq_obj), "matrix")
taxdf = as.data.frame(tax1)

#########################################################
#MERGE DF: Now we have a tax and otu dataframes... need to merge them into a single dataframe
merged_df <- cbind(OTUdf,taxdf)

colnames(merged_df)[colnames(merged_df)=="Rank1"] <- "Kingdom"
colnames(merged_df)[colnames(merged_df)=="Rank2"] <- "Phylum"
colnames(merged_df)[colnames(merged_df)=="Rank3"] <- "Class"
colnames(merged_df)[colnames(merged_df)=="Rank4"] <- "Order"
colnames(merged_df)[colnames(merged_df)=="Rank5"] <- "Family"
colnames(merged_df)[colnames(merged_df)=="Rank6"] <- "Genus"


merged_df
}

make_taxabundance_df(phyloseq_import)

lavatubes_ABC <- mothur_files_import(shared_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared", taxonomy_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")
KP_yc <- mothur_files_import(shared_file = "data/ETN_KP_yc.trim.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared", taxonomy_file = "data/ETN_KP_yc.trim.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy")
KP_wc <- mothur_files_import(shared_file = "data/ETN_KP_wc.trim.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared", taxonomy_file = "data/ETN_KP_wc.trim.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy")

lavatubes_ABC_df <- make_tax_df(lavatubes_ABC)
KP_wc_df <- make_tax_df(KP_wc)
KP_yc_df <- make_tax_df(KP_yc)

#make OTU tables for each phyloseq
lavatubes_otu_table <- make_otudf(lavatubes_ABC)
KP_wc_otu_table <- make_otudf(KP_wc)
KP_yc_otu_table <- make_otudf(KP_yc)

#Make tax tables using both otu and tax tables from the phyloseq
#to be fed into merging samples

make_tax_df <- function(phyloseq_obj){
  library(tidyverse)
  library(phyloseq)
  #x = phyloseq_obj made with mothur_files_import (phyloseq object class)
  #y = csv file with sample names, number of samples, metadata
  
  ##Convert otu_table into a dataframe + DON'T need to transpose
  OTU1 = as(otu_table(phyloseq_obj), "matrix")
  # Coerce to data.frame
  OTUdf = as.data.frame(OTU1)
  #contains the number of each OTU in each sample 
  ##convert tax_table to dataframe
  tax1 = as(tax_table(phyloseq_obj), "matrix")
  taxdf = as.data.frame(tax1)
  
  #########################################################
  #MERGE DF: Now we have a tax and otu dataframes... need to merge them into a single dataframe
  merged_df <- cbind(OTUdf,taxdf)
  
  colnames(merged_df)[colnames(merged_df)=="Rank1"] <- "Kingdom"
  colnames(merged_df)[colnames(merged_df)=="Rank2"] <- "Phylum"
  colnames(merged_df)[colnames(merged_df)=="Rank3"] <- "Class"
  colnames(merged_df)[colnames(merged_df)=="Rank4"] <- "Order"
  colnames(merged_df)[colnames(merged_df)=="Rank5"] <- "Family"
  colnames(merged_df)[colnames(merged_df)=="Rank6"] <- "Genus"
  merged_df
  
  ###########################################
  # tax1 = as(tax_table(phyloseq_import), "matrix")
#   taxdf = as.data.frame(tax1)
#   
#   
#   merged_df <- cbind(OTUdf,taxdf)
#   
#   colnames(merged_df)[colnames(merged_df)=="Rank1"] <- "Kingdom"
#   colnames(merged_df)[colnames(merged_df)=="Rank2"] <- "Phylum"
#   colnames(merged_df)[colnames(merged_df)=="Rank3"] <- "Class"
#   colnames(merged_df)[colnames(merged_df)=="Rank4"] <- "Order"
#   colnames(merged_df)[colnames(merged_df)=="Rank5"] <- "Family"
#   colnames(merged_df)[colnames(merged_df)=="Rank6"] <- "Genus"
#   
#   
#   merged_df
#   
#   
#   OTU1 = as(otu_table(phyloseq_import), "matrix")
#   # Coerce to data.frame
#   OTUdf = as.data.frame(OTU1)
#   OTUdf
  ##number of samples 
  count_samples <- ncol(OTUdf)

  #make three new columns calculating percentages for each sampled 1:count_samples
  merged_df_percent <- merged_df %>% 
    mutate_at(.funs = funs(percent_ab = (./sum(merged_df[1:count_samples]))*100), .vars = c(1:count_samples))
  merged_df_percent
}




##bin samples from a single mothur run into separate files. Will group by Genus
bin_by_genus <- function(taxabundance_df,sample_key) {
  key <- read.csv(sample_key)
  sample_count <- nrow(key) #retrieves number of samples
  sample_list <- key[,1] # retrieves names of samples
  
  listofdfs <- list()
  if(sample_count == 1) {
       sample_n <- sample_list[1]
  # sample2 <- sample_list[2]
  # sample3 <- sample_list[3]
  
    df <- taxabundance_df %>% select(sample_n, "Genus")
 }
     else {
  for (i in 1:sample_count) {
    sample_n <- sample_list[i]
    # sample2 <- sample_list[2]
    # sample3 <- sample_list[3]
    
    df <- taxabundance_df %>% select(sample_n, "Genus")
    # rating[i] <- mean(dat$Rating[dat$Season == taxabundance_df[i]])
    
  # ab_df2 <- lavatubes_ABC_df %>% select(sample1, "Genus")
  # ab_df3 <- lavatubes_ABC_df %>% select(sample1, "Genus")
     
    listofdfs[[i]] <-df
    
   } 
  return(listofdfs)
  
  }
}

############use function to create a list of the now-separated samples
lava_tubes_binned <- bin_by_genus(lavatubes_ABC_df,"data/lava_tubes_key.csv")
print(lava_tubes_binned)

#unlist binned samples
lavatubes_A <- lava_tubes_binned[[1]]
lavatubes_B <- lava_tubes_binned[[2]]
lavatubes_C <- lava_tubes_binned[[3]]

#testing if statement by inputing file with one sample. function is still useful to select Genus
KP_yc_binned <- bin_by_genus(KP_yc_df,"data/ETN_KP_yc_key.csv")
print(KP_yc_binned)
#When comparing multiple datasets with different OTU number assignments, we need to base abundance on identity only, 
#disregarding OTU#

###this dataframe will contain a certain number of columns depending on the number of samples. Ideally,
###the number and name of samples needs to be provided to customize the next steps to each dataset.

##we could bin each sample into its own dataframe to merge later? That would make samples from other datasets easier to compare?

########################################################
#FORMAT DF: re-order and re-name columns



merged_df$otu.number <- seq.int(nrow(merged_df))

merged_df1 <- merged_df[,c(10,7,8,9,1,2,3,4,5,6)]


merged_df_avg <- cbind(merged_df1,IVMean=rowMeans(merged_df1[2:4], na.rm=TRUE))

merged_df_rmsing <- merged_df_avg %>% filter(IVMean >= 2)

merged_df_rmsing
}

##########making the charts
bar <- ggplot() + 
  geom_bar(aes(x=A, y = A_percent_ab, fill = Phylum), 
           data=lavatubes_ABC_df, 
           stat = "identity")
print(bar)


pie <- plot_ly(merged_df_final(),labels = ~x, values = ~y, type = 'pie',
                  textposition = 'inside',
                  showlegend = FALSE,
                  textinfo = 'label+percent') %>%
            layout(title = 'Community Composition',
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
