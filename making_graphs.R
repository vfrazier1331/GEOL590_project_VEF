library(phyloseq)
    #convert files into a phyloseq class object
    phyloseq_import <- import_mothur(mothur_shared_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared", mothur_constaxonomy_file = "lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")
    #makes an otu_table file and a tax_table file
    ###################################################
    #MAKE DF
    ##Convert otu_table into a dataframe + DON'T need to transpose
    OTU1 = as(otu_table(phyloseq_import), "matrix")
    # Coerce to data.frame
    OTUdf = as.data.frame(OTU1)
    #contains the number of each OTU in each sample (A-C)
    
    ##convert tax_table to dataframe
    tax1 = as(tax_table(phyloseq_import), "matrix")
    taxdf = as.data.frame(tax1)
    
    #########################################################
    #MERGE DF: Now we have a tax and otu dataframes... need to merge them into a single dataframe to graph
    merged_df <- cbind(OTUdf,taxdf)
    
    ########################################################
    #FORMAT DF: re-order and re-name columns
    merged_df$otu.number <- seq.int(nrow(merged_df))
    
    # merged_df1 <- merged_df[,c(10,7,8,9,1,2,3,4,5,6)]
    
    colnames(merged_df)[colnames(merged_df)=="Rank1"] <- "Kingdom"
    colnames(merged_df)[colnames(merged_df)=="Rank2"] <- "Phylum"
    colnames(merged_df)[colnames(merged_df)=="Rank3"] <- "Class"
    colnames(merged_df)[colnames(merged_df)=="Rank4"] <- "Order"
    colnames(merged_df)[colnames(merged_df)=="Rank5"] <- "Family"
    colnames(merged_df)[colnames(merged_df)=="Rank6"] <- "Genus"
    
    merged_df
    #We should now have a merged and formatted dataframe (actually: function that is reactive to inputFile) that can be fed into the plot....
    
    merged_df_csv <- write.csv(merged_df, file = "data/merged_df_rmsing.csv")
  })
  
  # filt_otu <- reactive({
  #  low.otu <- input$otu.number[1]
  # hi.otu <- input$otu.number[2]
  
  #  merged_df_final %>%
  #   filter(otu.number >= low.otu) %>%
  #  filter(otu.number <= hi.otu)
  
  
  # Make the plot
  # eventReactive listens for a change in state of input$go, and only runs the code when the state changes
  # Note use of aes_string, since input$xvar returns a string, not a reference to the object ("carat" not carat)
  # Create a dynamic plot plot
  # I moved the ggplot into its own reactive context.
  # Note that creating a ggplot object is very fast - it is actually drawing the plot that is slow.
    
    #  ggplot(filt_otu(), aes_string(x = input$xvar, y = input$yvar, fill = input$color)) + # Note that you need () after filt_dia, since filt_dia() is a function to get the object you want, not the actual object
    #   geom_bar(position="stack", stat="identity") +
    #  theme_minimal() +
    
    # theme(text = element_text(size=14),
    #      axis.text.x = element_text(size=10,angle=45, hjust=1))
  )
  
}

# make_pie <- function(taxonomy_table, taxa, sample, plot_title){
# plot_ly(data = taxonomy_table, labels = ~taxa, values = ~sample, type = 'pie', 
#              textposition = 'inside', 
#              showlegend = FALSE,
#              textinfo = 'label+percent') %>%
#   layout(title = plot_title,
#          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
#          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
# 
# }
# 
# 
# make_pie(taxonomy_table = merged_df, taxa = ~Order, sample = ~A, plot_title = 'Sample A by Phylum')
count.data <- merged_df %>%
  arrange(desc(Phylum)) %>%
  mutate(lab.ypos = cumsum(A) - 0.5*A)
print(count.data)

library(reshape2)
# phylum_samples <- merged_df[,c(1:3,5)]
# 
# phylum_melt <- melt(phylum_samples, id.vars = "Phylum", variable.name = "sample")
# phylum_melt
# 
# 
# bp<- ggplot(phylum_melt, aes(x=sample, y=Phylum))+
#   geom_bar(width = 1, stat = "identity")
# bp

#NOT WORKING


make_ggplot_pie <- function(taxonomydf,sample,taxa){
  count.data <- merged_df %>%
    arrange(desc(Phylum)) %>%
    mutate(lab.ypos = cumsum(A) - 0.5*A)
  

  
  p <- ggplot(count.data, aes(x="", y=A, fill=Phylum)) + 
    geom_bar(stat="identity", width=1) + 
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    geom_text(aes(y = lab.ypos, label = Phylum), color = "white")+
    theme_void()

  print(p)
  }

make_ggplot_pie(taxonomydf = merged_df, sample = A, taxa = Phylum)

print(merged_df)

#convert files into a phyloseq class object

#makes an otu_table file and a tax_table file

###################################################
#MAKE DF
##Convert otu_table into a dataframe + DON'T need to transpose

#We should now have a merged and formatted dataframe (actually: function that is reactive to inputFile) that can be fed into the plot....

