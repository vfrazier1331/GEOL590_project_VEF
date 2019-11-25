########
# Shiny app to create a dynamically-filterable visualization of the diamonds app
########
#install.packages("plotly")
# These bits get run before any of the rest of the code
# Note: contrary to what I told you on Monday, the use of the global.R file is no longer recommended.
# At present, I'm not sure why.
# install.packages('ggthemes')
# install.packages('extrafont')
library(shiny)
library(tidyverse)
library(phyloseq)
library(plotly)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(viridis)
library(vegan)
options(shiny.maxRequestSize = 30*1024^2)


#phyloseq_import <- import_mothur(mothur_group_file = "data/lavatube_mat1.good.groups", mothur_shared_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared",
#             mothur_constaxonomy_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")
# Need a vector of axis variables as characters
# axis_vars <- names(merged_df_rmsing)

# Create a character vector of those columns of diamonds that are 
# factor.indices <- vapply(merged_df_rmsing, is.factor, TRUE) # vapply is a base R function like sapply & lapply that we haven't talked about
# It applies a function (is.factor) to every element of a list (diamonds, remember that a data frame is a list)
# The 'TRUE' is at the end to convey that we want vapply to return the results as a logical vector:
#   I could just as well have used FALSE or c(TRUE, FALSE, TRUE) or any other logical vector
# factor.columns <- axis_vars[factor.indices]

#min.otu <- min(merged_df_rmsing$otu.number)
#max.otu <- max(merged_df_rmsing$otu.number)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Taxonomy Charts"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("tax_file", "Choose *.taxonomy File",
                      multiple = FALSE,
                      accept = c(".taxonomy")),
            fileInput("shared_file", "Choose *.shared File",
                      multiple = FALSE,
                      accept = c(".shared")),
            fileInput("meta_file", "Choose *.csv File",
                      multiple = FALSE,
                      accept = c(".csv"))
            
            # uiOutput("sample"),
            # uiOutput("taxonomy")
            # uiOutput("color")
            
        ),
        
        # Show a plot of diamonds data frame. This output doesn't care what that plot is, only that it will be associated with output$diamonds_plot
        mainPanel(
            plotOutput("taxonomy_plot"),
            plotOutput("rarefaction_plot")
        )
    )
)


server <- function(input, output) {
    ###making dataframe out of two files: .taxonomy and .shared files
    abundance_data <- reactive({ 
        
        sharedfile <- input$shared_file
        taxfile <- input$tax_file
        meta.file <- input$meta_file
        
        metafile <- read.csv(meta.file$datapath)
        #assign rownames to be sample.id names
        rownames(metafile) <- metafile$sample.id
        
        phyloseq_import <- import_mothur(mothur_shared_file = sharedfile$datapath, mothur_constaxonomy_file = taxfile$datapath)
        #             mothur_constaxonomy_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")
        map_phy <- sample_data(metafile)
         
        
        #merged metadata into phyloseq object
        moth_merge <- merge_phyloseq(phyloseq_import,map_phy)
        
        colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
        
        moth_merge
        #rarefaction not needed if only one dataset is used....
        #mm.rarefied_phylum = rarefy_even_depth(moth_merge, rngseed=1, sample.size=0.9*min(sample_sums(moth_merge_KPLT)), replace=F)
        
        
        # ###group and rename phyla less than 1% of entire community
        # Remainders_phy <- (moth_phylum) %>%
        #     group_by(sample.id,new.name) %>% 
        #     summarise(Abundance = (100-sum(Abundance))) %>% 
        #     as.data.frame()
        # Remainders_phy$Order<-"_Phylums < 1%"
        # Remainders_phy$Phylum<-"Phylum < 1%"
        # Remainders_phy$Class<-"_Phylums < 1%"
        # Remainders_phy$Kingdom<-"Bacteria"
        ##Convert otu_table into a dataframe + DON'T need to transpose
        # OTU1 = as(otu_table(phyloseq_import), "matrix")
        # # Coerce to data.frame
        # OTUdf = as.data.frame(OTU1)
        # #contains the number of each OTU in each sample (A-C)
        # # class(OTUdf)
        # ##convert tax_table to dataframe
        # tax1 = as(tax_table(phyloseq_import), "matrix")
        # taxdf = as.data.frame(tax1)
        # 
        # #########################################################
        # #MERGE DF: Now we have a tax and otu dataframes... need to merge them into a single dataframe to graph
        # merged_df <- cbind(OTUdf, taxdf)
        # 
        # ########################################################
        # #FORMAT DF: re-order and re-name columns
        # merged_df$otu.number <- seq.int(nrow(merged_df))
        # 
        # 
        # colnames(merged_df)[colnames(merged_df)=="Rank1"] <- "Kingdom"
        # colnames(merged_df)[colnames(merged_df)=="Rank2"] <- "Phylum"
        # colnames(merged_df)[colnames(merged_df)=="Rank3"] <- "Class"
        # colnames(merged_df)[colnames(merged_df)=="Rank4"] <- "Order"
        # colnames(merged_df)[colnames(merged_df)=="Rank5"] <- "Family"
        # colnames(merged_df)[colnames(merged_df)=="Rank6"] <- "Genus"
        # 
        # merged_df_avg <- cbind(merged_df,IVMean=rowMeans(merged_df[1:3], na.rm=TRUE))
        # 
        # merged_df_rmsing <- merged_df_avg %>% filter(IVMean >= 2)
        # 
        # count_samples <- ncol(OTUdf)
        # 
        # #make three new columns calculating relative abundance/percentages for each sampled 1:count_samples
        # merged_df_percent <- merged_df_rmsing %>% 
        #     mutate_at(.funs = funs(percent_ab = (./sum(merged_df_rmsing[1:count_samples]))*100), .vars = c(1:count_samples))
        # 
        # # remainder <- merged_df_percent[merged_df_percent$]
        #  merged_df_percent
        #in order to use for stacked bar plot, we need to calculate percentage of each community member in each sample
    #     count_samples <- ncol(OTUdf)
    #     names_samples <- colnames(OTUdf)
    #     if (count_samples == 1) {
    #         merged_df_percent <- merged_df_rmsing %>% mutate(merged_df_rmsing, percent_ab = merged_df_rmsing[1]/sum(merged_df_rmsing[1]))
    #     } else {
    #     for(i in 1:count_samples){
    #     
    #         sample_n = names_samples[i]
    #         
    #         merged_df_percent <- merged_df_rmsing %>% mutate(merged_df_rmsing, percent_ab = sample_n/sum(sample_n))
    #     }
    # }
    #     return(merged_df_percent)
        
    })
    
    # output$sample = renderUI({
    #     df <- abundance_data()
    #     axis_vars = names(df)
    #     selectInput(inputId = "sample",
    #                 label = "Sample",
    #                 choices = axis_vars,
    #                 selected = "A")
    
    
    # moth_phylum <- abundance_data() %>%
    #     tax_glom(taxrank = input$taxonomy) %>% #agglomerate at phylum level
    #     transform_sample_counts(function(x) {x/sum(x)}) %>% #transform to rel abund
    #     psmelt() %>%
    #     filter(Abundance > 0.02) %>%
    #     arrange(input$taxonomy)
    # 
    # output$taxonomy = renderUI({
    #     df <- abundance_data()
    #     axis_vars = rank_names(df)
    #     selectInput(inputId = "taxonomy",
    #                 label = "Taxonomy",
    #                 choices = axis_vars,
    #                 selected = "Phylum")
    # })
    
    # output$color = renderUI({
    #     df <- abundance_data()
    #     axis_vars = sample_names(df)
    #     selectInput(inputId = "color",
    #                 label = "Color",
    #                 choices = axis_vars,
    #                 selected = "Order")
    # })
    
    # calculate_plot <- reactive({
    #     
    #     
    #     moth_tax <- abundance_data() %>%
    #         tax_glom(taxrank = "Phylum") %>% #agglomerate at phylum level
    #         transform_sample_counts(function(x) {x/sum(x)}) %>% #transform to rel abund
    #         psmelt() %>%
    #         filter(Abundance > 0.02) %>%
    #         arrange(Phylum)
    # 
    # moth_tax
    # 
    # })
    output$rarefaction_plot = renderPlot({
        
        
        col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
        lty <- c("solid", "dashed", "longdash", "dotdash")
        pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
        head(pars)
        rarefaction <- with(pars[1:7,],
                            rarecurve(t(otu_table(abundance_data())), col = col, lty = lty, step=50, cex=0.5))
        
        
        
    })
    output$taxonomy_plot = renderPlot({
       
        moth_phylum <- abundance_data() %>%
            tax_glom(taxrank = "Phylum") %>% #agglomerate at phylum level
            transform_sample_counts(function(x) {100 * x/sum(x)}) %>% #transform to rel abund
            psmelt() %>%
            filter(Abundance > 1.0) %>%
            arrange(desc(Phylum))
        
        ggplot(moth_phylum, aes(x = new.name, y = Abundance, fill = Phylum)) +
            geom_bar(stat = "identity") +
            scale_fill_viridis(discrete = TRUE, option="D") +
            theme_minimal() +
            geom_text(aes(label = Phylum), color = "white", size = 3, hjust = 0.5, vjust = 3, position =     "stack") +
            ggtitle(label = "Sub-aerial Biofilms: Relative Abundance of Taxa by Phylum") +
            theme(axis.title.x = element_blank(), legend.position = "none")
        
        # 
        # 
        # tax_rank <- input$taxonomy
        # 
        # moth_tax <- abundance_data() %>%
        #     tax_glom(taxrank = tax_rank) %>% #agglomerate at phylum level
        #     transform_sample_counts(function(x) {x/sum(x)}) %>% #transform to rel abund
        #     psmelt() %>%
        #     filter(Abundance > 0.02) %>%
        #     arrange(tax_table(tax_rank))
        # 
        # 
        # phylum_colors <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
        #                    "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
        #                    "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
        # ggplot(moth_tax, aes(x = sample.id, y = Abundance, fill = tax_rank)) +
        #     geom_bar(stat = "identity") +
        #     scale_fill_manual(values = phylum_colors) 
        # 
        
        # ggplot(abundance_data(), aes_string(x = input$taxonomy, y = input$sample, fill = input$color)) + # Note that you need () after filt_dia, since filt_dia() is a function to get the object you want, not the actual object
        #     geom_bar(position="stack", stat="identity") +
        #     theme_minimal() +
        #     theme(text = element_text(size=14),
        #           legend.text = element_text(size=6),
        #           axis.text.x = element_text(size=10,angle=45, hjust=1))
       #  count.data <- abundance_data() %>%
       #      arrange(desc(input$sample)) %>%
       #      mutate(lab.ypos = cumsum(input$taxonomy) - 0.5*input$taxonomy)
       #  
       #  
       #  
       # ggplot(count.data, aes(x=input$taxonomy, y=input$sample)) + 
       #      geom_bar(stat="identity", width=1) + 
       #      coord_polar("y", start=0) + 
       #      geom_text(aes(label = paste0(round(value*100), "%")), position = position_stack(vjust = 0.5)) + 
       #      scale_fill_manual(values=c("#55DDE0", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999")) +
       #      labs(x = NULL, y = NULL, fill = NULL, title = "Phones - Market Share") +
       #      theme_classic() + theme(axis.line = element_blank(),
       #                                      axis.text = element_blank(),
       #                                      axis.ticks = element_blank(),
       #                                      plot.title = element_text(hjust = 0.5, color = "#666666"))
       #  # output$taxonomy_plot <- renderPlotly({
        #   
        #   # 
        #   plot_ly(merged_df_final(),labels = ~x, values = ~y, type = 'pie',
        #           textposition = 'inside',
        #           showlegend = FALSE,
        #           textinfo = 'label+percent') %>%
        #     layout(title = 'Community Composition',
        #            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        #            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
        
        #  ggplot(filt_otu(), aes_string(x = input$xvar, y = input$yvar, fill = input$color)) + # Note that you need () after filt_dia, since filt_dia() is a function to get the object you want, not the actual object
        #   geom_bar(position="stack", stat="identity") +
        #  theme_minimal() +
        # 
        # theme(text = element_text(size=14),
        #      axis.text.x = element_text(size=10,angle=45, hjust=1))
    })
    
}
#making random changes to test commits. Heppy Herlloween
# Run the application 
shinyApp(ui = ui, server = server)