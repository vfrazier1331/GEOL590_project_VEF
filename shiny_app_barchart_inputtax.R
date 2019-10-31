########
# Shiny app to create a dynamically-filterable visualization of the diamonds app
########
#install.packages("plotly")
# These bits get run before any of the rest of the code
# Note: contrary to what I told you on Monday, the use of the global.R file is no longer recommended.
# At present, I'm not sure why.

# Slight edits by Warren
library(shiny)
library(tidyverse)

# source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
#        local = TRUE)

library(phyloseq)
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
            
            # This is a range slider (i.e. there's a max and min). It is set that way by "value" (the starting value), which is a 2-element vector
            #    sliderInput("otu.number",
            #               "OTUs",
            #              min = min.otu,
            #             max = max.otu,
            #            value = c(min.otu, max.otu)),
            actionButton(inputId = "upload",label="Upload"),
            
            uiOutput("sample"),
            uiOutput("taxonomy"),
            uiOutput("color")
            
        ),
        
        # Show a plot of diamonds data frame. This output doesn't care what that plot is, only that it will be associated with output$diamonds_plot
        mainPanel(
            plotOutput("taxonomy_plot")
        )
    )
)


server <- function(input, output) {
    ###making dataframe out of two files: .taxonomy and .shared files
    abundance_data <- reactive({ 
        
        sharedfile <- input$shared_file
        taxfile <- input$tax_file
        
        phyloseq_import <- import_mothur(mothur_shared_file = sharedfile$datapath, mothur_constaxonomy_file = taxfile$datapath)
        #             mothur_constaxonomy_file = "data/lava_tube_merged.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy")
        #makes an otu_table file and a tax_table file
        ###################################################mothur_shared_file = input$shared_file, mothur_constaxonomy_file = input$tax_file
        #MAKE DF
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
        merged_df <- cbind(taxdf, OTUdf)
        
        ########################################################
        #FORMAT DF: re-order and re-name columns
        merged_df$otu.number <- seq.int(nrow(merged_df))
        
        merged_df1 <- merged_df[,c(10,7,8,9,1,2,3,4,5,6)]
        
        names(merged_df1)[5] <- "Kingdom"
        names(merged_df1)[6] <- "Phylum"
        names(merged_df1)[7] <- "Class"
        names(merged_df1)[8] <- "Order"
        names(merged_df1)[9] <- "Family"
        names(merged_df1)[10] <- "Genus"
        merged_df_avg <- cbind(merged_df1,IVMean=rowMeans(merged_df1[2:4], na.rm=TRUE))
        
        merged_df_rmsing <- merged_df_avg %>% filter(IVMean >= 2)
        
        merged_df_rmsing
    })
    
    output$sample = renderUI({
        df <- abundance_data()
        axis_vars = names(df)
        selectInput(inputId = "sample",
                    label = "Sample",
                    choices = axis_vars,
                    selected = "A")
    })
    output$taxonomy = renderUI({
        df <- abundance_data()
        axis_vars = names(df)
        selectInput(inputId = "taxonomy",
                    label = "Taxonomy",
                    choices = axis_vars,
                    selected = "Phylum")
    })
    
    output$color = renderUI({
        df <- abundance_data()
        axis_vars = names(df)
        selectInput(inputId = "color",
                    label = "Color",
                    choices = axis_vars,
                    selected = "Order")
    })
    
    output$taxonomy_plot = renderPlot({
        ggplot(abundance_data(), aes_string(x = input$taxonomy, y = input$sample, fill = input$color)) + # Note that you need () after filt_dia, since filt_dia() is a function to get the object you want, not the actual object
            geom_bar(position="stack", stat="identity") +
            theme_minimal() +
            theme(text = element_text(size=14),
                  legend.text = element_text(size=6),
                  axis.text.x = element_text(size=10,angle=45, hjust=1))
        
        # output$taxonomy_plot <- renderPlotly({
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

# Run the application 
shinyApp(ui = ui, server = server)

