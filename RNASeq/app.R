library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(htmltools)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(SummarizedExperiment)
library('DESeq2')
#library('biomaRt')
#library('fgsea')
library(hrbrthemes) #theme_ipsum
library(pheatmap)
library(Cairo)
library(patchwork)
#library(profvis)
library(plotly)


# Increase the maximum upload size (in bytes)
options(shiny.maxRequestSize = 30 * 1024^2)  # Set to 30 MB

#profvis({
# Define UI for application that draws a histogram
ui <- fluidPage(tabsetPanel(
                  tabPanel("Samples", 
                           sidebarLayout(
                             sidebarPanel(
                              fileInput("Count_File", "Choose a txt file", accept = '.txt')
                            ),
                           mainPanel(
                             tabsetPanel(
                               tabPanel("Summary", tableOutput("Summary_table")),
                               tabPanel("DataTable", dataTableOutput("DataTable_table")),
                               tabPanel("Continous",
                                            HTML("Cell Type Options:<br>"),
                                            HTML("e = Explanted, dissociated adult cardiomyocytes<br>"),
                                            HTML("v = In vivo ventricular myocardium<br>"),
                                            HTML("i = In vitro isolated Cardiomyocytes<br><br>"),
                                            radioButtons(inputId = "subset_data", "Choose the Cell Type to plot", 
                                                        choices = c("e", "v", "i"), selected = 'e'),
                                            actionButton("Sample_btn", "Violin Plot", icon = icon("play"), class = "btn-block"),
                                          plotOutput("Continous_plot"), width = "100%", height =  "1500px")
                             )
                            )
                           )
                  ),
                  tabPanel("Counts",
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("Normalized_Count_File", "Choose a (normalized) txt file", accept = '.txt'),
                               HTML("All analysis on the Counts tab takes normalized filtered data based on the sliders as an input"),
                               sliderInput(inputId =  "percentile_variance", "Include genes that have at least X percentile of variance", min = 0, max = 1, value = .75),
                               sliderInput(inputId =  "min_non_zero", "Include genes that have at least X samples that are non-zero", min = 0, max = 36, value = 3),
                               actionButton(inputId = 'Counts_btn', 'Add Filter', icon = icon('plus'), class = 'btn-block')
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Filtered Summary", tableOutput("Filtered_Summary_table")),
                              tabPanel("Diagnostic ScatterPlot", plotOutput("Diagnostic_plot")),
                              tabPanel("Heatmap", shiny::div(style = "overflow-x: auto; overflow-y: auto; white-space: nowrap;", plotOutput("Heatmap_plot", height = '400px'))),
                              tabPanel("PCA", plotOutput("PCA_plot"))
                            )
                          )
                  ),
                  tabPanel("DE"),
                  tabPanel("Unsure")
            )
                # titlePanel("BF591: Assignment 7"), # app title
                # sidebarLayout( # Sidebar layout with input and output definitions
                #   sidebarPanel(# Sidebar panel for inputs
                #     fileInput("file", "Load differential expression results", accept = '.csv'), #load file
                #     HTML("A volcano plot can be generated with <b> log2 fold-change </b> on the x-axis and <b> p-adjusted </b> on the y-axis. <p>"), 
                #     radioButtons(inputId = "x_axis", "Choose the column for the x-axis", 
                #                  choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), selected = 'log2FoldChange'),
                #     radioButtons(inputId =  "y_axis", "Choose the column for the y-axis", 
                #                  choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), selected = 'padj'),
                #     colourInput("color_base", "Select Base Point Color", "#F033A4"),
                #     colourInput("color_highlight", "Select Highlight Point Color", "#F7D513"),
                #     tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: purple}")), #color slider purple
                #     sliderInput(inputId =  "padj_color", "Select the magnitude of the p adjusted coloring:", min = -300, max = 0, value = -100),
                #     actionButton("run_button", "Plot", icon = icon("play"), class = "btn-block")
                #   ),
                #   mainPanel( # Main panel for displaying outputs
                #     tabsetPanel(
                #       tabPanel("Volcano Plot", plotOutput("volcano")), #Volcano Plot tab
                #       tabPanel("Table", tableOutput("table"))  #Table tab
                #     )
                #   )
                # )
  )

)
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Load_Data: Sample information matrix in CSV format
  load_data <- reactive({
    
    #require an input file
    req(input$Count_File)
    # read the file
    input_file <- read_delim(file = input$Count_File$datapath, delim = "\t") #results of the file upload are nested
    
    return(input_file)
  })
  
  # Tab with a summary of the table that includes a summary of the type and values in each column, 
  # e.g.: Number of rows: X Number of columns: Y
  column_summary <- function(data){
    
    #find mean over column
    data_num_mean <- data %>% select_if(is.numeric) %>% summarise_each(funs(mean(., na.rm=TRUE))) %>% pivot_longer(., cols = everything())
    colnames(data_num_mean) <- c('Column Name','Mean')
    
    #find stdev over column
    data_num_stdev <- data %>% select_if(is.numeric) %>% summarise_each(funs(sd(., na.rm=TRUE))) %>% pivot_longer(., cols = everything())
    colnames(data_num_stdev) <- c('Column Name','sd')
    
    #merge mean and sd
    data_num_updated <- merge(data_num_mean, data_num_stdev, on = 'Column Name')
    
    #format like example
    data_num_updated['Mean (sd) or Distinct Values'] <- paste0(round(data_num_updated$Mean,2), ' +/- (', round(data_num_updated$sd,2), ')')
    data_num_updated <- data_num_updated[c('Column Name', 'Mean (sd) or Distinct Values')]
    
    # Replicate
    data_num_updated['Replicate'] <- str_sub(data_num_updated$'Column Name', -1) 
    
    # Timepoint
    data_num_updated <- data_num_updated %>% mutate("Timepoint" = case_when(
      str_starts(`Column Name`, "ex") ~  str_sub(data_num_updated$`Column Name`, 3,-3),
      (str_sub(`Column Name`, 2, 2) == "P") ~ str_sub(data_num_updated$`Column Name`, 2,3),
      (str_sub(`Column Name`, 2, 2) == "D") ~ str_sub(data_num_updated$`Column Name`, 2,4),
      (str_sub(`Column Name`, 2, 3) == "Ad") ~ "Ad",
      TRUE ~ NA_character_
    ))
    #Cell Type
    data_num_updated <- data_num_updated %>% mutate("Cell Type" = case_when(
      str_starts(`Column Name`, "e") ~ "Explanted, dissociated adult cardiomyocytes",
      str_starts(`Column Name`, "v") ~ "In vivo ventricular myocardium",
      str_starts(`Column Name`, "i") ~ "In vitro isolated Cardiomyocytes",
      TRUE ~ NA_character_
    ))
    
    # # Concatenate unique instances of the Value column into one cell
    data_cat_gene <- paste0(data$Gene[1], ", etc...")
    data_cat_geneID <- paste0(data$GeneID[1], ", etc...")
    data_cat_coord <- paste0(data$Coordinates[1], ", etc...")
    
    # Create a new dataframe with the result
    data_cat_updated <- data.frame(Gene = data_cat_gene, GeneID = data_cat_geneID, Coordinates = data_cat_coord)
    data_cat_updated <- pivot_longer(data_cat_updated, cols = everything())
    colnames(data_cat_updated) <- c('Column Name', 'Mean (sd) or Distinct Values')
    data_cat_updated['Replicate'] <- "N/A"
    data_cat_updated['Timepoint'] <- "N/A"
    data_cat_updated['Cell Type'] <- "N/A"
    
    # Row bind character and numeric values
    updated <- rbind(data_cat_updated, data_num_updated)
    
    # Find type of each column
    col_type <- data.frame(sapply(data, type))
    col_type <- rownames_to_column(col_type)
    colnames(col_type) <- c('Column Name','Type')
    
    updated <- merge(updated, col_type, on = 'Column Name')
    updated <- updated %>% select('Column Name', 'Mean (sd) or Distinct Values', 'Type', everything()) %>% dplyr::arrange('Timepoint')
    
    return(updated)
  }
  
  #' Tab with histograms, density plots, or violin plots of continuous variables.
  #'    If you want to make it fancy, allow the user to choose which column to plot and another column to group by!
  plot_continous_var <- function(data, cell_type_choosen){
    
    # subset data
    data <- data %>% select(starts_with(cell_type_choosen)) 
    # format data
    plot_df  <- pivot_longer(data, cols = everything(), values_to = "Count", names_to = "Sample") #keep columns that start with
    # timepoint
    plot_df <- plot_df %>% mutate("Timepoint" = case_when(
      str_starts(Sample, "e") ~  str_sub(plot_df$Sample, 1,-3),
      str_starts(Sample, "v") ~  str_sub(plot_df$Sample, 1,-3),
      str_starts(Sample, "i") ~ str_sub(plot_df$Sample, 1,-3),
      TRUE ~ NA_character_
    ))
    
    # histogram plot
    h_plt <- ggplot(plot_df, aes(x = Sample, y = -log10(Count), fill = Timepoint)) + 
      geom_violin(adjust=3, alpha=.4, binwidth = 50) +
      theme_ipsum() +
      geom_boxplot(width=0.1) +     
      theme(legend.position="bottom")
    
    return(h_plt)
  } 
  
  filtered_data <- function(data, percentile_variance, min_non_zero){
    
    #'Calculate Variance per row
    data <- data %>% rowwise() %>% dplyr::mutate(Variance = var(c_across(4:39))) #calculate variance across each row
    #'Calculate number of non-zeros per row
    data$Number_of_Non_Zeros <- rowSums(data[4:39] != 0) #sum number of non-zeros and make new column
    #'Calculate number of zeros per row
    data$Number_of_Zeros = 36 - data$Number_of_Non_Zeros
    #'Calculate Row Median across row
    row_median = apply(data[, 4:39], 1, median)
    data$Median_Count <- row_median
    
    #'Include genes with at least X percentile of variance and genes with at least X samples that are non-zero (0-35 columns)
    filtered_data_on_var <- data[data$Variance < quantile(data$Variance, percentile_variance), ] #keep rows when Variance is greater than certain percentile 
    filtered_data_on_var_non_zero <- filtered_data_on_var[filtered_data_on_var$Number_of_Non_Zeros >= min_non_zero, ] #count number of zeros per row 

    return(list(data, filtered_data_on_var_non_zero)) 
    
  }
  
  #' Tab with text or a table summarizing the effect of the filtering, including:
  #' number of samples, total number of genes, number and % of genes passing current filter, number and % of genes not passing current filter
  filtered_summary_table <- function(data, f_data){
    
    num_sample <- ncol(data) -2-3-2
    total_num_genes <- nrow(data)
    num_genes_pass <- nrow(f_data)
    perc_genes_pass <- round(num_genes_pass/total_num_genes, 4)*100
    num_genes_fail <- total_num_genes-num_genes_pass
    perc_genes_fail <- round(num_genes_fail/total_num_genes, 4)*100
    
    table_summarize <- data.frame(
      Metric = c("Number of Samples", 'Total Number of Genes', 'Number of Genes that Pass Filter', 'Percent of Genes that Pass Filter', 'Number of Genes that Fail Filter', 'Percent of Genes that Fail Filter'),
      Value = c(num_sample, total_num_genes, num_genes_pass, paste0(perc_genes_pass, '%'), num_genes_fail, paste0(perc_genes_fail, '%'))
    )
    
    return(table_summarize)
    
  }
  
  #' Tab with diagnostic scatter plots, where genes passing filters are marked in a darker color, and genes filtered out are lighter:
  diagnostic_scatter_plots <- function(data, f_data){

    #' median count vs variance (consider log scale for plot)
    median_v_var <- ggplot() + 
      geom_point(data = data, aes(x=-log10(Median_Count), y=-log10(Variance), color = "Fail_Filter"), alpha = 0.5) +
      geom_point(data = f_data, aes(x=-log10(Median_Count), y=-log10(Variance), color = "Pass_Filter"), alpha = 0.5) +
      scale_color_manual(values = c("Fail_Filter" = "lightblue", "Pass_Filter" = "darkblue")) +
      ggtitle('Median count vs variance') +
      coord_cartesian(clip = 'off')
    
    #' median count vs number of zeros
    median_v_zero <- ggplot() + 
      geom_point(data = data, aes(x=-log10(Median_Count), y=Number_of_Zeros, color = "Fail_Filter"), alpha = 0.5) +
      geom_point(data = f_data, aes(x=-log10(Median_Count), y=Number_of_Zeros, color = "Pass_Filter"), alpha = 0.5) +
      scale_color_manual(values = c("Fail_Filter" = "lightblue", "Pass_Filter" = "darkblue")) +
      ggtitle('Median count vs number of zeros') +
      coord_cartesian(clip = 'off')
    
    return(median_v_var | median_v_zero)
  }
  #' Tab with a clustered heatmap of counts remaining after filtering
  #' consider enabling log-transforming counts for visualization
  clustered_heatmap <- function(f_data) {
    
    #Subset data
    filtered_data_numeric <- f_data[, c(1, 4:39)]
    
    #rownames as Gene
    filtered_data_numeric <- filtered_data_numeric %>% column_to_rownames(var = "Gene")
    #heatmap takes in matrix
    filtered_data_matrix <- as.matrix(filtered_data_numeric)
    
    # Choose a color palette from RColorBrewer
    my_palette <- colorRampPalette(brewer.pal(9, "RdBu"))(100)
    
    #clustered heatmap with scale
    hetmap <- pheatmap(filtered_data_matrix, 
                       col = my_palette, 
                       main = "Clustered Heatmap of Counts",
                       height = 100, 
                       annotation_row_text_size = 8,
                       fontsize_row = 5,
                       cellheight = 5, 
                       cellwidth = 5)
  
    
    return(hetmap)
  }
  #' 
  # # Helper function to render plot
  # renderVolcanoPlot <- function() {
  #   req(load_data()) #require load_data()
  #   isolate({ #isolate the volcano plot function from changes in other reactive values.
  #     volcano_plot(load_data(), input$x_axis, input$y_axis, input$padj_color, input$color_base, input$color_highlight)
  #   })
  # }
  # 
  # # Helper function to render table
  # renderDataTable <- function() {
  #   req(load_data()) #require load_data()
  #   isolate({ #isolate the draw table function from changes in other reactive values.
  #     draw_table(load_data(), input$padj_color)
  #   })
  # }
  # 
  # # Render the volcano plot
  # output$volcano <- renderPlot({renderVolcanoPlot()})
  #
  # # Update data and re-render on actionButton click
  # observeEvent(input$run_button, {
  #   output$volcano <- renderPlot({renderVolcanoPlot()})
  #   output$table <- renderTable({renderDataTable()})
  # })
  # 
  # Render objects for Sample tab
  output$Summary_table <- renderTable({column_summary(load_data())})
  
  output$DataTable_table <- renderDataTable({load_data()})
  
  renderContinoutVarPlot <- function() {
    req(input$subset_data)
    isolate({ 
      plot_continous_var(load_data(), input$subset_data)
    })
  }
  
  output$Continous_plot <- renderPlot({renderContinoutVarPlot()}, height = 1500)
  
  observeEvent(input$Sample_btn, {
    output$Continous_plot <- renderPlot({renderContinoutVarPlot()}, height = 1500)})
  
  # Render objects for Counts tab
  renderFiltered_S_table <- function() {
    req(load_data())
    isolated_data <- isolate({filtered_data(load_data(), input$percentile_variance, input$min_non_zero)
    })
    filtered_summary_table(isolated_data[[1]], isolated_data[[2]])
  }
  output$Filtered_Summary_table <- renderTable({renderFiltered_S_table()})
  observeEvent(input$Counts_btn, {
    output$Filtered_Summary_table <- renderTable({renderFiltered_S_table()})})
  
  
  renderScatterePlot <- function() {
    req(load_data())
    isolated_data <- isolate({filtered_data(load_data(), input$percentile_variance, input$min_non_zero)
    })
    diagnostic_scatter_plots(isolated_data[[1]], isolated_data[[2]])
  }
  output$Diagnostic_plot <- renderPlot({renderScatterePlot()})
  observeEvent(input$Sample_btn, {
    output$Diagnostic_plot <- renderPlot({renderScatterePlot()})})
  
  renderHeatmapPlot <- function() {
    req(load_data())
    isolated_data <- isolate({filtered_data(load_data(), input$percentile_variance, input$min_non_zero)
    })
    clustered_heatmap(isolated_data[[2]])
  }
  output$Heatmap_plot <-  renderPlot({renderHeatmapPlot()}, height = 400)
  observeEvent(input$Sample_btn, {
    output$Heatmap_plot <- renderPlot({renderHeatmapPlot()}, height = 400)})
  
  

}
# Run the application
shinyApp(ui = ui, server = server)