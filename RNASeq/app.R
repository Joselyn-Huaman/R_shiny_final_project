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
library(plotly)
library(future)
plan(multisession)

# Increase the maximum upload size (in bytes)
options(shiny.maxRequestSize = 30 * 1024^2)  # Set to 30 MB

#profvis({
# Define UI for application that draws a histogram
ui <- fluidPage(theme = bs_theme(version = 5, bootswatch = "solar"), 
                titlePanel("BF591: Analysis of RNA-Seq Data"), 
                tabsetPanel(
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
                                            radioButtons(inputId = "subset_data", label = HTML("<strong>Choose the Cell Type to plot</strong>"), 
                                                        choices = c("e", "v", "i"), selected = 'e'),
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
                               sliderInput(inputId =  "percentile_variance", label = HTML("Include genes that have at least <em>X</em> percentile of variance"), min = 0, max = 1, value = .75),
                               sliderInput(inputId =  "min_non_zero", label = HTML("Include genes that have at least <em>X</em> samples that are non-zero"), min = 0, max = 36, value = 3),
                               actionButton(inputId = 'Counts_btn', 'Add Filter', icon = icon('plus'), class = 'btn-block')
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Filtered Summary", tableOutput("Filtered_Summary_table")),
                              tabPanel("Diagnostic ScatterPlot", plotOutput("Diagnostic_plot")),
                              tabPanel("Heatmap", shiny::div(style = "overflow-x: auto; overflow-y: auto; white-space: nowrap;", plotOutput("Heatmap_plot", height = '400px'))),
                              tabPanel("PCA",
                                       HTML("<strong>Choose Principal Components to Plot</strong> <br>"),
                                       fluidRow(
                                         column(6, numericInput("First_PC", "Select principal components to plot on the x-axis", 1, min = 1, max = 40)),
                                         column(6, numericInput("Second_PC", "Select principal components to plot on the y-axis", 2, min = 1, max = 40))
                                       ), 
                                       plotOutput("PCA_plot"))
                            )
                          )
                          )
                  ),
                  tabPanel("DE",
                           mainPanel(
                             HTML("<strong>The normalized count data has undergone different expression analysis via DESeq2</strong> <br> <br>"),
                             radioButtons(inputId = "cell_stage", label = HTML("<strong>Choose the Cell Stage to subset</strong>"), 
                                         choices = c("ex", "vP", "vD", "iP", "iD"), selected = 'vP'),
                             tabsetPanel(
                               tabPanel("Differential Expression Results", dataTableOutput("Diff_eq_table")),
                               tabPanel("Volcano Plot",
                                        sidebarLayout(
                                            sidebarPanel(
                                              HTML("A volcano plot can be generated with <b> log2 fold-change </b> on the x-axis and <b> p-adjusted </b> on the y-axis. <p>"),
                                              radioButtons(inputId = "x_axis", "Choose the column for the x-axis",
                                                           choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), selected = 'log2FoldChange'),
                                              radioButtons(inputId =  "y_axis", "Choose the column for the y-axis",
                                                           choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), selected = 'padj'),
                                              colourInput("color_base", "Select Base Point Color", "#F033A4"),
                                              colourInput("color_highlight", "Select Highlight Point Color", "#F7D513"),
                                              tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: purple}")), #color slider purple
                                              sliderInput(inputId =  "padj_color", "Select the magnitude of the p adjusted coloring:", min = -300, max = 0, value = -100),
                                              actionButton("run_button", "Plot", icon = icon("image"), class = "btn-block")
                                            ),
                                            mainPanel(
                                              plotOutput("volcano_plot")
                                            )
                                        )
                               )
                             )
                           )
                  ),
                  tabPanel("Unsure")
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
    
    #Cell Type
    data_num_updated <- data_num_updated %>% mutate("Cell Type" = case_when(
      str_starts(`Column Name`, "e") ~ "Explanted, dissociated adult cardiomyocytes",
      str_starts(`Column Name`, "v") ~ "In vivo ventricular myocardium",
      str_starts(`Column Name`, "i") ~ "In vitro isolated Cardiomyocytes",
      TRUE ~ NA_character_
    ))
    
    # Timepoint
    data_num_updated <- data_num_updated %>% mutate("Timepoint" = case_when(
      str_starts(`Column Name`, "ex") ~  str_sub(data_num_updated$`Column Name`, 3,-3),
      (str_sub(`Column Name`, 2, 2) == "P") ~ str_sub(data_num_updated$`Column Name`, 2,3),
      (str_sub(`Column Name`, 2, 2) == "D") ~ str_sub(data_num_updated$`Column Name`, 2,4),
      (str_sub(`Column Name`, 2, 3) == "Ad") ~ "Ad",
      TRUE ~ NA_character_
    ))
    
    # Replicate
    data_num_updated['Replicate'] <- str_sub(data_num_updated$'Column Name', -1) 
    
    #Cell Stage
    data_num_updated <- data_num_updated %>% mutate("Cell Stage" = case_when(
      str_starts(`Column Name`, "e") ~ "Adult CM Exmplant",
      str_starts(`Column Name`, "vP") ~ "In vivo Maturation using ventricular samples",
      str_starts(`Column Name`, "vD") ~ "In vivo Regeneration using ventricular samples",
      str_starts(`Column Name`, "iP") ~ "In vivo Maturation using iCM samples",
      str_starts(`Column Name`, "iD") ~ "In vivo Regeneration using iCM samples",
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
    data_cat_updated['Cell Stage'] <- "N/A"
    
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
    filtered_data_on_var <- data[data$Variance >= quantile(data$Variance, percentile_variance), ] #keep rows when Variance is greater than certain percentile 
    filtered_data_on_var_non_zero <- filtered_data_on_var[filtered_data_on_var$Number_of_Non_Zeros >= min_non_zero, ] #count number of zeros per row 

    return(list(data, filtered_data_on_var_non_zero)) 
    
  }
  #' Tab with text or a table summarizing the effect of the filtering, including:
  #' number of samples, total number of genes, number and % of genes passing current filter, number and % of genes not passing current filter
  filtered_summary_table <- function(data, f_data){
    
    num_sample <- ncol(data) -2-3-3
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
  #' Tab with a scatter plot of principal component analysis projections.
  PCA_plot <- function(data, PC_x, PC_y) {
    #Subset data
    data <- data[, c(4:39)]
    
    #convert to dataframe to complete PC analysis
    data <- as.data.frame(data)
    
    #PCA calculation
    pca_results <- prcomp(t(data))
    
    # % variance explained
    pc_variance_explained <- round((((pca_results$sdev)**2)/sum(((pca_results$sdev)**2)))*100, 0)
    #name pc_variance_explained (vector) to variance_explained (tibble)
    var_exp_tib <- tibble("variance_explained_percent" = pc_variance_explained)
    #extract principal components that required
    pca_results_subset <- data.frame(pca_results$x[,c(PC_x, PC_y)])
    #Get column_names
    first_column_name <- colnames(pca_results_subset)[1]
    second_column_name <- colnames(pca_results_subset)[2]
    
    #'allow the user to select which principal components to plot in a scatter plot (e.g. PC1 vs PC2)
    scatterplot_PCA <- pca_results_subset %>%
                        ggplot(aes(x = .data[[first_column_name]], y = .data[[second_column_name]])) +
                        geom_point() +
                        ggtitle('Principal Component Scatterplot') +
                        xlab(paste0('PC',PC_x, ": ", var_exp_tib[PC_x,], '% Variance')) +
                        ylab(paste0('PC', PC_y, ": ", var_exp_tib[PC_y,], '% Variance'))
    
    return(scatterplot_PCA)
  }
  
  #' Tab with sortable table displaying differential expression results
  diff_eq <- function(filtered_data, metadata, cell_stage){
    
    if (cell_stage == "vP") {
      cell_stage <- c("vP", "vA")
      factor_level <- c("P0", "P4", "P7", "Ad")
    }
    
    else if (cell_stage == "ex") {
      factor_level <- c("0hr", "24hr", "48hr", "72hr")
    } 
    else if (cell_stage == "iP") {
      factor_level <- c("P0", "P4")
    }
    else if (cell_stage == "iD") {
      factor_level <- c("D7S", "D7R")
    }
    else {
      factor_level <- c("D1S", "D1R", "D7S", "D7R")
    }
    
    #Subset data
    filtered_data <- filtered_data[, c(1,4:39)]
    #rownames as Gene, select certain columns, make all values integers
    filtered_data <- filtered_data %>% column_to_rownames(var = "Gene") %>% dplyr::select(starts_with(cell_stage)) %>% round(.)
    #Deseq2 takes in matrix
    filtered_data_matrix <- as.matrix(filtered_data)
    
    #column metadata - keep only important info
    metadata <- metadata %>% select("Column Name", "Cell Stage", "Cell Type", "Timepoint", "Replicate") %>% 
      subset(Replicate != "N/A") %>% 
      rename(Sample = "Column Name") %>% 
      filter(Sample %in% colnames(filtered_data_matrix))
    
    #tibble of column metadata
    metadata <- as_tibble(metadata)
    # Specify reference levels for Timepoint and Cell Type
    metadata$Timepoint <- factor(metadata$Timepoint, levels = factor_level)
    
    #store counts matrix and sample df in a SummarizedExperiments object
    se <- SummarizedExperiment(assays = list(counts = filtered_data_matrix), #subsetted counts matrix
                               colData = metadata) #store your sample dataframe as colData
    ddsSE <- DESeqDataSet(se, design = ~Timepoint)
    
    #results from DESeq2 as df
    dds <- DESeq(ddsSE)
    dds_results <- results(dds)   #DESeqDataSet object updated ???
    dds_results_df <- as.data.frame(dds_results)
    dds_results_df <- tibble::rownames_to_column(dds_results_df, "Genes") 
    return(dds_results_df)
  }
  
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
  
  output$Continous_plot <- renderPlot({
    req(input$subset_data)
    isolate({plot_continous_var(load_data(), input$subset_data)})
  }, height = 1500)
  
  # Render objects for Counts tab
  # Function to process data based on file input and filtering parameters
  processData <- function(file, percentile_variance, min_non_zero) {
    req(file)
    data <- filtered_data(file, percentile_variance, min_non_zero)
    list(
      summary_table = filtered_summary_table(data[[1]], data[[2]]),
      diagnostic_plot = diagnostic_scatter_plots(data[[1]], data[[2]]),
      heatmap_plot = clustered_heatmap(data[[2]]),
      pca_plot = PCA_plot(data[[2]], input$First_PC, input$Second_PC),
      diff_eq_table = diff_eq(data[[2]], column_summary(load_data()), input$cell_stage)
    )
  }
  
  # Reactive expression for file input
  data <- eventReactive(load_data(), {
    processData(load_data(), input$percentile_variance, input$min_non_zero)
  })
  
  # Render plots based on reactive expression
  output$Filtered_Summary_table <- renderTable({
    data()$summary_table
  })
  
  output$Diagnostic_plot <- renderPlot({
    data()$diagnostic_plot
  })
  
  output$Heatmap_plot <- renderPlot({
    data()$heatmap_plot
  }, height = 400)
  
  output$PCA_plot <- renderPlot({
    data()$pca_plot
  })
  
  output$Diff_eq_table <- renderDataTable({
    data()$diff_eq_table
  })
#   renderFiltered_S_table <- function() {
#     req(load_data())
#     isolated_data <- isolate({filtered_data(load_data(), input$percentile_variance, input$min_non_zero)
#     })
#     filtered_summary_table(isolated_data[[1]], isolated_data[[2]])
#   }
#   output$Filtered_Summary_table <- renderTable({renderFiltered_S_table()})
#   observeEvent(input$Counts_btn, {
#     output$Filtered_Summary_table <- renderTable({renderFiltered_S_table()})})
#   
#   
#   renderScatterePlot <- function() {
#     req(load_data())
#     isolated_data <- isolate({filtered_data(load_data(), input$percentile_variance, input$min_non_zero)
#     })
#     diagnostic_scatter_plots(isolated_data[[1]], isolated_data[[2]])
#   }
#   output$Diagnostic_plot <- renderPlot({renderScatterePlot()})
#   observeEvent(input$Counts_btn, {
#     output$Diagnostic_plot <- renderPlot({renderScatterePlot()})})
#   
#   renderHeatmapPlot <- function() {
#     req(load_data())
#     isolated_data <- isolate({filtered_data(load_data(), input$percentile_variance, input$min_non_zero)
#     })
#     clustered_heatmap(isolated_data[[2]])
#   }
#   output$Heatmap_plot <-  renderPlot({renderHeatmapPlot()}, height = 400)
#   observeEvent(input$Counts_btn, {
#     output$Heatmap_plot <- renderPlot({renderHeatmapPlot()}, height = 400)})
#   
#   renderPCPlot <- function() {
#     req(load_data())
#     isolated_data <- isolate({filtered_data(load_data(), input$percentile_variance, input$min_non_zero)
#     })
#     PCA_plot(isolated_data[[2]], input$First_PC, input$Second_PC)
#   }
#   output$PCA_plot <-  renderPlot({renderPCPlot()})
#   observeEvent(input$Counts_btn, {
#     output$Heatmap_plot <- renderPlot({renderPCPlot()})})
# 
}
# Run the application
shinyApp(ui = ui, server = server)