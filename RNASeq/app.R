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
library(heatmaply)
library(gplots)

# Increase the maximum upload size (in bytes)
options(shiny.maxRequestSize = 30 * 1024^2)  # Set to 30 MB

#profvis({
# Define UI for application that draws a histogram
ui <- fluidPage(theme = bs_theme(version = 5, bootswatch = "pulse"), 
                HTML("<h1>Analysis of RNA-Seq Data</h1>"), 
                HTML("<h5><em>Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration</em></h5>"),
                HTML("Joselyn Huaman Argandona"),
                tabsetPanel(
                  tabPanel("Sample Information Exploration", 
                           sidebarLayout(
                             sidebarPanel(
                              fileInput("Count_File", "Choose a txt file", accept = '.txt')
                            ),
                           mainPanel(
                             tabsetPanel(
                               tabPanel("Column Summary of Counts File", dataTableOutput("Summary_table")),
                               tabPanel("Counts DataTable", dataTableOutput("DataTable_table")),
                               tabPanel("Plotting Counts of Continous Variables",
                                        wellPanel(
                                            radioButtons(inputId = "subset_data", label = HTML("Choose the Cell Stage to plot"), 
                                                        choices = c("Adult CM Exmplant", "In vivo Maturation using ventricular samples", "In vivo Regeneration using ventricular samples",  "In vivo Maturation using iCM samples",  "In vivo Regeneration using iCM samples"), selected = "Adult CM Exmplant")
                                        ),
                                            plotOutput("Continous_plot"), width = "100%", height =  "1500px")
                             )
                            )
                           )
                  ),
                  tabPanel("Counts Exploration",
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("Normalized_Count_File", "Choose a txt file", accept = '.txt'),
                               HTML("Analysis on the Counts tab takes filtered data based on the below sliders <br>"),
                               sliderInput(inputId =  "percentile_variance", label = HTML("Include genes that have at least <em>X</em> percentile of variance"), min = 0, max = 1, value = .75),
                               sliderInput(inputId =  "min_non_zero", label = HTML("Include genes that have at least <em>X</em> samples that are non-zero"), min = 0, max = 36, value = 3),
                               downloadButton("download_filter_count", "Download")
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Filtered Summary", tableOutput("Filtered_Summary_table")),
                              tabPanel("Diagnostic ScatterPlot", 
                                       HTML("Rank is in ascending order <br>"),
                                       plotOutput("Diagnostic_plot")),
                              tabPanel("Heatmap", plotOutput("Heatmap_plot", height = "1500px")),
                              tabPanel("PCA",
                                       wellPanel(
                                         fluidRow(
                                           column(6, numericInput("First_PC", "Select principal components to plot on the x-axis", 1, min = 1, max = 40)),
                                           column(6, numericInput("Second_PC", "Select principal components to plot on the y-axis", 2, min = 1, max = 40))
                                       )), 
                                       plotOutput("PCA_plot"))
                            )
                          )
                          )
                  ),
                  tabPanel("DE",
                           mainPanel(
                             sidebarLayout(
                               sidebarPanel(
                                 fileInput("Normalized_Count_File_for_DE", "Choose a file to perform DESeq2", accept = '.csv'),
                                 radioButtons(inputId = "cell_stage", label = HTML("Choose the Cell Stage to subset"), 
                                             choices = c("ex", "vP", "vD", "iP", "iD"), selected = 'vP')),
                                 mainPanel(
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
                                                  colourInput("color_base", "Select Base Point Color", "grey"),
                                                  colourInput("color_highlight", "Select Highlight Point Color", "#3C0E8C"),
                                                  tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: purple}")), #color slider purple
                                                  sliderInput(inputId =  "padj_color", "Select the magnitude of the p adjusted coloring:", min = -30, max = 0, value = -5),
                                                  actionButton("run_button", "Plot", icon = icon("image"), class = "btn-block")
                                                ),
                                                mainPanel(
                                                  plotOutput("volcano_plot")
                                                )
                                        )))
                               )
                             )
                           )
                  ),
                  tabPanel("Gene Set Ennrichment Analysis",
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("FGSEA_file", "Choose a FGSEA txt file", accept = '.txt')
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Barplot of Top Pathways", 
                                          sidebarPanel(
                                            sliderInput(inputId =  "filter_by_adjp_3", "Adjust maxiumum included adjusted p-value (10^X)", min = -26, max = 0, value = -10)
                                          ), 
                                          
                                          mainPanel(
                                          plotOutput("barplot_fgsea"))
                                        , width = "100%", height =  "1000px"),
                                 tabPanel("DataTable of FGSEA Results",
                                          sidebarLayout(
                                            sidebarPanel(
                                              sliderInput(inputId =  "filter_by_adjp_1", "Adjust maxiumum included adjusted p-value (10^X)", min = -26, max = 0, value = -10),
                                              radioButtons(inputId =  "NES_type", "Select type of NES pathways",
                                                           choices = c("positive", "negative", "both"), selected = 'both'),
                                              downloadButton("download_fgsea_result", "Download")
                                            ),
                                          mainPanel(
                                            dataTableOutput("fgsea_table")
                                            )
                                          )
                                 ),
                                 tabPanel("Scatterplot of NES vs padj",
                                          sidebarLayout(
                                            sidebarPanel(
                                              sliderInput(inputId = "filter_by_adjp_2", "Adjust maxiumum included adjusted p-value (10^X)", min = -26, max = 0, value = -10)
                                            ),
                                          mainPanel(
                                            plotOutput("Scatterplot_NES")
                                            )
                                 )
                               )
                             )
                           )
                           )
                   )
                )

)
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  #Move slider for FGSEA
  observe({
    val <- input$filter_by_adjp_3
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "filter_by_adjp_2", value = val)
    updateSliderInput(session, "filter_by_adjp_1", value = val)
    
  })
  
  observe({
    val <- input$filter_by_adjp_2
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "filter_by_adjp_1", value = val)
    updateSliderInput(session, "filter_by_adjp_3", value = val)
    
  })
  
  observe({
    val <- input$filter_by_adjp_1
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "filter_by_adjp_2", value = val)
    updateSliderInput(session, "filter_by_adjp_3", value = val)
    
  })
  
  # Load_Data: Sample information matrix in CSV format
  load_data <- reactive({
    
    #require an input file
    req(input$Count_File)
    # read the file
    input_file <- read_delim(file = input$Count_File$datapath, delim = "\t") #results of the file upload are nested
    
    return(input_file)
  })
  
  # Tab 1: Sample
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
      str_starts(`Column Name`, "vA") ~ "In vivo Maturation using ventricular samples",
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
    data_cat_updated['Cell Stage'] <- "N/A"
    
    # Row bind character and numeric values
    updated <- rbind(data_cat_updated, data_num_updated)
    
    # Find type of each column
    col_type <- data.frame(sapply(data, type))
    col_type <- rownames_to_column(col_type)
    colnames(col_type) <- c('Column Name','Type')
    
    updated <- merge(updated, col_type, on = 'Column Name') 
    updated <- updated[order(updated$Timepoint), ]
    
    return(updated)
  }
  #' Tab with histograms, density plots, or violin plots of continuous variables.
  #'    If you want to make it fancy, allow the user to choose which column to plot and another column to group by!
  plot_continous_var <- function(data, cell_stage_choosen){
    
    if (cell_stage_choosen == "Adult CM Exmplant"){ 
      data_sub <- dplyr::select(data, starts_with("e")) } 
    else if (cell_stage_choosen == "In vivo Maturation using ventricular samples"){ 
      data_1 <- dplyr::select(data, starts_with("vP")) 
      data_2 <- dplyr::select(data, starts_with("vA")) 
      data_sub <- cbind(data_1, data_2) } 
    else if (cell_stage_choosen == "In vivo Regeneration using ventricular samples"){ 
      data_sub <- dplyr::select(data, starts_with("vD")) } 
    else if (cell_stage_choosen == "In vivo Maturation using iCM samples"){ 
      data_sub <- dplyr::select(data, starts_with("iP")) } 
    else if (cell_stage_choosen == "In vivo Regeneration using iCM samples"){ 
      data_sub <- dplyr::select(data, starts_with("iD")) }
    else {
      stop("Invalid value for cell_stage_choosen.")
    }
    # subset data
    #data <- data %>% select(starts_with(cell_type_choosen)) 
    # format data
    plot_df  <- pivot_longer(data_sub, cols = everything(), values_to = "Count", names_to = "Sample") #keep columns that start with
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
      theme(legend.position="bottom") +
      ggtitle(paste("Counts Distribution of cells in", cell_stage_choosen))
    
    return(h_plt)
  } 
  
  # Tab 2: Counts
  load_data_count <- reactive({
    #require an input file
    req(input$Normalized_Count_File)
    # read the file
    input_file <- read_delim(file = input$Normalized_Count_File$datapath, delim = "\t") #results of the file upload are nested
    
    return(input_file)
  })
  
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
    
    num_sample <- if (!is.null(data)) ncol(data) - 2 - 3 - 3 else 0
    total_num_genes <- if (!is.null(data)) nrow(data) else 0
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
      geom_point(data = data, aes(x=rank(Median_Count), y= Variance, color = "Fail_Filter"), alpha = 0.5) +
      geom_point(data = f_data, aes(x=rank(Median_Count), y= Variance, color = "Pass_Filter"), alpha = 1) +
      scale_color_manual(values = c("Fail_Filter" = "#00BFC4", "Pass_Filter" = "#7CAE00")) +
      ggtitle('Median count vs variance') +
      coord_cartesian(clip = 'off')
    
    #' median count vs number of zeros
    median_v_zero <- ggplot() + 
      geom_point(data = data, aes(x=rank(Median_Count), y= Number_of_Zeros, color = "Fail_Filter"), alpha = 0.5) +
      geom_point(data = f_data, aes(x=rank(Median_Count), y= Number_of_Zeros, color = "Pass_Filter"), alpha = 1) +
      scale_color_manual(values = c("Fail_Filter" = "#00BFC4", "Pass_Filter" = "#7CAE00")) +
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
    my_palette <- colorRampPalette(brewer.pal(9, "PRGn"))(100)
    
    # Log-transform the data
    log_transformed_data <- log2(filtered_data_matrix + 1)  # Adding 1 to avoid log(0)
    
    #clustered heatmap with scale
    heatmap.2(log_transformed_data, 
             col = my_palette, 
             main = "Clustered Heatmap of Filtered Counts",
             key = TRUE,
             key.title = "Log-Transformed Counts",
             key.xlab = "Log2(count + 1)",
             trace = "none",
             height = 1000)
  
    
   # return((hetmap))
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
                        geom_point(color = "#3C0E8C") +
                        ggtitle('Principal Component Scatterplot') +
                        xlab(paste0('PC',PC_x, ": ", var_exp_tib[PC_x,], '% Variance')) +
                        ylab(paste0('PC', PC_y, ": ", var_exp_tib[PC_y,], '% Variance'))
    
    return(scatterplot_PCA)
  }
  
  #' Tab 3: Differential Expression
  load_data_de <- reactive({
    #require an input file
    req(input$Normalized_Count_File_for_DE)
    # read the file
    input_file <- read_csv(file = input$Normalized_Count_File_for_DE$datapath) #results of the file upload are nested
    
    return(input_file)
  })
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
    filtered_data <- filtered_data %>% column_to_rownames(var = "Gene") %>% dplyr::select(starts_with(cell_stage)) %>% ceiling(.)
    #Deseq2 takes in matrix
    filtered_data_matrix <- as.matrix(filtered_data)
    
    #column metadata - keep only important info
    metadata <- metadata %>% dplyr::select("Column Name", "Cell Stage", "Timepoint", "Replicate") %>% 
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
  #' Tab with content similar to that described in [Assignment 7] (volcano plot)
  volcano_plot <-function(dataf, x_name, y_name, slider, color1, color2) {
    
    slider <- 1 * 10**slider #calculate slider to fit chart
    
    plot_v <- ggplot(data = dataf) + #load data
      geom_point(aes(x = .data[[x_name]], y = -log10(.data[[y_name]]), col = .data[[y_name]] < slider)) + #get data columns and color based on if statement 
      scale_colour_manual(name = paste(y_name, '>', slider), values = setNames(c(color2,color1),c(T, F))) + #set colors and name
      theme(legend.position="bottom") + #change legend position
      ggtitle(paste0(x_name, "vs -log10(", y_name, ")"))
    
    return(plot_v)
  }

  #' Tab 4: FGSEA
  # Load_Data: Sample information matrix in CSV format
  load_fgsea_data <- reactive({
    
    #require an input file
    req(input$FGSEA_file)
    # read the file
    input_file <- read_delim(file = input$FGSEA_file$datapath, delim = ",") #results of the file upload are nested
  
    return(input_file)
  })
  #' Barplot of fgsea NES for top pathways selected by slider
  barplot_fgsea <- function(file, num_paths){
    
    num_paths <- 10**num_paths
    
    #keep rows that have a padj or lower and arrange by NES
    fgsea_results_top <- file %>% filter(padj < num_paths) %>% arrange(desc(NES))
    
    stacked_bar <- fgsea_results_top %>% 
      mutate(fill_color = ifelse(NES > 0, "positive NES", "negative NES")) %>%
      ggplot(aes(x = reorder(pathway, NES), y = NES, fill = fill_color)) +
      geom_col() +
      scale_fill_manual(values = c("positive NES" = "#7CAE00", "negative NES" = "#00BFC4")) +
      coord_flip() +
      ggtitle("fgsea results for CP MSigDB genes") +
      ylab("Normalized Enrichment Score (NES)") +
      xlab("Pathway")
    
    return(stacked_bar)
  }
  #' Filtered data table displaying the FGSEA results
  filter_fgsea_res <- function(file, num_paths, sign){
    
    if (sign == "positive"){
      file <- file %>% filter(NES > 0)
    }
    else if (sign == "negative"){
      file <- file %>% filter(NES < 0)
    }
    else {
      file <- file
    }
    
    num_paths <- 10**num_paths
    
    #keep rows that have a padj or lower and arrange by NES
    file <- filter(file, padj < num_paths)
    
    return(file)
  }
  #' Scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis, with gene sets below threshold in grey color
  scatter_fgsea <- function(file, num_paths){
    
    num_paths <- 10**num_paths
    
    scatter <- file %>% 
      ggplot(aes(x = NES, y = -log10(padj))) +
      geom_point(aes(color = ifelse(-log10(padj) < -log10(num_paths), "below_threshold", "above_threshold"))) +
      scale_color_manual(values = c("below_threshold" = "grey", "above_threshold" = "#3C0E8C"), name = paste("padj threshold =", num_paths)) +
      ggtitle("Scatterplot of NES vs -log10(padj)") +
      xlab("Normalized Enrichment Score (NES)")
    
    return(scatter)
  }

  
  # Render objects for Sample tab
  output$Summary_table <- renderDataTable({column_summary(load_data())})
  output$DataTable_table <- renderDataTable({load_data()})
  output$Continous_plot <- renderPlot({
    req(input$subset_data)
    plot_continous_var(load_data(), input$subset_data)
  }, height = 1500)
  
  # Render objects for Counts
  # Function to process data based on file input and filtering parameters
  # Define reactive values to store filter values

  
  # Render plots based on reactive expression
  filtered_data_react <- reactive({
    new_data <- filtered_data(load_data_count(), input$percentile_variance, input$min_non_zero)
    return(new_data)
  })
  
  output$download_filter_count <- downloadHandler(
    filename = function() {
      paste("filtered_count_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      filtered_data <- filtered_data_react()[[2]]
      write.csv(filtered_data, file, row.names = FALSE)})
  
  output$Filtered_Summary_table <- renderTable({
    filtered_summary_table(filtered_data_react()[[1]], filtered_data_react()[[2]])
    })
  output$Diagnostic_plot <- renderPlot({
    diagnostic_scatter_plots(filtered_data_react()[[1]], filtered_data_react()[[2]])
  })
  output$Heatmap_plot <- renderPlot({
    clustered_heatmap(filtered_data_react()[[2]])
  })
  
  
  output$PCA_plot <- renderPlot({
    PCA_plot(filtered_data_react()[[2]], input$First_PC, input$Second_PC)
  })
  
  # Render objects for Sample tab
  output$Diff_eq_table <- renderDataTable({
    diff_eq(load_data_de(), column_summary(load_data()), input$cell_stage)
  })
  output$volcano_plot <- renderPlot({
    volcano_plot(diff_eq(load_data_de(), column_summary(load_data()), input$cell_stage), input$x_axis, input$y_axis, input$padj_color, input$color_base, input$color_highlight)
  })
  
  # Render objects for FGSEA tab
  output$barplot_fgsea <- renderPlot({
    req(load_fgsea_data()) 
    barplot_fgsea(load_fgsea_data(), input$filter_by_adjp_3)}, height = 1500)
  renderFGSEADataTable <- function() {
    req(load_fgsea_data())
    filter_fgsea_res(load_fgsea_data(), input$filter_by_adjp_1, input$NES_type)}
  output$fgsea_table <- renderDataTable({renderFGSEADataTable()})
  output$download_fgsea_result <- downloadHandler(
    filename = function() {
      paste("filtered_fgsea_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      filtered_data <- renderFGSEADataTable()
      write.csv(filtered_data, file, row.names = FALSE)})
  output$Scatterplot_NES <- renderPlot({
      req(load_fgsea_data()) 
      scatter_fgsea(load_fgsea_data(), input$filter_by_adjp_2)})
  
}
# Run the application
shinyApp(ui = ui, server = server)