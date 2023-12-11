library(tidyverse)
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('fgsea')
library(dplyr)
library('hrbrthemes') #theme_ipsum
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(patchwork)

#'Load Data
load_data <- function(file){
  df_txt <- read_delim(file, delim = "\t")
  
  return(df_txt)
}

#' Sample Info exploration
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
  
  # Timepoint
  data_num_updated <- data_num_updated %>% mutate("Timepoint" = case_when(
    str_starts(`Column Name`, "ex") ~  str_sub(data_num_updated$`Column Name`, 3,-3),
    (str_sub(`Column Name`, 2, 2) == "P") ~ str_sub(data_num_updated$`Column Name`, 2,3),
    (str_sub(`Column Name`, 2, 2) == "D") ~ str_sub(data_num_updated$`Column Name`, 2,4),
    (str_sub(`Column Name`, 2, 3) == "Ad") ~ "Ad",
    TRUE ~ NA_character_))
  
  # Replicate
  data_num_updated['Replicate'] <- str_sub(data_num_updated$'Column Name', -1) 
  
  
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

plot_continous_var <- function(data, column_name){
  
  #subset data
  data <- data %>% select(starts_with(column_name)) 
  
  # format data
  plot_df  <- pivot_longer(data, cols = everything()) #keep columns that start with
  # plot_df$name <- factor(plot_df$name) #make name column into factors to plot
  # 
  # # violin plot
  # v_plt <- ggplot(plot_df, aes(x = -log10(value), color = name, fill = name)) + 
  #           geom_density(adjust=3, alpha=.4) +
  #           theme_ipsum() +
  #           facet_wrap(~name) +
  #           theme(legend.position="bottom")
  
  return(plot_df)
}

#' Counts Matrix exploration 
#' Input controls that filter out genes based on their statistical properties:
#' Slider to include genes with at least X percentile of variance
#' Slider to include genes with at least X samples that are non-zero
counts_matrix_exploration <- function(data, filter_var, filter_non_zero){
  
  #'Include genes with at least X percentile of variance

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
  filtered_data_on_var <- data[data$Variance >= quantile(data$Variance, filter_var), ] #keep rows when Variance is greater than certain percentile 
  filtered_data_on_var_non_zero <- filtered_data_on_var[filtered_data_on_var$Number_of_Non_Zeros >= filter_non_zero, ] #count number of zeros per row 
  #filtered_data_on_var_non_zero$Passed_Filter = 'True'
  
  return(list(data, filtered_data_on_var_non_zero)) 
}

summary_table <- function(data, filtered_data){
  
  #' Tab with text or a table summarizing the effect of the filtering, including:
  #' number of samples, total number of genes, number and % of genes passing current filter, number and % of genes not passing current filter
  num_sample <- ncol(data) -2-3-2
  total_num_genes <- nrow(data)
  num_genes_pass <- nrow(filtered_data)
  perc_genes_pass <- round(num_genes_pass/total_num_genes, 4)*100
  num_genes_fail <- total_num_genes-num_genes_pass
  perc_genes_fail <- round(num_genes_fail/total_num_genes, 4)*100
  
  table_sumarize <- matrix(c(num_sample, total_num_genes, num_genes_pass, paste0(perc_genes_pass, '%'), num_genes_fail, paste0(perc_genes_fail, '%')))
  rownames(table_sumarize) <- c("Number of Samples", 'Total Number of Genes', 'Number of Genes that Pass Filter', 'Percent of Genes that Pass Filter', 'Number of Genes that Fail Filter', 'Percent of Genes that Fail Filter')
  
  return(table_sumarize)
  
}

diagnostic_scatter_plots <- function(data, filtered_data){
  #' Tab with diagnostic scatter plots, where genes passing filters are marked in a darker color, and genes filtered out are lighter:
  
  #' median count vs variance (consider log scale for plot)
  median_v_var <- ggplot() + 
          geom_point(data = data, aes(x=-log10(Median_Count), y=-log10(Variance), color = "Fail_Filter"), alpha = 0.5, clip = FALSE) +
          geom_point(data = filtered_data, aes(x=-log10(Median_Count), y=-log10(Variance), color = "Pass_Filter"), alpha = 0.5, clip = FALSE) +
          scale_color_manual(values = c("Fail_Filter" = "lightblue", "Pass_Filter" = "darkblue")) +
          ggtitle('Median count vs variance') +
          coord_cartesian(clip = 'off')
   #' median count vs number of zeros
  median_v_zero <- ggplot() + 
            geom_point(data = data, aes(x=-log10(Median_Count), y=Number_of_Zeros, color = "Fail_Filter"), alpha = 0.5) +
            geom_point(data = filtered_data, aes(x=-log10(Median_Count), y=Number_of_Zeros, color = "Pass_Filter"), alpha = 0.5) +
            scale_color_manual(values = c("Fail_Filter" = "lightblue", "Pass_Filter" = "darkblue")) +
            ggtitle('Median count vs number of zeros') +
            xlim(0,500)
  return(median_v_var | median_v_zero)
}

clustered_heatmap <- function(filtered_data) {
#' Tab with a clustered heatmap of counts remaining after filtering
#' consider enabling log-transforming counts for visualization
  
  #Subset data
  filtered_data_numeric <- filtered_data[, c(1, 4:39)]
  
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
          cellheight = 10)
  
  return(hetmap)
}

PCA_plot <- function(data, PC_x, PC_y) {
  #' Tab with a scatter plot of principal component analysis projections. You may either:
  
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

  #' #' allow the user to select which principal components to plot in a scatter plot (e.g. PC1 vs PC2)
  scatterplot_PCA <- pca_results_subset %>%
                      ggplot(aes(x = .data[[first_column_name]], y = .data[[second_column_name]])) +
                      geom_point() +
                      ggtitle('Principal Component Scatterplot') +
                      xlab(paste0('PC',PC_x, ": ", var_exp_tib[PC_x,], '% Variance')) + #' be sure to include the % variance explained by each component in the plot labels
                      ylab(paste0('PC', PC_y, ": ", var_exp_tib[PC_y,], '% Variance'))
  return(scatterplot_PCA)
}

#' Differential Expression
#' Results of a differential expression analysis in CSV format.
#' If results are already made available, you may use those
#' Otherwise perform a differential expression analysis using DESeq2

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
   metadata <- metadata %>% select("Column Name", "Cell Stage", "Cell Type", "Timepoint", "Replicate") %>% 
                        subset(Replicate != "N/A") %>% 
                        rename_with(~ "Sample", `Column Name`) %>%
                          filter(Sample %in% colnames(filtered_data_matrix))
  # 
  # #tibble of column metadata
  metadata <- as_tibble(metadata)
  # # Specify reference levels for Timepoint and Cell Type
  metadata$Timepoint <- factor(metadata$Timepoint, levels = factor_level)

  #store counts matrix and sample df in a SummarizedExperiments object
  se <- SummarizedExperiment(assays = list(counts = filtered_data_matrix), #subsetted counts matrix
                             colData = metadata) #store your sample dataframe as colData
  ddsSE <- DESeqDataSet(se, design = ~Timepoint)

  #results from DESeq2 as df
  dds <- DESeq(ddsSE)
  dds_results <- results(dds)
  dds_results <- as.data.frame(dds_results)
  dds_results <- tibble::rownames_to_column(dds_results, "Genes")

  # dds <- DESeqDataSetFromMatrix(countData= filtered_data_matrix,
  #                               colData=metadata,
  #                               design= ~Timepoint)
  # # compute normalization factors
  # dds <- estimateSizeFactors(dds)
  #
  # # extract the normalized counts
  # dds <- as_tibble(counts(dds, normalized=TRUE)) %>%
  #   mutate(gene = count_data$gene) %>%
  #   relocate(gene) # relocate changes column order
  # # Perform differential expression analysis
  # dds <- DESeq(dds)
  # results <- results(dds)
  
  return(dds_results)
}
#' Tab with content similar to that described in [Assignment 7] (R shiny App)
volcano_plot <-function(dataf, x_name, y_name, slider, color1, color2) {
    
    if (y_name == "padj" || x_name == "padj"){
    dataf <- dataf[!is.na(dataf$padj), ]
    }
  
    slider <- 1 * 10**slider #calculate slider to fit chart
    
    plot_v <- ggplot(data = dataf) + #load data
      geom_point(aes(x = .data[[x_name]], y = -log10(.data[[y_name]]), col = .data[[y_name]] < slider)) + #get data columns and color based on if statement 
      scale_colour_manual(name = paste(y_name, '>', slider), values = setNames(c(color2,color1),c(T, F))) + #set colors and name
      theme(legend.position="bottom") #change legend position
    
    return(plot_v)
}

fgsea_lst <- function(labeled_results, gmt_file) {
  
  # Connect to the appropriate BioMart database
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = "mmusculus_gene_ensembl")
  
  # Grab gene names
  ensembl_ids <- labeled_results$Genes
  
  # Retrieve gene symbols using biomaRt
  ensembl_to_symbol <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                              filters = "ensembl_gene_id", 
                              values = ensembl_ids, 
                              mart = ensembl)
   
  # Merge the original data frame with the gene symbols
  labeled_results_with_symbols <- merge(labeled_results, ensembl_to_symbol, by.x = "Genes", by.y = "ensembl_gene_id", all.x = TRUE)
  
  #log2FC values in descending order and select 2 columns
  labeled_results_with_symbols <- labeled_results_with_symbols %>% arrange(desc(log2FoldChange)) %>% #descending order
                                  drop_na(log2FoldChange) %>%
                                  drop_na(external_gene_name) %>%
                                  dplyr::select(external_gene_name, log2FoldChange)

  #generate a named vector of symbols and log2FC values
  labeled_results_vector <- deframe(labeled_results_with_symbols)
  
  #Run fgsea using a ranked list of descending log2FC against the M2 canonical pathways gene set
  cp_pathways <- fgsea::gmtPathways(gmt_file)
  
  fgseaRes <- fgsea(cp_pathways, 
                    labeled_results_vector,
                    minSize  = 15,
                    maxSize  = 500)
  
  fgseaRes <- fgseaRes %>% as_tibble()
  
  return(fgseaRes)
}

#' Barplot of fgsea NES for top pathways selected by slider
barplot_fgsea <- function(file, num_paths){
  
  fgsea_results <- read_csv(file)
  
  #keep rows that have a padj or lower and arrange by NES
  fgsea_results_top <- fgsea_results %>% filter(padj < num_paths) %>% arrange(desc(NES))
  
  #select necessary columns 
  fgsea_results_10 <- fgsea_results_top %>% dplyr::select(pathway,NES)
  
  #bar chart
   stacked_bar <- fgsea_results_10 %>% 
     ggplot(aes(x = reorder(pathway, NES), y = NES)) +
     geom_col(aes(fill = NES > 0)) +
     coord_flip() +
     ggtitle("fgsea results for CP MSigDB genes") +
     xlab("Normalized Enrichment Score (NES)")
  
  return(stacked_bar)
}
#' Filtered data table displaying the FGSEA results
filter_fgsea_res <- function(file, num_paths, sign){
  
  fgsea_results <- read_csv(file)
  
  #keep rows that have a padj or lower and arrange by NES
  fgsea_results <- fgsea_results %>% filter(padj < num_paths)
  
  if (sign == "positive"){
    fgsea_results <- fgsea_results %>% filter(NES > 0)
  }
  else if (sign == "negative"){
    fgsea_results <- fgsea_results %>% filter(NES < 0)
  }
  else {
    fgsea_results <- fgsea_results
  }
  
  return(fgsea_results)
}
#' Scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis, with gene sets below threshold in grey color
scatter_fgsea <- function(file, num_paths){
  
  fgsea_results <- read_csv(file)
  
  scatter <- fgsea_results %>% 
    ggplot(aes(x = NES, y = -log10(padj))) +
    geom_point(aes(color = ifelse(-log10(padj) < -log10(num_paths), "below_threshold", "above_threshold"))) +
    scale_color_manual(values = c("below_threshold" = "grey", "above_threshold" = "pink"), name = paste("padj threshold =", num_paths)) +
    ggtitle("Scatterplot of NES vs -log10(padj)") +
    xlab("Normalized Enrichment Score (NES)")
  
  return(scatter)
}
