library(tidyverse)
library('SummarizedExperiment')
library('DESeq2')
#library('biomaRt')
#library('testthat')
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
                          filter("Column Name" %in% colnames(filtered_data_matrix))
  #colnames(metadata)["Column Name"] <- "Sample"
  # 
  # #tibble of column metadata
  # metadata <- as_tibble(metadata)
  # # Specify reference levels for Timepoint and Cell Type
  # metadata$Timepoint <- factor(metadata$Timepoint, levels = factor_level)
  # 
  # #store counts matrix and sample df in a SummarizedExperiments object
  # se <- SummarizedExperiment(assays = list(counts = filtered_data_matrix), #subsetted counts matrix
  #                            colData = metadata) #store your sample dataframe as colData
  # ddsSE <- DESeqDataSet(se, design = ~Timepoint)
  # 
  # #results from DESeq2 as df
  # dds <- DESeq(ddsSE)
  # dds_results <- results(dds)
  # dds_results <- as.data.frame(dds_results)
  # dds_results <- tibble::rownames_to_column(dds_results, "Genes")

  
  
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
  
  return(metadata)
}
#' Tab with content similar to that described in [Assignment 7] (R shiny App)
volcano_plot <-function(dataf, x_name, y_name, slider, color1, color2) {
    
    slider <- 1 * 10**slider #calculate slider to fit chart
    
    plot_v <- ggplot(data = dataf) + #load data
      geom_point(aes(x = x_name, y = -log10(y_name), col = y_name < slider)) #+ #get data columns and color based on if statement 
     # scale_colour_manual(name = paste(y_name, '>', slider), values = setNames(c(color2,color1),c(T, F))) + #set colors and name
     # theme(legend.position="bottom") #change legend position
    
    return(plot_v)
  }

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  
  #make row names to column "gene"
  deseq2_res <- tibble::rownames_to_column(deseq2_res, "genes") 
  
  #convert to tibble
  deseq2_res <- as_tibble(deseq2_res)
  #add column volc_plot_status based on padj threshold
  deseq2_res <- deseq2_res %>% mutate(volc_plot_status = case_when(padj >= padj_threshold ~ 'NS',
                                                                   padj < padj_threshold & log2FoldChange > 0 ~ 'UP',
                                                                   padj < padj_threshold & log2FoldChange < 0 ~ 'DOWN'))
  
  return(deseq2_res)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  
  gg_hist <- labeled_results %>% 
    ggplot(aes(x=pvalue)) + 
    geom_histogram(color="black", fill="lightblue",  binwidth = .02) +
    ggtitle("Histogram of raw pvalues obtained from DE Analysis (vPO vs VAd)")
  
  return(gg_hist)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  
  #filter by padj threshold
  labeled_results <- labeled_results %>% filter(padj < padj_threshold)
  
  #plot histograam
  gg_hist <- labeled_results %>% 
    ggplot(aes(x=log2FoldChange)) + 
    geom_histogram(color="black", fill="lightblue", binwidth = 0.2) +
    ggtitle("Histogram of log2FoldChange obtained from DE Genes (vPO vs VAd)")
  
  return(gg_hist)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  
  #getting top 10 most significant genes
  dds_results <- labeled_results %>% arrange(padj) #ascending order
  dds_results <- head(dds_results, num_genes) #get top 10 genes
  top_genes <-dds_results$genes
  
  #DESeq2 normalized counts by size factors
  dds <- estimateSizeFactors(dds_obj)
  norm_counts <- counts(dds, normalized=TRUE) #normalized counts
  norm_counts_df <- data.frame(norm_counts)
  
  #keep only necessary genes in norm_counts
  filtered_norm_counts <- norm_counts_df %>% filter(row.names(.) %in% top_genes)
  #add genes as a column
  filtered_norm_counts$genes <- row.names(filtered_norm_counts)
  #pivot longer
  data <- filtered_norm_counts %>% pivot_longer(!genes, names_to = "samplenames", values_to = "norm_counts")
  
  #plot scatterplot
  scatterplot <- data  %>% #final tibble
    ggplot(aes(x=genes,y=log10(norm_counts), color = samplenames)) + #gene vs norm counts
    geom_point() + #scatterplot
    ggtitle("Plot of log10(normalized counts) for top 10 DE genes") +  #add title
    theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))
  
  return(scatterplot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  
  #remove na rows 
  labeled_results <- labeled_results %>% drop_na(volc_plot_status)
  
  #plot
  scatterplot <- labeled_results  %>% #final tibble
    ggplot(aes(x=log2FoldChange, y= -log10(padj), color = volc_plot_status)) + #log2FoldChange vs padj
    geom_point() + #scatterplot
    ggtitle("Volcano Plot of DESeq2 Differential Expresison Results") #add title
  
  return(scatterplot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  
  #first load in the id2gene.txt appropriately
  id <- read_delim(id2gene_path, col_names = FALSE, delim = "\t")
  colnames(id) <- c('genes', 'symbols')
  
  #add a new column in your labeled results that matches IDs to symbols
  labeled_results <- merge(x=labeled_results, y=id, by = "genes", all.x=TRUE)
  
  #log2FC values in descending order and select 2. columns
  labeled_results <- labeled_results %>% arrange(desc(log2FoldChange)) %>% #descending order
    drop_na(log2FoldChange) %>%
    dplyr::select(symbols, log2FoldChange)
  
  #generate a named vector of symbols and log2FC values 
  labeled_results_vector <- deframe(labeled_results)
  
  return(labeled_results_vector)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  
  #Run fgsea using a ranked list of descending log2FC against the C2 canonical pathways gene set
  c2_pathways <- fgsea::gmtPathways(gmt_file_path)
  
  fgseaRes <- fgsea(c2_pathways, 
                    rnk_list,
                    minSize  = min_size,
                    maxSize  = max_size)
  
  fgseaRes <- fgseaRes %>% as_tibble()
  
  return(fgseaRes)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  
  #arrange by NES
  fgsea_results <- fgsea_results %>% arrange(NES)
  
  #top X rows
  fgsea_results_top <- head(fgsea_results, num_paths)
  #bottom X rows
  fgsea_results_bottom <- tail(fgsea_results, num_paths)
  
  #rowbind top and bottom rows
  fgsea_results_10 <- rbind(fgsea_results_top, fgsea_results_bottom)
  
  #select necessary columns 
  fgsea_results_10 <- fgsea_results_10 %>% dplyr::select(pathway,NES)
  
  #bar chart
  stacked_bar <- fgsea_results_10 %>% 
    ggplot(aes(x = reorder(pathway, NES), y = NES)) +
    geom_col(aes(fill = NES > 0)) +
    coord_flip() +
    ggtitle("fgsea results for Hallmark MSigDB genes") +
    xlab("Normalized Enrichment Score (NES)")
  
  
  return(stacked_bar)
}

