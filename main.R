library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')
library('dplyr')



#'Load Data
load_data <- function(){
  
}

#' Sample Info exploration


column_summary <- function(){
  #summary of the type and mean (sd) or distinct values in each column, 
  #e.g.: Number of rows: X Number of columns: Y
  
}

plot_continous_var <- function(){
  #Tab with histograms, density plots, or violin plots of continuous variables.
  #If you want to make it fancy, allow the user to choose which column to plot and another column to group by!
}


#' Counts Matrix exploration 


#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  
  #read files
  count_data <- readr::read_tsv(counts_csv)
  metadata <- readr::read_csv(metafile_csv)
  
  #gene name
  gene <- count_data$gene
  
  #keep count columns that contain third parameter, select does not work with tibble
  count_df <- data.frame(count_data) %>% dplyr::select(contains(selected_times)) 
  count_matrix <- data.matrix(count_df)
  dimnames(count_matrix) <- list(gene) #make genes into row names; dimnames requires list
  colnames(count_matrix) <- colnames(count_df)
  
  #metadata needs only samplename, timepoint, filter to keep only rows in selected times
  metadata_df <- data.frame(metadata) %>% dplyr::select(c(samplename, timepoint)) %>% filter(timepoint %in% selected_times)
  
  #timepoint is a factor
  metadata_df$timepoint <- factor(metadata_df$timepoint, levels = selected_times)
  
  #store counts matrix and sample df in a SummarizedExperiments object
  se <- SummarizedExperiment(assays = list(counts = count_matrix), #subsetted counts matrix
                             colData = metadata_df) #store your sample dataframe as colData
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  
  ddsSE <- DESeqDataSet(se, design = design)
  
  #results from DESeq2 as df
  dds <- DESeq(ddsSE)
  dds_results <- results(dds)   #DESeqDataSet object updated ???
  dds_results_df <- as.data.frame(dds_results)
  
  #return list
  results_list <- list(Results = dds_results_df, DESeqDataSet = dds)
  
  return(results_list)
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

