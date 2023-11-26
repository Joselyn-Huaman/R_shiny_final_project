#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) # you might need to install this
library(htmltools)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library('SummarizedExperiment')
library('DESeq2')
#library('biomaRt')
#library('testthat')
#library('fgsea')
library(hrbrthemes) #theme_ipsum
#library(pheatmap)
library(Cairo)


# Increase the maximum upload size (in bytes)
options(shiny.maxRequestSize = 30 * 1024^2)  # Set to 30 MB

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
                                            actionButton("run_cont_plt_button", "Violin Plot", icon = icon("play"), class = "btn-block"),
                                          plotOutput("Continous_plot"), width = "100%", height =  "1500px")
                             )
                            )
                           )
                  ),
                  tabPanel("Counts"),
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
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Load_Data: Sample information matrix in CSV format
  load_data <- reactive({
    
    #require an input file
    req(input$Count_File)
    # read the file
    input_file <- read_delim(file = input$Count_File$datapath, delim = "\t") #results of the file upload are nested
    #colnames(input_file)[1] <- 'gene' #change column name
    
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
  
  
  #' Volcano plot
  #'
  #' @param dataf The loaded data frame.
  #' @param x_name The column name to plot on the x-axis
  #' @param y_name The column name to plot on the y-axis
  #' @param slider A negative integer value representing the magnitude of
  #' p-adjusted values to color. Most of our data will be between -1 and -300.
  #' @param color1 One of the colors for the points.
  #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
  #'
  #' @return A ggplot object of a volcano plot
  #' @details I bet you're tired of these plots by now. Me too, don't worry.
  #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
  #' Write a normal volcano plot using geom_point, and integrate all the above 
  #' values into it as shown in the example app. The testing script will treat 
  #' this as a normal function.
  #' 
  #' !!sym() may be required to access column names in ggplot aes().
  #'
  #' #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
  #' volcano_plot <-
  #'   function(dataf, x_name, y_name, slider, color1, color2) {
  #'     
  #'     slider <- 1 * 10**slider #calculate slider to fit chart
  #'     
  #'     plot_v <- ggplot(data = dataf) + #load data
  #'       geom_point(aes(x = .data[[x_name]], y = -log10(.data[[y_name]]), col = .data[[y_name]] < slider)) + #get data columns and color based on if statement 
  #'       scale_colour_manual(name = paste(y_name, '>', slider), values = setNames(c(color2,color1),c(T, F))) + #set colors and name
  #'       theme(legend.position="bottom") #change legend position
  #'     
  #'     return(plot_v)
  #'   }
  #' 
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
    
    #plot_df$name <- factor(plot_df$name) #make name column into factors to plot
    
    # histogram plot
    h_plt <- ggplot(plot_df, aes(x = Sample, y = -log10(Count), fill = Timepoint)) + 
      geom_violin(adjust=3, alpha=.4, binwidth = 50) +
      theme_ipsum() +
      geom_boxplot(width=0.1) +     
      theme(legend.position="bottom")
    
    return(h_plt)
  } 
  #' Draw and filter table
  #'
  #' @param dataf Data frame loaded by load_data()
  #' @param slider Negative number, typically from the slider input.
  #'
  #' @return Data frame filtered to p-adjusted values that are less than 
  #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
  #' displayed.
  #' @details Same as above, this function is a standard R function. Tests will 
  #' evaluate it normally. Not only does this function filter the data frame to 
  #' rows that are above the slider magnitude, it should also change the format 
  #' of the p-value columns to display more digits. This is so that it looks 
  #' better when displayed on the web page. I would suggest the function 
  #' `formatC()`
  #'
  #' @examples draw_table(deseq_df, -210)
  #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
  #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
  #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
  draw_table <- function(dataf, slider) {
    
    # Ensure that the slider is not NULL or NA
    if (is.null(slider) || is.na(slider)) { # Handle the case where slider is not valid
      return(dataf)
    }
    
    slider <- 1 * 10**slider #calculate slider to fit chart
    
    dataf_filter <- dataf %>% filter(padj < slider) %>% mutate(padj = formatC(padj, digits=6))  #filter based on slider value
    
    return(dataf_filter)
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
  # Render objects foor Sample tab
  output$Summary_table <- renderTable({column_summary(load_data())})
  
  output$DataTable_table <- renderDataTable({load_data()})
  
  renderContinoutVarPlot <- function() {
    req(input$subset_data)
    isolate({ 
      plot_continous_var(load_data(), input$subset_data)
    })
  }
  
  output$Continous_plot <- renderPlot({renderContinoutVarPlot()}, height = 1500)
  
  observeEvent(input$run_cont_plt_button, {
    output$Continous_plot <- renderPlot({renderContinoutVarPlot()}, height = 1500)})
  
  # # Update data and re-render on actionButton click
  # observeEvent(input$run_button, {
  #   output$volcano <- renderPlot({renderVolcanoPlot()})
  #   output$table <- renderTable({renderDataTable()})
  # })
}
# Run the application
shinyApp(ui = ui, server = server)
