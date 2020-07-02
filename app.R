library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(Seurat)
library(ggplot2)
library(plotly)
library(RColorBrewer)
# library(ggpubr)


## Only run examples in interactive R sessions
# if (interactive()) {
##################################################################################################################
##################################################################################################################
################################################## INTERFACE #####################################################
##################################################################################################################
##################################################################################################################
ui <- fluidPage(
  sidebarLayout(fluid = TRUE,
    
    #########################
    #### Parameter entry ####
    #########################
    sidebarPanel(width = 3,
      # Load Metadata
      fileInput(inputId = "metadata",
                label = "Metadata + 2D embedings (coord_x, coord_y)",
                multiple = FALSE),
      
      # Load Expression data
      fileInput(inputId = "data",
                label = "Expression matrix",
                multiple = FALSE),
      
      # Select genes to use as markers
      selectizeInput("gene_ls", 'Select marker genes:',
                     selected = NULL,
                     choices = NULL,
                     options = list(create = TRUE),
                     multiple = TRUE),
      actionButton(inputId = "apply_markers", label = "Update markers"),
      
      # Enter value to group by
      selectizeInput("groupby", 'Coloring feature:',
                     selected = NULL,
                     choices = NULL,
                     options = list(create = TRUE),
                     # multiple = TRUE,
                     multiple = FALSE),
      actionButton(inputId = "apply_groupby", label = "Update coloring"),
      
      # Which labels to add to interactive output
      selectizeInput("interactive_labels", 'Interactive labels:',
                     selected = NULL,
                     choices = NULL,
                     options = list(create = TRUE),
                     multiple = TRUE),
      
      # Slider to determine point size of UMAP plots
      sliderInput("size", "Dot size:",
                  min = 0,
                  max = 10,
                  value = 4,
                  step = 0.1),
      
      # Which labels to add to interactive output
      selectizeInput("filter_var", 'Interactive labels:',
                     selected = NULL,
                     choices = NULL,
                     options = list(create = TRUE),
                     multiple = TRUE),
      
      
    ),
    
    ############################
    #### Plot visualization ####
    ############################
    mainPanel(
      tabsetPanel(
        tabPanel(title = 'UMAP Plot',
                 plotlyOutput("dimPlot", height = "800")),
        tabPanel(title = 'Feature Plots',
                 plotlyOutput("FeaturePlot", height = "800")),
        tabPanel(title = 'Violin Plots',
                 plotlyOutput("ViolinPlot", height = "800"))
      )
    )
  )
)

##################################################################################################################
##################################################################################################################
################################################ SERVER ##########################################################
##################################################################################################################
##################################################################################################################

server <- function(input, output, session) {
  # Setting maximum file size to 4GB
  options(shiny.maxRequestSize=8000*1024^2)
  
  ##########################################
  ##### Defining environment variables #####
  ##########################################
  # x <- NULL; makeReactiveBinding("x")
  # marker_ds <- NULL
  # se_obj <- NULL
  
  ##############################
  ##### Update input fields ####
  ##############################
  observe({
    
    nrow_dictionary <<- list()
    nrow_dictionary[[1]] <<- 1
    nrow_dictionary[[2]] <<- 1
    nrow_dictionary[[3]] <<- 2
    nrow_dictionary[[4]] <<- 2
    nrow_dictionary[[5]] <<- 2
    nrow_dictionary[[6]] <<- 2
    nrow_dictionary[[7]] <<- 3
    nrow_dictionary[[8]] <<- 3
    nrow_dictionary[[9]] <<- 3
    nrow_dictionary[[10]] <<- 4
    
    # Load data
    # Read marker list
    file1 <- input$metadata
    if (is.null(file1) ) { return() }
    tmp1_ds <- readRDS(file1$datapath)
    metadata_df <<- tmp1_ds

    file2 <- input$data
    if( is.null( file2 ) ) { return() }
    tmp2_ds <- readRDS(file2$datapath)
    expr_mtrx <<- tmp2_ds
    
    # Update Cluster labels
    # updateTextInput(session,
    #                 "labels_vec",
    #                 value = paste0(levels(tmp1_ds),collapse = ','))
    
    # Update marker selection
    updateSelectizeInput(session,
                         inputId = "gene_ls",
                         choices = c("", rownames(expr_mtrx)),
                         selected = rownames(expr_mtrx)[1])
    
    # Update groupby selection
    updateSelectizeInput(session,
                         inputId = "groupby",
                         choices = c("", 
                                     colnames(metadata_df)),
                         selected = colnames(metadata_df)[1])
    
    # Update interactive_labels selection
    updateSelectizeInput(session,
                         inputId = "interactive_labels",
                         choices = c("", 
                                     colnames(metadata_df)),
                         selected = colnames(metadata_df)[1])
    
  })
  
  ########################################
  ######## Setting reactive events #######
  ########################################
  gene_list <- eventReactive(input$apply_markers, {
    input$gene_ls
  })
  
  groupby_var <- eventReactive(input$apply_groupby, {
    input$groupby
  })
  
  interactive_labels <- eventReactive(input$apply_labs, {
    input$interactive_labels
  })
  
  ########################################
  ######### Plot visualization ###########
  ########################################

  # UMAP clusters
  output$dimPlot <- renderPlotly({
    # tmp_dim <- DimPlot(se_obj, reduction = "tsne", label = F, pt.size = 0.5, group.by = groupby_var())
    # dim_plot <- HoverLocator(plot = tmp_dim, information = FetchData(se_obj, vars = input$interactive_labels))
    
    labs_dim <- lapply(input$interactive_labels, function(i){
      paste(sprintf('\n%s: ',i), metadata_df[,i], sep = '')
    }) %>% purrr::pmap_chr(., paste)
    
    dim_plot <- plot_ly(x = metadata_df[,"coord_x"],
                        y = metadata_df[,"coord_y"],
                        color = metadata_df[, groupby_var()],
                        # Hover text:
                        text = labs_dim,
                        marker = list(size = as.numeric(input$size))
    )
    
    return(dim_plot)
  })
  
  # Feature plot and Violin plot are dependent on apply_markers button to be pressed
  
  # Feature plot
  output$FeaturePlot <- renderPlotly({
    # tmp_feat <- FeaturePlot(se_obj, features = gene_list())
    labs_feat <- lapply(input$interactive_labels, function(i){
      paste(sprintf('\n%s: ',i), metadata_df[,i], sep = '')
    }) %>% purrr::pmap_chr(., paste)
    
    ## Plot all genes
    plt_ls <- lapply(gene_list(), function(gene) {
      
      if (gene %in% rownames(expr_mtrx)) {
        ## Title
        # font style
        f <- list(
          family = "Courier New, monospace",
          size = 18,
          color = "black")
        
        titl <- list(
          text = sprintf("Gene: %s", gene),
          font = f,
          xref = "paper",
          yref = "paper",
          yanchor = "bottom",
          xanchor = "center",
          align = "center",
          x = 0.5,
          y = 1,
          showarrow = FALSE,
          face = "bold"
        )
        
        ## Plot
        feat_plt <- plotly::plot_ly(x = metadata_df[, "coord_x"],
                                    y = metadata_df[, "coord_y"],
                                    color = expr_mtrx[gene, ],
                                    colors = "Blues",
                                    marker = list(size = as.numeric(input$size)),
                                    # Hover text:
                                    text = labs_feat) %>%
          plotly::layout(annotations = titl,
                         showlegend = FALSE)
      } else {
        warning(sprintf("Gene %s not found in the expression matrix.", gene))
        return(NULL)
      }
      
    })
    
    # Set the number of rows
    if (length(gene_list()) > length(nrow_dictionary)) {
      n_row <- 4
    } else {
      n_row <- nrow_dictionary[[length(gene_list())]]
    }
    
    # Arrange all the plots
    feat_arr <- plotly::subplot(plt_ls,
                                nrows = n_row,
                                shareX = FALSE,
                                shareY = FALSE,
                                margin = 0.05)
    
    return(feat_arr)
  })
  
  # Violin plots
  output$ViolinPlot <- renderPlotly({
    labs_ls <- lapply(input$interactive_labels, function(i){
      paste(sprintf("%s: ", i), metadata_df[, i], sep = "")
    }) %>% purrr::pmap_chr(., paste)
    
    ## Plot all genes
    vln_ls <- lapply(gene_list(), function(gene) {
      ## Title
      # font style
      f <- list(
        family = "Courier New, monospace",
        size = 18,
        color = "black")
      
      titl <- list(
        text = sprintf("Gene: %s", gene),
        font = f,
        xref = "paper",
        yref = "paper",
        yanchor = "bottom",
        xanchor = "center",
        align = "center",
        x = 0.5,
        y = 1,
        showarrow = FALSE,
        face = "bold")
      
      ## X axis text
      f1 <- list(
        family = "Courier New, monospace",
        size = 10,
        color = "black"
      )
      
      a <- list(
        title = groupby_var(),
        tickfont = f1,
        showticklabels = TRUE,
        tickangle = 90
      )
      
      ## Plots
      vln_plt <- plotly::plot_ly(
        x = stringr::str_wrap(
          string = metadata_df[, groupby_var()],
          width = 15,
          indent = 1, # let's add extra space from the margins
          exdent = 1  # let's add extra space from the margins
        ),
        y = expr_mtrx[gene, ],
        color = metadata_df[, groupby_var()],
        type = "violin",
        text = labs_ls,
        points = "all",
        jitter = 1,
        pointpos = 0,
        box = list(
          visible = FALSE
        ),
        meanline = list(
          visible = TRUE
        )
      ) %>%
        plotly::layout(annotations = titl,
                       xaxis = a,
                       showlegend = FALSE)
      
      return(vln_plt)
      })
    
    # Set the number of rows
    if (length(gene_list()) > length(nrow_dictionary)) {
      n_row <- 4
    } else {
      n_row <- nrow_dictionary[[length(gene_list())]]
    }
    
    # Arrange all the plots
    vln_arr <- plotly::subplot(vln_ls,
                                nrows = n_row,
                                shareX = FALSE,
                                shareY = FALSE,
                                margin = 0.075)
    
    
    return(vln_arr)
  })
  
}

shinyApp(ui, server)
# }


# p <- plot_ly(
#     # x = se_obj@meta.data[,groupby_var()],
#     # y = se_obj@assays$SCT@data[,gene_list()[1]],
#     # split = se_obj@meta.data[,groupby_var()],
#     x = se_obj@meta.data[,"maturation_stage"],
#     y = se_obj@assays$SCT@data['TRIM13',],
#     split = se_obj@meta.data[,"maturation_stage"],
#     type = 'violin',
#     points = 'all',
#     jitter = 5,
#     pointpos = 0,
#     box = list(
#       visible = F
#     ),
#     meanline = list(
#       visible = T
#     )
#   ) %>%
#   layout(
#     xaxis = list(
#       title = groupby_var()
#     ),
#     yaxis = list(
#       title = "Expression Level",
#       zeroline = F
#     )
#   )
