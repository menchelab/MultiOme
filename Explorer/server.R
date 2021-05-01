
# load required packages
#library(pacman)
#p_load(ggplot2, visNetwork, igraph, RColorBrewer, shiny, cowplot, tidyverse)

library(tidyverse)
library(visNetwork)
library(RColorBrewer)
library(shiny)
library(cowplot)
library(igraph)


# load required data

load("./data/precomputedData_server.RData")

# Define server logic 
shinyServer(function(input, output) {
  
  
  values <- reactiveValues()
  #############
  ## Tab 1 ##
  #############
  
  #network plot
  output$networkplot <- renderVisNetwork({


    disease_of_interest <- input$disease
    graph_of_interest <- network_info$network[network_info$subtype==input$network]

    disease_graph <- subgraphs[[disease_of_interest]][[graph_of_interest]]
     
    nodes_df <- tibble(id = V(disease_graph)$name, label = V(disease_graph)$name)
    edges_df <- igraph::as_data_frame(disease_graph)

    visNetwork(nodes = nodes_df, edges = edges_df, width = "200%") %>%
      visPhysics(stabilization=F) %>%
      visEdges(smooth=F, color="grey", width=0.3)  %>%
      visNodes(color = list(background="white", highlight="magenta", border="black")) %>%
      visIgraphLayout() %>%
      visOptions(nodesIdSelection = list(enabled=T, useLabels=T,
                                         style = 'width: 200px; height: 26px;
                                         background: #f8f8f8;
                                         color: black;
                                         border:none;
                                         outline:none;'),
                 highlightNearest = list(enabled =TRUE, degree = 1, hover = T))%>%
      visInteraction(multiselect = TRUE) %>%
      #    visGroups(groupname = "upSIG", color = "red") %>%
      #  visGroups(groupname = "downSIG", color = "blue") %>%
      visExport(type = "png", name = "Network",
                float = "left", label = "Save network (png)", background = "white", style= "") %>%
      visEvents(select = "function(nodes) {
            Shiny.onInputChange('current_node_id', nodes.nodes);
            ;}")



  })
  
  
  
  output$selectednodes <- renderTable({
 # tibble(nodes = input$networkplot_selected)
    tibble(nodes = input$current_node_id)
  })

  
  
  # image2 sends pre-rendered images
  output$DiffmodImage <- renderImage({
    if (is.null(input$diseases_DifMod))
      return(NULL)
    
    dis <- input$diseases_DifMod
    
    return(list(
      src = paste0("Figs/module/", dis ,".png"),
      contentType = "image/png",
      width = 800,
      alt =  dis
    ))
  }, deleteFile = FALSE)
  
  
  
  # tSNE plot landscape
  output$LandscapeImage <- renderImage({
    if (is.null(input$Network_landscape))
      return(NULL)
    
    net <- input$Network_landscape
    netid <- network_info$network[network_info$subtype==net]
    
    return(list(
      src = paste0("Figs/landscape/", netid, ".png"),
      contentType = "image/png",
      width = 800,
      alt =  net
    ))
  }, deleteFile = FALSE)
  

  # LCC table
  output$LCCtable <- renderTable({
    processed_result_df %>%
      dplyr::filter(name == input$diseases_DifMod, LCC.signif != "none") %>%
      arrange(desc(LCC.zscore)) %>% 
      dplyr::select(subtype, source, N_in_graph, LCC.size, LCC.mean, LCC.sd, LCC.zscore, LCC.signif)
  })
  
  output$LCCbarplot <- renderGirafe({
    girafe(ggobj = processed_result_df %>%
      dplyr::filter(name == input$diseases_DifMod, LCC.signif != "none") %>%
      arrange(desc(LCC.zscore)) %>% 
      ggplot(., aes(y=subtype, x=LCC.zscore, tooltip = paste0("p-value: ",correctedPval), data_id = subtype)) + 
        geom_col_interactive() + 
      facet_grid(main_type~., scales = "free", space = "free") +
      ylab("Network") + xlab("LCC z-score") + 
      theme_cowplot() 
    )
    })
  
  
})
