#library(pacman)
#p_load(shiny, visNetwork, ggiraph)

library(shiny)
library(visNetwork)
library(ggiraph)

load("./data/precomputedData_ui.RData")

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Rare Disease MultiOme explorer"),
  
  # Main panel has 4 tab panels
  mainPanel(
    tabsetPanel(
      #############
      ## Panel 1 ##
      #############
      tabPanel('Disease Differential Modularity',
               fluidRow(
                 column(width=3,
                        selectInput('diseases_DifMod', 'Select rare disease group', 
                                    diseases_below_1000,
                                    selected = "Ciliopathy"),
                        
                        radioButtons("pType", "Information to show:",
                                     c("Network" = "network",
                                       "Table" = "table",
                                       "Barplot" = "barplot"
                                     )),
                        

                 ), #column 1: plot options
                 column(width=8,
                         tags$p('Differential modularity of the selected disease group showing underlying causal mechanisms played by gene connectivity on different molecular organizations'),
                         tags$p('Click dropdown button on the left to change disease groups, and the radio button to change the information shown (Network, Table, or Bar chart)'),
                        

                        fluidRow(
                          conditionalPanel('input.pType=="network"', imageOutput("DiffmodImage")),
                          conditionalPanel('input.pType=="table"', tableOutput("LCCtable")),
                          conditionalPanel('input.pType=="barplot"', girafeOutput("LCCbarplot"))
                        )
                        # plotOutput
                 )
               ))
      ,
      #######
      # Tab 2
      ########
      tabPanel('Network Landscape',
               fluidRow(
                 column(width=3,
                        selectInput('Network_landscape', 'Select network', 
                                    network_name,
                                    selected = "Human Phenotype"))
                        
                 ), #column 1: plot options
                 column(width=8,
                         tags$p('Network-disease landscape (node2vec + t-SNE) shown for a particular embedded network. The white dots represent causal genes for the corresponding disease group'),
                        
                        
                        fluidRow(
                          imageOutput("LandscapeImage")
                        )
                        # plotOutput
                 )
               )
      ,
      ##########
      # Tab 3
      ##########
      tabPanel('Network-disease inspection', 
               fluidRow(
                 column(width=3, 
                        selectInput('disease', 'Select rare disease group', diseases, selected = "Ciliopathy"),
                        selectInput('network', 'Select network to display', network_name, selected = "Human Phenotype"),

                        
                 ), # /Column 1: plot options
                 column(width=8,
                        tags$p("Inspect gene connectivity for a particular disease on the selected network."),
                       visNetworkOutput("networkplot", width = "100%"),
                 tags$p('Zoom in a particular area by scrolling. Gene labels will appear and can be selected upon click.'),
                 tags$p('Select multiple genes by clicking and holding Cmd (Mac) or Ctrl (Windows/Linux). Selected genes can be enriched via:'),
                 tags$a(href="https://maayanlab.cloud/Enrichr/", "Enrichr"),
                 ),#, # /column 2: plot 
                 column(width=1,
                        tags$p('Selected genes:'),
                       # DTOutput("selectednodes", width = "15%")
                        tableOutput("selectednodes"),
                       
                   )# /column 3: download buttons
                 
               ) 
               
      ) # end tabPanel 3
    

      
    ) # /tabsetPanel
  ) # /mainPanel
))
