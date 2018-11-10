pkgVersion = function() as.character(read.dcf(system.file("DESCRIPTION", package="pcmp"))[,"Version"])

 uiMaker = function(sce) {

  rd = reducedDims(sce)
  nrd = names(rd)
  ncomps = vapply(rd, ncol, numeric(1))
  stopifnot(all(ncomps == ncomps[1]))
  ncomp <- ncomps[1]
  discv = discreteColdVars(sce)

  fluidPage(
   sidebarPanel(width=3,
    fluidRow(
     column(12,
       helpText(h3(sprintf("pcmp %s: crosstalk-based interactive graphics \
for dimension reduction in single-cell transcriptomics. \ 
See the 'about' tab for more information.", pkgVersion())))
       )
      ),
    fluidRow( 
     column(12,
       selectInput("pickedStrat", "stratby", discv, discv[1])
      )
     ),
    fluidRow( 
     column(12,
       helpText(h4("projection methods:"))
      )
     ),
    fluidRow(
      column(6,
       selectInput("meth1", "left", nrd, nrd[1])),
      column(6,
       selectInput("meth2", "right", nrd, nrd[2]))
     ),
    fluidRow( 
     column(12,
       helpText(h4("dimensions to use:"))
      )
     ),
    fluidRow(
      column(6,
       numericInput("topx", "top x", 1, min=1, max=ncomp-1, step=1)),
      column(6,
       numericInput("topy", "top y", 2, min=2, max=ncomp, step=1))
     ),
    fluidRow(
      column(6,
       numericInput("botx", "bot x", 2, min=1, max=ncomp-1, step=1)),
      column(6,
       numericInput("boty", "bot y", 3, min=2, max=ncomp, step=1))
     ),
    fluidRow(
      column(12,
       helpText(h4("downloads:")))
     ),
    fluidRow(
      column(6,
       downloadButton("downloadData", "DE genes")),
      column(6,
       downloadButton("downloadData2", "cellSets"))
     ),
    fluidRow(
      column(12,
       helpText(h4("to conclude:")))
     ),
    fluidRow(
      column(12,
       actionButton("btnSend", "Stop app"))
     )
     ),
   mainPanel(
    tabsetPanel(
    tabPanel("scatter",
     fluidRow(
       column(6, d3scatterOutput("scatter1")),
       column(6, d3scatterOutput("scatter2"))
       ),
      fluidRow(
       column(6, d3scatterOutput("scatter3")),
       column(6, d3scatterOutput("scatter4"))
       )
     ), # end panel
#    tabPanel("topr3d",
#       scatterplotThreeOutput("try3d")
#     ),
    tabPanel("selTable",
     DT::dataTableOutput("summary")
     ),
    tabPanel("accum",
       helpText("Method and dimensions taken from bottom left panel"),
       plotOutput("accum")
        ),
    tabPanel("about",
     helpText(h3("pcmp demonstrates crosstalk-based interactive graphics for surveying different dimension reduction procedures for data in SingleCellExperiment containers.  The reducedDims component must be populated with several reductions, each including at least 4 dimensions.   Different methods are used in the left and right columns, and different projection components can are used in the top and bottom rows, as selected using the method/top/bot controls.  The 'stratby' button will recolor points according to discrete covariates in the colData of the input object.")),
     helpText(h3("current input data structure:")),
     verbatimTextOutput("scedump"),
     helpText(h3("metadata strings:")),
     verbatimTextOutput("scedump2"),
     helpText(h3("pcmpApp can be demonstrated with the object pcmp::sce300xx, an extract from the Allen Brain Atlas RNA-seq data on anterior cingulate cortex (ACC) and primary visual cortex (VIS) brain regions.  Strata were formed using donor (3 levels) and region (2 levels) and 300 cells were sampled at random in each stratum.  The murLung3k app at vjcitn.shinyapps.io uses an extract from the Tabula Muris project data focused on a collection of cells from the mouse lung; the uncorrected data are available in the github repo vjcitn/pcmpshin in data/mouse3kf.rda."))
   )
  )
  )
  ) 
# end UI
} # end UImaker

library(SingleCellExperiment)
library(pcmp)
library(pcmpshin)
load("rscmou.rda")
uiMaker(rscmou)

