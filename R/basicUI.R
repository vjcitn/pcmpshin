pkgVersion = function() as.character(read.dcf(system.file("DESCRIPTION", package="pcmp"))[,"Version"])

 uiMaker = function(sce) {

  rd = reducedDims(sce)
  nrd = names(rd)
  ncomps = vapply(rd, ncol, numeric(1))
  stopifnot(all(ncomps == ncomps[1]))
  ncomp <- ncomps[1]
  discv = discreteColdVars(sce)

  fluidPage(
   sidebarPanel(width=2,
     helpText(sprintf("pcmp %s: crosstalk-based interactive graphics \
for dimension reduction in single-cell transcriptomics. \ 
See the 'about' tab for more information.", pkgVersion())),
     selectInput("pickedStrat", "stratby", discv, discv[1]),
     selectInput("meth1", "method left", nrd, nrd[1]),
     selectInput("meth2", "method right", nrd, nrd[2]),
     numericInput("topx", "top x", 1, min=1, max=ncomp-1, step=1),
     numericInput("topy", "top y", 2, min=2, max=ncomp, step=1),
     numericInput("botx", "bot x", 2, min=1, max=ncomp-1, step=1),
     numericInput("boty", "bot y", 3, min=2, max=ncomp, step=1),
     actionButton("btnSend", "Stop app")
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
     helpText("pcmpApp can be demonstrated with the object pcmp::sce300xx, an extract from the Allen Brain Atlas RNA-seq data on anterior cingulate cortex (ACC) and primary visual cortex (VIS) brain regions.  Strata were formed using donor (3 levels) and region (2 levels) and 300 cells were sampled at random in each stratum.  The murLung3k app at shinyapps.io uses an extract from the Tabula Muris project data focused on a collection of cells from the mouse lung; the data are available in the github repo vjcitn/pcmpshin in data/mouse3k.rda.")
   )
  )
  )
  ) 
}
#
