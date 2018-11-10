# server for pcmp app
# defines data flow for a pair of projection types, each shown
# in two view based on different choices of dimensions

basicServer <- function(sce) function(input, output, session) {
#SERVER <- function(input, output, session) {
  requireNamespace("limma")
#
# add colnames as a column in colData -- uses .cellid field silently
#
  sce$.cellid = colnames(sce)
#
# create a named data.frame combining the reducedDims with the colData
#
  rd = reducedDims(sce)
  nmeth = length(rd) # list of matrices of projected data
  methnames = names(rd)
  nrd = names(rd)
  ncomps = vapply(rd, ncol, numeric(1))
  stopifnot(all(ncomps == ncomps[1]))  # requires balanced representation 
    # of all projections
  ncomp <- ncomps[1]
  indf = data.frame(do.call(cbind, as.list(rd)))
  cn = paste0(nrd[1], 1:ncomp)
  if (length(nrd) > 1) {
      for (j in 2:length(nrd)) cn = c(cn, paste0(nrd[j], 1:ncomp))
      }
  colnames(indf) = cn
  indf <- as.data.frame(cbind(indf, colData(sce)))
  
#
# build the formulas needed for d3scatter
#
  fmlist = lapply(methnames, function(x) list())
  names(fmlist) = methnames
  for (i in 1:nmeth) {
   curtags = paste0(methnames[i], 1:ncomp)
   fmlist[[i]] = lapply(curtags, function(x) as.formula(c("~", x)))
   names(fmlist[[i]]) = curtags
  }
  
#
# build the shared data
#
  enhDf = reactive({
   indf$strat = colData(sce)[[input$pickedStrat]]
   indf$key = 1:nrow(indf)
   indf
   })  

  shared_dat <- SharedData$new(enhDf) #enhDf, key=~key)

#
# set up reactive download entities: table of limma results, table of selected cells with selection sequence number
#
    output$downloadData <- downloadHandler(
       filename = function() {
         paste('data-', Sys.Date(), '.csv', sep='')
       },
       content = function(con) {
         write.csv(.GlobalEnv$.pcmpTab, con)
       }
     )
    output$downloadData2 <- downloadHandler(
       filename = function() {
         paste('data-', Sys.Date(), '.csv', sep='')
       },
       content = function(con) {
         dat = .GlobalEnv$.pcmpSelCells
         nsel = length(dat)
         selind = rep(1:nsel,sapply(dat,length))
         ans = data.frame(group=selind, cellid=unlist(dat))
         write.csv(ans, con)
       }
     )

#
# for the 'about' tab, show the SCE in use and some metadata
#

    output$scedump = renderPrint({
        print(sce)
    })
    output$scedump2 = renderPrint({
        print(metadata(sce)[c("note", "origin")])
    })


#
# produce the panels
#
  output$scatter1 <- renderD3scatter({
    methx = paste0(input$meth1, input$topx)
    methy = paste0(input$meth1, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~strat, width = "100%")
  })
  output$scatter2 <- renderD3scatter({
    methx = paste0(input$meth2, input$topx)
    methy = paste0(input$meth2, input$topy)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~strat, width = "100%")
  })
  output$scatter3 <- renderD3scatter({
    methx = paste0(input$meth1, input$botx)
    methy = paste0(input$meth1, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth1]][[methx]], 
            fmlist[[input$meth1]][[methy]], ~strat, width = "100%")
  })
  output$scatter4 <- renderD3scatter({
    methx = paste0(input$meth2, input$botx)
    methy = paste0(input$meth2, input$boty)
    d3scatter(shared_dat, fmlist[[input$meth2]][[methx]], 
            fmlist[[input$meth2]][[methy]], ~strat, width = "100%")
  })
#  output$try3d <- renderScatterplotThree({
#    colors = palette(rainbow(30))[ as.numeric(
#                      factor(colData(sce)[[input$pickedStrat]])) ]
#    print(head(shared_dat))
#    scatterplot3js(PC1, PC2, PC3, crosstalk=shared_dat, brush=TRUE,
#       color=colors)
#    })

#
# collect the information on selections so far
#
 output$accum = renderPlot({
 ans = list(cells = .GlobalEnv$.pcmpSelCells, limmaTab=.GlobalEnv$.pcmpTab)
 tmp = new("PcmpSels", cellSets=ans$cells, geneTable=ans$limmaTab)
 replay(sce, tmp, input$meth1, input$botx, input$boty) 
 })
    
#
# very rudimentary approach to acquiring a signature of a selected group of cells
# use limma on log-transformed counts comparing selected to non-selected
# could do something to balance sample sizes ...
#

output$summary <- DT::renderDataTable({
    df <- shared_dat$data(withSelection = TRUE) %>%
      filter(selected_ | is.na(selected_)) %>%
      mutate(selected_ = NULL)
    sel=rep(0, ncol(sce))
    names(sel) = colnames(sce)
    sel[df$.cellid] = 1
    mm = stats::model.matrix(~sel, data=data.frame(sel=sel))
   showNotification(paste("starting table processing", date()), id="limnote")
    X = log(assay(sce)+1)
    f1 = lmFit(X, mm)
    ef1 = eBayes(f1)
    options(digits=3)

    tt = topTable(ef1, 2, n=20)
    if (!(".pcmpSelNum" %in% ls(.GlobalEnv, all.names=TRUE))) assign(".pcmpSelNum", 1, .GlobalEnv)
      else assign(".pcmpSelNum", .GlobalEnv$.pcmpSelNum + 1, .GlobalEnv)
    if (!(".pcmpSelCells" %in% ls(.GlobalEnv, all.names=TRUE))) assign(".pcmpSelCells", list(df$.cellid), .GlobalEnv)
      else assign(".pcmpSelCells", c(.GlobalEnv$.pcmpSelCells, list(df$.cellid)), .GlobalEnv)
    tt = cbind(tt, selnum=.GlobalEnv$.pcmpSelNum[1])
    if (!(".pcmpTab" %in% ls(.GlobalEnv, all.names=TRUE))) assign(".pcmpTab", tt, .GlobalEnv)
      else assign(".pcmpTab", rbind(.GlobalEnv$.pcmpTab, tt), .GlobalEnv)
    ans = DT::formatRound(DT::datatable(tt), 1:7, digits=3)
   removeNotification(id="limnote")
    ans
  })

#
# prepare stop button
#

   observe({
                    if(input$btnSend > 0)
                        isolate({
                           stopApp(returnValue=0)
                        })  
           })  
  } # end server

library(SingleCellExperiment)
library(pcmp)
library(pcmpshin)
library(limma)
load("rscmou.rda")
basicServer(rscmou)

