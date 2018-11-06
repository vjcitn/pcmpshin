
basicServer <- function(sce) function(input, output, session) {
  sce$.cellid = colnames(sce)
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
  
  fmlist = lapply(methnames, function(x) list())
  names(fmlist) = methnames
  for (i in 1:nmeth) {
   curtags = paste0(methnames[i], 1:ncomp)
   fmlist[[i]] = lapply(curtags, function(x) as.formula(c("~", x)))
   names(fmlist[[i]]) = curtags
  }
  
  #PCtags = paste0("PC", 1:ncomp)
  #UMtags = paste0("UM", 1:ncomp)
  #TStags = paste0("TS", 1:ncomp)
  #PCfmlas = lapply(PCtags, function(x) as.formula(c("~", x)))
  #UMfmlas = lapply(UMtags, function(x) as.formula(c("~", x)))
  #TSfmlas = lapply(TStags, function(x) as.formula(c("~", x)))
  #names(PCfmlas) = PCtags
  #names(UMfmlas) = UMtags
  #names(TSfmlas) = TStags
  #fmlist = list(PC=PCfmlas, UM=UMfmlas, TS=TSfmlas)

  enhDf = reactive({
   indf$strat = colData(sce)[[input$pickedStrat]]
   indf$key = 1:nrow(indf)
   indf
   })  

  shared_dat <- SharedData$new(enhDf) #enhDf, key=~key)

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
 output$accum = renderPlot({
 ans = list(cells = .GlobalEnv$.pcmpSelCells, limmaTab=.GlobalEnv$.pcmpTab)
 tmp = new("PcmpSels", cellSets=ans$cells, geneTable=ans$limmaTab)
 replay(sce, tmp, input$meth1, input$botx, input$boty) 
 })
    

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
    f1 = limma::lmFit(X, mm)
    ef1 = limma::eBayes(f1)
  print(paste0("finish lmFit", date()))
    options(digits=3)

    tt = limma::topTable(ef1, 2, n=20)
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
load("mouse3k.rda")
basicServer(mouse3k)

