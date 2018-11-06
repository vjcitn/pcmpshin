runner = function(sce=pcmp::sce300xx) {
stores = c(".pcmpTab", ".pcmpSelNum", ".pcmpSelCells")
 sapply(stores, function(x) if(x %in% ls(.GlobalEnv, all.names=TRUE))
    rm(list=x, envir=.GlobalEnv))

 theServer = basicServer(sce)
 theUI = uiMaker(sce)
 shiny::runApp(list(ui=theUI, server=theServer))
 pcmp::PcmpSels(list(cells = .GlobalEnv$.pcmpSelCells, limmaTab = .GlobalEnv$.pcmpTab))
}
