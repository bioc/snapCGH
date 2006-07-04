"IDProbes" <-
function (input, array = 1, naut = 22, 
    Y = FALSE, X = FALSE, status, values, pch, cex, col, chrominfo = chrominfo.Mb, 
    ylim = c(-2, 2), ylb = "Log2Ratio", chrom.to.plot = 1, xlim=c(0,NA)){

  genomePlot(input = input, array = array, naut = naut, 
    Y = Y, X = X, status, values, pch, cex, col, chrominfo = chrominfo, 
    ylim = ylim, ylb = ylb, chrom.to.plot = chrom.to.plot, xlim=xlim)

    ##MALists haven't been adjust with respect to which channel is the test.
  ##This is done here, but isn't neccessary for SegLists
  if(class(input) == "MAList"){
    if (is.null(input$design)) 
        stop("MA$design component is null")
    for(i in 1:length(input$design)){
      temp <- input$design[i]* input$M[,i]
      input$M[,i] <- temp
    }
  }
  else if(class(input) == "SegList"){} 
  else{
    stop("Class must be either MAList or SegList")
  }

  chrom = input[input$genes$Chr == chrom.to.plot, array]

  identify(x = chrom$genes$Position, y = log2ratios(chrom), labels = chrom$genes$Name, tolerance = 0.5)

}
