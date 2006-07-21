zoomGenome <- function(..., array = 1, colors = NULL, chrominfo = chrominfo.Mb){

  objects <- list(...)
  nobjects <- length(objects)
  
    #### check they are SegLists ####
  for (i in 1:nobjects){
    if(class(objects[[i]]) != "SegList"){
      stop("This function currently only takes objects of type SegList")}
  }
  
  dist <- vector()
  dist[1] <-  0
  for(i in 1:nrow(chrominfo)){
    dist[i+1] <- sum(chrominfo$length[1:i])
  }

  close.screen(all.screens = TRUE)
  split.screen(c(2,1))
  screen(1)
  par(mar = c(2, 4, 4, 2))
  plotSegmentedGenome(..., array = array, chrominfo = chrominfo, colors = colors)
  
  doIt = TRUE
  i = 0
  
  while(doIt){
    loc <- locator(1)
                                        # print(loc$x)
    close.screen(all.screens = TRUE)
    split.screen(c(2,1))

### THIS CODE IS AS UGLY AS SIN AND IS TEMPORARY.  I'LL MAKE IT NICER SOON. #####
    
  if(loc$x >= dist[1] & loc$x < dist[2]){
    i=1
  }
  else if(loc$x >= dist[2] & loc$x < dist[3]){
    i=2
  }
  else if(loc$x >= dist[3] & loc$x < dist[4]){
    i=3
  }
  else if(loc$x >= dist[4] & loc$x < dist[5]){
    i=4
  }
  else if(loc$x >= dist[5] & loc$x < dist[6]){
    i = 5
  }
  else if(loc$x >= dist[6] & loc$x < dist[7]){
    i = 6
  }
  else if(loc$x >= dist[7] & loc$x < dist[8]){
    i = 7
  }
  else if(loc$x >= dist[8] & loc$x < dist[9]){
    i = 8
  }
  else if(loc$x >= dist[9] & loc$x < dist[10]){ 
    i = 9
  }
  else if(loc$x >= dist[10] & loc$x < dist[11]){
    i = 10
  }
  else if(loc$x >= dist[11] & loc$x < dist[12]){
    i = 11
  }
  else if(loc$x >= dist[12] & loc$x < dist[13]){
    i = 12
  }
  else if(loc$x >= dist[13] & loc$x < dist[14]){
    i = 13
  }
  else if(loc$x >= dist[14] & loc$x < dist[15]){
    i = 14
  }
  else if(loc$x >= dist[15] & loc$x < dist[16]){
    i = 15
  }
  else if(loc$x >= dist[16] & loc$x < dist[17]){
    i = 16
  }
  else if(loc$x >= dist[17] & loc$x < dist[18]){
    i = 17
  }
  else if(loc$x >= dist[18] & loc$x < dist[19]){
    i = 18
  }
  else if(loc$x >= dist[19] & loc$x < dist[20]){
    i = 19
  }
  else if(loc$x >= dist[20] & loc$x < dist[21]){
    i = 20
  }
  else if(loc$x >= dist[21] & loc$x < dist[22]){
    i = 21
  }
  else if(loc$x >= dist[22] & loc$x < dist[23]){
    i = 22
  }
  else{
    doIt = FALSE
  }
    screen(2)
    plotSegmentedGenome(..., array = array, chrominfo = chrominfo, colors = colors, chrom.to.plot = i)
    screen(1)
    par(mar = c(2, 4, 4, 2))
    plotSegmentedGenome(..., array = array, chrominfo = chrominfo, colors = colors)
  }
}

zoomChromosome <- function(..., array = 1, chrom.to.plot, colors = NULL, chrominfo = chrominfo.Mb, ylim = c(-2,2)){

  objects <- list(...)
  close.screen(all.screens = TRUE)
  split.screen(c(2,1))
  par(mar = c(2, 4, 4, 2))
  screen(1)

  if(is.null(objects[[1]]$M.predicted)){
      genomePlot(objects[[1]], chrom.to.plot = chrom.to.plot, array = array, chrominfo = chrominfo, ylim = ylim)
      doIt = TRUE
      while(doIt){
        loc <- locator(2)
        close.screen(all.screens = TRUE)
        split.screen(c(2,1))
        screen(2)
        genomePlot(objects[[1]], array = array, xlim = c(min(loc$x), max(loc$x)), chrom.to.plot = chrom.to.plot, chrominfo = chrominfo, ylim = ylim)
        screen(1)
        par(mar = c(2, 4, 4, 2))
        genomePlot(objects[[1]], array = array, chrom.to.plot = chrom.to.plot, chrominfo = chrominfo, ylim = ylim)
        rect(xleft = 0, xright = min(loc$x), ybottom = (min(ylim)-2), ytop = (max(ylim)+2), density = 10)
        rect(xleft = max(loc$x), xright = max(objects[[1]]$genes$Position), ybottom = (min(ylim)-2), ytop = (max(ylim)+2), density = 10)
      }
    }
  else{
    plotSegmentedGenome(..., chrom.to.plot = chrom.to.plot, array = array, chrominfo = chrominfo, colors = colors, ylim = ylim)
    doIt = TRUE
    while(doIt){
      loc <- locator(2)
      close.screen(all.screens = TRUE)
      split.screen(c(2,1))
      screen(2)
      plotSegmentedGenome(..., array = array, xlim = c(min(loc$x), max(loc$x)), chrom.to.plot = chrom.to.plot, chrominfo = chrominfo, colors = colors, ylim = ylim)
      screen(1)
      par(mar = c(2, 4, 4, 2))
      plotSegmentedGenome(..., array = array, chrom.to.plot = chrom.to.plot, chrominfo = chrominfo, colors = colors, ylim = ylim)
      rect(xleft = 0, xright = min(loc$x), ybottom = (min(ylim)-2), ytop = (max(ylim)+2), density = 10)
      rect(xleft = max(loc$x), xright = max(objects[[1]]$genes$Position), ybottom = (min(ylim)-2), ytop = (max(ylim)+2), density = 10)
    }
  }
}
