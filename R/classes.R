setClass("LargeDataObject")

setClass("SegList",
representation("list")
)

setIs("SegList", "LargeDataObject")

dim.SegList <- function(x) if(is.null(x[[1]])) c(0,0) else dim(as.matrix(x[[1]]))
length.SegList <- function(x) prod(dim(x))
#dimnames.SegList <- function(x) dimnames(x[[1]])

#allows the subsetting of the SegList object.  
assign("[.SegList",
function(object, i, j, ...) {
  if (nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if (missing(i))
		if (missing(j))
			return (object)
		else {
			object$state <- object$state[,j,drop=FALSE]
			object$rpred <- object$rpred[,j,drop=FALSE]
			object$prob <- object$prob[,j,drop=FALSE]
			object$M.predicted <- object$M.predicted[,j,drop=FALSE]
			object$dispersion <- object$dispersion[,j,drop=FALSE]
                        object$variance <- object$variance[,j,drop=FALSE]
			object$M.observed <- object$M.observed[,j,drop=FALSE]
                        object$num.states <- object$num.states[,j,drop=FALSE]
                      }
	else
		if (missing(j)) {
			object$state <- object$state[i,,drop=FALSE]
			object$rpred <- object$rpred[i,,drop=FALSE]
			object$prob <- object$prob[i,,drop=FALSE]
			object$M.predicted <- object$M.predicted[i,,drop=FALSE]
			object$dispersion <- object$dispersion[i,,drop=FALSE]
                        object$variance <- object$variance[i,,drop=FALSE]
			object$M.observed <- object$M.observed[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
#                        object$num.states <- object$num.states[i,,drop=FALSE]
		} else {
			object$state <- object$state[i,j,drop=FALSE]
			object$rpred <- object$rpred[i,j,drop=FALSE]
			object$prob <- object$prob[i,j,drop=FALSE]
			object$M.predicted <- object$M.predicted[i,j,drop=FALSE]
			object$dispersion <- object$dispersion[i,j,drop=FALSE]
                        object$variance <- object$variance[i,j,drop=FALSE]
			object$M.observed <- object$M.observed[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
		}
	object
})


cbind.SegList <- function(..., deparse.level=1){
  objects <- list(...)
  nobjects <- length(objects)
  out <- objects[[1]]
  if (nobjects > 1) {
    for (i in 2:nobjects){
      out$M.predicted <- cbind(out$M.predicted, objects[[i]]$M.predicted)
      out$dispersion <- cbind(out$dispersion, objects[[i]]$dispersion)
      out$variance <- cbind(out$variance, objects[[i]]$variance)
      out$M.observed <- cbind(out$M.observed, objects[[i]]$M.observed)
      out$state <- cbind(out$state, objects[[i]]$state)
      out$num.states <- cbind(out$num.states, objects[[i]]$num.states)
    }
    out
  }
}

rbind.SegList <- function(..., deparse.level=1){
  objects <- list(...)
  nobjects <- length(objects)
  out <- objects[[1]]
  if(nobjects > 1){
    for(i in 2:nobjects){
      out$M.predicted <- rbind(out$M.predicted, objects[[i]]$M.predicted)
      out$dispersion <- rbind(out$dispersion, objects[[i]]$dispersion)
      out$variance <- rbind(out$variance, objects[[i]]$variance)
      out$M.observed <- rbind(out$M.observed, objects[[i]]$M.observed)
      out$state <- rbind(out$state, objects[[i]]$state)
      out$num.states <- rbind(out$num.states, objects[[i]]$num.states)
      out$genes <- rbind(out$genes, objects[[i]]$genes)
    }
    out
  }
}

dimnames.SegList <- function(segList){
  dimnames(segList$M.observed)
}

#Taken from Limma
#I can't get it to use the print.LargeDataObject method for a SegList
#when a NAMESPACE is used.  Hopefully this will be sorted in the future.
print.SegList <- function(x, ...) {
	cat("An object of class \"",class(x),"\"\n",sep="")
	for (what in names(x)) {
		y <- x[[what]]
		cat("$",what,"\n",sep="")
		printHead(y)
		cat("\n")
	}
	for (what in setdiff(slotNames(x),".Data")) {
		y <- slot(x,what)
		if(length(y) > 0) {
			cat("@",what,"\n",sep="")
			printHead(y)
			cat("\n")
		}
	}
}
