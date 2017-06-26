#' Basic methods for class demand
#' This document lists a series of basic methods for the class demand

#' @param object Object of class demand
#' @return returns summary information about the demand object
#' @rdname printDeMAND
#' @examples
#' data(inputExample)
#' dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork)
#' printDeMAND(dobj)
#' @export
printDeMAND <- function(x) {
        cat("\nAn object of class demand\n\n")
        cat("Slot exp:\n")
        cat("\tExpression data of", nrow(x@exp), "features by", ncol(x@exp), "samples\n\n")
        cat("Slot anno:\n")
        cat("\tAnnotation of", nrow(x@anno), "probes and", length(unique(x@anno[,2])), "genes\n\n")
        cat("Slot network:\n")
        cat("\tInteractome with", length(unique(x@network)), "nodes,", nrow(x@anno), "edges\n\n")
        cat("Slot moa (head):\n")
        if (length(x@moa)>0) print(head(x@moa))
        else cat("\tEmpty\n\n")
        }

