# This function creates a list with all graphical paramenters
# used by CNAnorm. It makes an object then loaded by CNAnorm. It not used in
# the package, but I leave it in case I want to change the object.


makeDefaultGraphParamteres <- function() {
    ### genomePlot ###
    # colors
    colors <- list()
    colors$neutral.dot <- 'grey60'
    colors$gain.dot <- 'red'
    colors$loss.dot <- 'cornflowerblue'
    colors$over.dot <- 'darkgreen'
    colors$below.dot <- 'darkgreen'
    colors$segLine <- 'black'
    colors$smoothLine <- 'cyan'
    colors$yAxes.left <- 'blue'
    colors$yAxes.right <- 'black'
    colors$grid <- 'gray'
    colors$chrVertLines <- 'black'
    colors$closestPeak <- 'gray20'
    colors$ratioMedian <- 'gray20'
    colors$centrm <- 'gray60'

    # cex
    cex <- list()
    cex$neutral.dot <- .2
    cex$gain.dot <- .2
    cex$loss.dot <- .2
    cex$over.dot <- .2
    cex$below.dot <- .2
    cex$yAxes.left <- 1
    cex$yAxes.right <- 1

    #lwd
    lwd <- list()
    lwd$segLine <- 5
    lwd$smoothLine <- 5
    lwd$grid <- 1
    lwd$closestPeak <- 2
    lwd$ratioMedian <- 2
    lwd$chrVertLines <- 1
    lwd$centrm <- 1

    # lty
    lty <- list()
    lty$grid <- 1 # solid line
    lty$chrVertLines <- 1
    lty$closestPeak <- 1
    lty$ratioMedian <- 3
    lty$centrm <- 3

    # pch
    pch <- list()
    pch$neutral.dot <- 19
    pch$gain.dot <- 19
    pch$loss.dot <- 19
    pch$over.dot <- 2
    pch$below.dot <- 6

    # lables
    lab <- list()
    lab$xAxes.low <- "Genomic location"
    lab$yAxes.left <- "Estimated ploidy"
    lab$yAxes.right <- "Ratio centered on most common"


    # generic
    mar <- c(5, 4, 4, 5)



    genome <- list(colors = colors, cex = cex, pch = pch, lwd = lwd, 
        lty = lty, lab = lab, mar = mar)

        
    gPar <- list(genome = genome)
    return (gPar)
}
