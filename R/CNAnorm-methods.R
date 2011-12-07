#### "helper" functions that will be called by function (called by methods) or
# other functions not exported in the environment, and inaccessible to user


.discreteNorm <- function(object, normBy = object) {
    toNorm = object
    if (length(valid.peaks(normBy)) == 0) {
        normBy <- validation(normBy)
    } 
    ploidy <- valid.ploidy(normBy)
    peaks <- valid.peaks(normBy)

    if (length(ploidy) != length(peaks)) {
        stop("Peaks and Ploidy vector have different length")
    } else if (length(ploidy) == 1){
        shift = 0
        scale = 1/peaks[1]
    } else {
        MQR <- getMQR(ploidy, peaks)
        interpDiploidXloc <- c(0,2) * MQR[1] + MQR[2] # y = mx + q
        shift = interpDiploidXloc[1]
        scale = 1 / ( interpDiploidXloc[2] - interpDiploidXloc[1] )
    }
    
    toNorm@DerivData@ratio.n <- (ratio(toNorm) - shift) * scale
    toNorm@Res@validated.ratioMedian <- (normBy@Res@ratioMedian - shift) * scale
    toNorm@Res@validated.closestPeak <- (normBy@Res@closestPeak - shift) * scale
    if ( length(ratio.s(toNorm)) > 0 ) {
        toNorm@DerivData@ratio.s.n <- (ratio.s(toNorm) - shift) * scale
    } 
    if ( length(segMean(toNorm)) > 0 ) {
        toNorm@DerivData@segMean.n <- ( 2^segMean(toNorm) - shift ) * scale
    }

    
    toNorm@Res@shifting = shift
    toNorm@Res@scaling = scale
    return (toNorm)
}

.validation <- function (object, peaks = sugg.peaks(object), ploidy = sugg.ploidy(object)) {
    if (length(peaks) != length(ploidy)) {
        stop ("Peaks and ploidy must have same length\n" )
    } else if (any(round(ploidy) != ploidy)) {
        stop ("Ploidy must be integer\n")
    } else {
        object@Res@validated.peaks = peaks
        object@Res@validated.ploidy = ploidy
        MQR <- getMQR(ploidy, peaks)
        interpDiploidXloc <- c(0,2) * MQR[1] + MQR[2]
        shift = interpDiploidXloc[1]
        scale = 1 / ( interpDiploidXloc[2] - interpDiploidXloc[1] )
        if (length(interpDiploidXloc) == 1) {
            content = "Unknown"
        } else if (length( interpDiploidXloc) == 2){
            content = as.character(
                signif(100*(1-interpDiploidXloc[1]/interpDiploidXloc[2]), 4))
        } else {
            stop("interpDiploidXloc has more than 2 elements!")
        }

        object@Res@validated.tumContent <- content
        object@Res@validated.ratioMedian <- 
            (object@Res@ratioMedian - shift) * scale
        object@Res@validated.closestPeak <- 
            (object@Res@closestPeak - shift) * scale
        object@Res@validated.R           <- signif(MQR[3], 4) 
    }

    return (object)
}


.plotGenome <- function (object, maxRatio = 8, minRatio = -1, 
    superimpose = character(0),  supLineColor = character(0), 
    supLineCex = character(0), numHorLables = 10, ...) {
    # superimpose can be one of the following: "DNAcopy" or "smooth"    

    if (length(ratio.n(object)) == 0){
        stop("Object does not contains normalized ratio\n")
    }

    defaultColors <- c("black", "cyan", 'grey60')
    defaultsLineCex <- c(.5, .5)
    par(mar=c(5, 4, 4, 5))
    
    bolT <- ratio.n(object) < maxRatio & ratio.n(object) > minRatio
    
    where <- which(bolT)
    outOfPlotPlus <- which(ratio.n(object) > maxRatio)
    outOfPlotMinus <- which(ratio.n(object) < minRatio)
    xAxes <- 1:length(ratio.n(object))
    xAxesoutOfPlotPlus <- pos(object)[outOfPlotPlus] 
    xAxesoutOfPlotMinus <- pos(object)[outOfPlotMinus] 
    # 1:length(outOfPlotPlus)
    numDots <- length(xAxes)

    yAxesoutOfPlotPlus <- rep(maxRatio, length(outOfPlotPlus))
    yAxesoutOfPlotMinus <- rep(minRatio, length(outOfPlotMinus))
    xAxisLab <- "Genomic position:"
    xRange <- range(xAxes[where], na.rm = TRUE)
    yRange <- range(ratio.n(object)[where], na.rm = TRUE)
    yRange[1] <- yRange[1] - 0.1
    yRange[2] <- yRange[2] + 0.1

    plot(xRange, yRange, type = 'n', xlab = "", xaxt="n", yaxt="n", ylab="", ...) 
    lines(xAxes[where], ratio.n(object)[where], type = 'p', cex = .5, pch=19, 
        col = defaultColors[3])
    if (length(outOfPlotPlus) > 0) {
        lines(xAxes[outOfPlotPlus], yAxesoutOfPlotPlus + 0.1, pch = 2, 
            cex = .5, type = 'p', col = 'red')
        lines(xAxes[outOfPlotMinus], yAxesoutOfPlotMinus - 0.1, pch = 6, 
            cex = .5, type = 'p', col = 'red')
    }
    # numDots <- length(where)
    xTickWidth <- numDots/(numHorLables + 1)
    breakPoints <- seq(from = 1, to = numDots, length.out = numHorLables + 1)
    middlePoints = round(breakPoints[1:length(breakPoints)-1] + 

    breakPoints[2:length(breakPoints)])/2
    lableNames <- chrs(object)
    tickNames <- lableNames[middlePoints]
    tickLocation = rep(NA, length(tickNames))
    for (i in 1: length(tickNames)) {
        tickLocation[i] = round(mean(which(chrs(object) %in% tickNames[i])) )
    }
    axis(1, at = tickLocation, labels = tickNames)

    yMax = floor(max(ratio.n(object)[where], na.rm = TRUE))
    axis(2, at = seq(0, yMax, by = .5), labels=2*(seq(0, yMax, by = .5)), 
        col.axis = 'blue')
    mtext("Estimated ploidy", side = 2, line = 3, col= 'blue')
    roundPloidy <- signif(object@Res@validated.closestPeak, digits=2)
    axis(4, at = seq(0, yMax, by = roundPloidy), 
        labels=signif(seq(0, yMax, by = roundPloidy)/roundPloidy, digits = 2))
    mtext("ratio centered on most common", side=4, line=3)

    abline(h = seq(0, yMax, by = .5), col = 'gray')
    
    
    # allChr <- Res$data$Chr[where]
    allChr <- chrs(object)
    chrChange <- c(1, which(allChr[1:length(allChr)-1] != 
        allChr[2:length(allChr)]), length(allChr))
    abline(v = chrChange, col = 'red')
    xAxisLab <- paste(xAxisLab, "Chromosomes") 
    
    abline(h = object@Res@validated.closestPeak, col = 'gray20', lwd = 2)
    abline(h = object@Res@validated.ratioMedian, col = 'gray20', lwd = 2, lty = 3 )

    # superimpose DNAcopy and/or smooth
    if (any(superimpose == 'smooth', na.rm = TRUE)){
        thisCol <- userOrDefault(supLineColor, defaultColors, which(superimpose == 'smooth')) 
        thisCex <- userOrDefault(supLineCex, defaultsLineCex, 
            which(superimpose == 'smooth'))
        lines(xAxes[where], ratio.s.n(object)[where], type = 'p', cex = thisCex, 
            pch=19, col = thisCol) #, type = 'l', lw = thisLw, col = thisCol)
    }
    if (any(superimpose == 'DNACopy', na.rm = TRUE)){
        if (length(segMean.n(object)) != length(object)) {
            object <- addDNACopy(object)
        } 
        thisCol <- userOrDefault(supLineColor, defaultColors, which(superimpose == 'DNACopy')) 
        thisCex <- userOrDefault(supLineCex, defaultsLineCex, 
            which(superimpose == 'DNACopy'))
        lines(xAxes[where], segMean.n(object)[where], type = 'p', cex = thisCex, 
            pch=19, col = thisCol) # type = 'l', lw = thisLw, col = thisCol)
    }
    if (!all(superimpose %in% c('smooth', 'DNACopy'))) {
        warning("'superimpose' must be either 'smooth' or 'DNACopy'. skipping...\n")
    }
}

.peakPloidy <- function(object, method = 'mixture', exclude = character(0),
    ploidyToTest = 12, sd = 5, dThresh = 0.01, n = 2048, adjust = 0.9, 
    # force.smooth = TRUE, reg=FALSE, ds=1.5, zero.cont=FALSE, ...) {
    force.smooth = TRUE, reg=FALSE, ds=1.5, zero.cont=FALSE, ...) {
   
    whatM <- checkMethod(method, c('mixture', 'density'))
    # check if we have smoothed signal...
#     if (length(ratio.s(object)) == 0 & force.smooth){
#         message("smoothing ratio...")
#         object <- addSmooth(object)
#     }
    
    chrName <- chrs(object)
    ratio <- ratio2use(object)
    isOutlier <- findOutliers(ratio, sd, dThresh, n = n, adjust = adjust)
    ok4density <- (! chrName %in% exclude) & (! isOutlier) & (! is.na(ratio))
    if (whatM == 'density'){
        object <- .guessPeaksAndPloidy(object, ok4density = ok4density, ...)
        object@Params@method = 'density'
    } else if (whatM == 'mixture') {
        options(show.error.messages = FALSE)
        sk <- try(suggest.k(ratio[ok4density], chrName[ok4density], np = seq(3, ploidyToTest, by = 1), 
            n = n, adjust = adjust, reg=reg, ds=ds, zero.cont=zero.cont, ...))
        if (class(sk) == "try-error"){
            message ("method 'mixture' failed with the following message:")
            message (sk[1])
            message ("falling back to method 'density'...")
            object <- .peakPloidy(object, method = 'density', exclude = exclude,
                ploidyToTest = ploidyToTest, sd = sd, dThresh = dThresh, n = n,
                adjust = adjust, reg = reg, ds = ds, zero.cont = zero.cont, ...)
            return(object)
        }
        options(show.error.messages = TRUE)
#         sk <- suggest.k(ratio[ok4density], chrName[ok4density], np = seq(3, ploidyToTest, by = 1), 
#             n = n, adjust = adjust, reg=reg, ds=ds, zero.cont=zero.cont, ...)   
        gn <- global.norm(ratio[ok4density], chrName[ok4density], k = sk[[1]][1], reg=reg,
	    ds=ds, zero.cont=zero.cont, ...)
        object@Params@method <- 'mixture'
        pl <- gn$pl[[1]]
        peaks <- gn$muest1[[1]]
        object@Res@suggested.ploidy <- pl
        object@Res@suggested.peaks <- peaks
        MQR <- getMQR(pl, peaks)
        interpDiploidXloc <- c(0,2) * MQR[1] + MQR[2]
        shift = interpDiploidXloc[1]
        scale = 1 / ( interpDiploidXloc[2] - interpDiploidXloc[1] )

        if (length(interpDiploidXloc) == 1) {
            content = "Unknown"
        } else if (length( interpDiploidXloc) == 2){
            content = as.character(
                signif(100*(1-interpDiploidXloc[1]/interpDiploidXloc[2]), 4))
        } else {
            stop("interpDiploidXloc has more than 2 elements!")
        }
        object@Res@suggested.tumContent <- content
        object@Res@suggested.R <- signif(MQR[3], 4)
        ratioMedian <- median(ratio[ok4density], na.rm = TRUE)
        object@Res@ratioMedian <- ratioMedian
        medianPeaksIndex <- which(abs(peaks-ratioMedian)  ==
            min(abs(peaks-ratioMedian), na.rm = TRUE))
        object@Res@closestPeak <- peaks[medianPeaksIndex]
        

    } else {
        stop('method must be either "mixture" or "density"\n')
    }
    
    object@Params@gp.excludeFromDensity <- exclude
    object@DerivData@ok4density         <- ok4density
    object@Params@density.n             <- n
    object@Params@density.adjust        <- adjust

    return(object) 


}

findOutliers <- function(ratio, sd, dThresh, ...){
    s <- median(ratio, na.rm = TRUE) + sd * sd(ratio, na.rm = TRUE)
    
    d <- density(ratio, na.rm = TRUE, ...)
    selected <- which(d$y >= dThresh*max(d$y))
    dT <- max(d$x[selected], ra.rm=T)

    isOutlier <- ratio > dT & ratio > s
    isOutlier[is.na(isOutlier)] <- FALSE
    return(isOutlier)
}


checkMethod <- function(method, avaiables){
    if (length(method) != 1 | !is.character(method)){
        stop ("method must be one element character vector\n")
    }
    where <- grep(method, avaiables, ignore.case = TRUE)
    if (length(where) == 0){
        stop (paste("method '", method, "' not recognized\n", sep = ''))
    } else if (length(where) > 1) {
        stop (paste("method '", method, "' matches multiple formal options\n", sep = ''))
    } else {
        return(avaiables[where])
    }
}




userOrDefault <- function(user, default, elNum) {
    if (elNum > length(user) | is.na(user[elNum])) {
        return (default[elNum])
    } else {
        return (user[elNum])
    }
}



smoothSignal <- function (signal, chrName, lambda = 7, ...) {
    ratio.s <- rep (NA, length(signal))
    uniqChrNames <- unique(chrName)
    for (n in 1:length(uniqChrNames)){
        isThisChr <- chrName %in% uniqChrNames[n]
        if (length(which(!is.na(signal[isThisChr]))) < 5) {
            ratio.s[isThisChr] <- rep(mean(signal[isThisChr], na.rm = TRUE), length(which(isThisChr)))     
        } else {
            ratio.s[isThisChr] <- smooth.it(signal[isThisChr], chrName[isThisChr], lambda = lambda, ...)
        }   
    } 
    return(ratio.s)
}

smooth.it <- function(a, ch, lambda = 7, ...){
    if(sum(is.na(a))> 0 | sum(is.infinite(a))> 0 | sum(is.nan(a))>0){
        a[is.infinite(a)] <- NA
        a[is.na(a)] <- NA
        a[a==0] <- 0.0001
        temp <- smoothseg(a,chrom=ch, lambda = lambda, impute.only=TRUE)
        y <- smoothseg(temp, chrom=ch, lambda = lambda, ...)
        return(y)
    } else {
        a[a==0] <- 0.0001
        y <- smoothseg(a, chrom=ch, lambda = lambda, ...)
        return(y)
    }
} # end of function

gcNormalize <- function (gcC, ratio, points2use = NULL, maxNumPoints = 10000) {
    # perform GC content normalization. Returns normalizing vector
    # gcC: vector with gc content (x)
    # ratio: vector with number of reads (y)
    # points2use: vector (TRUE/FALSE) with gcC and ratio to use. If NULL, all
    #   gcC and ratio will be used
    # maxNumPoints: maximal number of points to use to get the vector. Owing to
    #   computational time, and possible huge size of dataset, only maxNumPoints
    #   randomly selected are used to estimate the normalization curve
    # 
    # Returns a vector of the value of Loess Curve for all points in gcC. It
    #   keeps the same order as gcC. NA will be present where gcC is NA and
    #   where points2use is FALSE
    validX <- which(!is.na(gcC) & ratio > 0 & !is.na(ratio))
    # check if all points are good or only some of them
    if (is.null(points2use)) {
        validGc <- validX
    } else {
        validGc <- which(!is.na(gcC) & ratio > 0 & !is.na(ratio) & points2use)
    }   
    # if there are too many points, randomly select (maxNumPoints)
    if (length(validGc) > maxNumPoints ) { 
#         TandF <- c(rep(TRUE, maxNumPoints), rep(FALSE, length(validGc) - maxNumPoints))
#         randomNums <- runif(length(TandF))
#         o <- order(randomNums)
#         shuffle <- TandF[o]
        validGc <- sample(validGc, maxNumPoints)
    }   
    # make a data fram with relevant names (important!)
    gcDist <- data.frame(gc = gcC[validGc], ratio = ratio[validGc])
    # fit loess with ratio as function of gc (names must correspond to those in
    # the data frame)
    gcDist.lo <- loess(ratio ~ gc, gcDist)
    # initialize vector with prediciton
    normVector <- rep(NA, length(ratio))
    # predict the expected number of reads as function of the observed GC
    # content
    normVector[validX] <- predict(gcDist.lo, data.frame(gc = gcC[validX]))

    medianNorm <- median(normVector, na.rm = TRUE)
    vectorRatio <- normVector/medianNorm
    ratio.gcC <- ratio/vectorRatio
    return(ratio.gcC)
} # end of function

specialDensity <- function (object, special, adjust, n) {
    if (length(special) != 0){
        specialPlot <- chrs(object) %in% special
        if (length(which(specialPlot)) < 2 | length(which(!is.na(ratio(object)[specialPlot]))) < 2) {
            message(paste("    ", special,  "does not have enough values to estimate density. Skipping."))
            X <- NA
            Y <- NA
            Dn <- list(x = X, y = Y)
            return (Dn)
        }
        Dn <- myDensity(ratio(object)[specialPlot], adjust = adjust)
        X <- Dn$x
        Y <- Dn$y/length(object)*length(which(specialPlot))
    } else {
        X <- NA
        Y <- NA
    }
    Dn <- list(x = X, y = Y)
    return (Dn)
}

selectGoodGuess <- function (Density, MQRDS, spacingTolerance, interceptRatio, 
    possiblePloidy, peakRatio, peakSpan){

    
    isApeak <- myPeaks(Density$y, span = peakSpan)
    peakPos <- which(isApeak)
    peakValues <- Density$y[peakPos]
    highPeaks <- isApeak & (Density$y >= max(peakValues/peakRatio))

    # get the highest of the first three peaks
    first3validPeaks <- peakValues[peakValues %in% Density$y[highPeaks]][seq(1,3)]
    highestPeakIndex <- which(first3validPeaks == max(first3validPeaks, na.rm = TRUE))

    # select max R, positive intercept, highest peak is on ploidy == 2
    # which possible ploidy has max R
    maxR <- MQRDS[,1] >= max(MQRDS[,1])*spacingTolerance
    # which possible ploidy would have a positive intercept?
    positiveIntercept <- MQRDS[,2] >= interceptRatio
    # which possible ploidy has the highest peak either on 2?
    highestPeakOn2 <- possiblePloidy[,highestPeakIndex]==2
    highest1orMore <- possiblePloidy[,highestPeakIndex]>=1
    highest2orMore <- possiblePloidy[,highestPeakIndex]>=2

    # first get only possible solutions (no negative counts for P = 0 =>
    # intercept must be >= 0)
    goodGuess <- positiveIntercept
    if (length(which(positiveIntercept)) > 1) {
    # next look for the best/very-good alignments of points
        goodGuess <- maxR & positiveIntercept
    }

    # next look for the soltion that have a peak at ploidy = 1 or more
    if (length(which(goodGuess)) > 1) {
        goodGuess <- goodGuess & highest2orMore
    }

    # here, in same rare cases, it happens that no solutions are
    # available. We need to re-run using a smaller threshold to allow some
    # "elasticity". Recursive function.
    if (length(which(goodGuess)) == 0) {
        return (logical(0))
#             warning("No possible solutions. deacreasing parameter 'spacingTolerance'")
#             spacingTolerance = spacingTolerance * .95
#             cnaList <- guessPeaksAndPloidy(cnaList, smooth = smooth, 
#                 excludeFromDensity = excludeFromDensity, 
#                 peakRatio = peakRatio, ploidyToTest = ploidyToTest, 
#                 peakSpan = peakSpan  , adjust = adjust,
#                 spacingTolerance = spacingTolerance , interceptRatio = interceptRatio) 
#             return (cnaList)
    }



    # if still more than one option, check wich one has at least a peak for
    # ploidy >=1 (avoid genome with 1 and 0 copy)
    if ( length(which(goodGuess)) > 1 ) {
        okPloidy <- possiblePloidy[,dim(possiblePloidy)[2]] >= 2
        goodGuess <- goodGuess & okPloidy
    }

    # if still more than one option get the one(s) that require smallest
    # range of ploidy (less gaps)
    if ((length(which(goodGuess))) > 1 ) {
        ploidyDiff <- possiblePloidy[goodGuess,dim(possiblePloidy)[2]] -  possiblePloidy[goodGuess,1]
        lessGap <- which(ploidyDiff==min(ploidyDiff))
        goodGuessTmp <- rep(FALSE, length(goodGuess))
        goodGuessTmp[which(goodGuess)[lessGap]] <- TRUE
        goodGuess <- goodGuessTmp
    }

    # if still there is more than than one options, chose the one with the
    # smallest ploidy
    if (length(which(goodGuess)) > 1 ) {
        maxPloidy <- apply(possiblePloidy, 1, max)
        minMaxPloidy <- which(goodGuess)[which(maxPloidy[goodGuess] == min(maxPloidy[goodGuess]))]
        goodGuess[] <- FALSE
        goodGuess[minMaxPloidy] <- TRUE
    }

    # if still more than a solution (?) warn and pick the first.
    if (length(which(goodGuess)) > 1 ) {
        warning("Two or more solution equally possible. Picking one first available.")
        stillGood <- which(goodGuess)
        goodGuess[] <- FALSE
        goodGuess[stillGood[1]] <- TRUE
    }

    return (goodGuess)
    
}

MQRDS <- function(MQR){
    MQRDS <- data.frame(matrix(NA, length(MQR[,2]), 3))        
    colnames(MQRDS) <- c("R", "shift", "scale")

    for (n in 1:length(MQR[,1])) {
        M <- MQR[n,1]
        Q <- MQR[n,2]
        R <- MQR[n,3]
        y0 <- Q
        y2 <- 2 * M + Q
        shift <- MQR[n,2]
        scale <- 1 / (y2 - y0)
        MQRDS[n,] <- c(R, shift, scale)
    }
    return(MQRDS)
}

MQR <- function(possibleMatrix, ploidyToTest, peakToUse){
    MQR <- data.frame(matrix(NA, length(possibleMatrix[,2]), 3))
    colnames(MQR) <- c("M", "Q", "R")
    for (n in 1:length(possibleMatrix[,1])) {
        valToUse <- which(possibleMatrix[n,]==1)
        MQR[n,] <- getMQR(seq(0,ploidyToTest)[valToUse], peakToUse)
    }
    return(MQR)
}

possibleBinNumbers <- function (numDigits, numOnes) {
# produce a binary matrix of all combination on how to sort numOnes of 1s with
# numDigits digits (0 or 1)
    numOfnums = 2^numDigits 
    m <- matrix(nrow=numOfnums, ncol=numDigits)
    for (n in 1:numDigits) {
        # cat (paste("Cicle", n, "of", numDigits, "\n"))
        m[,n] <- c(rep(1, numOfnums/(2^n)), rep(0, numOfnums/(2^n)))

    }
    sumByRow <- apply(m,1,sum)
    m <- m[sumByRow == numOnes,]
    return(m)
} # end of function

getMQR <- function (xs, ys) {
# fit line across x and y, return coefficient (m) and correlation coefficient
# (R)
    mod.summary <- summary(lm(ys ~ xs))
    M <- mod.summary$coefficients[2]
    Q <- mod.summary$coefficients[1]
    R2 <- mod.summary$r.squared
    return (c(M, Q, R2))
}

myPeaks <- function (series, span = 3) {
# found here 
# http://finzi.psych.upenn.edu/R/Rhelp02a/archive/33097.html
# and adapted
    z <- embed(series, span)
    s <- span%/%2
    v<- max.col(z, "first") == 1 + s
    result <- c(rep(FALSE,s),v,rep(FALSE,s))
    return(result)
}

myDensity <- function (signal, span = 3, adjust = 1) {
    # we don't want outliers. Get only values that are smaller than the max of:
    # span times the median (up to pentaploid?)
    # Median plus 3 times the SD (if it is really aneuploid with many "peaks")
    Median <- median(signal, na.rm = TRUE)
    SD <- sd(signal, na.rm = TRUE)
    bigY <- max(c(span*Median, Median + span * SD))
    OKsignal <- signif(signal[signal <= bigY & !is.na(signal)], 5)
    Density <- density(OKsignal, adjust = adjust, na.rm = TRUE)
    return(Density)
}

ratio2use <- function(object){
    # is there is smoothed ratio, use it, otherwise use ratio
    if (length(ratio.s(object)) == length(object)){
        ratio <- ratio.s(object)
    } else {
        ratio <- ratio(object)
    }
    return (ratio)
}


# findOutliers <- function(ratio, quantile, sd){
#     q <- quantile(ratio, quantile/100, na.rm = TRUE)
#     s <- median(ratio, na.rm = TRUE) + sd * sd(ratio, na.rm = TRUE)
#     isOutlier <- ratio > q & ratio > s
#     return(isOutlier)
# }

# .peaksAndPloidy <- function (object, method = "mixture", exclude = character(0),
#     peakRatio = 50, ploidyToTest = 14, spacingTolerance = .999,
#     interceptRatio = -0.1, quantile = 99.5, sd = 5, ...){
#     
#     }
# 

## function called by the Method (setMethod)

.guessPeaksAndPloidy <- function (object, ok4density = rep(TRUE, length(object)),
    peakRatio = 50, ploidyToTest = 14, spacingTolerance = .999,         
    interceptRatio = -0.1, adjust = 0.9, n = 2048) {
  
    # a variable that should not be changed
    peakSpan = 3

    # get data in handy variables
    ratio <- ratio2use(object)
    chrName <- chrs(object)
    
    interpDiploidXloc = NA
    # isOutlier <- findOutliers(ratio, quantile, sd)
    # ok4density <- (! chrName %in% exclude) & (! isOutlier) & (! is.na(ratio))
    


#     Density <- myDensity(ratio[ok4density], adjust = adjust),
    Density <- density(ratio[ok4density], adjust = adjust, n = n)
    isApeak <- myPeaks(Density$y, span = peakSpan)
    peakPos <- which(isApeak)
    peakValues <- Density$y[peakPos]
    highPeaks <- isApeak & (Density$y >= max(peakValues/peakRatio))

    if (length(which(highPeaks)) == 1) { # only one peak, scale but don't shift (impossible to calculate )
        scaling <- 1/Density$x[highPeaks]
        interpDiploidXloc <- Density$x[highPeaks]
        MQR <- data.frame(M = Density$x[highPeaks]/2, Q = 0, R = 1)
    } else {
        possibleMatrix <- possibleBinNumbers(ploidyToTest+1, length(which(highPeaks)))
        # make table with ploidy
        possiblePloidy <- matrix(NA, length(possibleMatrix[,1]), length(which(highPeaks)))
        for (n in 1:length(possibleMatrix[,1]) ) {
            possiblePloidy[n,] <- which(possibleMatrix[n,]==1) - 1
        }
        peakToUse <- Density$x[highPeaks]
        MQR <- MQR(possibleMatrix, ploidyToTest, peakToUse)
        
        # MQRDS actually contains R shift and scale
        MQRDS <- MQRDS(MQR)
        goodGuess <- selectGoodGuess(Density, MQRDS, spacingTolerance, 
            interceptRatio, possiblePloidy, peakRatio, peakSpan)

        # here, in same rare cases, it happens that no solutions are
        # available. Rerun using a smaller threshold to allow some
        # "elasticity". Recursive function.
        if (length(goodGuess) == 0) {
            warning("No possible solutions. deacreasing parameter 'spacingTolerance'")
            spacingTolerance = spacingTolerance * .95
            object <- .guessPeaksAndPloidy(object, ok4density = ok4density, 
                peakRatio = peakRatio, ploidyToTest = ploidyToTest,
                spacingTolerance = spacingTolerance , 
                interceptRatio = interceptRatio, adjust = adjust, n = n)
                
            return(object)
        }
 
        solutions <- MQRDS[goodGuess,]
        # which one is the highest valid peak?
        scaling <- solutions$scale
        shifting <- solutions$shift
        
        thisMQR <- MQR[goodGuess,]
        interpDiploidXloc <- c(0,2) * thisMQR$M + thisMQR$Q # y = mx + q
    }

    if (exists("possiblePloidy")) {
        PP <- possiblePloidy[goodGuess,]
    } else {
        PP <- 2
    }
    if (length(interpDiploidXloc) == 1) {
        content = "Unknown"
    } else if (length( interpDiploidXloc) == 2){
        content = as.character(
            signif(100*(1-interpDiploidXloc[1]/interpDiploidXloc[2]), 3))
    } else {
        stop("interpDiploidXloc has more than 2 elements!")
    }

    # peak closest to the median
    ratioMedian <- median(ratio, na.rm = TRUE)
    peakX <- Density$x[highPeaks]
    peakBol <- (Density$x %in% peakX) & highPeaks
    medianPeaksIndex <- which(abs(Density$x[peakBol]-ratioMedian)  ==  
        min(abs(Density$x[peakBol]-ratioMedian), na.rm = TRUE))
    peakClosestToMedian <- Density$x[peakBol][medianPeaksIndex] 

    if (PP[Density$x[highPeaks] == peakClosestToMedian] == 0 
            & PP[length(PP)] == ploidyToTest){
        PP <- PP + 2
    } else if (PP[Density$x[highPeaks] == peakClosestToMedian] == 1
            & PP[length(PP)] == ploidyToTest) {
        PP <- PP + 1
    }


    object@Res@suggested.ploidy        <- PP
    object@Res@suggested.peaks         <- Density$x[highPeaks]
    object@Res@suggested.tumContent    <- content
    object@Res@interpDiploidXloc       <- interpDiploidXloc
    object@Res@notExcluded.density     <- Density
    object@Res@notExcluded.isAPeak     <- highPeaks
    object@Res@ratioMedian             <- ratioMedian
    object@Res@closestPeak             <- peakClosestToMedian
    if (exists ("thisMQR")) {
        object@Res@MQR                 <- thisMQR
        object@Res@suggested.R         <- signif(thisMQR$R, 4)
    } else {
        object@Res@MQR <- MQR
    }
    
    # object@Params@gp.excludeFromDensity <- exclude
    return(object)
}

.addDNACopy <- function (object) {
    # object is CNAnorm  
    toSegment <- object@DerivData@ratio

    toSegment[toSegment <= 0] <- .05
    # using DNA copy, segment the data
    CNA.object <- CNA(log2(toSegment), chrs(object), pos(object), data.type = 'logratio')
    smoothed.CNA.object <- smooth.CNA(CNA.object)
    segObj <- segment(smoothed.CNA.object, verbose = 0)
    segID <- rep(NA, length(object))
    segMean <- rep(NA, length(object))

    for (segNum in 1:length(segObj$output$chrom)) {
        thisChr <- as.character(segObj$output$chrom[segNum])
        thisStart <- segObj$output$loc.start[segNum]
        thisEnd <- segObj$output$loc.end[segNum]
        thisMean <- segObj$output$seg.mean[segNum]
        # pointInThisSeg <- cnaList$data$Chr == thisChr & cnaList$data$Pos >= thisStart & cnaList$data$Pos < thisEnd
        pointInThisSeg <- (chrs(object) == thisChr & pos(object) >= thisStart &
            pos(object) <= thisEnd)
        segID[pointInThisSeg] <- segNum
        segMean[pointInThisSeg] <- thisMean
    }
    object@DerivData@segID <- segID
    object@DerivData@segMean <- segMean
    return(object)
}

# version that requires DNA copy
# .addSmooth <- function (object, ...) {
#     # obj is CNAnorm  
#     if (length(segID(object)) != length(object)){
#         object <- .addDNACopy(object)
#     }
#     chrNames <- segID(object)
#     notNAIndex <- which(!is.na(chrNames))
#     tmpChr <- chrNames[notNAIndex]
#     tmpRatio <- ratio(object)[notNAIndex]
#     tmpRatio.s <- smoothSignal(tmpRatio, tmpChr, ...)
#     ratio.s <- rep(NA, length(object))
#     ratio.s[notNAIndex] <- tmpRatio.s 
#     object@DerivData@ratio.s <- ratio.s
#     return (object)
# }


# smooth chromosomes, intemediate vlaues across "jumps"
.addSmooth <- function (object, lambda = 7, ...) {
    # obj is CNAnorm  
    ratio <- ratio(object)
    chrNames <- as.numeric(as.factor(chrs(object)))
    
    # notNAIndex <- which(!is.na(chrNames))
    notNAIndex <- which(!is.na(ratio))
    tmpChr <- chrNames[notNAIndex]
    tmpRatio <- ratio[notNAIndex]
    tmpRatio.s <- smoothSignal(tmpRatio, tmpChr, lambda = lambda, ...)
    ratio.s <- rep(NA, length(object))
    ratio.s[notNAIndex] <- tmpRatio.s 
    object@DerivData@ratio.s <- ratio.s
    return (object)
}




.gcNorm <- function(object, exclude = character(0), maxNumPoints = 10000){
    if (length(object@InData@GC) == 0){
        stop("Impossible to perform GC content normalization: no GC content available\n")
    } 
    if (length(exclude) == 0){
        points2use = NULL
    } else {
        points2use <- ! object@InData@Chr %in% exclude
        if (length(which(points2use)) < 2) {
            warning ("Not enough points to perform GC normalization.\n")
            return (object)
        }
    }
    object@DerivData@ratio <- gcNormalize(gcC(object), object@InData@ratio,
        points2use = points2use, maxNumPoints = maxNumPoints)
    # now save parameters used
    object@Params@gc.excludeFromGCNorm = exclude
    object@Params@gc.maxNumPoints = maxNumPoints

    return (object)    
}
makeTextLabels <- function (xLoc, Ploidy) {
    if (length(xLoc) != length(Ploidy)) {
        stop ("xLoc and Ploidy must have the same length")
    }
    textLab <- rep(NA, length(xLoc))
    for (i in 1:length(textLab)) {
#         textLab[i] <- paste("x = ", as.character(signif(xLoc[i],2)), 
#         ", Pl = ", as.character(Ploidy[i]), sep = '')
        textLab[i] <- paste(as.character(signif(xLoc[i],2)), " -> ", as.character(Ploidy[i]), sep = '')
    }
    return(textLab)
}

plotPeaksMixture <- function (object, special1 = character(0), special2 = character(0),
    show ='suggested', adjust = get.adjust(object), n = get.n(object), bins = 200, ...) {
    
    ratioMedian <- median(ratio(object), na.rm = TRUE)
#     D1 <- specialDensity(object, special1, adjust, n)
#     D2 <- specialDensity(object, special2, adjust, n)
    
    # xRange <- c(0, max(1.2*ratio.s(object), na.rm = TRUE))
    R.s <- ratio.s(object)
    hist(R.s, n = bins, xlab = 'Ratio Test/Control' )

    # initialize legend variables
    if (length(object@Params@gp.excludeFromDensity) > 0) {
        notExclLab <- paste("Everything but", 
            paste(object@Params@gp.excludeFromDensity, collapse = ', '))
    } else {
        notExclLab <- "Everything"
    }
    textLegend <- as.character(c(notExclLab))
    legCol = c('black')
    legLty = c(1)
    legLwd = c(3)
 
    if (show == 'suggested') {
        peakX <- object@Res@suggested.peaks
        ploidy <- object@Res@suggested.ploidy/2
        tumContent <- object@Res@suggested.tumContent
    } else if (show == 'validated') {
        peakX <- object@Res@validated.peaks
        ploidy <- object@Res@validated.ploidy/2
        tumContent <- object@Res@validated.tumContent
    } else {
        stop ("`show` must be either `suggested` or `validated`")
    }

    abline(v = ratioMedian, col = 'red', lty = 2, lwd = 2)
    abline(v = peakX, col = 'cyan', lwd = 2)


}

.plotPeaks <- function (object, special1 = character(0), special2 = character(0),
    show ='suggested', adjust = get.adjust(object), n = get.n(object), ...) {
    # as input provide the output of guessPeaksAndPloidy
    # Normalize by number of reads in ok4density
    bins = 200
    DD <- density(ratio2use(object)[object@DerivData@ok4density], adjust = adjust, n = n)
    X <- DD$x
    Y <- DD$y/length(object)*length(which(object@DerivData@ok4density))
    ratioMedian <- median(ratio(object), na.rm = TRUE)
    D1 <- specialDensity(object, special1, adjust, n)
    D2 <- specialDensity(object, special2, adjust, n)
   
    xRange <- c(0, max(X*1.2, na.rm = TRUE))
    yRange <- c(0, max(Y*1.2, na.rm = TRUE))

#     xRange <- c(0, max(c(X, D1$x, D2$x)*1.2, na.rm = TRUE))
#     yRange <- c(0, max(c(Y, D1$y, D2$y)*1.2, na.rm = TRUE))
    vShift = yRange[2]/20 
    # create an empty plot space
    plot(xRange, yRange, type = 'n', xlab = 'Ratio Test/Control', 
        ylab = 'Density', ...)

    # plot density of "non special" genome
    lines (X, Y, type = 'l', lw = 3)

    # initialize legend variables
    if (length(object@Params@gp.excludeFromDensity) > 0) {
        notExclLab <- paste("Everything but", 
            paste(object@Params@gp.excludeFromDensity, collapse = ', '))
    } else {
        notExclLab <- "Everything"
    }
    textLegend <- as.character(c(notExclLab))
    legCol = c('black')
    legLty = c(1)
    legLwd = c(3)
    

   
    # plot special 1
    if (!all(is.na(D1$x))) {
        lines (D1$x, D1$y, type = 'l', lty = 1, lw = 3, col = 'orchid4')
        specialLab1 <- paste(special1, collapse = ', ')
        textLegend <- c(textLegend, specialLab1)
        legCol = c(legCol, 'orchid4')
        legLty = c(legLty, 1)
        legLwd = c(legLwd, 3)
    }
    # plot special 2
    if (!all(is.na(D2$x))) {
        lines (D2$x, D2$y, type = 'l', lw = 3, lty = 1, col = 'lightblue')
        specialLab2 <- paste(special2, collapse = ', ')
        textLegend <- c(textLegend, specialLab2)
        legCol = c(legCol, 'lightblue')
        legLty = c(legLty, 1)
        legLwd = c(legLwd, 3)
    }
    
    # plot annotation of "non special" genome
    if (show == 'suggested') {
        peakX <- object@Res@suggested.peaks
        ploidy <- object@Res@suggested.ploidy/2
        tumContent <- object@Res@suggested.tumContent
    } else if (show == 'validated') {
        peakX <- object@Res@validated.peaks
        ploidy <- object@Res@validated.ploidy/2
        tumContent <- object@Res@validated.tumContent
    } else {
        stop ("`show` must be either `suggested` or `validated`")
    }


    abline(v = ratioMedian, col = 'red', lty = 2)
    
    mth <- object@Params@method
    ##### THIS IS WHERE THE FUNCTION DOES NOT WORK ANYMORE FOR MIXTURE MODEL...
    ### check into peakBol
    if (mth == 'density'){
        peakBol <- (X %in% peakX) & object@Res@notExcluded.isAPeak
        ww <- which(peakBol)    
    } else {
        isAPeak <- mapPeaks(X, peakX)
        ww <- which(isAPeak)
        for (i in 1:length(ww)){
            lines ( c(X[ww[i]], X[ww[i]]), c(0, Y[ww[i]]))
        }
    }
    
#     medianPeaksIndex <- which(abs(X[peakBol]-ratioMedian)  ==  min(abs(X[peakBol]-ratioMedian), na.rm = TRUE))
#     # actual plotting
#     lines (X[peakBol], Y[peakBol], type = 'p', pch=19, col='blue')
#     lines (X[peakBol][medianPeaksIndex], Y[peakBol][medianPeaksIndex], type = 'p', cex = 1.8, col='red')
#     textLab <- makeTextLabels(X[peakBol], ploidy)
#     text(X[peakBol], Y[peakBol] + vShift, textLab, cex = .8, col = 'blue')
#     
    medianPeaksIndex <- which(abs(X[ww]-ratioMedian)  ==  min(abs(X[ww]-ratioMedian), na.rm = TRUE))
    # actual plotting
    lines (X[ww], Y[ww], type = 'p', pch=19, col='blue')
    lines (X[ww][medianPeaksIndex], Y[ww][medianPeaksIndex], type = 'p', cex = 1.8, col='red')
    textLab <- makeTextLabels(X[ww], ploidy)
    text(X[ww], Y[ww] + vShift, textLab, cex = .8, col = 'blue')
   
    # add tumour content label
    legCol = c(legCol, 'black')
    legLty = c(legLty, 0)
    legLwd = c(legLwd, 0)
    textLegend <- c(textLegend,  paste("Est Tumour Content: ", tumContent, "%", sep = '') )
    

    legend("topright", inset=0.01, textLegend, lty = legLty, 
        lwd = legLwd, col = legCol, cex = .8)
}

mapPeaks <- function(map, peaks) {
    # map is a vector (X) where the peak will be mapped.
    # returns A logic vector of values of X closest to peaks
    isAPeak <- rep(FALSE, length(map))
    for (i in 1: length(peaks)){
        thisPeak <- peaks[i]
        ww <- which(abs(map - thisPeak) == (min(abs(map - thisPeak), na.rm = TRUE)))
        if (length(ww) > 1){
            ww <- ww[1]
        }
        isAPeak[ww] <- TRUE
    }
    return (isAPeak)
}



.plotPeaks_old <- function (object, special1 = character(0), special2 = character(0),
    show ='suggested', adjust = .9, n = 2048, ...) {
    # as input provide the output of guessPeaksAndPloidy
    # Normalize by number of reads in ok4density
    DD <- object@Res@notExcluded.density
    X <- DD$x
    Y <- DD$y/length(object)*length(which(object@DerivData@ok4density))
    ratioMedian <- median(ratio(object), na.rm = TRUE)
    D1 <- specialDensity(object, special1, adjust, n)
    D2 <- specialDensity(object, special2, adjust, n)
   
    xRange <- c(0, max(X*1.2, na.rm = TRUE))
    yRange <- c(0, max(Y*1.2, na.rm = TRUE))

#     xRange <- c(0, max(c(X, D1$x, D2$x)*1.2, na.rm = TRUE))
#     yRange <- c(0, max(c(Y, D1$y, D2$y)*1.2, na.rm = TRUE))
    vShift = yRange[2]/20 
    # create an empty plot space
    plot(xRange, yRange, type = 'n', xlab = 'Ratio Test/Control', 
        ylab = 'Density', ...)

    # plot density of "non special" genome
    lines (X, Y, type = 'l', lw = 3)

    # initialize legend variables
    if (length(object@Params@gp.excludeFromDensity) > 0) {
        notExclLab <- paste("Everything but", 
            paste(object@Params@gp.excludeFromDensity, collapse = ', '))
    } else {
        notExclLab <- "Everything"
    }
    textLegend <- as.character(c(notExclLab))
    legCol = c('black')
    legLty = c(1)
    legLwd = c(3)
    

   
    # plot special 1
    if (!all(is.na(D1$x))) {
        lines (D1$x, D1$y, type = 'l', lty = 1, lw = 3, col = 'orchid4')
        specialLab1 <- paste(special1, collapse = ', ')
        textLegend <- c(textLegend, specialLab1)
        legCol = c(legCol, 'orchid4')
        legLty = c(legLty, 1)
        legLwd = c(legLwd, 3)
    }
    # plot special 2
    if (!all(is.na(D2$x))) {
        lines (D2$x, D2$y, type = 'l', lw = 3, lty = 1, col = 'lightblue')
        specialLab2 <- paste(special2, collapse = ', ')
        textLegend <- c(textLegend, specialLab2)
        legCol = c(legCol, 'lightblue')
        legLty = c(legLty, 1)
        legLwd = c(legLwd, 3)
    }
    
    # plot annotation of "non special" genome
    if (show == 'suggested') {
        peakX <- object@Res@suggested.peaks
        ploidy <- object@Res@suggested.ploidy/2
        tumContent <- object@Res@suggested.tumContent
    } else if (show == 'validated') {
        peakX <- object@Res@validated.peaks
        ploidy <- object@Res@validated.ploidy/2
        tumContent <- object@Res@validated.tumContent
    } else {
        stop ("`show` must be either `suggested` or `validated`")
    }
    peakBol <- (X %in% peakX) & object@Res@notExcluded.isAPeak
    medianPeaksIndex <- which(abs(X[peakBol]-ratioMedian)  ==  min(abs(X[peakBol]-ratioMedian), na.rm = TRUE))
    # actual plotting
    abline(v = ratioMedian, col = 'red', lty = 2)
    lines (X[peakBol], Y[peakBol], type = 'p', pch=19, col='blue')
    lines (X[peakBol][medianPeaksIndex], Y[peakBol][medianPeaksIndex], type = 'p', cex = 1.8, col='red')
    textLab <- makeTextLabels(X[peakBol], ploidy)
    text(X[peakBol], Y[peakBol] + vShift, textLab, cex = .8, col = 'blue')
    
    # add tumour content label
    legCol = c(legCol, 'black')
    legLty = c(legLty, 0)
    legLwd = c(legLwd, 0)
    textLegend <- c(textLegend,  paste("Est Tumour Content: ", tumContent, "%", sep = '') )
    

    legend("topright", inset=0.01, textLegend, lty = legLty, 
        lwd = legLwd, col = legCol, cex = .8)
}
.subSet <- function (x, i, j, drop) {
    if (!missing(j)) {
        stop("Error: incorrect number of dimensions\n")
    } else {
        Chr <- chrs(x)[i]
        Pos <- pos(x)[i]
        Test <- if(length(x@InData@Test) > 0){x@InData@Test[i]}else{numeric()}
        Norm <- if(length(x@InData@Norm) > 0){x@InData@Norm[i]}else{numeric()}
        ratio <- if(length(x@InData@ratio) > 0){x@InData@ratio[i]}else{numeric()}
        GC <- if(length(x@InData@GC) > 0){x@InData@GC[i]}else{numeric()}
        newInData <- new("InData", Chr = Chr, Pos = Pos, Test = Test, Norm = Norm, 
            ratio = ratio, GC = GC)
        
        ratio <- if(length(ratio(x)) > 0){ratio(x)[i]}else{numeric()}
        ratio.n <- if(length(ratio.n(x)) > 0){ratio.n(x)[i]}else{numeric()}
        ratio.s <- if(length(ratio.s(x)) > 0){ratio.s(x)[i]}else{numeric()}
        ratio.s.n <- if(length(ratio.s.n(x)) > 0){ratio.s.n(x)[i]}else{numeric()}
        segID <- if(length(x@DerivData@segID[i]) > 0){x@DerivData@segID[i]}else{numeric()}
        segMean <- if(length(segMean(x)) > 0){segMean(x)[i]}else{numeric()}
        segMean.n <- if(length(segMean.n(x)) > 0){segMean.n(x)[i]}else{numeric()}
        ok4density <- if(length(x@DerivData@ok4density) > 0){x@DerivData@ok4density[i]}else{logical()}
        
        newDerivData <- new("DerivData", ratio = ratio, ratio.n = ratio.n, ratio.s = ratio.s,
            ratio.s.n = ratio.s.n, segID = segID, segMean = segMean, segMean.n = segMean.n,
            ok4density = ok4density)
        newRes <- x@Res
        newParams <- x@Params

        
        newObject <- new("CNAnorm", InData = newInData, DerivData = newDerivData, 
            Res = newRes, Params = newParams)
    }
    return(newObject)
}

.exportTable <- function (object, file = "CNAnorm_table.tab", show = 'ratio', sep = "\t", row.names = FALSE, ...){
    roundPloidy <- object@Res@validated.closestPeak
    regExp <- paste("^", show, sep = '')
    showOptions <- c('ratio', 'ploidy', 'center')
    ww <- grep(regExp, showOptions, ignore.case = TRUE)
    if (length(ww) > 1){
        warning("ambigous value for 'show', using 'ratio'...\n")
        ww <- 1
    } else if (length(ww) == 0){
        warning("'show' not matching available options, using 'ratio'...\n")
        ww <- 1
    }
    if (ww == 1){
        mulP <- 1
    } else if (ww == 2){
        mulP <- 2
    } else {
        mulP <- 1/roundPloidy
    }


    numEl <- length(chrs(object))
    if (length(segMean(object)) == numEl){
        SegMean <- segMean(object) * mulP
        SegMean.n <- segMean.n(object) * mulP
    } else {
        SegMean <- rep(NA, length(ratio.s.n(object)))
        SegMean.n <- rep(NA, length(ratio.s.n(object)))
    }
    if (length(ratio.s.n(object)) == numEl){
        Ratio.s.n = ratio.s.n(object) * mulP
    } else {
        Ratio.s.n = rep(NA, numEl)
    }

    df <- data.frame(
        Chr = chrs(object),
        Pos = pos(object),
        Ratio = ratio(object),
        Ratio.n = ratio.n(object) * mulP,
        Ratio.s.n = Ratio.s.n,
        SegMean = SegMean,
        SegMean.n = SegMean.n
    )
#     return (df)
    write.table(df, file = file, sep = sep, row.names = row.names, ...)
}

setMethod(f = "gcNorm", signature = "CNAnorm", definition = .gcNorm)
setMethod(f = "addSmooth", signature = "CNAnorm", definition = .addSmooth)
setMethod(f = "addDNACopy", signature = "CNAnorm", definition = .addDNACopy)
setMethod(f = "peakPloidy", signature = "CNAnorm", definition = .peakPloidy)
setMethod(f = "validation", signature = "CNAnorm", definition = .validation)
setMethod(f = "discreteNorm", signature = "CNAnorm", definition = .discreteNorm)
setMethod(f = "plotGenome", signature = "CNAnorm", definition = .plotGenome)
setMethod(f = "plotPeaks", signature = "CNAnorm", definition = .plotPeaks)
setMethod(f = "[", signature = "CNAnorm", definition = .subSet)
setMethod(f = "exportTable", signature = "CNAnorm", definition = .exportTable)
