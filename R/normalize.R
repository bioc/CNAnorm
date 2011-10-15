.manualValidation <- function (object, peaks = object@Res@suggested.peaks,
    ploidy = object@Res@suggested.ploidy) {
    peaks <- peaks
    ploidy <- ploidy
    if (length(peaks) != length(ploidy)) {
        stop ("Peaks and ploidy must have same length\n" )
    } else if (any(round(ploidy) != ploidy)) {
        stop ("Ploidy must be integer\n")
    } else {
        object@Res@suggested.peaks = peaks
        object@Res@suggested.ploidy = ploidy
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

        object@Res@validated.tumContent = content
        object@Res@validated.ratioMedian <- 
            (object@Res@ratioMedian - shift) * scale
        object@Res@validated.closestPeak <- 
            (object@Res@closestPeak - shift) * scale
        object@Res@validated.R           <- signif(MQR[3], 4) 
    }

    return (object)
}


