.InData.validity <-  function(object) {
    # check if the user passed ratio
    dataLength <- length(object@ratio)
    if (dataLength != 0) {
        # we have ratio
        if (length(object@Test) != 0 | length(object@Norm) != 0) {
            # we have both Test/Norm and ratio, let's check they match
            if (length(object@Test) != dataLength | length(object@Norm) != dataLength){
                stop ("Inconsistency detected: input ratio, Test and Normal have different length\n")
            } else if (any(object@ratio != justRatio(object), na.rm = TRUE)){
                # the ratio is not the ratio
                stop ("Provided ratio does not match Test/Control")
            } else {}
        }
    } else {
        # the user did not pass ratio
        if (!.sameLength(object@Test, object@Norm)){
            stop ("Test and Norm must have same length\n")
        } else {
            dataLength <- length(object@Test)
        }
    }
    # check if data was provided
    if (dataLength == 0) {
        stop ("data provided is zero Length\n")
    }
    
    # check length of data match length of Chr and Pos
    allLengths <- c(length(object@Chr), length(object@Pos), dataLength)
    if (!all(allLengths[1] == allLengths, na.rm = TRUE)){
        stop("Chr, Pos, Test and Norm (or ratio) must have same length\n")
    }
    if (length(object@GC) != 0 & length(object@GC) != dataLength) {
        stop ("GC must have same length than provided data\n")
    }
    return (TRUE)
}

# Input data. only for first calculation and reference 
setClass (
    Class = "InData", 

    #### SLOTS ####
    representation = representation (
        Chr = "character", # name of contig/chromosome
        Pos = "numeric", # starting position of window
        Test = "numeric", # number of reads in test
        Norm = "numeric", # number of reads in control/normal
        # either provide Test & Norm OR ratio
        ratio = "numeric",
        GC = "numeric" # GC content, optional
    ), 
    # check validity
    validity =.InData.validity
)

# Derived data. Alway work on this
setClass (
    Class = "DerivData",

    #### SLOTS ####
    representation = representation (
        # numRow = "numeric",
        ratio = "numeric",
        ratio.n = "numeric",
        ratio.s = "numeric", 
        ratio.s.n = "numeric",
        segID = "numeric",
        segMean = "numeric",
        segMean.n = "numeric",
        ok4density = "logical"
    )
)


setClass (
    Class = "Res",

    #### SLOTS ####
    representation = representation (
        suggested.ploidy = "numeric",
        suggested.peaks = "numeric",
        suggested.tumContent = "character",
        suggested.R = "numeric",
        validated.ploidy = "numeric",
        validated.peaks = "numeric",
        validated.tumContent = "character",
        validated.R = "numeric",
        validated.ratioMedian = "numeric",
        validated.closestPeak = "numeric",
        interpDiploidXloc = "numeric",
        notExcluded.density = "density",
        notExcluded.isAPeak = "logical",
        ratioMedian = "numeric",
        closestPeak = "numeric",
        MQR = "data.frame",
        shifting = "numeric",
        scaling = "numeric",
        content = "character"
    )

)

setClass (
    Class = "Params",

    #### SLOTS ####
    representation = representation (
        method = "character",
        density.n = "numeric",
        density.adjust = "numeric",
        gc.excludeFromGCNorm = "character",
        gc.maxNumPoints = "numeric",
        gp.excludeFromDensity = "character"
    )
)

setClass (
    Class = "CNAnorm", 

    #### SLOTS ####
    representation = representation (
        InData = "InData", 
        DerivData = "DerivData",
        Params = "Params",
        Res = "Res"
        )
)
