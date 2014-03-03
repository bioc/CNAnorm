# this function is a mega wrapper to use for a fully automated CNAnorm workflow where
# interactivity is not required. It contains MOST possible paramenters. Defaults
# are set to to run a standard and conservative workflow.


CNAnormWorkflow <- function(dataFrame, 
    gc.do = FALSE, gc.exclude = character(0), gc.maxNumPoints = 10000,
    
    smooth.do = TRUE, smooth.lambda = 7, smooth.other = list(), 
    
    peak.method = 'closest', 
    peak.exclude = character(0), peak.ploidyToTest = 12, peak.sd = 5, peak.dThresh = 0.01,
    peak.n = 2048, peak.adjust = .9, peak.force.smooth = TRUE, 
    peak.reg = FALSE, peak.ds = 1.5, peak.zero.count = FALSE, peak.other = list(), 

    DNAcopy.do = TRUE, DNAcopy.independent.arms = FALSE, DNAcopy.ideograms = NULL,
    DNAcopy.smooth = list(), DNAcopy.segment = list(),
    DNAcopy.weight=character(),

    dNorm.normBy = NULL
){

    # df is the input dataframe as loaded by data(LS041)

    ##### some basic error checking before we start
    ## non-independent arguments
    if (DNAcopy.do & DNAcopy.independent.arms){
        if (is.null(DNAcopy.ideograms)){
            warning("'DNAcopy.independent.arms' is set to TRUE, but DNAcopy.ideograms is NULL.\n", 
                call. = FALSE)
            warning("Re-setting 'DNAcopy.independent.arms' to FALSE. See ?hg19_hs_ideogr for more information",
                call. = FALSE)
            DNAcopy.independent.arms <- FALSE
        }
    }
    
    ## make sure all the *.other arguments are key/values pairs

    # convert to object
    CN <- dataFrame2object(dataFrame)

    # perform (optional) gc normalisation
    if(gc.do){
        CN <- gcNorm(CN, exclude = gc.exclude)
    } else {
        # nothing, we already have CN
    }

    # perform (optional) data smoothing
    if(smooth.do){
        # as we allowed optional arguments, we need to call with do.call
        CN <- do.call(addSmooth, c(list(CN, lambda = smooth.lambda), smooth.other))
    } else {
        # nothing, we already have CN
    }

    # Let's peak the peak(s)
    # as we allowed optional arguments, we need to call with do.call
    CN <- do.call(peakPloidy, c(list(CN, exclude = peak.exclude, 
        method = peak.method, ploidyToTest = peak.ploidyToTest, sd = peak.sd, 
        adjust = peak.adjust, n = peak.n, force.smooth = peak.force.smooth, 
        reg = peak.reg, ds = peak.ds, zero.count = peak.zero.count, 
        dThresh = peak.dThresh), peak.other))

    CN <- validation(CN)
    # Let's add segmentation 
    if (DNAcopy.do){
        CN <- addDNACopy(CN, independent.arms = DNAcopy.independent.arms,
            ideograms = DNAcopy.ideograms, DNAcopy.smooth = DNAcopy.smooth,
            DNAcopy.segment = DNAcopy.segment, DNAcopy.weight = DNAcopy.weight)
    } else {
        # # nothing, we already have CN
    }
    if (is.null(dNorm.normBy)){
        nCN <- CN
    } else {
        nCN <- dNorm.normBy
    }
    CN <- discreteNorm(CN, normBy = nCN)
    
    return(CN)

}



