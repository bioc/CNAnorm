## initialize

justRatio <- function(obj) {
    # get the ratio where Norm != 0
    # obj is InData
    normNoZ <- (obj@Norm != 0) & (!is.na(obj@Norm))
    ratio <- rep(NA, length(obj@Norm))
    ratio[normNoZ] <- obj@Test[normNoZ]/obj@Norm[normNoZ]
    return (ratio)
}



setMethod ("initialize", "CNAnorm", 
    # function (.Object, InData, DerivData){
    function (.Object, InData, DerivData = new("DerivData"), Res = new("Res"), 
        Params = new("Params")){
        .Object@InData <- InData
        # .Object@DerivData@numRow <- length(InData@Chr)
        # derivData is empty, but we already have InData ratio
        if (length(DerivData@ratio) == 0 & length(InData@ratio) == length(InData@Chr)){
            .Object@DerivData@ratio <- InData@ratio 
        # derivData is empty, and don't have a ratio yet
        } else if (length(DerivData@ratio) == 0 & length(InData@ratio) == 0) {
            .Object@DerivData@ratio <- justRatio(InData)
        # derivData is not empty, but coherent with InData (as when subSetting)
        } else if (length(DerivData@ratio) == length(InData@Chr)) {
            .Object@DerivData <- DerivData
        } else {
            stop ("length of slot ratio in DerivData different from length of ratio in InData\n")
        }

        
        if(!missing(Res)){
            .Object@Res <- Res
        }
        if (!missing(Params)){
            .Object@Params <- Params
        }
        return (.Object)
   }
)

setMethod ("initialize", "InData", 
    function (.Object, Chr, Pos, Test = numeric(0), Norm = numeric(0),
        ratio = numeric(0), GC = numeric(0)) {
        .Object@Chr <- Chr
        .Object@Pos <- Pos
        .Object@Test <- Test
        .Object@Norm <- Norm
        .Object@ratio <- ratio
        if (length(ratio) == 0) {
            .Object@ratio <- justRatio(.Object)
        }
        if (length(unique(GC)) == 1) {
#             warning("All GC content is constant. Ignoring GC")
            GC = numeric(0)
        }
        .Object@GC <- GC
        return (.Object)
    }
)


## validation
.sameLength <- function (a,b){
    if (length(a) == length(b)) {
        return (TRUE)
    } else {
        return (FALSE)
    }
}




