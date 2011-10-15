## Getters
setMethod ("gcC", "CNAnorm", function(object){object@InData@GC})
setMethod ("chrs", "CNAnorm", function(object){object@InData@Chr})
setMethod ("pos", "CNAnorm", function(object){object@InData@Pos})
setMethod ("segID", "CNAnorm", function(object){object@DerivData@segID})
setMethod ("ratio", "CNAnorm", function(object){object@DerivData@ratio})
setMethod ("ratio.n", "CNAnorm", function(object){object@DerivData@ratio.n})
setMethod ("ratio.s", "CNAnorm", function(object){object@DerivData@ratio.s})
setMethod ("ratio.s.n", "CNAnorm", function(object){object@DerivData@ratio.s.n})
setMethod ("segMean", "CNAnorm", function(object){object@DerivData@segMean})
setMethod ("segMean.n", "CNAnorm", function(object){object@DerivData@segMean.n})
setMethod ("sugg.peaks", "CNAnorm", function(object){object@Res@suggested.peaks})
setMethod ("sugg.ploidy", "CNAnorm", function(object){object@Res@suggested.ploidy})
setMethod ("valid.peaks", "CNAnorm", function(object){object@Res@validated.peaks})
setMethod ("valid.ploidy", "CNAnorm", function(object){object@Res@validated.ploidy})
setMethod ("get.adjust", "CNAnorm", function(object){object@Params@density.adjust})
setMethod ("get.n", "CNAnorm", function(object){object@Params@density.n})

## Setters
setGeneric("chrs<-", function(object, value){standardGeneric("chrs<-")})
setReplaceMethod("chrs", "CNAnorm", function(object, value){
    object@InData@Chr <- value
    validObject(object@InData)
    return(object) })
setGeneric("pos<-", function(object, value){standardGeneric("pos<-")})
setReplaceMethod("pos", "CNAnorm", function(object, value){
    object@InData@Pos <- value
    validObject(object@InData)
    return(object)})


