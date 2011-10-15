.buildWarn <- function(field) {
    msg <- paste("Data frame does not have column \"", field, "\".", sep = "")
    return(msg)
}

dataFrame2object <- function(dataFrame){
    N <- names(dataFrame)
    requiredFields <- c("Chr", "Pos", "Test", "Norm")
    for (i in 1:length(requiredFields)){ 
        if(! requiredFields[i] %in% N){
            thisM <- .buildWarn(requiredFields)
            stop(thisM)
        }
    }
    if ("GC" %in% N){
        InD <- new("InData", Chr = as.character(dataFrame$Chr), Pos = dataFrame$Pos, 
            Test = dataFrame$Test, Norm = dataFrame$Norm, GC = dataFrame$GC)
        
    } else {
        InD <- new("InData", Chr = as.character(dataFrame$Chr), Pos = dataFrame$Pos, 
            Test = dataFrame$Test, Norm = dataFrame$Norm)
    }

    CN <- new("CNAnorm", InData = InD)
    return(CN)
}

