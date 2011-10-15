# summary
inData.summary <- function(object, ...){
    DF <- makeDataFrameFromSlots(object)
    numSlots <- length(slotNames(object))
    rowNum <- length(DF[,1])
    colNum <- length(DF[1,])
    cat(paste("Object with ", numSlots, " slots. \n", sep = ""))
    if(numSlots == 0){
        cat("All slots are empty\n")
    } else if (colNum == 1) {
        cat(paste("One slot has ", rowNum, " elements\n", sep = ""))
        if (length(object) < 8){
            print (DF)
        } else {
            totNumRow <- length(DF[,1])
            skippedRows <- totNumRow - 7 
            smallDF <- as.data.frame(DF[1:7,])
            colnames(smallDF) <- colnames(DF)
            print (smallDF)
            cat(paste("... and other", skippedRows, "rows\n", sep = " "))
        } 
    } else {
        
        cat(paste(colNum, " slots have ", rowNum, " elements\n", sep = ""))
        if (length(object) < 8){
            print (DF)
        } else {
            totNumRow <- length(DF[,1])
            skippedRows <- totNumRow - 7 
            print (DF[1:7, ])
            cat(paste("... and other", skippedRows, "rows\n", sep = " "))
        }
    }
    return (invisible())
}

Res.summary <- function(object, ...) {
    SN <- slotNames(object)
    cat(paste( "Object with ", length(SN), " slots.\n", sep = ""))
    nonEmptySlots <- vector()
    elementNumber <- vector()
    for (i in 1:length(SN)) {
        thisSlotContent <- slot(object, SN[i])
        if (length(thisSlotContent) > 0) {
            nonEmptySlots <- c(nonEmptySlots, SN[i])
            elementNumber <- c(elementNumber, length(thisSlotContent))
        } 
    }
    if (length(nonEmptySlots) > 0) {
        cat("Non empty slots (num of elements):\n")
        for (i in 1:length(nonEmptySlots)) {
            cat(paste(nonEmptySlots[i], " (", elementNumber[i], ")\n", sep = ""))
        }
        cat ("\n")
    } else {
        cat("All slots are empty\n")
    }
    return (invisible())
}
setMethod ("summary", "InData", definition = inData.summary)
setMethod ("summary", "DerivData", definition = inData.summary)
setMethod ("summary", "Res", definition = Res.summary)
setMethod ("summary", "Params", definition = Res.summary)
setMethod ("summary", "CNAnorm", definition = function(object, ...){
        show(object)
        return(invisible())
    }
)

## show
showCNAnorm <- function (object) {
    className <- as.character(class(object))
    slots <- slotNames(object)

    cat (paste("An object of class", className, "\n", sep = " "))

    if (length(object) == 0) {
        cat ("Empty object\n")
    } else {
        cat (paste("with", as.numeric(length(slots)), "slots:\n", sep = " "))
        cat (paste(slots, collapse = ", "))
        cat ("\n")

        cat ("\nSummary of slot InData:\n")
        summary(object@InData)
        cat ("\nSummary of slot DerivData:\n")
        summary(object@DerivData)
        cat ("\nSummary of slot Res:\n")
        summary(object@Res)
        cat ("\nSummary of slot Params:\n")
        summary(object@Params)
    }
}
setMethod("show", "CNAnorm", definition = showCNAnorm)

