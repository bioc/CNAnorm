


makeDataFrameFromSlots <- function(x, row.names = NULL, optional = FALSE, ...) {
    SN <- slotNames(x)
    vectorDF <- vector()
    listDF <- list()
    thisObjectSize <- length(x)
    successfullyConverted = FALSE
    for (i in 1:length(SN)){
        thisSlotContent <- slot(x, SN[i])
        
        if (length(thisSlotContent) == length(x) & length(x) > 1) {
            listDF[[ SN[i] ]] <- thisSlotContent
            successfullyConverted = TRUE
        }
    }
    if (successfullyConverted){
        DF <- as.data.frame(listDF, ...)
        return (DF)
    } else {
        return(NULL)
    }
}

###########################################
## Already declared
# as.data.frame
setMethod ("as.data.frame", "InData", definition = makeDataFrameFromSlots)
setMethod ("as.data.frame", "DerivData", definition = makeDataFrameFromSlots)


## Replace

## validity


