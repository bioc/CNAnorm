# length
setMethod ("length", "CNAnorm", function (x) {length(x@InData)})
setMethod ("length", "InData", function (x) {length(x@Chr)})
setMethod ("length", "DerivData", function (x) {length(x@ratio)})

