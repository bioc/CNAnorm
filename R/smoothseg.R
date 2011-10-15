`smoothseg` <-
function(xdat, position=NULL, chrom=1, lambda=1, impute.only=FALSE,
show.plot=FALSE, display.sample=1, display.chrom=0, bw=1, maxiter=3, k=2)
{ # Function to perform smooth segmentation of CGH data
  # Arief.Gusnanto@mrc-bsu.cam.ac.uk and j.huang@ucc.ie
  
  if(is.null(xdat)) stop("Argument xdat is missing.")
  if(!(class(xdat)=="matrix" | class(xdat)=="data.frame")){
  X <- as.matrix(xdat)
  } else {
    X <- xdat }
  
  ncolX <- ncol(X)
  nrowX <- nrow(X)
  if(!(length(chrom)==1 | length(chrom)==nrowX)) stop("The length of argument chrom does not match. You can set it to 1 for genome-wide smoothing.")
  if(length(chrom)==1){
  if(chrom!=1) stop("Set the argument chrom to 1 for genomewide processing, or vectors of length nrow(xdat) for smoothing per chromosome.")
  chrom=rep(1,nrowX)
  }
  
  if(is.null(position)) {
  x.axis <- 1:nrowX
  } else {
    if(length(c(position)) != nrowX) stop("Argument position has different length to nrow(xdat).")
    x.axis <- position
    }

  if(impute.only){
  if(!any(is.na(X))) stop("Argument impute.only is set TRUE, but no missing values (NA's). Use default setting.")
  }
  
  if(any(is.na(X))){
        id.col.missing <- which(apply(X,2,function(a) any(is.na(a))))
          for(j in id.col.missing) {
              X[is.na(X[,j]), j] <- mean(X[!is.na(X[,j]), j])
           } # end of iteration in imputing each column
  } # End of if(any(is.na(X)))
  
  
  if(!impute.only){
  uni.chrom <- unique(chrom)
  smoothX <- matrix(NA, nrow=nrowX, ncol=ncolX)
  for(i in 1:ncolX) {
      for(j in uni.chrom) {
          smoothX[chrom==j,i] <- rseg(c(x.axis)[chrom==j], X[chrom==j,i], lambda=lambda,
                                 bw=bw, maxiter=maxiter, k=k)$y
      } # End of iteration for chromosomes in each sample
  } # End of iteration for samples
  } # End of !impute.only
  
  if(show.plot) {
    if(display.chrom==0){
        plot(x.axis, X[,display.sample], type="p", main=paste("Sample ",display.sample),
          xlab="Position", ylab="Intensity", pch=19, cex=0.4)
        if(!impute.only) lines(x.axis, smoothX[,display.sample], col=2)
      } else {
        plot(c(x.axis)[chrom==display.chrom], X[chrom==display.chrom,display.sample],
          type="p", main=paste("Sample ",display.sample),
          xlab="Position", ylab="Intensity", pch=19, cex=0.4)
        if(!impute.only) lines(c(x.axis)[chrom==display.chrom], smoothX[chrom==display.chrom,display.sample], col=2)
      }
  } # End of if(show.plot)
  
  if(impute.only){
    return(X)
    } else {
      return(smoothX)
      }
} # End of function smoothseg()

