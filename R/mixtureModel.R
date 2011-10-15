# is.ch <- function(x){
#     x <- c(x)
#     if(length(x)>1000){
#         x <- sample(x,1000)
#     }
#     return(length(unique(x[is.finite(x)])) <= 24)
# }
# 
suggest.k <- function(x, ch = NULL, np = NULL, smooth = FALSE, plot.it = FALSE, reg=FALSE,
ds=1.5, zero.cont=FALSE,  ...){
# This is the function to identify
# optimal number of mixture components (peaks)
# global normalisation for each sample pair
# Note: called in global.norm() function
# when opt=T
    x <- as.matrix(x)
    ncx <- ncol(x)
    nrx <- nrow(x)

    if (is.null(ch)){
        ch <- x[,1]
        x <- as.matrix(x[,-1])
        ncx <- ncx - 1
    }
    if(ncx==0){
        stop("No data on ratio. The information you input is considered only as chromosome.")
    }
    if(smooth & is.null(ch)){
    message("You require smooth-segmentation (smooth=T), but")
    message("did not provide chromosomal information.")
    stop("You can set ch=1 for genome-wide smooth segmentation.")
    }

    if(smooth){
        x <- smooth.it(x,ch)
    }

    if(is.null(np)){
        np <- seq(3, 10, by = 1)
    }

    lk <- length(np)
    if(lk <= 4){
        stop("Not enough number of values for np. Try more than four.")
    }
   
    original.ds <- ds
    suggested.k <- list()
    aic.list <- list()
    maxgr.list <- list()
    for(counter.col in 1:ncx){
        message("Iterating over k...")
        aic.vector <- vector()
        maxgr.vector <- vector()
        for(counter.k in np){
            message("k = ",counter.k)
             if(counter.k<=7) {
			ds=original.ds
                 } else {
                      if(counter.k>=8 & counter.k<=11){
			ds=0.8*original.ds 
                      } else {
			ds=c(2/3)*original.ds
                      }
                 }
            result <- pdetect(x[,counter.col], k=counter.k, reg=reg, ds=ds, zero.cont=zero.cont, ...)
            aic.vector <- c(aic.vector, result$aic)
            maxgr.vector <- c(maxgr.vector, result$maxgroup)
        } # end iteration on k
        if(plot.it){
            plot(np, aic.vector, type="l", xlab="No of peaks",
                ylab="AIC")
        }
        if(reg){
          result.k = get.optimal.k.reg(k=np, aic=aic.vector,maxgr=maxgr.vector)
        } else {
          result.k = get.optimal.k(k=np, aic=aic.vector,maxgr=maxgr.vector)
        }
        suggested.k <- c(suggested.k, list(result.k))
        aic.list <- c(aic.list, list(aic.vector))
        maxgr.list <- c(maxgr.list, list(maxgr.vector))
        if(length(result.k)==0){
            message("Sample in column ",counter.col," cannot find optimal number of peaks.")
            message("Try to increase argument fac to 0.01 to exclude extreme observations ")
            message("or test more number of peaks in np (default: 3-10). ")
            message("In addition, check the data for any suspicious pattern, e.g. outliers.")
        }
    }# end of multiple samples, iteration counter.col
    attr(suggested.k, "AIC") <- aic.list
    attr(suggested.k, "maxgr") <- maxgr.list
    return(suggested.k)
} # end of funtion suggest.k




pdetect <- function(x, k=NULL, init.p = NULL, init.mu=NULL, init.var=NULL,
niter = 50, fixed.mu = FALSE, n=2048, fac=0.001, adjust=0.9,
ds = 1.5, equal.var = FALSE, c1 = 0.001, reg = FALSE, zero.cont=FALSE, ...){

# note: k will overide the init.p, init.mu and init.var if they are NULL
# or their length is not the same as k
# the difference between this function and normal.EM is that
# in this function, we contrain the mu for estimation of
# contamination. note an extra argument "ds" as the minimum distance 
# in unit standard deviation between the components, with
# the second component as the reference

    old.x <- x

    id <- 1:length(x)
    id2 <- id[!is.na(x) & x != 0]
    x <- x[!is.na(x) & x != 0]

    if(zero.cont){
    zero.region <- x > -0.01 & x < 0.01
    zero.n <- sum(zero.region)
    err <- rnorm(zero.n,mean=0,sd=0.005)
    err[err < -0.02] <- -0.02
    err[err > 0.02] <- 0.02
    err <- scale(err,scale=F)
    x[zero.region] <- x[zero.region]+err
    } # end if(zero.cont)

    d <- density(x,n=n,adjust=adjust, ...)
    selected <- which(d$y >= fac*max(d$y))
    id3 <- id2[x >= min(d$x[selected], ra.rm=T) & x <= max(d$x[selected], ra.rm=T)]
    x <- x[x >= min(d$x[selected], ra.rm=T) & x <= max(d$x[selected], ra.rm=T)]

    # If k = 1, we need to stop this
    if(is.null(k)) k <- 7
    if(k<=1) stop("Value of k is invalid, k must at least be equal to 2.")

    mx <- mean(x, na.rm=T)
    vx <- var(x, na.rm=T)
    nx <- length(c(x))

    if(length(init.p) != k){
        if(!is.null(init.p)){
            warning("init.p is overriden (ignored), because its length differs from k.")
        }
        init.p <- rep(1/k,k)
    } # end if(length(init.p) != k)

    if(fixed.mu & length(init.mu) != k){
        stop("Since fixed.mu=T, the length of init.mu must be equal to k.")
    }
    if(length(init.mu) != k){
        if(!is.null(init.mu)){
            warning("init.mu is overriden (ignored), because its length differs from k.")
        }
         if(reg){
         res.init <- get.initial.values(x, k=k, fixed.mu=fixed.mu, ds=ds, equal.var=equal.var, c1=c1)
         init.mu <- res.init$init.mu
         } else {
         init.mu <- as.numeric(quantile(x, prob=c(seq(1/(k+1),1,by=1/(k+1)))[1:k]))
         }
    } # end if(length(init.mu) != k)

    if(length(init.var) != k){
        if(!is.null(init.var)){
            warning("init.var is overriden (ignored), because its length differs from k.")
        }
        init.var <- rep(vx/k,k)
        if(reg){
          init.var <- res.init$init.var
        } #end if(reg)
    } # end if(length(init.var) != k)

    # start iteration
    dmat <- matrix(1, nrow=nx, ncol=k)
    ymat <- matrix(1, nrow=nx, ncol=k)
    i <- 1
    while(i <= niter){
        #cat("iter: ",i," prop: ",round(init.p,2)," mu: ",round(init.mu,2),"\n")
        for(j in 1:k){
            dmat[,j] <- init.p[j] * dnorm(x,mean=init.mu[j], sd=sqrt(init.var[j]))
        } # end for j
        sdmat <- apply(dmat,1,sum)
        pij <- dmat/sdmat
        sum.pij <- apply(pij,2,sum)
        sum.pij.x <- apply(pij*x,2,sum)
        if(k<=5){
            init.p <- sum.pij/nx
        } else {
            init.p[1:5] <- c1+(1-5*c1)*c(sum.pij/nx)[1:5]
            init.p[6:k] <- (1-5*c1)*c(sum.pij/nx)[6:k]
        } # end if(k<=4)
        if(!fixed.mu){
            init.mu <- sum.pij.x/sum.pij
        }
        ord.init.mu <- order(init.mu)
        init.p <- init.p[ord.init.mu]
        init.mu <- init.mu[ord.init.mu]
        pij <- pij[,ord.init.mu]
        for(j in 1:k){
            ymat[,j] <- (x-init.mu[j])^2
        }
        init.var <- apply(pij*ymat,2,sum)/sum.pij
        if(equal.var){
            init.var <- rep(median(init.var),k)
        }
        med.sd <- median(sqrt(init.var))
        maxgr <- which.max(init.p)
        if(maxgr>5) {
            maxgr <- which.max(init.p[1:5])
        }
        if(maxgr==1){
            for(mm in c(maxgr+1):k){
                if(c(init.mu[mm]-init.mu[mm-1]) < ds*med.sd){
                   init.mu[mm] <- init.mu[mm-1] + ds*med.sd
                }
            } # for mm
        } # if maxgr==1
        if(k<=5) {
            if(maxgr==k){
                for(mm in c(maxgr-1):1){
                    if(c(init.mu[mm+1]-init.mu[mm])< ds*med.sd){
                        init.mu[mm] <- init.mu[mm+1] - ds*med.sd
                    }

                } # for mm
            } # if maxgr==k
        } # if(k<=5)
        if(maxgr>1 & maxgr<k){
            for(mm in c(maxgr-1):1){
                if(c(init.mu[mm+1]-init.mu[mm]) < ds*med.sd){
                    init.mu[mm] <- init.mu[mm+1] - ds*med.sd
                }
            } # for mm
            for(mm in c(maxgr+1):k){
                if(c(init.mu[mm]-init.mu[mm-1]) < ds*med.sd){
                   init.mu[mm] <- init.mu[mm-1] + ds*med.sd
                }
            } # for mm
        } # if maxgr>1 & maxgr<k

        i=i+1
    } # end iteration of EM

    ddmat <- matrix(nrow=length(x), ncol=k)
    for(counter in 1:k){
        ddmat[,counter] <- dnorm(x,mean=init.mu[counter],sd=sqrt(init.var[counter]))
    }
    ddmat2 <- t(t(ddmat)*init.p)
    like <- apply(ddmat2,1,sum)
    ll <- sum(log(like))
    if(equal.var){
        aic <- -2*ll+2*k 
    } else {
        aic <- -2*ll+4*k 
    }

    pij2 <- matrix(NA, nrow=length(old.x), ncol=k)  
    pij2[id3,] <- pij

    return(list(prop=init.p, mu=init.mu, var=init.var,
        weight=pij2, weight2= pij, maxgroup=which.max(init.p), 
        maxmu=init.mu[which.max(init.p)], fixed.mu=fixed.mu,
        c1=c1,aic=aic) )

} # end of function





get.optimal.k = function(k,aic,maxgr=NULL){
d1 = diff(aic)
d2 = d1/abs(d1)
d3 = diff(d2)
range.aic <- max(aic)-min(aic)
drop.aic <- c(10,d1/c(range.aic))

if(sum(d3==2)==0){
id.loc.min=which.min(aic)
} else {
id.loc.min=which(d3==2)+1
}
if(!is.null(maxgr)){
   maxgr.min= maxgr[id.loc.min]
   id.loc.min = id.loc.min[maxgr.min %in% c(3,5)]
 }
if(length(id.loc.min)>1)
  id.loc.min <- id.loc.min[which.min(drop.aic[id.loc.min])]
selected.k = k[id.loc.min]

if(length(selected.k)==0 & !is.null(maxgr)){
 id.loc.min = which(maxgr %in% c(3,5))
 if(length(id.loc.min)==1){
    selected.k <- k[id.loc.min]
  } else {
       if(length(id.loc.min)==2){
          selected.k <- which.min(aic[id.loc.min])
       } else {
            temp.loc.min <- id.loc.min[which.min(drop.aic[id.loc.min])]
            relative.aic <- (aic[id.loc.min]-aic[temp.loc.min])/range.aic
            temp.loc.min <- id.loc.min[abs(relative.aic)<0.1]
            selected.k <- k[temp.loc.min]
       }
   }
}# end of if selected.k
return(selected.k)
} # end of function

get.optimal.k.reg <- function(k,aic,maxgr=NULL){

d1 = diff(aic)
d2 = d1/abs(d1)
d3 = diff(d2)
range.aic <- max(aic)-min(aic)
drop.aic <- c(10,d1/c(range.aic))

if(is.null(maxgr)){
id.loc.min <- which.min(drop.aic)
selected.k = k[id.loc.min]
} else {
 id.loc.min = which(maxgr %in% c(3,5))
 if(length(id.loc.min)==1){
    selected.k <- k[id.loc.min]
  } else {
       if(length(id.loc.min)==2){
          selected.k <- which.min(aic[id.loc.min])
       } else {
	    temp.loc.min <- id.loc.min[which.min(aic[id.loc.min])]
            relative.aic <- (aic[id.loc.min]-aic[temp.loc.min])/range.aic
            temp.loc.min <- id.loc.min[abs(relative.aic)<0.1]
            selected.k <- k[temp.loc.min]
       }
   }
}# end of if selected.k
return(selected.k)
} # end of function



Rcheck <- function(a,addition=2) {
    # difference with Rcalc: 1 to 5 fixed
    if("mu" %in% names(a)) {
        a <- a$mu
    }
    la <- length(a)
    lb <- length(a)+addition
    if(la<4 | addition<2){
       stop("length of means should be >=4 and addition should be >=2.")
       }
    pa <- c(0:c(la-1))
    pb <- c(0:c(lb-1))
    Rvec <- vector()
    config.mat <- vector()
    for(j in 4:(lb-1)){
        for(k in (j+1):lb){
            x <- pb[-c(j,k)]
            config.mat <- rbind(config.mat,x)
            lm.temp <- lm(a ~ matrix(x))
            Rvec <- c(Rvec, summary(lm.temp)$r.squared)
        }  # end for k
    } # end for j
    mmm <- which.max(Rvec)
    return(list(r.squared=Rvec[mmm], ploidy=config.mat[mmm,],Rvec=Rvec, config.mat=config.mat))
} # end of function



global.norm <- function(x, ch=NULL, smooth = FALSE, k=NULL, Rt=0.95, maxgr=NULL, reg=FALSE,
ds=1.5, zero.cont=FALSE,  ...){
# function to perform global normalisation
# x should be the ratio data (not segmented ratio data)
# as we have the option to perform segmentation in this
# function
# this function can depend on smoothseg package
# if the argument smooth is set TRUE (which is the default)
# When this is set TRUE, the argument ch
# should be provided with a vector of chromosome id's

    x <- as.matrix(x)
    ncx <- ncol(x)
    nrx <- nrow(x)
    if(is.null(ch)){
        ch <- x[,1]
        x <- as.matrix(x[,-1])
        ncx <- ncx - 1
    }
    if(ncx==0){
        stop("No data on ratio. The information you input is considered only as chromosome.")
    }
    if(smooth & is.null(ch)){
        message("You require smooth-segmentation (smooth=T), but")
        message("did not provide chromosomal information.")
        stop("You can set ch=1 for genome-wide smooth segmentation.")
    }

    y <- matrix(nrow=nrx, ncol=ncx)
    z <- matrix(nrow=nrx, ncol=ncx)
    old.x <- x
    original.ds <- ds

    if(smooth){
        x <- smooth.it(x, ch=ch)
    } # end of if smooth

    if(is.null(k)){
      sk = suggest.k(x=x, ch=ch, np = NULL, smooth=smooth, plot.it=F, reg=reg, ds=original.ds,
        zero.cont=zero.cont, ...)
        k <- vector()
        for(j in 1:ncx){
            tempk = sk[[j]]
            if(length(tempk)==1){
                k[j]=tempk
            } else {
                tempk2 = abs(tempk-11)
                k[j]=tempk[which.min(tempk2)]
            } # End if tempk
        } # end for
    } else {
        if(length(k) != ncx){
            stop("Length of k is not the same as the number of columns in x.")
        }
    } # end if null k

    if(is.null(maxgr)){
        maxgr <- vector() # max group, for each sample
        update.maxgr <- TRUE
    } else {
        if(length(maxgr)==1){
            maxgr <- rep(maxgr,ncx)
        }
        if(length(maxgr)!=ncx) {
            stop("Length of maxgr is not the same as number of column in x.")
        }
        update.maxgr <- FALSE
    } # end if maxgr is null

    alpha <- vector() # contamination level, for each sample
    norm.factor <- vector() # mean of mixture component
                        # corresponding to maxgr-th group
    muest1 <- list()  # mean of mixture components (raw)
    muest2 <- list()  # mean of mixture components (after ratio set to 1)
    pl <- list() # ploidy


    for(i in 1:ncx){  # start iteration on each column

        # peak detection
        tempx <- x[,i]
        if(k[i]<=7){
          ds=original.ds
        } else {
	  if(k[i]>=8 & k[i]<=11) {
             ds=0.8*original.ds
             }  else {
		  ds=c(2/3)*original.ds 
                    }
        }
        res.em <- pdetect(tempx, k=k[i], reg=reg, ds=ds, zero.cont=zero.cont, ...)
        if(update.maxgr){
            maxgr <- c(maxgr, res.em$maxgroup)
        }
        muest1 <- c(muest1, list(res.em$mu))

        # ploidy check and normalisation: scaling to have the mean at one
        ploidy <- c(0:c(k[i]-1))
        lm.mu <- lm(c(res.em$mu) ~ ploidy)
        Rsq <- summary(lm.mu)$r.squared
        if(Rsq < Rt) {
            res.Rcheck <- Rcheck(res.em)
            ploidy <- res.Rcheck$ploidy
        } # end if
        pl <- c(pl, list(ploidy))
        lm.mu <- lm(c(res.em$mu) ~ ploidy)
        ppmu <- predict(lm.mu)
        if(update.maxgr){
            nf <- ppmu[res.em$maxgroup]
        } else {
        nf <- ppmu[c(maxgr[i])]
        }
        norm.factor <- c(norm.factor, nf)
  
        y[,i] <- old.x[,i]/nf # scaling to have ratio one
        norm.mu <- res.em$mu/nf
        muest2 <- c(muest2, list(norm.mu))

        #estimating contamination
        if(update.maxgr){
            temp.fac <- ploidy/ploidy[res.em$maxgroup]
        } else {
            temp.fac <- ploidy/ploidy[c(maxgr[i])]
        }
        dist.fac <- temp.fac[2]-temp.fac[1]
        dist.fac2 <- 0.5
        dist.mu <- abs(norm.mu[3]-norm.mu)
        dist.ploidy <- abs(ploidy - ploidy[3])
        contamination <- mean(1 - c(dist.mu/c(dist.fac2*dist.ploidy*norm.mu[3])), na.rm=T)
        alpha <- c(alpha, contamination)

        #  Normalising based on the contamination
        lm.mu <- lm(c(norm.mu) ~ ploidy)
        z[,i] <- (y[,i]-1)*(dist.fac/lm.mu$coef[2]) + 1

    } # end of global normalisation for each column i

    fin.res <- list(norm.data=z, scaled.data =y, alpha=alpha, maxgr=maxgr, norm.factor=norm.factor, k=k, smooth=smooth, muest1=muest1, muest2=muest2, pl=pl)
    class(fin.res) <- c("list", "gnorm")
    return(fin.res)

}# end of global normalisation function


get.initial.values <- function(x, k=NULL, init.p = NULL, init.mu=NULL, init.var=NULL,
tole =1e-5, niter=10, fixed.mu=FALSE, ds=1.5, equal.var=F, c1=0.001)
{
# This function generates initial values
# when reg=T in the pdetect().
# This function is called from pdetect()
# and is calling the function pdetect.iter()
# We check different initial values and
# select the best initial values after min. AIC
# in the first 10 iteration

#cat("In get.initial.values() \n")
range.x <- range(x, na.rm=T)
s2.x <- var(x, na.rm=T)
mean.x <- mean(x, na.rm=T)

final.p <- vector()
final.mu <- vector()
final.var <- vector()


aic.value <- 1E20
# Finch initial values, equal proportion, choice 1
         init.p <- rep(1/k,k)
         init.mu <- as.numeric(quantile(x, prob=c(seq(1/(2*k),by=1/k, length=k))))
         ss.init.mu <- sum((init.mu-mean.x)^2)/k
         init.var <- rep(s2.x-ss.init.mu,k)
         if(any(init.var < 0.01*s2.x/k)) init.var <- rep(s2.x/k,k)
         aic.temp <- pdetect.iter(x, k=k, init.p = init.p, init.mu=init.mu, init.var=init.var,
         niter=10, fixed.mu=fixed.mu, ds=ds, equal.var=equal.var, c1=c1)$aic
#         cat("In get.initial.values(). after Finch \n")
         if(aic.temp<aic.value){
		final.p <- init.p
		final.mu <- init.mu
		final.var <- init.var
		aic.value <- aic.temp
		choice <- "Finch"
         }
# our standard initial values ,choice 2
         init.mu <- as.numeric(quantile(x, prob=c(seq(1/(k+1),1,by=1/(k+1)))[1:k]))
         init.var <- rep(s2.x/k,k)
         aic.temp <- pdetect.iter(x, k=k, init.p = init.p, init.mu=init.mu, init.var=init.var,
         niter=10, fixed.mu=fixed.mu, ds=ds, equal.var=equal.var, c1=c1)$aic
#         cat("In get.initial.values(). after our standard \n")
         if(aic.temp<aic.value){
		final.p <- init.p
		final.mu <- init.mu
		final.var <- init.var
		aic.value <- aic.temp
		choice <- "Ourstandard"
         }
# moment matching method
	 if(k%%2==0){
                k2 <- c(-rev(1:(k/2)),(1:(k/2)))
	 } else {
                k2 <- c(-rev(1:((k-1)/2)),0,(1:((k-1)/2)))
	 }
	 init.mu <- mean.x+k2*sqrt(s2.x)/2
         ss.init.mu <- sum((init.mu-mean.x)^2)/k
         init.var <- rep(s2.x-ss.init.mu,k)
         if(any(init.var<0.01*s2.x/k)) init.var <- rep(s2.x/k,k)
#         cat("init.var=",init.var,"\n")
         aic.temp <- pdetect.iter(x, k=k, init.p = init.p, init.mu=init.mu, init.var=init.var,
         niter=10, fixed.mu=fixed.mu, ds=ds, equal.var=equal.var, c1=c1)$aic
#         cat("In get.initial.values(). after moment matching \n")
         if(aic.temp<aic.value){
		final.p <- init.p
		final.mu <- init.mu
		final.var <- init.var
		aic.value <- aic.temp
		choice <- "Momentmatching"
         }

#equal range
         init.mu <- as.numeric(quantile(range.x, prob=c(seq(1/(2*k),by=1/k, length=k))))
         init.var <- rep(s2.x/k,k)
         aic.temp <- pdetect.iter(x, k=k, init.p = init.p, init.mu=init.mu, init.var=init.var,
         niter=10, fixed.mu=fixed.mu, ds=ds, equal.var=equal.var, c1=c1)$aic
#         cat("In get.initial.values(). after equal range \n")
         if(aic.temp<aic.value){
		final.p <- init.p
		final.mu <- init.mu
		final.var <- init.var
		aic.value <- aic.temp
		choice <- "Equalrange"
         }
#message(choice)
return(list(init.p=final.p, init.mu=final.mu, init.var=final.var))

} # end of get.initial.values() function


pdetect.iter <- function(x, k=NULL, init.p = NULL, init.mu=NULL, init.var=NULL,
niter=10, fixed.mu=FALSE, ds=1.5, equal.var=F, c1=0.001)
{

# for iteration checks

range.x <- range(x, na.rm=T)
mx <- mean(x, na.rm=T)
vx <- var(x, na.rm=T)
nx <- length(c(x))


# start iteration
dmat <- matrix(1, nrow=nx, ncol=k)
ymat <- matrix(1, nrow=nx, ncol=k)
i <- 1
while(i <= niter){
  #cat("iter: ",i," prop: ",round(init.p,2)," mu: ",round(init.mu,2),"\n")
  for(j in 1:k){
  dmat[,j] <- init.p[j] * dnorm(x,mean=init.mu[j], sd=sqrt(init.var[j]))
  } # end for j
  sdmat <- apply(dmat,1,sum)
  pij <- dmat/sdmat
  sum.pij <- apply(pij,2,sum)
  sum.pij.x <- apply(pij*x,2,sum)
  if(k<=5){
     init.p <- sum.pij/nx
     } else {
        init.p[1:5]=c1+(1-5*c1)*c(sum.pij/nx)[1:5]
        init.p[6:k]=(1-5*c1)*c(sum.pij/nx)[6:k]
     } # end if(k<=4)
  if(!fixed.mu) init.mu <- sum.pij.x/sum.pij
  ord.init.mu <- order(init.mu)
  init.p <- init.p[ord.init.mu]
  init.mu <- init.mu[ord.init.mu]
  pij <- pij[,ord.init.mu]
  for(j in 1:k){
   ymat[,j] <- (x-init.mu[j])^2
  }
  init.var <- apply(pij*ymat,2,sum)/sum.pij
  if(equal.var) init.var=rep(median(init.var),k)
  med.sd <- median(sqrt(init.var))
  maxgr <- which.max(init.p)
  if(maxgr>5) maxgr <- which.max(init.p[1:5])

    if(maxgr==1){
            for(mm in c(maxgr+1):k){
                   if(c(init.mu[mm]-init.mu[mm-1]) < ds*med.sd)
                   init.mu[mm] <- init.mu[mm-1] + ds*med.sd
                } # for mm
           } # if maxgr==1
    if(k<=5) {
      if(maxgr==k){
            for(mm in c(maxgr-1):1){
                   if(c(init.mu[mm+1]-init.mu[mm])< ds*med.sd)
                   init.mu[mm] <- init.mu[mm+1] - ds*med.sd
                } # for mm
           } # if maxgr==k
           } # if(k<=4)
    if(maxgr>1 & maxgr<k){
            for(mm in c(maxgr-1):1){
                   if(c(init.mu[mm+1]-init.mu[mm])< ds*med.sd)
                   init.mu[mm] <- init.mu[mm+1] - ds*med.sd
                } # for mm
            for(mm in c(maxgr+1):k){
                   if(c(init.mu[mm]-init.mu[mm-1]) < ds*med.sd)
                   init.mu[mm] <- init.mu[mm-1] + ds*med.sd
                } # for mm
           } # if maxgr>1 & maxgr<k

  i=i+1
} # end iteration of EM

ddmat = matrix(nrow=length(x), ncol=k)
for(counter in 1:k){
ddmat[,counter] = dnorm(x,mean=init.mu[counter],sd=sqrt(init.var[counter]))
}
ddmat2 = t(t(ddmat)*init.p)
like = apply(ddmat2,1,sum)
ll = sum(log(like))
if(equal.var){
  aic = -2*ll+2*k 
   } else {
      aic = -2*ll+4*k 
      }

return(list(prop=init.p, mu=init.mu, var=init.var,
maxgroup=which.max(init.p), 
maxmu=init.mu[which.max(init.p)], fixed.mu=fixed.mu,
c1=c1,aic=aic))

} # end of function



