#  robust 'segmentation', based on Cauchy structure for the 
#  underlying pattern, and allowing for outliers in the measurements
#
# x,y are paired data
#
rseg = function(x,y,bw=1, maxiter=3, k=2, lambda = 1)
{
 nx = length(x)
 np = floor(nx/bw) # number of function evaluations
 x= x[1:(np*bw)]
 y = y[1:(np*bw)]
 # mean of x inside window
 mx = matrix(x,ncol=bw, byrow=T) %*% rep(1/bw, bw)  
 m = median(y)

 old= srsmooth(y) -m   # starting value using a simple smoother
                       # max number of points = 20000
 sig = median(abs(y-m-old))   # scale param

ERR =NULL
lambda.n = lambda*nx^(.75)  # make lambda depend on sample size
for (i in 1:maxiter){
   d = (y-m-old)/sig
   wy =(k+1)/(k+ d^2)    # robust weight see Pawitan page 354
  ZWZ = c(matrix(wy,ncol=bw, byrow=T) %*% rep(1,bw))
  avew = mean(ZWZ)
  lambda.w= avew * lambda.n   # avew makes it scale invariant
  ZWy  = c(matrix((y-m)*wy,ncol=bw, byrow=T) %*% rep(1,bw))
  new = bandseg(old, 1, ZWZ, ZWy,lambda.w)
  err = mean((old-new)^2)
  ERR = c(ERR,err)
  old = new
}
return(list(x=mx, y= m+new, err = ERR))
}

## banded-matrix inversion
## b     = current estimate, used for weights
## sigma = scale parameter, set to 1
bandseg = function(b, sigma, ZWZ,ZWy, lambda){
n = length(ZWZ)
#  ..............................  band-limited form
d = 2
a = diff(b, diff=d)            # difference pattern
D = 1/(a^2 + sigma^2)
### diagonal and off-diagonal terms terms: see band-trial.r 
#rb0 = convolve(c(1,4,1),rev(D),type='o')   
#rb1 = convolve(c(-2,-2),rev(D),type='o')  
a = c(filter(D, c(1,4,1)))
  nd = length(D)
  na = length(a)
  rb0= c(D[1], 4*D[1]+D[2], a[-c(1,na)], D[nd-1]+4*D[nd], D[nd])
a = c(filter(D,c(-2,-2)))
  rb1= c(-2*D[1], a[-length(a)], -2*D[length(D)])
rb2 = D

# ...............  band storage for Rinv    
ml = 2
mu = 2
m = ml + mu + 1   # diagonal position
R = matrix(0, nrow=2*ml+mu+1, ncol=n)      
  R[m,] = rb0
  R[7,1:(n-2)] = R[3,3:n] = rb2
  R[6,1:(n-1)] = R[4,2:n] = rb1

# working matrices  
ZWZR = lambda* R
  ZWZR[m,] = ZWZR[m,] + ZWZ

lda = 2*ml + mu + 1
ipvt = rep(0,n)
job=0
info=0
call1 = .Fortran("dgbfa", 
            abd = as.double(ZWZR),
            as.integer(lda),
            as.integer(n),
            as.integer(ml),
            as.integer(mu),
            ipvt=as.integer(ipvt),
            as.integer(info))
            # ,
            # PACKAGE="CNAnorm")
call2 = .Fortran("dgbsl", 
            abd = as.double(call1$abd),
            as.integer(lda),
            as.integer(n),
            as.integer(ml),
            as.integer(mu),
            ipvt=as.integer(call1$ipvt),
            b=as.double(ZWy),
            as.integer(job))
            # ,
            # PACKAGE="CNAnorm")

return(call2$b)
}


brute.seg = function(x,y,bw=1, maxiter=5, k=2, nx = length(x), lambda = sqrt(nx))
{
 nx = length(x)
 np = floor(nx/bw) # number of function evaluations
 x= x[1:(np*bw)]
 y = y[1:(np*bw)]
 # mean of x inside window
 mx = matrix(x,ncol=bw, byrow=T) %*% rep(1/bw, bw)  
 m = median(y)

 old= srsmooth(y) -m   # starting value using a simple smoother
                       # max number of points = 20000
 sig = median(abs(y-m-old))   # scale param
ERR =NULL
for (i in 1:5){
   d = (y-m-old)/sig
   wy =(k+1)/(k+ d^2)    # robust weight see Pawitan page 354
  ZWZ = c(matrix(wy,ncol=bw, byrow=T) %*% rep(1,bw))
  avew = mean(ZWZ)
  lambda= avew * sqrt(nx)   # avew makes it scale invariant
  ZWy  = c(matrix((y-m)*wy,ncol=bw, byrow=T) %*% rep(1,bw))
  new = brute.solve(old, 1, ZWZ, ZWy,lambda)
  err = mean((old-new)^2)
  ERR = c(ERR,err)
  old = new
}
return(list(x=mx, y= m+new, err = err))
}


# segmentation 
brute.solve = function(b, sigma, ZWZ,ZWy, lambda,d=2){
 n = length(ZWZ)
 ZWZR = diag(ZWZ) + lambda* Rb(b, sigma, d=d) 
  b = solve(ZWZR, ZWy)
  return(b)
}



# R inverse matrix; see Pawitan Chapter 18, page 481
R = function(n, d=2){
d0 = diag(n); 
  d1 = diff(d0, diff=d)
  R = t(d1) %*% d1
return(R)
}

# Weighted R inverse matrix; see Pawitan, page 466 and 501
# b is the current estimate
# sigma is the scale parameter of the differences
Rb = function(b, sigma, d=2){
n = length(b)
d0 = diag(n); 
  d1 = diff(d0, diff=d)
a = diff(b, diff=d)            # difference pattern
D = 1/(a^2 + sigma^2)
R = t(d1 * D) %*% d1
  return(R)
}

bsolve00 = function(ZWZ,ZWy, lambda){
 n = length(ZWZ)
 ZWZR = diag(ZWZ) + lambda* R(n) 
  b = solve(ZWZR, ZWy)
  return(b)
}

# ...................................  simple robust smooth
srsmooth = function(x, maxiter=2){
 n = length(x)
 xs = rep(0,n)
 isplit= rep(0,n)
 diff = rep(0,n)
 call1 = .Fortran('smooth',
           as.integer(n),
           x = as.single(x),
           as.integer(maxiter),
           as.single(xs),
           as.integer(isplit),
           as.single(diff))
           # , 
           # PACKAGE="CNAnorm")
 return(call1$x)
}



