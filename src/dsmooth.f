c	parameter (n=15)
c	dimension x(n)
c	data x /4,7,9,3,4,11,12,99,10,15,12,13,17,20,24/
c	write(6,'(15f5.1)')x
c	call smooth(n,x)
c	write(6,'(15f5.1)')x
c	end


c
c simple robust smoothing of a sequence
c maximum number of points = 20000
c
	subroutine smooth(n,x, maxiter, xs,isplit,diff)
	dimension x(n), xs(n), isplit(n), diff(n) 
	do 10 i=1, maxiter
	  call med3r(n,x,10,xs)
	  call split(n,x, xs,isplit,diff)
	  call split(n,x, xs,isplit,diff)
	  call hann(n,x, xs)
	  call hann(n,x,xs)
	  call hann(n,x,xs)
10	continue
	call med3r(n,x,1,xs)
	return
	end
	
	
	subroutine med3r(n,x,maxiter, oldx)
	dimension x(n),oldx(n)
	iter =0
1	continue
	iter = iter+1
	do 5 i=1,n
	  oldx(i) = x(i)
5 	continue
	do 10 i=2,n-1
	  x(i) = amed3(oldx(i-1),oldx(i),oldx(i+1))
10	continue
	err=0.
	do 20 i=1,n
	  err = err + abs(x(i) - oldx(i))
20 	continue
	if (err/n .gt. .01 .and. iter .lt. maxiter) goto 1
	x(1) = amed3(x(1),x(2),(3*x(2) - 2*x(3)))
	x(n) = amed3(x(n),x(n-1),(3*x(n-1) - 2*x(n-1)))
	return
	end

	function amed3(a,b,c)
	 if (a .le. b .and. b .le. c) amed3=b
	 if (a .le. c .and. c .le. b) amed3=c
	 if (b .le. c .and. c .le. a) amed3=c
	 if (b .le. a .and. a .le. c) amed3=a
	 if (c .le. a .and. a .le. b) amed3=a
	 if (c .le. b .and. b .le. a) amed3=b
	return
	end
	
	subroutine split(n,x, xs,isplit,diff)
	dimension x(n),xs(n),isplit(n),diff(n)
	do 10 i=1,n-1	
	  diff(i) = x(i+1) - x(i)
10	continue
	do 15 i=1,n-1
	  isplit(i)=0
15	continue
	do 20 i=3,n-3
	  if (abs(diff(i)) .lt. .000001 .and.
     &	      diff(i-1)*diff(i+1) .lt. 0.) isplit(i)=1
20	continue
	do 30 i=1,n
	  xs(i) = x(i)
30	continue
	do 40 i=1,n-1
	 if (isplit(i) .eq. 1) then
	   xs(i)   =  amed3((3*x(i-1) - 2*x(i-2)),x(i),x(i-1))
	   xs(i+1) =  amed3((3*x(i+2) - 2*x(i+3)),x(i+1),x(i+2))
	 endif
40	continue
	do 50 i=1,n
	  x(i) = xs(i)
50	continue
	return
	end

	subroutine hann(n,x, xs)
	dimension x(n),xs(n)
	do 10 i=2,n-1
	  xs(i) = (x(i-1) + 2*x(i) + x(i+1))/4.
10	continue
	do 20 i=2,n-1
	  x(i) = xs(i)
20	continue
	return
	end

