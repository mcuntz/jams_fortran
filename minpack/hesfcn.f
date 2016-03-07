      subroutine hesfcn(n,x,h,ldh,nprob)
      integer n,ldh,nprob
      double precision x(n),h(ldh,n)
c     **********
c
c     subroutine hesfcn
c
c     this subroutine defines the hessian matrices of eighteen
c     nonlinear unconstrained minimization problems. the problem
c     dimensions are as described in the prologue comments of objfcn.
c     the upper triangle of the (symmetric) hessian matrix is
c     computed columnwise. storage locations below the diagonal
c     are not disturbed until the final step, which reflects the
c     upper triangle to fill the square matrix.
c
c     the subroutine statement is
c
c       subroutine hesfcn(n,x,h,ldh,nprob)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c       h is an n by n array. on output h contains the hessian
c         matrix of the nprob objective function evaluated at x.
c
c       ldh is a positive integer input variable not less than n
c         which specifies the leading dimension of the array h.
c
c       nprob is a positive integer input variable which defines the
c         number of the problem. nprob must not exceed 18.
c
c     subprograms called
c
c       fortran-supplied ... dabs,datan,dcos,dexp,dlog,dsign,dsin,
c                            dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,iev,ivar,j,k,n2
      double precision ap,arg,cp0001,cp1,cp25,cp5,c1p5,c2p25,
     *                 c2p625,c3p5,c19p8,c25,c29,c100,c200,c10000,d1,
     *                 d2,eight,fifty,five,four,one,r,s1,s2,s3,t,t1,
     *                 t2,t3,ten,th,three,tpi,twenty,two,zero
      double precision d3,r1,r2,r3,u1,u2,v,v1,v2
      double precision fvec(50),fvec1(50),y(15)
      double precision dfloat
      double precision six,xnine,twelve,c120,c200p2,c202,c220p2,c360,
     *                 c400,c1200
      data six,xnine,twelve,c120,c200p2,c202,c220p2,c360,c400,c1200
     *     /6.0d0,9.0d0,1.2d1,1.2d2,2.002d2,2.02d2,2.202d2,3.6d2,
     *      4.0d2,1.2d3/
      data zero,one,two,three,four,five,eight,ten,twenty,fifty
     *     /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,8.0d0,1.0d1,2.0d1,
     *      5.0d1/
      data cp0001,cp1,cp25,cp5,c1p5,c2p25,c2p625,c3p5,c19p8,c25,c29,
     *     c100,c200,c10000
     *     /1.0d-4,1.0d-1,2.5d-1,5.0d-1,1.5d0,2.25d0,2.625d0,3.5d0,
     *      1.98d1,2.5d1,2.9d1,1.0d2,2.0d2,1.0d4/
      data ap /1.0d-5/
      data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),
     *     y(12),y(13),y(14),y(15)
     *     /9.0d-4,4.4d-3,1.75d-2,5.4d-2,1.295d-1,2.42d-1,3.521d-1,
     *      3.989d-1,3.521d-1,2.42d-1,1.295d-1,5.4d-2,1.75d-2,4.4d-3,
     *      9.0d-4/
      dfloat(ivar) = ivar
c
c     hessian routine selector.
c
      go to (10,20,60,100,110,170,210,290,330,380,390,450,490,580,620,
     *       660,670,680), nprob
c
c     helical valley function.
c
   10 continue
      tpi = eight*datan(one)
      th = dsign(cp25,x(2))
      if (x(1) .gt. zero) th = datan(x(2)/x(1))/tpi
      if (x(1) .lt. zero) th = datan(x(2)/x(1))/tpi + cp5
      arg = x(1)**2 + x(2)**2
      r = dsqrt(arg)
      t = x(3) - ten*th
      s1 = ten*t/(tpi*arg)
      t1 = ten/tpi
      t2 = t1/arg
      t3 = (x(1)/r - t1*t2*x(1) - two*x(2)*s1)/arg
      h(1,1) = c200
     *         *(one - x(2)/arg*(x(2)/r - t1*t2*x(2) + two*x(1)*s1))
      h(1,2) = c200*(s1 + x(2)*t3)
      h(2,2) = c200*(one - x(1)*t3)
      h(1,3) = c200*t2*x(2)
      h(2,3) = -c200*t2*x(1)
      h(3,3) = c202
      go to 800
c
c     biggs exp6 function.
c
   20 continue
      do 40 j = 1, 6
         do 30 i = 1, j
            h(i,j) = zero
   30       continue
   40    continue
      do 50 i = 1, 13
         d1 = dfloat(i)/ten
         d2 = dexp(-d1) - five*dexp(-ten*d1) + three*dexp(-four*d1)
         s1 = dexp(-d1*x(1))
         s2 = dexp(-d1*x(2))
         s3 = dexp(-d1*x(5))
         t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2
         th = d1*t
         r1 = d1*s1
         r2 = d1*s2
         r3 = d1*s3
         h(1,1) = h(1,1) + r1*(th + x(3)*r1)
         h(1,2) = h(1,2) - r1*r2
         h(2,2) = h(2,2) - r2*(th - x(4)*r2)
         h(1,3) = h(1,3) - s1*(th + x(3)*r1)
         h(3,3) = h(3,3) + s1**2
         h(1,4) = h(1,4) + r1*s2
         h(2,4) = h(2,4) + s2*(th - x(4)*r2)
         h(3,4) = h(3,4) - s1*s2
         h(4,4) = h(4,4) + s2**2
         h(1,5) = h(1,5) + r1*r3
         h(2,5) = h(2,5) - r2*r3
         h(5,5) = h(5,5) + r3*(th + x(6)*r3)
         h(1,6) = h(1,6) - r1*s3
         h(2,6) = h(2,6) + r2*s3
         h(3,6) = h(3,6) + s1*s3
         h(4,6) = h(4,6) - s2*s3
         h(5,6) = h(5,6) - s3*(th + x(6)*r3)
         h(6,6) = h(6,6) + s3**2
   50    continue
      h(1,1) = two*x(3)*h(1,1)
      h(1,2) = two*x(3)*x(4)*h(1,2)
      h(2,2) = two*x(4)*h(2,2)
      h(1,3) = two*h(1,3)
      h(2,3) = two*x(4)*h(1,4)
      h(3,3) = two*h(3,3)
      h(1,4) = two*x(3)*h(1,4)
      h(2,4) = two*h(2,4)
      h(3,4) = two*h(3,4)
      h(4,4) = two*h(4,4)
      h(1,5) = two*x(3)*x(6)*h(1,5)
      h(2,5) = two*x(4)*x(6)*h(2,5)
      h(3,5) = two*x(6)*h(1,6)
      h(4,5) = two*x(6)*h(2,6)
      h(5,5) = two*x(6)*h(5,5)
      h(1,6) = two*x(3)*h(1,6)
      h(2,6) = two*x(4)*h(2,6)
      h(3,6) = two*h(3,6)
      h(4,6) = two*h(4,6)
      h(5,6) = two*h(5,6)
      h(6,6) = two*h(6,6)
      go to 800
c
c     gaussian function.
c
   60 continue
      do 80 j = 1, 3
         do 70 i = 1, j
            h(i,j) = zero
   70       continue
   80    continue
      do 90 i = 1, 15
         d1 = cp5*dfloat(i-1)
         d2 = c3p5 - d1 - x(3)
         arg = -cp5*x(2)*d2**2
         r = dexp(arg)
         t = x(1)*r - y(i)
         s1 = r*t
         s2 = d2*s1
         t1 = s2 + d2*x(1)*r**2
         t2 = d2*t1
         h(1,1) = h(1,1) + r**2
         h(1,2) = h(1,2) - t2
         h(2,2) = h(2,2) + d2**2*t2
         h(1,3) = h(1,3) + t1
         h(2,3) = h(2,3) + two*s2 - d2*x(2)*t2
         h(3,3) = h(3,3) + x(2)*t2 - s1
   90    continue
      h(1,1) = two*h(1,1)
      h(1,2) = h(1,2)
      h(2,2) = cp5*x(1)*h(2,2)
      h(1,3) = two*x(2)*h(1,3)
      h(2,3) = x(1)*h(2,3)
      h(3,3) = two*x(1)*x(2)*h(3,3)
      go to 800
c
c     powell badly scaled function.
c
  100 continue
      t1 = c10000*x(1)*x(2) - one
      s1 = dexp(-x(1))
      s2 = dexp(-x(2))
      t2 = s1 + s2 - one - cp0001
      h(1,1) = two*((c10000*x(2))**2 + s1*(s1 + t2))
      h(1,2) = two*(c10000*(one + two*t1) + s1*s2)
      h(2,2) = two*((c10000*x(1))**2 + s2*(s2 + t2))
      go to 800
c
c     box 3-dimensional function.
c
  110 continue
      do 130 j = 1, 3
         do 120 i = 1, j
            h(i,j) = zero
  120       continue
  130    continue
      do 140 i = 1, 10
         d1 = dfloat(i)
         d2 = d1/ten
         s1 = dexp(-d2*x(1))
         s2 = dexp(-d2*x(2))
         s3 = dexp(-d2) - dexp(-d1)
         t = s1 - s2 - s3*x(3)
         th = d2*t
         r1 = d2*s1
         r2 = d2*s2
         h(1,1) = h(1,1) + r1*(th + r1)
         h(1,2) = h(1,2) - r1*r2
         h(2,2) = h(2,2) - r2*(th - r2)
         h(1,3) = h(1,3) + r1*s3
         h(2,3) = h(2,3) - r2*s3
         h(3,3) = h(3,3) + s3**2
  140    continue
      do 160 j = 1, 3
         do 150 i = 1, j
            h(i,j) = two*h(i,j)
  150       continue
  160    continue
      go to 800
c
c     variably dimensioned function.
c
  170 continue
      t1 = zero
      do 180 j = 1, n
         t1 = t1 + dfloat(j)*(x(j) - one)
  180    continue
c     t = t1*(one + two*t1**2)
      t2 = two + twelve*t1**2
      do 200 j = 1, n
         do 190 i = 1, j
            h(i,j) = dfloat(i*j)*t2
  190       continue
         h(j,j) = h(j,j) + two
  200    continue
      go to 800
c
c     watson function.
c
  210 continue
      do 230 j = 1, n
         do 220 k = 1, j
            h(k,j) = zero
  220       continue
  230    continue
      do 280 i = 1, 29
         d1 = dfloat(i)/c29
         s1 = zero
         d2 = one
         do 240 j = 2, n
            s1 = s1 + dfloat(j-1)*d2*x(j)
            d2 = d1*d2
  240       continue
         s2 = zero
         d2 = one
         do 250 j = 1, n
            s2 = s2 + d2*x(j)
            d2 = d1*d2
  250       continue
         t = s1 - s2**2 - one
         s3 = two*d1*s2
         d2 = two/d1
         th = two*d1**2*t
         do 270 j = 1, n
            v = dfloat(j-1) - s3
            d3 = one/d1
            do 260 k = 1, j
               h(k,j) = h(k,j) + d2*d3*(v*(dfloat(k-1) - s3) - th)
               d3 = d1*d3
  260          continue
            d2 = d1*d2
  270       continue
  280    continue
      t1 = x(2) - x(1)**2 - one
      h(1,1) = h(1,1) + eight*x(1)**2 + two - four*t1
      h(1,2) = h(1,2) - four*x(1)
      h(2,2) = h(2,2) + two
      go to 800
c
c     penalty function i.
c
  290 continue
      t1 = -cp25
      do 300 j = 1, n
         t1 = t1 + x(j)**2
  300    continue
      d1 = two*ap
      th = four*t1
      do 320 j = 1, n
         t2 = eight*x(j)
         do 310 i = 1, j
            h(i,j) = x(i)*t2
  310       continue
         h(j,j) = h(j,j) + d1 + th
  320    continue
      go to 800
c
c     penalty function ii.
c
  330 continue
      t1 = -one
      do 340 j = 1, n
         t1 = t1 + dfloat(n-j+1)*x(j)**2
  340    continue
      d1 = dexp(cp1)
      d2 = one
      th = four*t1
      do 370 j = 1, n
         t2 = eight*dfloat(n-j+1)*x(j)
         do 350 i = 1, j
            h(i,j) = dfloat(n-i+1)*x(i)*t2
  350       continue
         h(j,j) = h(j,j) + dfloat(n-j+1)*th
         s1 = dexp(x(j)/ten)
         if (j .eq. 1) go to 360
         s3 = s1 + s2 - d2*(d1 + one)
         h(j,j) = h(j,j) + ap*s1*(s3 + three*s1 - one/d1)/fifty
         h(j-1,j) = h(j-1,j) + ap*s1*s2/fifty
         h(j-1,j-1) = h(j-1,j-1) + ap*s2*(s2 + s3)/fifty
  360    continue
         s2 = s1
         d2 = d1*d2
  370    continue
      h(1,1) = h(1,1) + two
      go to 800
c
c     brown badly scaled function.
c
  380 continue
c     t1 = x(1) - c1pd6
c     t2 = x(2) - c2pdm6
      t3 = x(1)*x(2) - two
      h(1,1) = two*(one + x(2)**2)
      h(1,2) = four*(one + t3)
      h(2,2) = two*(one + x(1)**2)
      go to 800
c
c     brown and dennis function.
c
  390 continue
      do 410 j = 1, 4
         do 400 i = 1, j
            h(i,j) = zero
  400       continue
  410    continue
      do 420 i = 1, 20
         d1 = dfloat(i)/five
         d2 = dsin(d1)
         t1 = x(1) + d1*x(2) - dexp(d1)
         t2 = x(3) + d2*x(4) - dcos(d1)
         t = t1**2 + t2**2
c        s1 = t1*t
c        s2 = t2*t
         s3 = two*t1*t2
         r1 = t + two*t1**2
         r2 = t + two*t2**2
         h(1,1) = h(1,1) + r1
         h(1,2) = h(1,2) + d1*r1
         h(2,2) = h(2,2) + d1**2*r1
         h(1,3) = h(1,3) + s3
         h(2,3) = h(2,3) + d1*s3
         h(3,3) = h(3,3) + r2
         h(1,4) = h(1,4) + d2*s3
         h(2,4) = h(2,4) + d1*d2*s3
         h(3,4) = h(3,4) + d2*r2
         h(4,4) = h(4,4) + d2**2*r2
  420    continue
      do 440 j = 1, 4
         do 430 i = 1, j
            h(i,j) = four*h(i,j)
  430       continue
  440    continue
      go to 800
c
c     gulf research and development function.
c
  450 continue
      do 470 j = 1, 3
         do 460 i = 1, j
            h(i,j) = zero
  460       continue
  470    continue
      d1 = two/three
      do 480 i = 1, 99
         arg = dfloat(i)/c100
         r = (-fifty*dlog(arg))**d1 + c25 - x(2)
         t1 = dabs(r)**x(3)/x(1)
         t2 = dexp(-t1)
         t = t2 - arg
         s1 = t1*t2*t
         s2 = t1*(s1 + t2*(t1*t2 - t))
         r1 = dlog(dabs(r))
         r2 = r1*s2
         h(1,1) = h(1,1) + s2 - s1
         h(1,2) = h(1,2) + s2/r
         h(2,2) = h(2,2) + (s1 + x(3)*s2)/r**2
         h(1,3) = h(1,3) - r2
         h(2,3) = h(2,3) + (s1 - x(3)*r2)/r
         h(3,3) = h(3,3) + r1*r2
  480    continue
      h(1,1) = two*h(1,1)/x(1)**2
      h(1,2) = two*x(3)*h(1,2)/x(1)
      h(2,2) = two*x(3)*h(2,2)
      h(1,3) = two*h(1,3)/x(1)
      h(2,3) = two*h(2,3)
      h(3,3) = two*h(3,3)
      go to 800
c
c     trigonometric function.
c
  490 continue
      u2 = dcos(x(n))
      s1 = u2
      if (n .eq. 1) go to 510
      u1 = dcos(x(n-1))
      s1 = s1 + u1
      if (n .eq. 2) go to 510
      n2 = n - 2
      do 500 j = 1, n2
         h(j,n-1) = dcos(x(j))
         s1 = s1 + h(j,n-1)
  500    continue
  510 continue
      v2 = dsin(x(n))
      s2 = dfloat(2*n) - v2 - s1 - dfloat(n)*u2
      r2 = dfloat(2*n)*v2 - u2
      if (n .eq. 1) go to 570
      v1 = dsin(x(n-1))
      s2 = s2 + dfloat(2*n-1) - v1 - s1 - dfloat(n-1)*u1
      r1 = dfloat(2*n-1)*v1 - u1
      if (n .eq. 2) go to 560
      do 520 j = 1, n2
         h(j,n) = dsin(x(j))
         t = dfloat(n+j) - h(j,n) - s1 - dfloat(j)*h(j,n-1)
         s2 = s2 + t
  520    continue
      do 540 j = 1, n2
         v = dfloat(j)*h(j,n-1) + h(j,n)
         t = dfloat(n+j) - s1 - v
         t1 = dfloat(n+j)*h(j,n) - h(j,n-1)
         do 530 i = 1, j
            th = dfloat(i)*h(i,n) - h(i,n-1)
            h(i,j) = two*(h(i,n)*t1 + h(j,n)*th)
  530       continue
         h(j,j) = h(j,j) + two*(h(j,n-1)*s2 + v*t + th**2)
  540    continue
      do 550 i = 1, n2
         th = dfloat(i)*h(i,n) - h(i,n-1)
         h(i,n-1) = two*(h(i,n)*r1 + v1*th)
         h(i,n) = two*(h(i,n)*r2 + v2*th)
  550    continue
  560 continue
      v = dfloat(n-1)*u1 + v1
      t = dfloat(2*n-1) - s1 - v
      th = dfloat(n-1)*v1 - u1
      h(n-1,n-1) = two*(v1*(r1 + th) + u1*s2 + v*t + th**2)
      h(n-1,n) = two*(v1*r2 + v2*th)
  570 continue
      v = dfloat(n)*u2 + v2
      t = dfloat(2*n) - s1 - v
      th = dfloat(n)*v2 - u2
      h(n,n) = two*(v2*(r2 + th) + u2*s2 + v*t + th**2)
      go to 800
c
c     extended rosenbrock function.
c
  580 continue
      do 600 j = 1, n
         do 590 i = 1, j
            h(i,j) = zero
  590       continue
  600    continue
      do 610 j = 1, n, 2
c        t1 = one - x(j)
         h(j,j) = c1200*x(j)**2 - c400*x(j+1) + two
         h(j,j+1) = -c400*x(j)
         h(j+1,j+1) = c200
  610    continue
      go to 800
c
c     extended powell function.
c
  620 continue
      do 640 j = 1, n
         do 630 i = 1, j
            h(i,j) = zero
  630       continue
  640    continue
      do 650 j = 1, n, 4
c        t = x(j) + ten*x(j+1)
c        t1 = x(j+2) - x(j+3)
c        s1 = five*t1
         t2 = x(j+1) - two*x(j+2)
c        s2 = four*t2**3
         t3 = x(j) - x(j+3)
c        s3 = twenty*t3**3
         r2 = twelve*t2**2
         r3 = c120*t3**2
         h(j,j) = two + r3
         h(j,j+1) = twenty
         h(j+1,j+1) = c200 + r2
         h(j+1,j+2) = -two*r2
         h(j+2,j+2) = ten + four*r2
         h(j,j+3) = -r3
         h(j+2,j+3) = -ten
         h(j+3,j+3) = ten + r3
  650    continue
      go to 800
c
c     beale function.
c
  660 continue
      s1 = one - x(2)
      t1 = c1p5 - x(1)*s1
      s2 = one - x(2)**2
      t2 = c2p25 - x(1)*s2
      s3 = one - x(2)**3
      t3 = c2p625 - x(1)*s3
      h(1,1) = two*(s1**2 + s2**2 + s3**2)
      h(1,2) = two
     *         *(t1 + x(2)*(two*t2 + three*x(2)*t3)
     *           - x(1)*(s1 + x(2)*(two*s2 + three*x(2)*s3)))
      h(2,2) = two*x(1)
     *         *(x(1) + two*t2
     *           + x(2)*(six*t3 + x(1)*x(2)*(four + xnine*x(2)**2)))
      go to 800
c
c     wood function.
c
  670 continue
      s1 = x(2) - x(1)**2
c     s2 = one - x(1)
c     s3 = x(2) - one
      t1 = x(4) - x(3)**2
c     t2 = one - x(3)
c     t3 = x(4) - one
      h(1,1) = c400*(two*x(1)**2 - s1) + two
      h(1,2) = -c400*x(1)
      h(2,2) = c220p2
      h(1,3) = zero
      h(2,3) = zero
      h(3,3) = c360*(two*x(3)**2 - t1) + two
      h(1,4) = zero
      h(2,4) = c19p8
      h(3,4) = -c360*x(3)
      h(4,4) = c200p2
      go to 800
c
c     chebyquad function.
c
  680 continue
      do 690 i = 1, n
         fvec(i) = zero
  690    continue
      do 710 j = 1, n
         t1 = one
         t2 = two*x(j) - one
         t = two*t2
         do 700 i = 1, n
            fvec(i) = fvec(i) + t2
            th = t*t2 - t1
            t1 = t2
            t2 = th
  700       continue
  710    continue
      d1 = one/dfloat(n)
      iev = -1
      do 720 i = 1, n
         fvec(i) = d1*fvec(i)
         if (iev .gt. 0) fvec(i) = fvec(i) + one/(dfloat(i)**2 - one)
         iev = -iev
  720    continue
      do 770 j = 1, n
         do 730 k = 1, j
            h(k,j) = zero
  730       continue
         t1 = one
         t2 = two*x(j) - one
         t = two*t2
         s1 = zero
         s2 = two
         r1 = zero
         r2 = zero
         do 740 i = 1, n
            h(j,j) = h(j,j) + fvec(i)*r2
            th = eight*s2 + t*r2 - r1
            r1 = r2
            r2 = th
            fvec1(i) = d1*s2
            th = four*t2 + t*s2 - s1
            s1 = s2
            s2 = th
            th = t*t2 - t1
            t1 = t2
            t2 = th
  740       continue
         do 760 k = 1, j
            v1 = one
            v2 = two*x(k) - one
            v = two*v2
            u1 = zero
            u2 = two
            do 750 i = 1, n
               h(k,j) = h(k,j) + fvec1(i)*u2
               th = four*v2 + v*u2 - u1
               u1 = u2
               u2 = th
               th = v*v2 - v1
               v1 = v2
               v2 = th
  750          continue
  760       continue
  770    continue
      d2 = two*d1
      do 790 j = 1, n
         do 780 k = 1, j
            h(k,j) = d2*h(k,j)
  780       continue
  790    continue
c
c     reflect the upper triangle to fill the square matrix.
c
  800 continue
      do 820 j = 1, n
         do 810 i = j, n
            h(i,j) = h(j,i)
  810       continue
  820    continue
      return
c
c     last card of subroutine hesfcn.
c
      end

