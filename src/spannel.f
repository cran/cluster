C Compute the MVE -- Minimum Volume Ellipsoid --  
C --------------- for clusplot.default(*, span = TRUE)
C 
      subroutine spannel(ncas,ndep,dat,eps,dstopt,cov,
     +     varsum,varss,prob,work,ierr)
c
      implicit none
      integer ncas, ndep, ierr
C    ncas = number of objects.
C    ndep = number of variables.
C    dstopt = squared distances 
      double precision dat(ncas,0:ndep), eps, dstopt(ncas),
     2     cov(0:ndep,0:ndep),
     3     varsum(ndep), varss(ndep),
     4     prob(ncas), work(0:ndep)

C Local Variables
      double precision p, baswet, scal, scalfc, 
     +     aver, deter, dist, dmax, tempo
      integer i, j, k, loop

      do 50 i=1,ndep
         varsum(i)=0.d0
         varss(i)=0.d0
 50   continue

      ierr=0
      p=ndep

      baswet = 1. / dble(ncas)

      do 1 j = 1, ncas
         do 1 i = 1, ndep
            varsum(i) = varsum(i) + dat(j,i)
            varss(i) = varss(i) + dat(j,i) ** 2
 1    continue
      scalfc = 1
      do 10 j = 1, ndep
         aver = varsum(j) / ncas
         scal = dsqrt(varss(j) / ncas - aver * aver)
         scalfc = scalfc * scal
         do 10 i = 1, ncas
            dat(i,j) = (dat(i,j) - aver) / scal
 10   continue
      do 165 i = 1, ncas
         prob(i) = baswet
 165  continue

c ---- Repeat { ... up to 3000 times ]
      loop = 0
 160  continue
      loop = loop + 1
      do 300 j = 0,ndep
         do 300 i = 0,j
            cov(i,j) = 0.
 300  continue
      do 200 i = 1, ncas
         do 205 j = 0, ndep
            work(j) = dat(i,j)
            tempo = prob(i) * work(j)
            do 205 k = 0,j
               cov(k,j) = cov(k,j) + tempo * work(k)
 205        continue
 200  continue
      do 325 j = 0,ndep
         do 325 i = 0,j
            cov(j,i) = cov(i,j)
 325  continue
      deter = 1
    
      do 210 i = 0, ndep
         if (deter .le. 0.) then
            ierr=2
            return
         endif
         call sweep(cov,ndep,0,i,deter)
 210  continue
      
      dmax = 0.
      do 215 i = 1, ncas
         dist = - 1.
         do 220 j = 0, ndep
            work(j) = 0
            do 225 k = 0, ndep
               work(j) = work(j) - cov(j,k) * dat(i,k)
 225        continue
            dist = dist + work(j) * dat(i,j)
 220     continue
         dstopt(i) = dist
         if (dist .gt. dmax) dmax = dist
 215  continue
c     dmax = max{ dstopt[i] }

      if (dmax .gt. p+eps) then
c not yet converged
         do 230 i = 1, ncas
            prob(i) = prob(i) * dstopt(i) / p
 230     continue

        if (loop .lt. 5000) go to 160
      endif

      return
      end


      subroutine sweep (cov,nord,ixlo,nel,deter)

      implicit none
      integer nord,ixlo,nel
      double precision cov(0:nord,0:nord), deter
c
      double precision temp
      integer i,j
c
      temp=cov(nel,nel)
      deter=deter*temp
      if (nord.le.1) then
         cov(1,1)=1./temp
         return
      endif
c else : nord > 1      
      do 30 i=ixlo,nord
        if (i.eq.nel) go to 30
        do 20 j=ixlo,i
          if (j.eq.nel) go to 20
          cov(j,i)=cov(i,j)-cov(i,nel)*cov(nel,j)/temp
          cov(i,j)=cov(j,i)
20      continue
30    continue
      cov(nel,nel)=1.
      do 40 i=ixlo,nord
        cov(nel,i)=(-cov(i,nel))/temp
        cov(i,nel)=cov(nel,i)
40    continue
      return
      end
