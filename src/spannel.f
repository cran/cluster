C Compute the MVE -- Minimum Volume Ellipsoid --
C --------------- for clusplot.default(*, span = TRUE)
C
      subroutine spannel(ncas,ndep,dat,dstopt,cov,
     +     varsum,varss,prob,work, eps,maxit, ierr)
c
      implicit none
      integer ncas, ndep, maxit, ierr
C    ncas = number of objects.
C    ndep = number of variables.
C    maxit= maximal # iterations (and returns #{iter.})
C    dstopt = squared distances
      double precision dat(ncas,0:ndep), dstopt(ncas),
     2     cov(0:ndep,0:ndep), varsum(ndep), varss(ndep),
     4     prob(ncas), work(0:ndep), eps

C Local Variables
      double precision p, scal, aver, deter, dist, dmax, tempo
      integer i, j, k, loop

c (in clusplot's call,)  dat[i,0] are all  == 1

c Scale Data dat[i,j] to mean = 0 and var{1/n} = 1  -- for j= 1:ndep (not j=0!)
      do 1 j=1,ndep
         varsum(j)= 0.d0
         varss (j)= 0.d0
 1    continue
      do 5 i = 1, ncas
         do 5 j = 1, ndep
            varsum(j)= varsum(j) + dat(i,j)
            varss (j)= varss (j) + dat(i,j) ** 2
 5    continue

      do 10 j = 1, ndep
         aver = varsum(j) / ncas
         scal = dsqrt(varss(j) / ncas - aver * aver)
         do 10 i = 1, ncas
            dat(i,j) = (dat(i,j) - aver) / scal
 10   continue

      p = 1. / dble(ncas)
      do 165 i = 1, ncas
         prob(i) = p
 165  continue

      ierr=0
      p = ndep
c ---- Repeat { ... up to `maxit' times ]
      loop = 0
 160  continue
      loop = loop + 1
c     Cov[,] = weighted covariance of dat[,]  {weights = prob[]}
      do 300 j = 0,ndep
         do 300 k = 0,j
            cov(k,j) = 0.
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
         do 325 k = 0,j
            cov(j,k) = cov(k,j)
 325  continue

      deter = 1.
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
c           work(j) = - sum_{k=0}^p  dat(i,k) * cov(k,j) { = cov(j,k) },
c     i.e., work_j = - X[i,] %*% COV[,j]
            dist = dist + work(j) * dat(i,j)
 220     continue
         dstopt(i) = dist
c        Dist{opt}_i = -1 - t(X[i,]) %*% COV %*% X[i,]
         if (dist .gt. dmax) dmax = dist
 215  continue
c     dmax = max{ dstopt[i] }

      if (dmax .gt. p+eps) then
c     not yet converged
         do 230 i = 1, ncas
            prob(i) = prob(i) * dstopt(i) / p
 230     continue

        if (loop .lt. maxit) go to 160
      endif
      maxit = loop

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
