      subroutine spannel(ncas,ndep,data,eps,dstopt,cov,
     +     varsum,varss,prob,work,nstop)
c
c
      implicit double precision (a-h,o-z)
      dimension data(ncas,0:ndep),dstopt(ncas),
     2 cov(0:ndep,0:ndep),
     4 varsum(ndep),varss(ndep),
     5 prob(ncas),work(0:ndep)
C
C    ncas = number of objects.
C    ndep = number of variables.
C    dstopt = squared distances 
C
      zero=0.d0
      one=1.d0

      do 50 i=1,ndep
         varsum(i)=0.d0
         varss(i)=0.d0
 50   continue

      nstop=0
      flp=ndep

      baswet = one / dble(ncas)

      do 1 j = 1, ncas
      do 1 i = 1, ndep
      varsum(i) = varsum(i) + data(j,i)
      varss(i) = varss(i) + data(j,i) ** 2
 1    continue
      scalfc = 1
      do 10 j = 1, ndep
      aver = varsum(j) / ncas
      scal = dsqrt(varss(j) / ncas - aver * aver)
      scalfc = scalfc * scal
      do 10 i = 1, ncas
      data(i,j) = (data(i,j) - aver) / scal
 10   continue
      do 165 i = 1, ncas
      prob(i) = baswet
 165  continue
      loop = 0
 160  continue
      loop = loop + 1
      do 300 j = 0,ndep
      do 300 i = 0,j
 300  cov(i,j) = zero
      do 200 i = 1, ncas
      indexa =i
      do 205 j = 0, ndep
      work(j) = data(indexa,j)
      tempo = prob(i) * work(j)
      do 205 k = 0,j
      cov(k,j) = cov(k,j) + tempo * work(k)
 205  continue
 200  continue
      do 325 j = 0,ndep
      do 325 i = 0,j
 325  cov(j,i) = cov(i,j)
      one = 1
      deter = 1
    
      do 210 i = 0, ndep
      if (deter .le. zero) then
        nstop=2
        return
        endif
 210  call sweep(cov,ndep,0,i,deter)
      
      dmax = zero
      do 215 i = 1, ncas
      dist = - 1
      indexa = i
      do 220 j = 0, ndep
      work(j) = 0
      do 225 k = 0, ndep
 225  work(j) = work(j) - cov(j,k) * data(indexa,k)
 220  dist = dist + work(j) * data(indexa,j)

      dstopt(i) = dist
      if (dist .gt. dmax) then
        dmax = dist
        endif
 215  continue
      if (dmax .gt. flp+eps) then
        do 230 i = 1, ncas
        prob(i) = prob(i) * dstopt(i) / flp
 230    continue
      

        if (loop .lt. 5000) go to 160
        endif
       
      return
      end

      SUBROUTINE SWEEP (COV,NORD,IXLO,NEL,DETER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COV(0:NORD,0:NORD)
      ONE=1.D0
      TEMP=COV(NEL,NEL)
      DETER=DETER*TEMP
      IF (NORD.GT.1) GO TO 10
      COV(1,1)=ONE/TEMP
      RETURN
10    DO 30 I=IXLO,NORD
        IF (I.EQ.NEL) GO TO 30
        DO 20 J=IXLO,I
          IF (J.EQ.NEL) GO TO 20
          COV(J,I)=COV(I,J)-COV(I,NEL)*COV(NEL,J)/TEMP
          COV(I,J)=COV(J,I)
20      CONTINUE
30    CONTINUE
      COV(NEL,NEL)=ONE
      DO 40 I=IXLO,NORD
        COV(NEL,I)=(-COV(I,NEL))/TEMP
        COV(I,NEL)=COV(NEL,I)
40    CONTINUE
      RETURN
      END



