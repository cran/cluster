      subroutine daisy(nn,jpp,x,valmd,jtmd,jdat,vtype,ndyst,disv)
cc
cc  Calculating dissimilarities between objects or variables
cc
      implicit double precision(a-h,o-z)

      integer nn, jpp
cc          nn  = number of objects
cc          jpp = number of variables used for the calculations

cc  The following vectors and matrices must be dimensioned in the
cc  main program :
      double precision x(nn,jpp), disv(1+nn*(nn-1)/2)
      double precision valmd(jpp)
      integer jtmd(jpp), jdat
      integer vtype(jpp)
cc  vtype was character

cc    vtype(j) is the type of variable j:
cc        = 1 (A) for an Asymmetric binary variable
cc        = 2 (S) for a  Symmetric  binary variable
cc        = 3 (N) for a  Nominal  variable
cc        = 4 (O) for an Ordinal  variable
cc        = 5 (I) for an Interval variable (additive)
cc        = 6 (T) for a  raTio    variable
cc
cc    vector jtmd is only read if there are missing values
cc jtmd(j) =  0 if variable j is binary
cc         = -1 if variable j is not binary and has missing values
cc         = +1 if variable j is not binary and has no missing values

cc 
cc
cc   calculation of the dissimilarities
cc

      if(jdat .eq. 1) then

         nlk=1
         nbad=0
         disv(1)=0.0
         do 450 l=2,nn
            la=l-1
            do 440 k=1,la
               nlk=nlk+1
               ppa=0.
               dlk=0.
               do 420 j=1,jpp
                  if(vtype(j).eq.1.or.vtype(j).eq.2)go to 410
                  if(jtmd(j).ge.0)go to 405
                  if(x(l,j).eq.valmd(j))go to 420
                  if(x(k,j).eq.valmd(j))go to 420
 405              continue
                  ppa=ppa+1.
                  if(vtype(j).eq.3.and.x(l,j).ne.x(k,j))dlk=dlk+1.
                  if(vtype(j).eq.1.or.vtype(j).eq.2.or.vtype(j).eq.3)
     *                 go to 420
                  dlk=dlk+dabs(x(l,j)-x(k,j))
                  go to 420
 410              continue
                  if(x(l,j).ne.0..and.x(l,j).ne.1.)go to 420
                  if(x(k,j).ne.0..and.x(k,j).ne.1.)go to 420
                  if(vtype(j).eq.2.or.x(l,j).ne.0.or.x(k,j).ne.0)
     *                 ppa=ppa+1.
                  if(x(l,j).ne.x(k,j))dlk=dlk+1.
 420           continue
               if(ppa.gt.0.5)go to 430
               nbad=nbad+1
               disv(nlk)=-1
               go to 440
 430           continue
               disv(nlk)=dlk/ppa
 440        continue
 450     continue
         
      else
         
 500     pp=jpp
         nlk=1
         disv(1)=0.0
         do 600 l=2,nn
            lsubt=l-1
            do 520 k=1,lsubt
               clk=0.0
               nlk=nlk+1
               npres=0
               do 530 j=1,jpp
                  if(jtmd(j).ge.0)goto 540
                  if(x(l,j).eq.valmd(j))goto 530
                  if(x(k,j).eq.valmd(j))goto 530
 540              continue
                  npres=npres+1
                  if(ndyst.ne.1)goto 550
                  clk=clk+(x(l,j)-x(k,j))*(x(l,j)-x(k,j))
                  goto 530
 550              continue
                  clk=clk+dabs(x(l,j)-x(k,j))
 530           continue
               rpres=npres
               if(npres.ne.0)goto 560
               disv(nlk)=-1.0
               goto 520
 560           if(ndyst.ne.1)goto 570
               disv(nlk)=dsqrt(clk*(pp/rpres))
               goto 520
 570           continue
               disv(nlk)=clk*(pp/rpres)
 520        continue
 600     continue
      endif
      
      end

