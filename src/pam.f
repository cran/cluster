c $Id: pam.f,v 1.13 2002/10/28 21:25:37 maechler Exp $
c
c PAM := Partitioning Around Medoids
c
      subroutine pam(nn,jpp,kk,x,dys,jdyss,valmd,jtmd,ndyst,nsend,
     f	   nrepr,nelem,radus,damer,ttd,separ,ttsyl,med,obj,ncluv,
     f	   clusinf,sylinf,nisol)
c
c     carries out a clustering using the k-medoid approach.
c
      integer nn,jpp,kk, jdyss,ndyst

      double precision x(nn,jpp), dys(1+nn*(nn-1)/2), valmd(jpp)
      integer jtmd(jpp), nsend(nn), nrepr(nn), nelem(nn),
     +	   ncluv(nn), nisol(kk), med(kk)
      double precision radus(nn),damer(nn),ttd(nn),separ(nn),ttsyl
      double precision clusinf(kk,5),sylinf(nn,4),obj(2)
c
      integer jhalt, nhalf, k,l
      double precision sky,s

c jdyss = 0 : compute distances from x
c	= 1 : distances provided  in x
      if(jdyss.ne.1) then
	 jhalt=0
	 call dysta(nn,jpp,x,dys,ndyst,jtmd,valmd,jhalt)
	 if (jhalt.ne.0) then
	    jdyss=-1
	    return
	 endif
      endif
c     nhalf := #{distances} = length(dys)
      nhalf=nn*(nn-1)/2+1
c     s := max( dys[.] ), the largest distance
      s=0.0
      do 10 l=2,nhalf
	 if(s .lt. dys(l)) s=dys(l)
 10   continue

c     Build + Swap :
      call bswap(kk,nn,	     nrepr,radus,damer,ttd, nhalf,dys,sky,s,obj)
c     Compute STATs :
      call cstat(kk,nn,nsend,nrepr,radus,damer,ttd, separ,   sky,s,
     f	   nhalf,dys,ncluv,nelem,med,nisol)
      do 135 k=1,kk
	 clusinf(k,1)=nrepr(k)
	 clusinf(k,2)=radus(k)
	 clusinf(k,3)=ttd(k)
	 clusinf(k,4)=damer(k)
	 clusinf(k,5)=separ(k)
 135  continue
      if(1 .lt. kk .and. kk .lt. nn) then
c	 Compute Silhouette info :
	 call dark(kk,nn,nhalf,ncluv,nsend,nelem,nrepr,
     f	      radus,damer,ttd,ttsyl,dys,s,sylinf)
      endif
      end
c     -----------------------------------------------------------

c     Compute Distances from X matrix
c
      subroutine dysta(nn,jpp,x,dys,ndyst,jtmd,valmd,jhalt)

      integer nn, jpp, ndyst, jtmd(jpp), jhalt
      double precision x(nn,jpp), dys(1+nn*(nn-1)/2), valmd(jpp)
c VARs
      integer nlk,j,l,k, lsubt, npres
      double precision pp, clk, rpres

      pp=jpp
      nlk=1
      dys(1)=0.0
      do 100 l=2,nn
	 lsubt=l-1
	 do 20 k=1,lsubt
	    clk=0.0
	    nlk=nlk+1
	    npres=0
	    do 30 j=1,jpp
	       if(jtmd(j).lt.0) then
		  if(x(l,j).eq.valmd(j))goto 30
		  if(x(k,j).eq.valmd(j))goto 30
	       endif
	       npres=npres+1
	       if(ndyst.eq.1) then
		  clk=clk+ (x(l,j)-x(k,j))*(x(l,j)-x(k,j))
	       else
		  clk=clk+ dabs(x(l,j)-x(k,j))
	       endif
 30	    continue
	    rpres=npres
	    if(npres.eq.0) then
	       jhalt=1
	       dys(nlk)=-1.0
	    else
	       if(ndyst.eq.1) then
		  dys(nlk)= dsqrt(clk*(pp/rpres))
	       else
		  dys(nlk)= clk*(pp/rpres)
	       endif
	    endif
 20	 continue
 100  continue
      end
c     -----------------------------------------------------------
c
c     bswap(): the clustering algorithm in 2 parts:  I. build,	II. swap
c
      subroutine bswap(kk,nn,nrepr,dysma,dysmb,beter,hh,dys,sky,s,obj)

      integer kk,nn, nrepr(nn), hh
c     nrepr[]: here is boolean (0/1): 1 = "is representative object"
      double precision dysma(nn),dysmb(nn),beter(nn),dys(hh), sky,s,
     +	   obj(2)
c     function called
      integer meet
      external meet
c
      integer i,j,k,kj, njn,ij,nmax, nbest,kbest
      double precision cmd, ammax, dz,dzsky,small
c     -Wall:
      nbest=-1
      kbest=-1
c
c     first algorithm: build.
c
      do 17 i=1,nn
	 nrepr(i)=0
	 dysma(i)=1.1*s+1.0
 17   continue

      do 20 k=1,kk
	 do 22 i=1,nn
	    if(nrepr(i) .eq. 0) then
	       beter(i)=0.
	       do 21 j=1,nn
		  ij=meet(i,j)
		  cmd=dysma(j)-dys(ij)
		  if(cmd.gt.0.0) beter(i)=beter(i)+cmd
 21	       continue
	    endif
 22	 continue
	 ammax=0.
	 do 31 i=1,nn
	    if(nrepr(i).eq.0 .and. ammax .le. beter(i)) then
c		    does .lt. (instead of .le.) work too? -- NO!
	       ammax=beter(i)
	       nmax=i
	    endif
 31	 continue
	 nrepr(nmax)= 1 ! = .true. : *is* a representative
	 do 41 j=1,nn
	    njn=meet(nmax,j)
	    if(dysma(j).gt.dys(njn)) dysma(j)=dys(njn)
 41	 continue
 20   continue

      sky=0.
      do 51 j=1,nn
	 sky=sky+dysma(j)
 51   continue
      obj(1)=sky/nn

      if(kk .gt. 1) then
c
c     second algorithm: swap.
c
C--   Loop :
 60	 do 63 j=1,nn
	    dysma(j)=1.1*s+1.0
	    dysmb(j)=1.1*s+1.0
	    do 62 i=1,nn
	       if(nrepr(i) .eq. 1) then
		  ij=meet(i,j)
		  if(dys(ij).lt.dysma(j)) then
		     dysmb(j)=dysma(j)
		     dysma(j)=dys(ij)
		  else
		     if(dys(ij).lt.dysmb(j)) dysmb(j)=dys(ij)
		  endif
	       endif
 62	    continue
 63	 continue
	 dzsky=1.0
	 do 73 k=1,nn
	    if(nrepr(k) .eq. 0) then
	       do 72 i=1,nn
		  if(nrepr(i) .eq. 1) then
		     dz=0.
		     do 71 j=1,nn
			ij=meet(i,j)
			kj=meet(k,j)
			if(dys(ij).eq.dysma(j))then
			   small=dysmb(j)
			   if(small.gt.dys(kj)) small=dys(kj)
			   dz=dz-dysma(j)+small
			else
 70			   if(dys(kj).lt.dysma(j))
     +				dz=dz-dysma(j)+dys(kj)
			endif
 71		     continue
		     if(dz.lt.dzsky) then
			dzsky=dz
			kbest=k
			nbest=i
		     endif
		  endif
 72	       continue
	    endif
 73	 continue
	 if(dzsky .lt. 0.)then
	    nrepr(kbest)=1
	    nrepr(nbest)=0
	    sky=sky+dzsky
	    go to 60
	 endif
      endif

      obj(2)=sky/nn
      end
c     -----------------------------------------------------------
c
c cstat():
c Compute STATistics (numerical output) concerning each partition
c
      subroutine cstat(kk,nn,nsend,nrepr,radus,damer,ttd,separ,z,s,
     f	   hh,dys,ncluv,nelem,med,nisol)

      integer kk,nn, hh, nsend(nn),nrepr(nn)
      integer ncluv(nn),nelem(nn), nisol(kk), med(kk)
      double precision radus(nn),damer(nn),ttd(nn),separ(nn),z,s,dys(hh)

c function called
      integer meet
      external meet
c
      logical kand
      integer j,k,m,ja,jb,jk, njaj, ksmal,nplac,numcl,numl,nel,njm,
     +	   nvn,ntt,mevj,nvna,jndz
      double precision dsmal, ttt,rtt,rnn,dam,sep,aja,ajb

      do 130 j=1,nn
	 if(nrepr(j).eq.0) then
	    dsmal=1.1*s+1.0
	    do 110 k=1,nn
	       if(nrepr(k) .eq. 1) then
		  njaj=meet(k,j)
		  if(dys(njaj).lt.dsmal) then
		     dsmal=dys(njaj)
		     ksmal=k
		  endif
	       endif
 110	    continue
	    nsend(j)=ksmal
	 else
	    nsend(j)=j
	 endif
 130  continue
      jk=1
      nplac=nsend(1)
      do 135 j=1,nn
	 ncluv(j)=0
	 if(nsend(j).eq.nplac)ncluv(j)=1
 135  continue
      do 145 ja=2,nn
	 nplac=nsend(ja)
	 if(ncluv(nplac).eq.0) then
	    jk=jk+1
	    do 140 j=2,nn
	       if(nsend(j).eq.nplac) ncluv(j)=jk
 140	    continue
	    if(jk.eq.kk)go to 148
	 endif
 145  continue
c
c     analysis of the clustering.
c
 148  do 160 numcl=1,kk
	 ntt=0
	 radus(numcl)=-1.0
	 ttt=0.0
	 do 150 j=1,nn
	    if(ncluv(j).eq.numcl) then
	       ntt=ntt+1
	       m=nsend(j)
	       nelem(ntt)=j
	       njm=meet(j,m)
	       ttt=ttt+dys(njm)
	       if(dys(njm).gt.radus(numcl)) radus(numcl)=dys(njm)
	    endif
 150	 continue
	 rtt=ntt
	 ttd(numcl)=ttt/rtt
	 med(numcl)=m
 160  continue
 230  rnn=nn

      if(kk.eq.1) then
	 damer(1)=s
	 nrepr(1)=nn
	 return
      endif
c  ELSE	  kk > 1 :
c
c     numl = number of l-clusters.
c
 240  numl=0

      do 40 k=1,kk
c
c     identification of cluster k:
c     nel  = number of objects
c     nelem= vector of objects
c
	 nel=0
	 do 23 j=1,nn
	    if(ncluv(j).eq.k) then
	       nel=nel+1
	       nelem(nel)=j
	    endif
 23	 continue
	 nrepr(k)=nel
	 if(nel.eq.1) then
	    nvn=nelem(1)
	    damer(k)=0.
	    separ(k)=1.1*s+1.0
	    do 250 j=1,nn
	       if(j.ne.nvn) then
		  mevj=meet(nvn,j)
		  if(separ(k).gt.dys(mevj)) separ(k)=dys(mevj)
	       endif
 250	    continue
c
c Is cluster k
c	1) an L-cluster	 or
c	2) an L*-cluster ?
	    if(separ(k).eq. 0.) numl=numl+1

	 else
c	       nel != 1 :
	    dam=-1.
	    sep=1.1*s+1.0
	    kand=.TRUE.
	    do 26 ja=1,nel
	       nvna=nelem(ja)
	       aja=-1.
	       ajb=1.1*s+1.0
	       do 25 jb=1,nn
		  jndz=meet(nvna,jb)
		  if(ncluv(jb).eq.k) then
		     if(dys(jndz).gt.aja) aja=dys(jndz)
		  else
		     if(dys(jndz).lt.ajb) ajb=dys(jndz)
		  endif
 25	       continue
	       if(kand .and. aja.ge.ajb) kand=.FALSE.
	       if(dam.lt.aja) dam=aja
	       if(sep.gt.ajb) sep=ajb
 26	    continue
	    separ(k)=sep
	    damer(k)=dam
	    if(kand) then
	       numl=numl+1
	       if(dam.ge.sep) then
c		  L-cluster
		  nisol(k)=1
	       else
c		  L*-cluster
 27		  nisol(k)=2
	       endif
	       go to 40
	    endif
	 endif
	 nisol(k)=0
 40   continue
 300  end
c     -----------------------------------------------------------
c
c     Compute Silhouette Information :
c
      subroutine dark(kk,nn,hh,ncluv,nsend,nelem,negbr,
     f	   syl,srank,avsyl,ttsyl,dys,s,sylinf)

      integer kk,nn,hh
      integer ncluv(nn),nsend(nn),nelem(nn),negbr(nn)
      double precision syl(nn),srank(nn),avsyl(nn),dys(hh),sylinf(nn,4)
      double precision ttsyl, s
c     function called
      integer meet
      external meet
c
      integer j,l, lang,lplac, nsylr, nclu,numcl,ntt,nj,njl,nl,nbb,mjl
      double precision dysa,dysb,db,btt,rtt,rnn,symax

      nsylr=0
      ttsyl=0.0
      do 100 numcl=1,kk
	 ntt=0
	 do 30 j=1,nn
	    if(ncluv(j).eq.numcl)then
	       ntt=ntt+1
	       nelem(ntt)=j
	    endif
 30	 continue
	 do 40 j=1,ntt
	    nj=nelem(j)
	    dysb=1.1*s+1.0
	    negbr(j)=-1
	    do 41 nclu=1,kk
	       if(nclu.ne.numcl) then
		  nbb=0
		  db=0.0
		  do 43 l=1,nn
		     if(ncluv(l).eq.nclu)then
			nbb=nbb+1
			if(l .ne. nj) then
			   mjl=meet(nj,l)
			   db=db+dys(mjl)
			endif
		     endif
 43		  continue
		  btt=nbb
		  db=db/btt
		  if(db.lt.dysb)then
		     dysb=db
		     negbr(j)=nclu
		  endif
	       endif
 41	    continue

	    if(ntt.gt.1) then
	       dysa=0.0
	       do 45 l=1,ntt
		  nl=nelem(l)
		  if(nj .ne. nl) then
		     njl=meet(nj,nl)
		     dysa=dysa+dys(njl)
		  endif
 45	       continue
	       dysa=dysa/(ntt - 1)
	       if(dysa.gt.0.0)then
		  if(dysb.gt.0.0) then
		     if(dysb.gt.dysa) then
			syl(j)=1.0-dysa/dysb
		     else if(dysb.lt.dysa) then
			syl(j)=dysb/dysa-1.0
		     else
c     dysb == dysa:
			syl(j)=0.0
		     endif
		     if(syl(j).le. -1.0) syl(j)=-1.0
		     if(syl(j).ge.  1.0) syl(j)= 1.0
		  else
		     syl(j)=-1.0
		  endif
	       else if(dysb.gt.0.0) then
		  syl(j)=1.0
	       else
		  syl(j)=0.0
	       endif
	    else
c     ntt == 1:
	       syl(j)=0.0
	    endif
 40	 continue
	 avsyl(numcl)=0.0
	 do 60 j=1,ntt
	    symax=-2.0
	    do 70 l=1,ntt
	       if(symax .lt. syl(l)) then
		  symax=syl(l)
		  lang=l
	       endif
 70	    continue
	    nsend(j)=lang
	    srank(j)=syl(lang)
	    avsyl(numcl)=avsyl(numcl)+srank(j)
	    syl(lang)=-3.0
 60	 continue
	 ttsyl=ttsyl+avsyl(numcl)
	 rtt=ntt
	 avsyl(numcl)=avsyl(numcl)/rtt

	 if(ntt.lt.2) then
	    nsylr=nsylr+1
	    sylinf(nsylr,1)=numcl
	    sylinf(nsylr,2)=negbr(1)
	    sylinf(nsylr,3)=0.0
	    sylinf(nsylr,4)=nelem(1)
	 else
	    do 80 l=1,ntt
	       nsylr=nsylr+1
	       lplac=nsend(l)
	       sylinf(nsylr,1)=numcl
	       sylinf(nsylr,2)=negbr(lplac)
	       sylinf(nsylr,3)=srank(l)
	       sylinf(nsylr,4)=nelem(lplac)
 80	    continue
	 endif
 100  continue
      rnn=nn
      ttsyl=ttsyl/rnn
      end
