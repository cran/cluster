      subroutine pam(nn,jpp,kk,x,dys,jdyss,valmd,jtmd,ndyst,nsend,
     f     nrepr,nelem,radus,damer,ttd,separ,ttsyl,med,obj,ncluv,
     f     clusinf,sylinf,nisol)
c     
c     partitioning around medoids
c     
c     carries out a clustering using the k-medoid approach.
c     
      integer nn,jpp,kk, jdyss,ndyst

      double precision x(nn,jpp), dys(1+nn*(nn-1)/2), valmd(jpp)
      integer jtmd(jpp), nsend(nn),nrepr(nn),nelem(nn), 
     +     ncluv(nn), nisol(kk), med(kk)
      double precision radus(nn),damer(nn),ttd(nn),separ(nn),ttsyl
      double precision clusinf(kk,5),sylinf(nn,4),obj(2)
c
      integer jhalt, nhalf, k,l
      double precision sky,s

c jdyss = 0 : compute distances from x 
c       = 1 : distances provided  in x      
      if(jdyss.ne.1) then
         jhalt=0
         call dysta(nn,jpp,x,dys,ndyst,jtmd,valmd,jhalt)
         if (jhalt.ne.0) then
            jdyss=-1
            return
         endif
      endif
c     
      s=0.0
      nhalf=nn*(nn-1)/2+1
      l=1
 130  l=l+1
      if(dys(l).gt.s)s=dys(l)
      if(l.lt.nhalf)go to 130

      call bswap(kk,nn,      nrepr,radus,damer,ttd, nhalf,dys,sky,s,obj)
      call cstat(kk,nn,nsend,nrepr,radus,damer,ttd, separ,    sky,s,
     f     nhalf,dys,ncluv,nelem,med,nisol)
      do 135 k=1,kk
         clusinf(k,1)=nrepr(k)
         clusinf(k,2)=radus(k)
         clusinf(k,3)=ttd(k)
         clusinf(k,4)=damer(k)
         clusinf(k,5)=separ(k)
 135  continue
      if(1 .lt. kk .and. kk .lt. nn) then
         call dark(kk,nn,nhalf,ncluv,nsend,nelem,nrepr,
     f        radus,damer,ttd,ttsyl,dys,s,sylinf)
      endif
      end
c     -----------------------------------------------------------
c     
      subroutine dysta(nn,jpp,x,dys,ndyst,jtmd,valmd,jhalt)

      integer nn, jpp, ndyst, jtmd(jpp), jhalt
      double precision x(nn,jpp), dys(1+nn*(nn-1)/2), valmd(jpp)

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
               if(jtmd(j).ge.0)goto 40 
               if(x(l,j).eq.valmd(j))goto 30 
               if(x(k,j).eq.valmd(j))goto 30 
 40            npres=npres+1 
               if(ndyst.ne.1)goto 50 
               clk=clk+(x(l,j)-x(k,j))*(x(l,j)-x(k,j)) 
               goto 30 
 50            clk=clk+dabs(x(l,j)-x(k,j))
 30         continue
            rpres=npres 
            if(npres.ne.0)goto 60 
            jhalt=1 
            dys(nlk)=-1.0
            goto 20 
 60         if(ndyst.ne.1)goto 70 
            dys(nlk)=dsqrt(clk*(pp/rpres)) 
            goto 20 
 70         dys(nlk)=clk*(pp/rpres) 
 20      continue
 100  continue
      end 
c     -----------------------------------------------------------
c     
      subroutine bswap(kk,nn,nrepr,dysma,dysmb,beter,hh,dys,sky,s,obj)

      integer kk,nn,nrepr(nn),hh
      double precision dysma(nn),dysmb(nn),beter(nn),dys(hh), sky,s,
     +     obj(2)
c function called
      integer meet
      external meet
c
      integer nny, j,ja,k,nkj, njn,njaj,nmax, nbest,kbest
      double precision cmd, ammax, rnn,dz,dzsky,small
c     
c     first algorithm: build.
c     
      nny=0
      do 17 j=1,nn
         nrepr(j)=0
         dysma(j)=1.1*s+1.0
 17   continue

C-- Loop :
 20   do 22 ja=1,nn
         if(nrepr(ja).ne.0)go to 22
         beter(ja)=0.
         do 21 j=1,nn
            njaj=meet(ja,j)
            cmd=dysma(j)-dys(njaj)
            if(cmd.gt.0.0)beter(ja)=beter(ja)+cmd
 21      continue
 22   continue
      ammax=0.
      do 31 ja=1,nn
         if(nrepr(ja).ne.0)go to 31
         if(beter(ja).lt.ammax)go to 31
         ammax=beter(ja)
         nmax=ja
 31   continue
      nrepr(nmax)=1
      nny=nny+1
      do 41 j=1,nn
         njn=meet(nmax,j)
         if(dys(njn).lt.dysma(j))dysma(j)=dys(njn)
 41   continue
      if(nny.ne.kk)go to 20

      sky=0.
      do 51 j=1,nn
         sky=sky+dysma(j)
 51   continue
      rnn=nn
      obj(1)=sky/rnn
      if(kk.eq.1)goto 75
c     
c     second algorithm: swap.
c     
C-- Loop :
 60   do 63 j=1,nn
         dysma(j)=1.1*s+1.0
         dysmb(j)=1.1*s+1.0
         do 62 ja=1,nn
            if(nrepr(ja).eq.0)go to 62
            njaj=meet(ja,j)
            if(dys(njaj).ge.dysma(j))go to 61
            dysmb(j)=dysma(j)
            dysma(j)=dys(njaj)
            go to 62
 61         if(dys(njaj).ge.dysmb(j))go to 62
            dysmb(j)=dys(njaj)
 62      continue
 63   continue
      dzsky=1.0
      do 73 k=1,nn
         if(nrepr(k).eq.1)go to 73
         do 72 ja=1,nn
            if(nrepr(ja).eq.0)go to 72
            dz=0.
            do 71 j=1,nn
               njaj=meet(ja,j)
               nkj=meet(k,j)
               if(dys(njaj).ne.dysma(j))go to 70
               small=dysmb(j)
               if(dys(nkj).lt.small)small=dys(nkj)
               dz=dz-dysma(j)+small
               go to 71
 70            if(dys(nkj).lt.dysma(j))dz=dz-dysma(j)+dys(nkj)
 71         continue
            if(dz.ge.dzsky)go to 72
            dzsky=dz
            kbest=k
            nbest=ja
 72      continue
 73   continue
      if(dzsky.ge.0.)go to 75
      nrepr(kbest)=1
      nrepr(nbest)=0
      sky=sky+dzsky
      go to 60

 75   rnn=nn
      obj(2)=sky/rnn
      end
c     -----------------------------------------------------------
c     
      subroutine cstat(kk,nn,nsend,nrepr,radus,damer,ttd,separ,z,s,
     f     hh,dys,ncluv,nelem,med,nisol)

      integer kk,nn, hh, nsend(nn),nrepr(nn)
      integer ncluv(nn),nelem(nn), nisol(kk), med(kk)
      double precision radus(nn),damer(nn),ttd(nn),separ(nn),z,s,dys(hh)

c function called
      integer meet
      external meet
c
      integer j,k,m,ja,jb,jk, njaj, ksmal,nplac,numcl,numl,nel,njm,
     +     nvn,ntt,mevj,kand,nvna,jndz
      double precision dsmal, ttt,rtt,rnn,dam,sep,aja,ajb

      do 130 j=1,nn
         if(nrepr(j).eq.1)go to 120
         dsmal=1.1*s+1.0
         do 110 k=1,nn
            if(nrepr(k).eq.0)go to 110
            njaj=meet(k,j)
            if(dys(njaj).ge.dsmal)go to 110
            dsmal=dys(njaj)
            ksmal=k
 110     continue
         nsend(j)=ksmal
         go to 130
 120     nsend(j)=j
 130  continue
      jk=1
      nplac=nsend(1)
      do 135 j=1,nn
         ncluv(j)=0
         if(nsend(j).eq.nplac)ncluv(j)=1
 135  continue
      do 145 ja=2,nn
         nplac=nsend(ja)
         if(ncluv(nplac).ne.0)go to 145
         jk=jk+1
         do 140 j=2,nn
            if(nsend(j).eq.nplac)ncluv(j)=jk
 140     continue
         if(jk.eq.kk)go to 148
 145  continue
c     
c     analysis of the clustering.
c     
 148  do 160 numcl=1,kk
         ntt=0
         radus(numcl)=-1.0
         ttt=0.0
         do 150 j=1,nn
            if(ncluv(j).ne.numcl)go to 150
            ntt=ntt+1
            m=nsend(j)
            nelem(ntt)=j
            njm=meet(j,m)
            ttt=ttt+dys(njm)
            if(dys(njm).gt.radus(numcl))radus(numcl)=dys(njm)
 150     continue
         rtt=ntt
         ttd(numcl)=ttt/rtt
         med(numcl)=m
 160  continue
 230  rnn=nn
      if(kk.ne.1)go to 240
      damer(1)=s
      nrepr(1)=nn
      go to 300
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
            if(ncluv(j).ne.k)go to 23
            nel=nel+1
            nelem(nel)=j
 23      continue
         nrepr(k)=nel
         if(nel.ne.1)go to 24
         nvn=nelem(1)
         damer(k)=0.
         separ(k)=1.1*s+1.0
         do 250 j=1,nn
            if(j.eq.nvn)go to 250
            mevj=meet(nvn,j)
            if(separ(k).gt.dys(mevj)) separ(k)=dys(mevj)
 250     continue
c     
c Is cluster k     
c	1) an l-cluster  or
c       2) an l*-cluster ?
c     
         if(separ(k).eq.0.)go to 400
         numl=numl+1
 400     go to 35
 24      dam=-1.
         sep=1.1*s+1.0
         kand=1
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
 25         continue
            if(aja.ge.ajb) kand=0
            if(dam.lt.aja) dam=aja
            if(sep.gt.ajb) sep=ajb
 26      continue
         separ(k)=sep
         damer(k)=dam
         if(kand.eq.0)go to 35
         numl=numl+1
         if(dam.lt.sep)go to 27
c     l-cluster
         nisol(k)=1
         go to 40
c     l*-cluster
 27      nisol(k)=2 
         go to 40
 35      nisol(k)=0
 40   continue
 300  end
c     -----------------------------------------------------------
c     
      subroutine dark(kk,nn,hh,ncluv,nsend,nelem,negbr,
     f     syl,srank,avsyl,ttsyl,dys,s,sylinf)

      integer kk,nn,hh
      integer ncluv(nn),nsend(nn),nelem(nn),negbr(nn)
      double precision syl(nn),srank(nn),avsyl(nn),dys(hh),sylinf(nn,4)
      double precision ttsyl, s
c function called
      integer meet
      external meet
c
      integer j,l, lang,lplac, nsylr, nclu,numcl,ntt,nj,njl,nl,nbb,mjl
      double precision dysa,dysb,db,att,btt,rtt,rnn,symax

      nsylr=0
      ttsyl=0.0 
      do 100 numcl=1,kk 
         ntt=0 
         do 30 j=1,nn
            if(ncluv(j).ne.numcl)go to 30
            ntt=ntt+1
            nelem(ntt)=j
 30      continue
         do 40 j=1,ntt 
            nj=nelem(j) 
            dysb=1.1*s+1.0
            negbr(j)=-1
            do 41 nclu=1,kk
               if(nclu.eq.numcl)go to 41
               nbb=0 
               db=0.0
               do 43 l=1,nn
                  if(ncluv(l).ne.nclu)go to 43
                  nbb=nbb+1 
                  mjl=meet(nj,l)
                  db=db+dys(mjl)
 43            continue
               btt=nbb 
               db=db/btt 
               if(db.ge.dysb)go to 41
               dysb=db 
               negbr(j)=nclu 
 41         continue
            if(ntt.eq.1)go to 50
            dysa=0.0
            do 45 l=1,ntt 
               nl=nelem(l) 
               njl=meet(nj,nl) 
               dysa=dysa+dys(njl)
 45         continue
            att=ntt-1
            dysa=dysa/att
            if(dysa.gt.0.0)go to 51 
            if(dysb.gt.0.0)go to 52
 50         syl(j)=0.0
            go to 40
 52         syl(j)=1.0
            go to 40
 51         if(dysb.le.0.0)go to 53 
            if(dysb.gt.dysa)syl(j)=1.0-dysa/dysb
            if(dysb.lt.dysa)syl(j)=dysb/dysa-1.0
            if(dysb.eq.dysa)syl(j)=0.0
            go to 54
 53         syl(j)=-1.0 
 54         if(syl(j).le.(-1.0))syl(j)=-1.0 
            if(syl(j).ge.1.0)syl(j)=1.0 
 40      continue
         avsyl(numcl)=0.0
         do 60 j=1,ntt 
            symax=-2.0
            do 70 l=1,ntt
               if(syl(l).le.symax)go to 70 
               symax=syl(l)
               lang=l
 70         continue
            nsend(j)=lang 
            srank(j)=syl(lang)
            avsyl(numcl)=avsyl(numcl)+srank(j)
            syl(lang)=-3.0
 60      continue
         ttsyl=ttsyl+avsyl(numcl)
         rtt=ntt 
         avsyl(numcl)=avsyl(numcl)/rtt 

         if(ntt.ge.2)goto 75
         nsylr=nsylr+1
         sylinf(nsylr,1)=numcl
         sylinf(nsylr,2)=negbr(1)
         sylinf(nsylr,3)=0.0
         sylinf(nsylr,4)=nelem(1)
         goto 100
 75      do 80 l=1,ntt 
            nsylr=nsylr+1
            lplac=nsend(l)
            sylinf(nsylr,1)=numcl
            sylinf(nsylr,2)=negbr(lplac)
            sylinf(nsylr,3)=srank(l)
            sylinf(nsylr,4)=nelem(lplac)
 80      continue
 100  continue
      rnn=nn
      ttsyl=ttsyl/rnn 
 96   end
