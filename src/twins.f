c
c     This program performs agglomerative nesting (AGNES) using the
c     group average method of Sokal and Michener (1958), as well as
c     divisive analysis (DIANA) using the method of
c     Mcnaughton-Smith, Williams, Dale, and Mockett (1964).
c
      subroutine twins(nn,jpp,x,dys,dys2,jdyss,valmd,jtmd,ndyst,jalg,
     f     method,kwan,ner,ban,coef,merge)

c Arguments
      integer nn, jpp
c     	nn  = maximal number of objects
c       jpp = maximal number of variables used in the analysis
c
      integer jdyss, jtmd(jpp), ndyst,jalg,method, kwan(nn), ner(nn)
      double precision x(nn,jpp), valmd(jpp)
      double precision dys((nn*(nn-1))/2 + 1), dys2((nn*(nn-1))/2 + 1)
      double precision ban(nn), coef
      integer merge(nn-1,2)
C VARs
      integer i,jhalt

      if(jdyss.ne.0) then
         jpp=1
      else
         jhalt=0
         call dysta(nn,jpp,x,dys,ndyst,jtmd,valmd,jhalt)
         if(jhalt.ne.0) then
            jdyss=-1
            return
         endif
      endif

      do 10 i=1,(nn*(nn-1))/2 + 1
         dys2(i)=dys(i)
 10   continue

      if(jalg.ne.2) then
         call averl(nn,kwan,ner,ban,dys,method,merge)
         call banag(nn,ban,ner,coef)
      else
         call splyt(nn,kwan,ner,ban,dys,merge)
         call bandy(nn,ban,ner,coef)
      endif
      end
c     -----------------------------------------------------------
c
      subroutine averl(nn,kwan,ner,ban,dys,method,merge)

      integer nn, kwan(nn), ner(nn), method, merge(nn-1,2)
      double precision ban(nn), dys(1+nn*(nn-1)/2)
c Function (defined in ./meet.f ):
      integer meet
c VARs
      integer j,l,la,lb,l1,l2,lq, lka, lenda, lendb
      integer lfyrs, llast, lmuch, lnext, lput, lnum
      integer naq, nbq, nab, nej,nlj, nns, nclu,nmerge
      double precision akb, d, dnew, fa,fb,fc, smald, ta,tb,tq
c
c
      nclu=nn-1
c     initialization
      do 10 l=1,nn
         kwan(l)=1
         ner(l)=l
 10   continue
c
c     find closest clusters
c
      nmerge=1
 100  j=1
 80   j=j+1
      if(kwan(j).eq.0)goto 80
      nej=meet(1,j)
      smald=dys(nej)*1.1+1.0
      nns=nn-1
      do 120 l=1,nns
         if(kwan(l).eq.0)go to 120
         lmuch=l+1
         do 110 j=lmuch,nn
            if(kwan(j).eq.0)go to 110
            nlj=meet(l,j)
            if(dys(nlj).gt.smald)go to 110
            smald=dys(nlj)
            la=l
            lb=j
 110     continue
 120  continue
c
c     merge-structure for plotting tree in S
c
      l1=-la
      l2=-lb
      if(nmerge.eq.1)go to 121
      do 122 j=1,(nmerge-1)
         if((merge(j,1).eq.l1).or.(merge(j,2).eq.l1))l1=j
         if((merge(j,1).eq.l2).or.(merge(j,2).eq.l2))l2=j
 122  continue
 121  merge(nmerge,1)=l1
      merge(nmerge,2)=l2
      nmerge=nmerge+1
c
c     determine lfyrs and llast
c
      do 200 l=1,nn
         if(ner(l).eq.la)lfyrs=l
         if(ner(l).eq.lb)llast=l
 200  continue
      ban(llast)=smald
c
c     if the two clusters are next to each other,
c     ner must not be changed
c
      lnext=lfyrs+kwan(la)
      if(lnext.eq.llast)goto 230
c
c     updating ner and ban
c
      lput=lfyrs+kwan(la)
      lnum=llast-lput
      do 220 l=1,lnum
         lka=ner(lput)
         akb=ban(lput)
         lenda=llast+kwan(lb)-2
         lendb=lenda+1
         do 210 j=lput,lenda
            ner(j)=ner(j+1)
            ban(j)=ban(j+1)
 210     continue
         ner(lendb)=lka
         ban(lendb)=akb
 220  continue
c
c     calculate new dissimilarities
c
 230  do 240 lq=1,nn
         if(lq.eq.la.or.lq.eq.lb)go to 240
         if(kwan(lq).eq.0)go to 240
         naq=meet(la,lq)
         nbq=meet(lb,lq)
         if(method.eq.2)go to 300
         if(method.eq.3)go to 310
         if(method.eq.4)go to 320
         if(method.eq.5)go to 330
c     group average method
         ta=kwan(la)
         tb=kwan(lb)
         fa=ta/(ta+tb)
         fb=tb/(ta+tb)
         dys(naq)=fa*dys(naq)+fb*dys(nbq)
         go to 240
c     single linkage
 300     dnew=dys(naq)
         if(dys(nbq).lt.dnew)dnew=dys(nbq)
         dys(naq)=dnew
         go to 240
c     complete linkage
 310     dnew=dys(naq)
         if(dnew.lt.dys(nbq))dnew=dys(nbq)
         dys(naq)=dnew
         go to 240
c     ward's method
 320     ta=kwan(la)
         tb=kwan(lb)
         tq=kwan(lq)
         fa=(ta+tq)/(ta+tb+tq)
         fb=(tb+tq)/(ta+tb+tq)
         fc=-tq/(ta+tb+tq)
         nab=meet(la,lb)
         d=fa*dys(naq)*dys(naq)+fb*dys(nbq)*dys(nbq)
         d=d+fc*dys(nab)*dys(nab)
         dys(naq)=sqrt(d)
         go to 240
c     weighted average linkage
 330     dys(naq)=(dys(naq)+dys(nbq))/2.d0
 240  continue
 250  kwan(la)=kwan(la)+kwan(lb)
      kwan(lb)=0
      nclu=nclu-1
      if(nclu.gt.0)goto 100
      end
c     -----------------------------------------------------------
c
      subroutine banag(nn,ban,ner,ac)

      integer nn, ner(nn)
      double precision ban(nn), ac
c VARs
      integer k,kearl,kafte
      double precision sup,syze,rnn

      sup=0.0
      do 70 k=2,nn
         if(ban(k).gt.sup)sup=ban(k)
 70   continue
      ac=0.0
      do 80 k=1,nn
         kearl=k
         if(k.eq.1)kearl=2
         kafte=k+1
         if(k.eq.nn)kafte=nn
         syze=ban(kearl)
         if(ban(kafte).lt.syze)syze=ban(kafte)
         ac=ac+1.0-(syze/sup)
 80   continue
      rnn=nn
      ac=ac/rnn
      end
c     -----------------------------------------------------------
c
      subroutine splyt(nn,kwan,ner,ban,dys,merge)

      integer nn, kwan(nn), ner(nn)
      double precision ban(nn), dys(1+nn*(nn-1)/2)
      integer merge(nn-1,2)
c Function (defined in ./meet.f ):
      integer meet
c VARs
      double precision arest, dyff, bdyff, bygsd, cs, da,db, dmin
      double precision rest, splyn, sd
      integer j,k,l, ja,jb, jan,jbn,jab, jma, jmb, jaway, jner
      integer lndsd, lchan, lner, lmm, lmma, lmmb
      integer l1, l2, lgrb, llq, lmz, lxf, lxg, lxx,lxy, lxxa, lxxp
      integer nj, nlj, nclu, nhalf, nmerge

c
c     initialization
c
      nclu=1
      nhalf=nn*(nn-1)/2+1
      do 10 l=1,nn
         kwan(l)=0
         ban(l)=0.
         ner(l)=l
 10   continue
      kwan(1)=nn
      ja=1
c
c     computation of diameter of data set
c
      cs=0.0
      k=0
 20   k=k+1
      if(dys(k).gt.cs)cs=dys(k)
      if(k.lt.nhalf)go to 20
c
c     prepare for splitting
c
 30   jb=ja+kwan(ja)-1
      jma=jb
c
c     special case of a pair of objects
c
      if(kwan(ja).ne.2)go to 50
      kwan(ja)=1
      kwan(jb)=1
      jan=ner(ja)
      jbn=ner(jb)
      jab=meet(jan,jbn)
      ban(jb)=dys(jab)
      go to 400
c
c     finding first object to be shifted
c
 50   bygsd=-1.
      do 110 l=ja,jb
         lner=ner(l)
         sd=0.
         do 100 j=ja,jb
            jner=ner(j)
            nlj=meet(lner,jner)
            sd=sd+dys(nlj)
 100     continue
         if(sd.le.bygsd)go to 110
         bygsd=sd
         lndsd=l
 110  continue
c
c     shifting the first object
c
      kwan(ja)=kwan(ja)-1
      kwan(jb)=1
      if(jb.eq.lndsd)go to 115
      lchan=ner(lndsd)
      lmm=jb-1
      do 112 lmma=lndsd,lmm
         lmmb=lmma+1
         ner(lmma)=ner(lmmb)
 112  continue
      ner(jb)=lchan
 115  splyn=0.
      jma=jb-1
c
c     finding the next object to be shifted
c
 120  splyn=splyn+1.
      rest=jma-ja
      bdyff=-1.
      do 150 l=ja,jma
         lner=ner(l)
         da=0.
         do 130 j=ja,jma
            jner=ner(j)
            nlj=meet(lner,jner)
            da=da+dys(nlj)
 130     continue
         da=da/rest
         db=0.
         jmb=jma+1
         do 140 j=jmb,jb
            jner=ner(j)
            nlj=meet(lner,jner)
            db=db+dys(nlj)
 140     continue
         db=db/splyn
         dyff=da-db
         if(dyff.le.bdyff)go to 150
         bdyff=dyff
         jaway=l
 150  continue
      jmb=jma+1
c
c     shifting the next object when necessary
c
      if(bdyff.le.0.)go to 200
      if(jma.eq.jaway)go to 165
      lchan=ner(jaway)
      lmz=jma-1
      do 160 lxx=jaway,lmz
         lxxp=lxx+1
         ner(lxx)=ner(lxxp)
 160  continue
      ner(jma)=lchan
 165  do 170 lxx=jmb,jb
         lxy=lxx-1
         if(ner(lxy).lt.ner(lxx))go to 180
         lchan=ner(lxy)
         ner(lxy)=ner(lxx)
         ner(lxx)=lchan
 170  continue
 180  kwan(ja)=kwan(ja)-1
      kwan(jma)=kwan(jmb)+1
      kwan(jmb)=0
      jma=jma-1
      jmb=jma+1
      if(jma.ne.ja)go to 120
c
c     switch the two parts when necessary
c
 200  if(ner(ja).lt.ner(jmb))go to 300
      lxxa=ja
      do 220 lgrb=jmb,jb
         lxxa=lxxa+1
         lchan=ner(lgrb)
         do 210 lxy=lxxa,lgrb
            lxf=lgrb-lxy+lxxa
            lxg=lxf-1
            ner(lxf)=ner(lxg)
 210     continue
         ner(lxg)=lchan
 220  continue
      llq=kwan(jmb)
      kwan(jmb)=0
      jma=ja+jb-jma-1
      jmb=jma+1
      kwan(jmb)=kwan(ja)
      kwan(ja)=llq
c
c     compute level for banner
c
 300  if(nclu.eq.1)ban(jmb)=cs
      if(nclu.eq.1)go to 400
      call supcl(dys,ja,jb,arest,nn,ner)
      ban(jmb)=arest
 400  nclu=nclu+1
      if(nclu.eq.nn)goto 500
c
c     continue splitting until all objects are separated
c
      if(jb.eq.nn)go to 430
 420  ja=ja+kwan(ja)
      if(ja.gt.nn)go to 430
      if(kwan(ja).le.1)go to 420
      go to 30
 430  ja=1
      if(kwan(ja).eq.1)go to 420
      go to 30
c
c     merge-structure for plotting tree in S
c
 500  do 550 nmerge=1,(nn-1)
         dmin=cs
         do 560 j=2,nn
            if ((kwan(j).lt.0).or.(ban(j).gt.dmin))goto 560
            dmin=ban(j)
            nj=j
 560     continue
         kwan(nj)=-1
         l1=-ner(nj-1)
         l2=-ner(nj)
         if(nmerge.eq.1)go to 570
         do 580 j=1,(nmerge-1)
            if((merge(j,1).eq.l1).or.(merge(j,2).eq.l1))l1=j
            if((merge(j,1).eq.l2).or.(merge(j,2).eq.l2))l2=j
 580     continue
 570     merge(nmerge,1)=l1
         merge(nmerge,2)=l2
 550  continue
      end
c     -----------------------------------------------------------
c
      subroutine supcl(dys,kka,kkb,arest,nn,ner)

      integer kka,kkb, nn, ner(nn)
      double precision dys(1+nn*(nn-1)/2), arest
c Function (defined in ./meet.f ):
      integer meet
c VARs
      integer j,l, kkc,kkd, jner,lner,mlj

      kkc=kkb-1
      arest=0.
      do 20 l=kka,kkc
         lner=ner(l)
         kkd=l+1
         do 10 j=kkd,kkb
            jner=ner(j)
            mlj=meet(lner,jner)
            if(dys(mlj).gt.arest)arest=dys(mlj)
 10      continue
 20   continue
      return
      end
c     -----------------------------------------------------------
c
      subroutine bandy(nn,ban,ner,dc)
c Arguments
      integer nn,ner(nn)
      double precision ban(nn), dc
c VARs
      integer k, kearl, kafte
      double precision sup, syze, rnn

      sup=0.0
      do 70 k=2,nn
         if(ban(k).gt.sup)sup=ban(k)
 70   continue
      dc=0.0
      do 80 k=1,nn
         kearl=k
         if(k.eq.1)kearl=2
         kafte=k+1
         if(k.eq.nn)kafte=nn
         syze=ban(kearl)
         if(ban(kafte).lt.syze)syze=ban(kafte)
         dc=dc+1.0-(syze/sup)
 80   continue
      rnn=nn
      dc=dc/rnn
      end

