      subroutine mona(nn,jpp, kx,jerr, nban,ner,kwan,lava, jlack)
cc
cc   MONothetic Analysis
cc
cc   Program for divisive hierarchical clustering of binary data,
cc   using association analysis.
cc
cc   list of functions and subroutines:
cc       function kab
cc
      implicit double precision(a-h,o-z)

      integer nn, jpp
cc          nn = number of objects 
cc          jpp = number of variables

cc The following vectors and matrices must be dimensioned in the main program:
      integer nban(nn),ner(nn),kwan(nn),lava(nn)
      integer jlack(jpp)
      integer kx(nn,jpp)
      integer jerr
cc jerr : error return code in {1,2,3,4}

      integer nzf

      lack=0
      jhalt=0 
      nnhal=(nn+1)/2
      jptwe=(jpp+4)/5 
      myst=0
      do 70 l=1,nn
         mysca=0 
         do 60 j=1,jpp 
            if(kx(l,j).eq.0)go to 60
            if(kx(l,j).eq.1)go to 60
            mysca=mysca+1 
 60      continue
         myst=myst+mysca 
         if(mysca .eq. jpp) then
c     all variables missing for this object
            jhalt=1 
            jerr=1
         endif
 70   continue
      if(jhalt.eq.1)return
      if(myst.eq.0)go to 290
      do 100 j=1,jpp
         jnul=0
         jeen=0
         do 80 l=1,nn
            if(kx(l,j).eq.0)jnul=jnul+1
            if(kx(l,j).eq.1)jeen=jeen+1
 80      continue
         jlack(j)=nn-jnul-jeen 
         if(jlack(j).ne.0)lack=lack+1
         if(jlack(j).ge.nnhal) then
c     at least 50% of the objects have missing values for this variable
            jhalt=1
            jerr=2
         endif
 90      if(jnul.eq.0)go to 95 
         if(jeen.eq.0)go to 95 
         go to 100 
cc   all non missing values are identical for this variable
 95      jhalt=1
         jerr=3
 100  continue
      if(jhalt.ne.0) return

      jpres=jpp-lack
      if(jpres.eq.0) then
c     all variables have missing values       
         jerr=4
         return
      endif
cc
cc   filling in missing values
cc
      do 260 j=1,jpp
         if(jlack(j).eq.0)go to 260
         lama=-1 
         nsyn=1
         do 240 ja=1,jpp 
            if(jlack(ja).ne.0)go to 240 
            jva=0 
            jvb=0 
            jvc=0 
            jvd=0 
            do 230 k=1,nn 
               if(kx(k,j).eq.1)go to 220
               if(kx(k,ja).eq.0)jva=jva+1
               if(kx(k,ja).eq.1)jvb=jvb+1
               go to 230
 220           if(kx(k,ja).eq.0)jvc=jvc+1
               if(kx(k,ja).eq.1)jvd=jvd+1
 230        continue
            kal=jva*jvd-jvb*jvc
            kalf=kab(kal)
            if(kalf.ge.lama)then
               lama=kalf 
               jma=ja
               if(kal.lt.0) nsyn=-1
            endif
 240     continue
         do 250 l=1,nn 
            if(kx(l,j).eq.0)go to 250
            if(kx(l,j).eq.1)go to 250
            if(nsyn.eq.1)then
               kx(l,j)=kx(l,jma)
            else
               if(kx(l,jma).eq.1)kx(l,j)=0
               if(kx(l,jma).eq.0)kx(l,j)=1
            endif
 250     continue
 260  continue
cc
cc   initialization 
cc
 290  do 300 k=1,nn 
         kwan(k)=0 
         ner(k)=k
         lava(k)=0 
 300  continue
      npass=1 
      kwan(1)=nn
cc
cc   algorithm
cc
      nclu=1
      ka=1
C --- Loop ---
 310  kb=ka+kwan(ka)-1
      lama=-1 
      jnat=jpp
      do 370 j=1,jpp
         if(nclu.eq.1)go to 330
         jnul=0
         jeen=0
         do 325 l=ka,kb
            nel=ner(l)
            if(kx(nel,j).eq.0)jnul=jnul+1
            if(kx(nel,j).eq.1)jeen=jeen+1
 325     continue
         if(jeen.eq.0)go to 370
         if(jnul.eq.0)go to 370
 330     jnat=jnat-1 
         lams=0
         do 360 jb=1,jpp 
            if(jb.eq.j)go to 360
            kva=0 
            kvb=0 
            kvc=0 
            kvd=0 
            do 350 l=ka,kb
               nel=ner(l)
               if(kx(nel,j).eq.1)go to 340
               if(kx(nel,jb).eq.0)kva=kva+1
               if(kx(nel,jb).eq.1)kvb=kvb+1
               go to 350
 340           if(kx(nel,jb).eq.0)kvc=kvc+1
               if(kx(nel,jb).eq.1)kvd=kvd+1
 350        continue
            lams=lams+kab(kva*kvd-kvb*kvc)
 360     continue
         if(lams.le.lama)go to 370 
         jtel=kvc+kvd
         jtelz=kva+kvb 
         lama=lams 
         jma=j 
 370  continue
      if(jnat.lt.jpp)go to 375
      kwan(ka)=-kwan(ka)
      go to 400 
cc
cc    splitting 
cc
 375  nel=ner(ka)
      if(kx(nel,jma).eq.1)then
         nzf=0
         jtel2=jtel
      else
         nzf=1
         jtel2=jtelz
      endif
      jres=kb-ka+1-jtel2
      km=ka+jtel2
      l=ka
c  -- inner loop --
 378  nel=ner(l)
      if(kx(nel,jma).eq.nzf)go to 380
      l=l+1 
      if(l.lt.km)go to 378
      go to 390 

 380  do 381 lbb=l,kb 
         nelbb=ner(lbb)
         if(kx(nelbb,jma).eq.nzf)go to 381
         lcc=lbb-1 
         go to 382 
 381  continue
 382  do 383 laa=l,lcc
         ldd=lcc+l-laa 
         lee=ldd+1 
         ner(lee)=ner(ldd) 
 383  continue
      ner(l)=nelbb
      go to 378 

 390  nclu=nclu+1 
      nban(km)=npass
      kwan(ka)=jtel2
      kwan(km)=jres
      lava(km)=jma
      ka=ka+kwan(ka)
 400  if(kb.eq.nn)go to 500 
 410  ka=ka+kab(kwan(ka)) 
      if(ka.gt.nn)go to 500 
      if(kwan(ka).lt.2)go to 410
      go to 310 

 500  npass=npass+1 
      do 510 ka=1,nn
         if(kwan(ka).ge.2)go to 310
 510  continue
      end 
cc  
cc kab(j) = |j| 
cc
      integer function kab(j) 

      implicit none
      integer j
      kab=j 
      if(j.lt.0) kab=-j
      return
      end
