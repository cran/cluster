c Utility called from CLARA, PAM & TWINS
c original had MEET(), MEET2(), and MEET3() in the corresponding 3 source files

      integer function meet(l,j)
      integer l,j

      if(l.gt.j) then
c     			l > j
         meet= (l-2)*(l-1)/2 + j+1
      else if(l.eq.j) then
         meet= 1
      else
c   			l < j
         meet= (j-2)*(j-1)/2 + l+1
      endif
      return
      end
