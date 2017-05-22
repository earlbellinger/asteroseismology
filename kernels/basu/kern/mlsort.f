      subroutine mlsort(a,nn,ia,icol,ncol,ind)
c
c  returns index ind(n), n = 1, nn, such that
c  a(icol(k), ind(n)), k = 1, ncol is sorted increasingly
c  with a varying more rapidly with increasing k
c
      logical resort
      parameter(nnmax = 10000)
      dimension a(ia,nn), ind(nn), icol(ncol), w(nnmax), iw(nnmax),
     *   iw1(nnmax)
c
c  initial sort
c
      call sort(a(icol(1),1),nn,ia,ind)
c
      if(ncol.eq.1) return
c
c  prepare for subsequent sorts
c
      ic=icol(1)
      do 20 nc=2,ncol
      icp=ic
      ic=icol(nc)
      n1=1
      ns=ind(1)
      ap=a(icp,ns)
c
      do 20 n=2,nn+1
      if(n.le.nn) then
        ns=ind(n)
c..        write(10,*) n, ns, ap, a(icp,ns)
        resort=a(icp,ns).ne.ap
      else
        resort=.true.
      end if
c
      if(resort) then
        if(n-n1.ge.2) then
          do 10 k=1,n-n1
          k1=n1+k-1
          ks=ind(k1)
          iw1(k)=ks
   10     w(k)=a(ic,ks)
c..          write(10,1000) (k, w(k), iw1(k), k=1,n-n1)
 1000     format(//' k, w(k), iw1(k)'//(i4,1pe13.5,i4))
          call sort(w,n-n1,1,iw)
c
          do 12 k=1,n-n1
          ks=iw(k)
   12     ind(n1+k-1)=iw1(ks)
c..          write(10,1010) (k, w(iw(k)), ind(n1+k-1), a(ic, ind(n1+k-1)),
c..     *      k=1,n-n1)
 1010     format(//' k, sorted w, resorted ind, resorted a(ic,.)'//
     *      (i4,1pe13.5,i4,1pe13.5))
        end if
        if(n.lt.nn) then
          n1=n
          ap=a(icp,ns)
        end if
      end if
c
   20 continue
c
      return
      end
      subroutine sort(a,nn,ia,ind)
c
c  returns index ind(n), n = 1, nn, such that a(1,ind(n)) is sorted
c  in increasing order.
c  note: a is not changed
c
      parameter(nnmax = 10000)
      dimension a(ia,nn),ind(nn),card(2,nnmax)
c
      do 10 n=1,nn
      card(1,n)=a(1,n)
   10 card(2,n)=n
c
      call hsort(card,nn)
c
      do 20 n=1,nn
      n1=nn+1-n
   20 ind(n1)=int(card(2,n))
c
      return
      end
c  ------------------------------------------------------------
      subroutine sift(card,nn,l,n,ln)
      dimension card(2,nn)
      i=l
      j=2*i
      x=card(1,i+ln)
100   if(j.gt.n) return
      if(j.lt.n.and.card(1,j+ln).gt.card(1,j+1+ln)) j=j+1
      if(x.gt.card(1,j+ln)) then
      do 200 k=1,2
      t1=card(k,i+ln)
      card(k,i+ln)=card(k,j+ln)
      card(k,j+ln)=t1
  200 continue
      i=j
      j=2*i
      else
      return
      endif
      go to 100
      end
c ----------------------------------------
      subroutine hsort(card,nn)
      dimension card(2,nn)
c     constructing the heap

      l=nn/2+1
      do 100 i=l,1,-1
      call sift(card,nn,i,nn,0)
100   continue

      do 200 m=nn,2,-1
      do 150 k=1,2
      x=card(k,1)
      card(k,1)=card(k,m)
      card(k,m)=x
  150 continue
      n1=m-1
      call sift(card,nn,1,n1,0)
200   continue
      end
