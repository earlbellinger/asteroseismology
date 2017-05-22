      subroutine symm(a,nd,nr,ndt)
c  symm completes the completely symmetric array a from those elements
c  a(i,j,k,l,...) for which i.le.j.le.k.le.l.le.... .
c   nd: range of i,j,k,.....
c   nr: rank of a (i. e. number of indices)
c   ndt: dimension of a (must be the same for all indices)
      real   da
      logical ord,ordnr2,ordit
      dimension a(1)
      dimension i(30),ord(30),iam(30)
c
      save
c
    5 if(nd.eq.1.or.nr.eq.1) return
      nr1=nr-1
      nr2=nr-2
      nd1=nd-1
      ndt1=ndt-1
      do 10 k=1,nr1
      ord(k)=.true.
   10 i(k)=1
      ndr=ndt**nr
      ndr1=ndr/ndt
      do 15 k=2,nd
   15 iam(k)=ndr
      iam(1)=ndr1
      nf=ndr1
      nfc=ndt
      ngc0=ndr1/ndt
      ng=1+ngc0
      ngc=ngc0
      it=nr1
      iit=2
      ordit=.true.
      inr2=1
      ordnr2=.true.
   25 im=iit-1
   26 do 30 l=1,im
      if(l.eq.1) go to 28
      l1=l-1
      iaml=iam(l1)/ndt
      nf=nf+iaml
      ng=ng+ndr1
      iam(l1)=iaml
   28 a(ng)=a(1+nf)
   30 continue
      im1=im-1
      ng=ng-im1*ndr1
      if(im.eq.nd) go to 38
      nfc=nfc/ndt
c  reset nf because of change of im to 1
      im2=im-2
      if(im2) 38,34,32
   32 do 33 k=1,im2
   33 nf=nf+k*(iam(k+1)-iam(k))
   34 nf=nf-im1*iam(im1)
      do 35 k=1,im1
   35 iam(k)=iam(k)*ndt
   38 if(iit.lt.nd) go to 40
      if(it.le.1) return
      it=it-1
      iit=i(it)
      nfc=nfc*ndt
      ngc=ngc/ndt
      go to 38
c  reset nf because of change of nd's to 1
   40 go to (44),nfc
      do 42 k=1,nd1
   42 iam(k)=nfc*iam(k)
   44 iaml=iam(iit)/ndt
      nf=nfc*nf+iaml-ndr*nd1*(nfc-1)/ndt1
      iam(iit)=iaml
      nfc=ndt
c  reset ng
      ng=ng-nd1*(ndr1-ngc*ndt)/ndt1+ngc
c  reset i and ord
      iit=iit+1
      i(it)=iit
      if(nr2) 25,25,45
   45 if(it.eq.nr1) go to 50
      ngc=ngc0
      it1=it+1
      do 46 k=it1,nr1
      i(k)=1
   46 ord(k)=.false.
      im=nd
      go to(48),it
      if(ord(it-1)) ord(it)=iit.ge.i(it-1)
   48 it=nr1
      iit=1
      inr2=i(nr2)
      ordnr2=ord(nr2)
      ordit=.false.
   50 if(ordnr2) ordit=iit.ge.inr2
      if(ordit) go to 25
      go to 26
      end
