      logical diag
      dimension cs(50),  t(2,1000), rk(2,1000)
      dimension data(8), xi0(2,1000)
      dimension aa(5,1000),x(1000)
      dimension aax(3,1000), xi0x(2,1000)
      common/consts/ fourpi,gconst
      common/work/dummy(10000)
      nk = 0
      diag = .false.
      ilseek = 0
      ilinc  = 1
      ichxi = 14
      ichxir= 13
      ichxis= 9
      ichxit= 8
      ichaa = 10
      ichkr1= 11
      ichkr2= 12
      call ofiles
      call openf(ichkr1,'n','u')
      call openf(ichkr2,'n','u')
      call openf(ichaa,'o','u')
      call openf(ichxi,'o','u')
      gconst = 6.6732e-8
      fourpi = 8.*acos(0.0)
      twopi  = fourpi/2.0
c
c     read model
c
      call rdaa5(ichaa,data,x,aa,nn,iend)
c
c     naughty model-dependent shuffle to eliminate x=0 and x>1
c
      nn  = nn-1
	 do 100 l1=1,nn
c        if (x(l1+1).gt.1.) goto 100
	 x(l1) = x(l1+1)
	    do 105 l2=1,5
	    aa(l2,l1) = aa(l2,l1+1)
  105       continue
  100    continue
      facp = data(1)/data(2)/data(2)/data(2)/data(2)
     +           *data(1)*gconst
      fac  = sqrt(facp/data(1)*data(2))
      gconst = 1.0
      write (6,1005) data
      if (diag) write(6,1023)(x(l1),(aa(l2,l1),l2=1,5),l1=1,nn)
 1023 format(' x,aa  ',e15.7,2x,5e15.7)
c
c     set up array of rho, c**2, d(rho)/dr
c
	 do 110 l1=1,nn
	 aax(1,l1) = aa(1,l1)*aa(5,l1)/fourpi
	 aax(2,l1) = aa(1,l1)/aa(2,l1)*x(l1)*x(l1)
         aax(3,l1) = -(aa(2,l1)+aa(4,l1))*aax(1,l1)/x(l1)
  110    continue
      write(ichkr1) nn,(x(l1),l1=1,nn)
c
c     read in next mode
      read(ichxi,end=900) nnw,(dummy(l1),l1=1,nnw)
      write(6,*) ' nnw = ',nnw
  200 read(ichxi,end=900) (cs(l1),l1=1,50),(xi0(1,l1),
     +  xi0(2,l1),l1=1,nnw)
      el     = cs(18)
      iel    = ifix(el+0.5)
c     if (iel.gt.150) goto 900
c 201 if (iel.lt.ilseek) goto 200
c     if (iel.eq.ilseek) goto 203
c     if (iel.gt.50) ilinc = 2
c     if (iel.gt.80) ilinc = 5
c     if (iel.gt.135) ilinc = 10
c     if (iel.gt.145) ilinc = 5
c     ilseek = ilseek+ilinc
c     goto 201
  203 continue
      iord   = ifix(abs(cs(19))+0.1)
      if (cs(19).lt.0.) iord = -iord
      omega0 = sqrt(cs(20))
      rnu    = cs(27)*1.e3
      om2    = cs(20)
      omv2 = cs(27)*1.e-3*twopi
      omv2 = omv2*omv2/6.6732e-8/data(1)*data(2)*data(2)*data(2)
      write(6,1006) iord,el,omv2
c     if (el.lt.2.0.or.iord.ne.2) goto 200
c
c     have to shuffle eigenfunction also
c
      do 202 l1=1,nn
      xi0(1,l1) = xi0(1,l1+1)
      xi0(2,l1) = xi0(2,l1+1)
  202 continue
         do 256 l1=1,nn
         facx = sqrt(aax(1,l1)*x(l1))*x(l1)
         xi0(1,l1)=xi0(1,l1)/facx
         xi0(2,l1)=xi0(2,l1)/facx/sqrt(el*(el+1.0))
  256    continue
      if (diag) write(6,8475)(x(l1),xi0(1,l1),xi0(2,l1),l1=1,nn)
 8475 format(' r,xi,eta '/(3e15.7))
c
c     calculate integral for denominator
      call calci(x,aa,xi0,nn,el,denom)
c
c     form d(xi)/dr
c
      ell1 = el*(el+1.)
      call derivk(x,xi0,xi0x,1,nn,2,2,2)
c     use Cowling approximation
c     do 115 l1=1,nn
c     xi0x(1,l1) =  (aa(2,l1)-2.)*xi0(1,l1)/x(l1)
c    +              + (ell1/x(l1)-om2/aax(2,l1)*x(l1))*xi0(2,l1)
c     write(6,1865) l1,x(l1),xi0x(1,l1),xi0x(2,l1)
  115 continue
 1865 format(i5,1p3e15.7)
c     ijk = 1
c     if (ijk.eq.1) stop
c
c     form divergence of xi
c
	 do 120 l1=1,nn
         xi0x(2,l1) = xi0x(1,l1) + 2.*xi0(1,l1)/x(l1)
     +                   - ell1*xi0(2,l1)/x(l1)
c     use Cowling approximation
c     xi0x(2,l1) = aa(2,l1)*xi0(1,l1)/x(l1)
c    +              - om2/aax(2,l1)*x(l1)*xi0(2,l1)
  120    continue
c
c
      call kercro(x,nn,aa,aax,xi0,xi0x,rk,2,omv2,
     +                            el,ell1,denom)
      write(ichkr1) iel,iord,rnu,(rk(1,l1),l1=1,nn)
      write(ichkr2)(rk(2,l1),l1=1,nn)
      nk = nk+1
      goto 200
  900 write(6,*) ' Running nk = ',nk
c     if (ichxi.eq.ichxit) stop
c     if (ichxi.eq.ichxis) ichxi = ichxit
c     if (ichxi.eq.ichxir) ichxi = ichxis
c     if (ichxi.ne.ichxit.and.ichxi.ne.ichxis) ichxi = ichxir
c     read(ichxi) nnew,(dummy(l1),l1=1,nnew)
c     goto 200
      stop
 1000 format(5e15.7)
 1001 format(i5/(5e15.7))
 1002 format(i5,f10.3,e15.7)
 1005 format(' array data: '/8e14.6/)
 1006 format(/'    ****************************'///
     +        ' new mode   n=',i5,
     1        '    l=',f10.1,'   omega0=',e15.4)
      end
      subroutine here(i)
      write(6,1) i
      return
    1 format(' here ',i5)
      end
      subroutine kercro(x,nn,aa,aax,xi0,xi0x,rk,ik,
     +                         om2,el,ell1,denom)
      dimension x(nn),aax(3,nn),xi0(2,nn),xi0x(2,nn)
      dimension aa(5,1000),rk(ik,nn)
      common/work/ dummy(10000)
      common/consts/ fourpi,gconst
      iopt = 0
         do 100 l1=1,nn
         dummy(nn-l1+1) = xi0(1,l1) * 
     +        ( 2.*aax(1,l1)*xi0x(2,l1)+aax(3,l1)*xi0(1,l1) )
         dummy(1001+nn-l1) = aax(1,l1)/x(l1)**el*(xi0(1,l1)
     1                            -el*xi0(2,l1))
         dummy(2000+nn-l1+1) = x(l1)
  100    continue
      call zero(dummy(4001),nn)
      call vintk(dummy(2001),dummy,dummy(4001),1,nn,1,1,2)
      call zero(dummy(3001),nn)
      call vintk(dummy(2001),dummy(1001),dummy(3001),1,nn,1,1,2)
c
c     the two minus signs in the next loop were inserted 2/9/89
c
         do 110 l1=1,nn
         dummy(l1) = -dummy(4000+nn-l1+1)
	 dummy(1000+l1) = -dummy(3000+nn-l1+1)
  110    continue
	 do 115 l1=1,nn
	 dummy(4000+l1) = x(l1)**(el+1.)*aax(1,l1)*
     +             (xi0(1,l1)+(el+1.)*xi0(2,l1))
  115    continue
      call zero(dummy(3001),nn)
      call vintk(x,dummy(4001),dummy(3001),1,nn,1,1,2)
	 do 116 l1=1,nn
	 dummy(1000+l1) = (el+1.)*dummy(1000+l1)
     +        -aax(1,l1)*xi0(1,l1)/x(l1)**(el-1.)
     1        +aax(1,nn)*xi0(1,nn)/x(nn)**(el-1.)
	 dummy(3000+l1) = -el*dummy(3000+l1)
     +        +aax(1,l1)*xi0(1,l1)*x(l1)**(el+2.)
  116    continue
c
c
c
         do 140 l1=1,nn
         irow = 1
         rx  = x(l1)
         rx2 = rx*rx
         rho = aax(1,l1)
         c2  = aax(2,l1)
         c   = sqrt(c2)
         chi = xi0x(2,l1)
         xi  = xi0(1,l1)
         eta = xi0(2,l1)
         rm  = aa(1,l1)*rx*rx2
         if (iopt.eq.2) goto 130
c       
c        sound speed kernel
         rk(irow,l1) = chi*chi*rho*c2*rx2/om2/denom
         if (iopt.eq.1) goto 140
         irow = 2
c
c        density kernel
  130    r1 = -(xi*xi+ell1*eta*eta)*om2
         r2 = chi*chi*c2
         r3 = -2.*gconst*rm*xi*
     +         (2.*xi-ell1*eta)/rx
         r4 = gconst*xi*xi*fourpi*rho
         r5 = -fourpi*gconst*dummy(l1)
         r1 = r1/om2/denom*rho*rx2
         r2 = r2/om2/denom*rho*rx2
         r12 = r1+r2
         r3 = r3/om2/denom*rho
         r4 = r4/om2/denom*rho*rx2
         r5 = r5/om2/denom*rho*rx2
	 r6 = el*rx**(el+1.)*(xi+(el+1.)*eta)*dummy(1000+l1)
     +        -(el+1.)/rx**el*(xi-el*eta)*dummy(3000+l1)
	 r6 = r6*fourpi*gconst*rho/denom/om2/(2.*el+1.)
c        write(6,*) r1,r2,r12,r3,r4,r5
         rk(irow,l1) = (r12+r3+r4+r5)/2.+r6
  140    continue
c     do 150 l1=1,nn
c     write(6,*) x(l1),rk(1,l1),rk(2,l1)
c 150 continue
      return
      end
C     FILED AS  MJT12.MGRTPDS:CALCI
C
      SUBROUTINE CALCI(X,AA,XI0,NN,EL,DENOM)
C
C     CALCULATE INTEGRAL I (DENOM)
      LOGICAL DIAG
      DIMENSION X(NN),AA(5,NN),XI0(2,NN)
      COMMON/WORK/ T1(1000), T2(1000)
      DIAG = .FALSE.
      ELL1 = EL*(EL+1.)
      FOURPI = 12.56637061
         DO 100 L1=1,NN
         RHO = AA(1,L1)*AA(5,L1)/FOURPI
         T1(L1) = RHO*(XI0(1,L1)*XI0(1,L1)+
     +             ELL1*XI0(2,L1)*XI0(2,L1))*X(L1)*X(L1)
  100    CONTINUE
      call zero(t2,nn)
      call vintk(x,t1,t2,1,nn,1,1,2)
      DENOM = T2(NN)
      DIAG = .TRUE.
      IF(DIAG) WRITE(6,1001) DENOM
      DIAG = .FALSE.
 1001 FORMAT(' INTEGRAL I = ',E15.7)
      RETURN
      END
