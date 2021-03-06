c                                                           ************
c                                                           *          *
c                                                           *   main   *
c                                                           *          *
c                                                           ************
c
c
c     first calculate exact Eulerian kernels ln(c^2) and ln(rho)
c     (assuming hydrostatic equilibrium), then if requested proceed
c     to calculate kernels for different pairs of variables
c
c     parameter nnmax must be large enough to accommodate the model
c
      implicit double precision(a-h,o-z)
      logical diag,zero1
      parameter (nnmax=1500,nwork=6*nnmax)
c
c     model variables
      dimension data(8), xold(nnmax), aaold(5,nnmax)
      dimension data(8), x(nnmax), aa(5,nnmax)
      dimension aax(3,nnmax),dq(nnmax)
c
c     mode variables
      dimension cs(50), xi0(2,nnmax)
      dimension xi0x(2,nnmax)
c
c     kernels
      dimension rk(2,nnmax),rka(2,nnmax)
c
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz
      common/work/dummy(nwork)
c     
      diag = .false.
c
c     define constants (gconst is set to unity for the calculation,
c     which is nondimensional with G=M=R=1)
c
c     the dimensional value of G is kept in gdimen
c
      gdimen = 6.6732d-8
      gconst = 1.d0
      fourpi = 8.d0*acos(0.d0)
      twopi  = fourpi/2.0
      tiny   = 1.d-10
c
c     counter for kernels
      nk = 0
c
c     define degrees for which kernels required
c
      write(*,*) 'Input lmin, lmax and lstep'
      write(*,*) '(min and max degrees for which kernels required, ',
     +           'and increment)'
      read(*,*) lmin,lmax,lstep
      write(*,*) 'Choose which kernels to calculate:'
      write(*,*) '   1.  Eulerian ln(c^2) and ln(rho)'
      write(*,*) '   2.  Eulerian ln(Gamma1) and ln(rho)'
      write(*,*) '   3.  Eulerian ln(c^2) and ln(Gamma1)'
      write(*,*) '   4.  Eulerian ln(u) and ln(Gamma1)'
      write(*,*) '  11.  Lagrangian ln(c^2) and ln(rho)'
      write(*,*) '  12.  Lagrangian ln(Gamma1) and ln(rho)'
      read(*,*) icase
      if ( (icase.lt.1.or.icase.gt.4).and.
     +     (icase.lt.11.or.icase.gt.12) ) then
         write(*,*) 'Inappropriate choice of case number'
         stop
      endif
c
c     unit numbers for files:
c     input  -- model (ichaa),  eigensolutions (ichxi)
c     output -- kernels (ichkr1, ichkr2)
      ichaa = 10
      ichxi = 14
      ichkr1= 11
      ichkr2= 12
c
c     open files
c
      call ofiles
      call openf(ichaa,'o','u')
      call openf(ichxi,'o','u')
      call openf(ichkr1,'u','u')
      call openf(ichkr2,'u','u')
c
c     read in model
c
      call rdaa5(ichaa,data,xold,aaold,nnold,iend)
      write (6,1005) data
      if (nnold.gt.nnmax) then
         write(*,*) 'nnold = ',nnold,' exceeds nnmax'
         stop
      endif
c
c     read in mode mesh
c
      read(ichxi,end=900) nn,(x(l1),l1=1,nn)
      write(6,*) ' nn = ',nn
      if (nn.gt.nnmax) then
         write(*,*) 'nn = ',nn,' exceeds nnmax'
         stop
      endif
         do 100 l1=1,nn
         call lir(x(l1),xold,aa(1,l1),aaold,5,5,nnold,l1,inter)
  100    continue
c
c     flag whether x(1) = 0
c
      if (abs(x(1)).lt.tiny) then
         zero1 = .true.
         inonz = 2
         x(1) = 0.
      else
         zero1 = .false.
         inonz = 1
      endif
         
c
c
c     set up dimensional factors facp = GM^2/R^4 (units of pressure)
c     and fac = (GM/R^3)^(1/2) (units of speed)
c
      facp = data(1)/data(2)/data(2)/data(2)/data(2)
     +           *data(1)*gdimen
      fac  = sqrt(facp/data(1)*data(2))
c     
c +++ diagnostic
      if (diag) write(6,701)(x(l1),(aa(l2,l1),l2=1,5),l1=1,nn)
  701 format(' x,aa  ',e15.7,2x,5e15.7)
c +++ end diagnostic
c
c     set up array of ancillary model variables
c           aax(1,.) = rho
c           aax(2,.) = c^2
c           aax(3,.) = d(rho)/dr
	 do 110 l1=inonz,nn
	 aax(1,l1) = aa(1,l1)*aa(5,l1)/fourpi
         aax(2,l1) = aa(1,l1)/aa(2,l1)*x(l1)**2
         aax(3,l1) = -(aa(2,l1)+aa(4,l1))*aax(1,l1)/x(l1)
  110    continue
      if (zero1) then
         x22 = x(2)*x(2)
         x32 = x(3)*x(3)
         aax(1,1) = aa(1,1)*aa(5,1)/fourpi
	 aax(2,1) = ( aax(2,2)*x32 - aax(2,3)*x22 ) / (x32-x22)
	 aax(3,1) = 0.
      endif
c
c    set up d(ln c^2)/dr or d(ln gamma)/dr if necessary
c
      if (icase.eq.11.or.icase.eq.12) then
c        first compute d(ln gamma)/dr
         call derivk(x,aa(3,1),dq,1,nn,5,1,2)
            do 255 l1=1,nn
            dq(l1) = dq(l1)/aa(3,l1)
  255       continue
c        then compute d(ln c^2)/dr
         if (icase.eq.11) then
            do 253 l1=1,nn
            dq(l1) = dq(l1)+(1.-aa(3,l1))*aa(2,l1)/x(l1)
     +                    + aa(4,l1)/x(l1)
  253       continue
         endif
      endif

c
c     write mesh to kernel files
      write(ichkr1) nn,(x(l1),l1=1,nn)
      write(ichkr2) nn,(x(l1),l1=1,nn)
c
c     loop to read in modes and (if suitable degree) compute kernels
c
      lseek = lmin
  200 continue
      read(ichxi,end=900)
     +   (cs(l1),l1=1,50),(xi0(1,l1),xi0(2,l1),l1=1,nn)
      el     = cs(18)
      iel    = ifix(el+0.5)
c     
  201 if (iel.gt.lseek) then
         lseek = lseek+lstep
	 if (lseek.gt.lmax) goto 900
         goto 201
      endif
c     
      if (iel.lt.lseek) goto 200
c
  203 continue
      iord   = ifix(abs(cs(19))+0.1)
      if (cs(19).lt.0.) iord = -iord
      omega0 = sqrt(cs(20))
      rnu    = cs(27)*1.e3
      om2    = cs(20)
      omv2 = cs(27)*1.e-3*twopi
      omv2 = omv2**2/gdimen/data(1)*data(2)*data(2)*data(2)
      write(6,1006) iord,el,omv2
c
c     read in Joergen's functions z1 and z2 (mode energy is proportional to
c     integral of z1^2 + z2^2 d(ln r)
c
c     convert this to xi and eta and store in xi0(1,.) and xi0(2,.)
c
         do 256 l1=inonz,nn
         facx = sqrt(aax(1,l1)*x(l1))*x(l1)
         xi0(1,l1)=xi0(1,l1)/facx
         xi0(2,l1)=xi0(2,l1)/facx/sqrt(el*(el+1.0))
  256    continue
      if (zero1) then
         xi0(1,1) = 0.
c  l=1 case might be improved upon
         if (iel.eq.1) then
            xi0(1,1) = xi0(1,2)
         endif
      xi0(2,1) = xi0(1,1)
      endif
c +++ diagnostic
      if (diag) write(6,702)(x(l1),xi0(1,l1),xi0(2,l1),l1=1,nn)
  702 format(' r,xi,eta '/(3e15.7))
c +++ end diagnostic
c
c     calculate integral for denominator
      call calci(x,aa,xi0,nn,el,denom)
c
c     form d(xi)/dr
c
      ell1 = el*(el+1.)
      call derivk(x,xi0,xi0x,1,nn,2,2,2)
c
c     form divergence of xi
c
	 do 120 l1=inonz,nn
         xi0x(2,l1) = xi0x(1,l1) + 2.*xi0(1,l1)/x(l1)
     +                   - ell1*xi0(2,l1)/x(l1)
  120    continue
c  l=0,1 cases should be thought about more for r=0
      if (zero1) then
         xi0x(2,1) = 0.
         if (iel.eq.0.or.iel.eq.1) xi0x(2,1) = xi0x(2,2)
      endif
c
c
      call kercro(x,nn,aa,aax,xi0,xi0x,rk,2,omv2,
     +                            el,ell1,denom)


      if (icase.eq.1) then
         write(ichkr1) iel,iord,rnu,(rk(1,l1),l1=1,nn)
         write(ichkr2) iel,iord,rnu,(rk(2,l1),l1=1,nn)
      endif

      if (icase.eq.11) then
         call rholag(x,nn,aax(2,1),3,dq,1,aax,rk)
         write(ichkr1) iel,iord,rnu,(rk(1,l1),l1=1,nn)
         write(ichkr2) iel,iord,rnu,(rk(2,l1),l1=1,nn)
      endif

      if (icase.eq.2) then
         call kergro(x,nn,aa,aax,rk,rka,2)
         write(ichkr1) iel,iord,rnu,(rka(1,l1),l1=1,nn)
         write(ichkr2) iel,iord,rnu,(rka(2,l1),l1=1,nn)
      endif

      if (icase.eq.12) then
         call kergro(x,nn,aa,aax,rk,rka,2)
         call rholag(x,nn,aa(3,1),5,dq,1,aax,rka)
         write(ichkr1) iel,iord,rnu,(rka(1,l1),l1=1,nn)
         write(ichkr2) iel,iord,rnu,(rka(2,l1),l1=1,nn)
      endif

      if(icase.eq.3.or.icase.eq.4) then
         call kercg(x,nn,aa,aax,rk,rka,2)
         write(ichkr1) iel,iord,rnu,(rka(1,l1),l1=1,nn)
         if (icase.eq.3)
     +     write(ichkr2) iel,iord,rnu,(rka(2,l1),l1=1,nn)
         if (icase.eq.4)
     +     write(ichkr2) iel,iord,rnu,(rka(1,l1)+rka(2,l1),
     1     l1=1,nn)
      endif

c     return to start of loop, to read in another mode
  250 nk = nk+1
      goto 200

c
c     finished computing kernels
c
  900 write(6,*) 'Number of modes for which kernels stored = ',nk
      stop

c
 1000 format(5e15.7)
 1001 format(i5/(5e15.7))
 1002 format(i5,f10.3,e15.7)
 1005 format(' array data: '/8e14.6/)
 1006 format(/'    ****************************'///
     +        ' new mode   n=',i5,
     1        '    l=',f10.1,'   omega0=',e15.4)
      end
c
c
c
c                                                           ************
c                                                           *          *
c                                                           *  kercro  *
c                                                           *          *
c                                                           ************
      subroutine kercro(x,nn,aa,aax,xi0,xi0x,rk,ik,
     +                         om2,el,ell1,denom)
c
c     calculates kernels for ln(c^2) and ln(rho)
c
      implicit double precision(a-h,o-z)
      logical zero1
      parameter(nnmax=1500)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension xi0(2,nn), xi0x(2,nn)
      dimension rk(ik,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz
      common/work/dummy(nnmax,6)
      iopt = 0
         do 100 l1=inonz,nn
         lback=nn-l1+1
         dummy(lback,1) = xi0(1,l1) * 
     +        ( 2.*aax(1,l1)*xi0x(2,l1)+aax(3,l1)*xi0(1,l1) )
         dummy(lback,2) = aax(1,l1)/x(l1)**el*(xi0(1,l1)
     1                            -el*xi0(2,l1))
         dummy(lback,3) = x(l1)
  100    continue
c   linear extrapolation of dummy(.,2) might be improved upon
      if (zero1) then
         dummy(nn,1) = xi0(1,1) * 
     +        ( 2.*aax(1,1)*xi0x(2,1)+aax(3,1)*xi0(1,1) )
         dummy(nn,2) = (dummy(nn-2,2)*dummy(nn-1,3) 
     +      - dummy(nn-1,2)*dummy(nn-2,3))/(dummy(nn-1,3)-dummy(nn-2,3))
         dummy(nn,3) = x(1)
      endif
      call zero(dummy(1,5),nn)
      call vintk(dummy(1,3),dummy(1,1),dummy(1,5),1,nn,1,1,2)
      call zero(dummy(1,4),nn)
      call vintk(dummy(1,3),dummy(1,2),dummy(1,4),1,nn,1,1,2)
c
         do 110 l1=1,nn
         dummy(l1,1) = -dummy(nn-l1+1,5)
	 dummy(l1,2) = -dummy(nn-l1+1,4)
  110    continue
c
	 do 115 l1=inonz,nn
	 dummy(l1,5) = x(l1)**(el+1.)*aax(1,l1)*
     +             (xi0(1,l1)+(el+1.)*xi0(2,l1))
  115    continue
      if (zero1) dummy(1,5) = 0.
      call zero(dummy(1,4),nn)
c
      call vintk(x,dummy(1,5),dummy(1,4),1,nn,1,1,2)
c
	 do 116 l1=inonz,nn
	 dummy(l1,2) = (el+1.)*dummy(l1,2)
     +        -aax(1,l1)*xi0(1,l1)/x(l1)**(el-1.)
     1        +aax(1,nn)*xi0(1,nn)/x(nn)**(el-1.)
	 dummy(l1,4) = -el*dummy(l1,4)
     +        +aax(1,l1)*xi0(1,l1)*x(l1)**(el+2.)
  116    continue
c   linear extrapolation below might be improved upon
      if (zero1) then
         call interp(x(2),dummy(2,2),1,x(1),dummy(1,2))
         call interp(x(4),dummy(2,4),1,x(1),dummy(1,4))
      endif
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
c           final factor of 2 inserted 6/11/92 to make kernel for c^2
c           rather than for c
         rk(irow,l1) = chi*chi*rho*c2*rx2/om2/denom/2.
         if (iopt.eq.1) goto 140
         irow = 2
c
c        density kernel
  130    r1 = -(xi**2+ell1*eta**2)*om2
         r2 = chi*chi*c2
         if (zero1.and.l1.eq.1) then
            r3 = 0.
         else
            r3 = -2.*gconst*rm*xi*
     +            (2.*xi-ell1*eta)/rx
         endif
         r4 = gconst*xi*xi*fourpi*rho
         r5 = -fourpi*gconst*dummy(l1,1)
         r1 = r1/om2/denom*rho*rx2
         r2 = r2/om2/denom*rho*rx2
         r12 = r1+r2
         r3 = r3/om2/denom*rho
         r4 = r4/om2/denom*rho*rx2
         r5 = r5/om2/denom*rho*rx2
         if (zero1.and.l1.eq.1) then
            r6 = 0.
         else
	    r6 = el*rx**(el+1.)*(xi+(el+1.)*eta)*dummy(l1,2)
     +           -(el+1.)/rx**el*(xi-el*eta)*dummy(l1,4)
         endif
	 r6 = r6*fourpi*gconst*rho/denom/om2/(2.*el+1.)
c        write(6,*) r1,r2,r12,r3,r4,r5
         rk(irow,l1) = (r12+r3+r4+r5)/2.+r6
  140    continue
      return
      end

      subroutine interp(xi,yi,iorder,x,y)
c
c     given xi(j) and yi(j), j=1,...,iorder+1 and a value x
c        fits a polynomial of order iorder through xi,yi and evaluate it
c        at x  -  returns value in y
c
      implicit double precision (a-h,o-z)
      dimension xi(1), yi(1)
c
c     currently only implemented for iorder=0,1
c
      if (iorder.lt.0) then
         write(*,*) 'In interp - iorder must be non-negative!!!'
         stop
      endif
      if (iorder.ge.2) then
         write(*,*) 'In interp - iorder>2 but this not implemented'
         stop
      endif
      if (iorder.eq.0) then
         y = yi(1)
         return
      endif
      if (iorder.eq.1) then
         y = ( yi(1)*(x-xi(2)) + yi(2)*(xi(1)-x) ) / ( xi(1)-xi(2) )
         return
      endif 
      end

c                                                           ************
c                                                           *          *
c                                                           *  kergro  *
c                                                           *          *
c                                                           ************
      subroutine kergro(x,nn,aa,aax,rkcro,rkgro,ik)
c
c     calculates kernels for ln(gamma1) and ln(rho)
c       given equilibrium state and kernels for ln(c^2) and ln(rho)
c
      implicit double precision(a-h,o-z)
      logical zero1
      parameter(nnmax=1500)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rkcro(ik,nn),rkgro(ik,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz
      common/work/dummy(nnmax,5),press(nnmax)
c
c     compute pressure
c
         do 100 l1=1,nn
	 press(l1) = aax(1,l1)*aax(2,l1)/aa(3,l1)
  100    continue
c
c     compute Ktilde == int_0^r K_{c^2,rho}/p 
c       and store it in dummy(.,2)
c
         do 110 l1=1,nn
	 dummy(l1,1) = rkcro(1,l1)/press(l1)
  110    continue
      call zero(dummy(1,2),nn)
      call vintk(x,dummy(1,1),dummy(1,2),1,nn,1,1,2)
c
c     compute 4 pi G rho r^2 * int_r^R Ktilde rho / s^2 ds
c       and store it in dummy(.,1)
c
         do 120 l1=inonz,nn
	 lback = nn-l1+1
	 dummy(lback,1) = x(l1)
	 dummy(lback,3) = dummy(l1,2)*aax(1,l1)/x(l1)**2
  120    continue
      if (zero1) then
         npt = nn-1
      else
         npt = nn
      endif
      call zero(dummy(1,4),nn)
      call vintk(dummy(1,1),dummy(1,3),dummy(1,4),1,npt,1,1,2)
         do 130 l1=inonz,nn
	 lback = nn-l1+1
	 dummy(l1,1) = 
     +       -dummy(lback,4)*fourpi*gconst*aax(1,l1)*x(l1)**2
  130    continue
      if (zero1) dummy(1,1) = 0.
c
c     finally assemble the kernels 
c
         do 150 l1=1,nn
	 rkgro(1,l1) = rkcro(1,l1)
	 rkgro(2,l1) = rkcro(2,l1) - rkcro(1,l1) + dummy(l1,1)
     +        + gconst*dummy(l1,2)*aa(1,l1)*aax(1,l1)*x(l1)
  150    continue
c
      return
      end
c                                                           ************
c                                                           *          *
c                                                           *  kercg   *
c                                                           *          *
c                                                           ************
      subroutine kercg(x,nn,aa,aax,rkcro,rkcg,ik)
c
c     calculates kernels for ln(c^2) and ln(gamma1) 
c       given equilibrium state and kernels for ln(c^2) and ln(rho)
c
c     solves d(psi)/dr  -4 pi G rho r^2 int_r^R (rho /(r^2 p) psi) dr  
c                = K_{rho,c^2} + some multiple of (rho r^2)
c
c     (for more details see s/r kercg1)
c
c     and then K_{gamma,c^2} = p d(p^-1 psi)/dr
c              K_{c^2,gamma} = K_{c^2,rho}-K_{gamma,c^2}
c
      implicit double precision(a-h,o-z)
      logical zero1
      parameter (nnmax=1500)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rkcro(ik,nn),rkcg(ik,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz
      common/work/dummy(nnmax,6)

      
c     kercg1 returns psi and (rho r^2)^-1 d(psi)/dr
      call kercg1(x,nn,aa,aax,rkcro,rkcg,ik)
      
         do 100 l1=inonz,nn
         dummy(l1,1) = rkcg(1,l1)
         dummy(l1,2) = rkcg(2,l1)*aax(1,l1)*x(l1)**2
  100    continue
         
         do 110 l1=inonz,nn
         rkcg(2,l1) = dummy(l1,2) + aa(2,l1)*aa(3,l1)/x(l1)*dummy(l1,1)
  110    continue
      if (zero1) rkcg(1,l1) = 0.

         do 120 l1=1,nn
         rkcg(1,l1) = rkcro(1,l1)-rkcg(2,l1)
  120    continue

      return
      end
c                                                           ************
c                                                           *          *
c                                                           *  kercg1  *
c                                                           *          *
c                                                           ************
      subroutine kercg1(x,nn,aa,aax,rkcro,rkcg,ik)
c
c     solves d(psi)/dr  -4 pi G rho r^2 int_r^R (rho /(r^2 p) psi) dr  
c                = K_{rho,c^2} + some multiple of (rho r^2)
c
c     The arbitrary multiple of (rho r^2) arises principally because of
c     numerical problems near the surface, where the density varies on
c     a subgrid scale. This adds a constant to the integral (from
c     near-surface errors). But K_{rho,c^2} is arbitrary up to an 
c     additive multiple of (rho r^2), so it doesn't seem worth worrying
c     further. Because of this freedom, the second boundary condition
c     (to the solde into which the original integral equation is cast) is 
c     chosen rather arbitrarily to be psi(R)=1. This of course adds some
c     complementary function into the solution, but this also just contributes
c     some multiple of (rho r^2) into the integral equation.
c
c     the dependent variables used are 
c          z1 = psi,     z2 = (rho r^2)^-1 d(psi)/dr
c
c     these are returned in rkcg(1,.) and rkcg(2,.) respectively
c
      implicit double precision(a-h,o-z)
      double precision intlhs,intrhs,intres
      integer*4 v
      logical zero1
      parameter (nnmax=1500)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rkcro(ik,nn),rkcg(ik,nn)
      dimension ea(2,3),v(2)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23
      common/work/dummy(nnmax,6)
      external rhs,bc

      if (zero1) then
         nn1 = nn - 1
      else
         nn1 = nn
      endif

      call zero(b1,nn)
      call zero(b2,nn)
         do 100 l1=inonz,nn
	 press = aax(1,l1)*aax(2,l1)/aa(3,l1)
         a11(l1) = 0.
         a12(l1) = aax(1,l1)*x(l1)**2
	 a21(l1) = -fourpi*gconst*aax(1,l1)/(press*x(l1)**2)
	 a22(l1) = 0.
  100    continue
c     put d(Kroc/(rho r^2))/dr into b2
         do 105 l1=inonz,nn
         dummy(l1,6) = rkcro(2,l1)/aax(1,l1)/x(l1)**2
  105    continue
      call derivk(x(inonz),dummy(inonz,6),b2(inonz),1,nn1,1,1,2)
      write(*,*) 'After call to derive'

      g23 = dummy(nn,6)

c
c     solve differential equation for psi
c
      ii = 2
      kk = 0
      ka = 1
      kb = 1
      ki = 0
      id = 2
      v(1) = 1
      v(2) = 2
      ucy = 1.

      iter = 0
  150 continue 

      
      write(*,*) 'Before call to nrk'
      call nrk(x(inonz),rkcg(1,inonz),zk,ap,aq,rhs,bc,ii,kk,ka,kb,ki,
     +         nn1,id,ucy,ea,det,v)
      write(*,*) 'After call to nrk'
      write(*,*) 'ea:'
      write(*,*) ea
      iter = iter+1
      if (iter.lt.2) goto 150

      if (zero1) then
         rkcg(1,1) = 0.
         rkcg(2,1) = 0.
      endif

c
c     checking answer
c
         do 200 l1=1,nn
         dummy(l1,1) = rkcg(2,l1)*aax(1,l1)*x(l1)**2
  200    continue
         do 210 l1=inonz,nn
         lback = nn-l1+1
         dummy(lback,3) = x(l1)
         dummy(lback,4) = rkcg(1,l1)*aa(3,l1)/(x(l1)**2*aax(2,l1))
  210    continue
      call zero(dummy(1,5),nn)
      call vintk(dummy(1,3),dummy(1,4),dummy(1,5),1,nn1,1,1,2)
         do 220 l1=inonz,nn
         lback = nn-l1+1
         dummy(l1,3) = -dummy(lback,5)
  220    continue
      if (zero1) then
         dummy(1,3) = dummy(2,3)
      endif
         do 230 l1=1,nn
         dummy(l1,2) = dummy(l1,1)
     +      - fourpi*gconst*aax(1,l1)*x(l1)*x(l1)*dummy(l1,3)
         dummy(l1,3) = dummy(l1,2)-rkcro(2,l1)
         dummy(l1,4) = aax(1,l1)*x(l1)**2
         dummy(l1,3) = dummy(l1,3)*dummy(l1,4)
         dummy(l1,4) = dummy(l1,4)**2
  230    continue
      call zero(dummy(1,5),nn)
      call zero(dummy(1,6),nn)
      call vintk(x,dummy(1,3),dummy(1,5),1,nn,1,1,2)
      call vintk(x,dummy(1,4),dummy(1,6),1,nn,1,1,2)
      alpha = dummy(nn,5)/dummy(nn,6)
      write(*,*) 'alpha = ',alpha
         do 250 l1=1,nn
         dummy(l1,3) = dummy(l1,2)**2
  250    continue
      call zero(dummy(1,5),nn)
      call vintk(x,dummy(1,3),dummy(1,5),1,nn,1,1,2) 
      intlhs = dummy(nn,5)
         do 260 l1=1,nn
         dummy(l1,3) = rkcro(2,l1)**2
  260    continue
      call zero(dummy(1,5),nn)
      call vintk(x,dummy(1,3),dummy(1,5),1,nn,1,1,2) 
      intrhs = dummy(nn,5)
         do 270 l1=1,nn
         dummy(l1,3) = dummy(l1,2)-rkcro(2,l1)
     +                    - alpha*aax(1,l1)*x(l1)**2
         dummy(l1,3) = dummy(l1,3)**2
  270    continue
      call zero(dummy(1,5),nn)
      call vintk(x,dummy(1,3),dummy(1,5),1,nn,1,1,2) 
      intres = dummy(nn,5)
      write(*,*) 'Integral of squared residue = ',intres
      write(*,*) 'compared with integrals of squared'
      write(*,*) 'LHS, RHS respectively = ',intlhs,intrhs

 1001 format(5e15.7)

      return
      end

     
      subroutine rhs(x,y,zk,ap,aq,f,fd,h,hd,ifd,n)
      implicit real*8 (a-h,o-z)
      logical zero1
      parameter (nnmax=1500)
      dimension y(1),f(1),fd(ifd,1),h(1),hd(ifd,1),zk(1),ap(1),aq(1)
      common/centre/zero1,inonz
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23
      if (zero1) then 
         n1 = n+1
      else
         n1 = n
      endif
      fd(1,1) = a11(n1)
      fd(1,2) = a12(n1)
      f(1) = fd(1,1)*y(1) + fd(1,2)*y(2)
      fd(2,1) = a21(n1)
      fd(2,2) = a22(n1)
      f(2) = fd(2,1)*y(1) + fd(2,2)*y(2) + b2(n1)
      return
      end
 
      subroutine bc(x1,x2,y1,y2,zk,ap,aq,g,gd,ig,id,n)
      implicit real*8(a-h,o-z)
      logical zero1
      parameter (nnmax=1500)
      dimension y1(1),y2(1),g(1),gd(ig,1),zk(1),ap(1),aq(1)
      common/centre/zero1,inonz
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23

      gd(1,1) = -3.
      gd(1,2) = x1*a12(inonz)
      g(1) = gd(1,1)*y1(1) + gd(1,2)*y1(2)
      gd(2,1) = 1.0
      gd(2,2) = 0.
      g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) - 0.0
c     gd(2,1) = 0.
c     gd(2,2) = 1.
c     g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) - g23
      return
      end

      subroutine rholag(x,nn,q,idq,dq,iddq,aax,rk)
c
c     convert Eulerian kernels for (q,rho) into Lagrangian kernels for (q,rho)
c
c     dq contains d(ln q)/dr
c
      implicit double precision (a-h,o-z)
      logical diag,zero1
      parameter (nnmax=1500)
      dimension x(nn),q(idq,nn),dq(iddq,nn),aax(3,nn),rk(2,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz
      common/work/dummy(nnmax,6)
  
      diag = .false.

      if (zero1) dummy(1,1) = 0.
         do 100 l1=inonz,nn
         dummy(l1,1) = rk(1,l1)*dq(1,l1)
     +                   + rk(2,l1)*aax(3,l1)/aax(1,l1)
         dummy(l1,1) = dummy(l1,1) / x(l1)**2
  100    continue
      
         do 105 l1=1,nn
         ipt = nn-l1+1
         dummy(l1,2) = x(ipt)
         dummy(l1,3) = dummy(ipt,1)
  105    continue
      
      dummy(1,4) = 0.
      call vintk(dummy(1,2),dummy(1,3),dummy(1,4),1,nn,1,1,2)
         do 107 l1=1,nn
         ipt = nn-l1+1
         dummy(l1,2) = dummy(ipt,4)
  107    continue
c
c     if x(1)=0 just integrate from x(2)  --  introduces a multiple of
c     r^2 into Lagrangian kernel for rho, which is OK
c
c     call zero(dummy(1,2),nn)
c     nn1 = nn - (inonz-1)
c     call vintk(x(inonz),dummy(inonz,1),dummy(inonz,2),1,nn1,1,1,2)
c     if (zero1) dummy(1,2) = 0.

         do 110 l1=1,nn
         rk(2,l1) = rk(2,l1) - dummy(l1,2)*x(l1)**2
  110    continue

      return
      end
          
