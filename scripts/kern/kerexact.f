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
      parameter (nnmax=10000,nwork=6*nnmax)
c
c     model variables
      dimension data(8), xold(nnmax), aaold(5,nnmax)
      dimension x(nnmax), aa(5,nnmax), xefn(nnmax), psi(nnmax)
      dimension aax(3,nnmax),dq(nnmax),dgamma(3,nnmax),
     +          xdum(nnmax),dgammad(3,nnmax)
c
c     mode variables
      dimension cs(50), xi0(2,nnmax)
      dimension xi0x(2,nnmax)
      real xr(nnmax),xi0r(2,nnmax),csr(50)
c
c     kernels
      dimension rk(2,nnmax),rka(2,nnmax),rkb(2,nnmax)
c
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
      common/comple/rkcomp(2,nnmax),rkcomi,psic(2,nnmax),icompl
      common/work/dummy(nwork)
c     
      diag = .false.
c.....I=IEEE_HANDLER('see','common',SIGFPE_ABORT)
c
c     define constants (gconst is set to unity for the calculation,
c     which is nondimensional with G=M=R=1)
c
c     the dimensional value of G is kept in gdimen
c
c      gdimen = 6.6732d-8
      gdimen = 6.67408d-8
c      gdimen = 6.67428d-8
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
      write(*,*) 'Input precision of eigenfunctions'
      write(*,*) '1 = single,   2 = double precision'
c     read(*,*) iprec1
      read(*,1000) iprec1
      write(*,*) 'Input lmin, lmax and lstep'
      write(*,*) '(min and max degrees for which kernels required, ',
     +           'and increment)'
      read(*,*) lmin,lmax,lstep
      write(*,*) 'Choose which kernels to calculate:'
      write(*,*) '  -1.  Energy density, normalized to be unimodular'
      write(*,*) '   1.  Eulerian ln(c^2) and ln(rho)'
      write(*,*) '   2.  Eulerian ln(Gamma1) and ln(rho)'
      write(*,*) '   3.  Eulerian ln(c^2) and ln(Gamma1)'
      write(*,*) '   4.  Eulerian ln(u) and ln(Gamma1)'
      write(*,*) '   5.  Eulerian ln(u) and Y'
      write(*,*) '   6.  Eulerian ln(rho) and Y'
      write(*,*) '   7.  Eulerian ln(c) and ln(Gamma1/c)'
      write(*,*) '  11.  Lagrangian ln(c^2) and ln(rho)'
      write(*,*) '  12.  Lagrangian ln(Gamma1) and ln(rho)'
      write(*,*) '  13.  Lagrangian ln(c^2) and ln(Gamma1)'
      write(*,*) '  14.  Lagrangian ln(u) and ln(Gamma1)'
      write(*,*) '  17.  Lagrangian ln(c) and ln(Gamma1/c)'
      read(*,*) icase
      if ( (icase.ne.-1).and.(icase.lt.1.or.icase.gt.7).and.
     +     (icase.lt.11.or.icase.gt.14) .and. icase.ne.17) then
         write(*,*) 'Inappropriate choice of case number'
         stop
      endif
      write(*,*) 'Input fractional radius at which model should be ',
     +   'truncated'
      write(*,*) '   (negative for no truncation)'
      read(*,*) rmax
c
c     unit numbers for files:
c     input  -- model (ichaa),  eigensolutions (ichxi)
c     output -- kernels (ichkr1, ichkr2)
c     inserted by Earl Bellinger on May 11 2017 
c     adds output for psi and p d(p^-1 psi)/dr
      ichaa = 10
      ichxi = 14
      ichkr1= 11
      ichkr2= 12
      ichpsi= 13
c
c     open files
c
      call ofiles
      call openf(ichaa,'o','u')
      call openf(ichxi,'o','u')
      call openf(ichkr1,'u','u')
      call openf(ichkr2,'u','u')
      call openf(ichpsi,'u','u')
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
      if (iprec1.eq.2) then
         read(ichxi,end=900) nn,(x(l1),l1=1,nn)
         write(*,*) "x(1), x(nn)", x(1), x(nn)
      else
         read(ichxi,end=900) nn,(xr(l1),l1=1,nn)
         do 811 l1=1,nn
         x(l1)=xr(l1)
  811    continue
      endif
      
      write(*,*) "read in nn", nn
c
c     truncate model if necessary
c
         do 81 l1=1,nn
         xefn(l1) = x(l1)
   81    continue
      nnefn = nn
      if (rmax.le.x(1)) goto 84
      iflag = 0
         do 82 l1=2,nn
         if (x(l1-1).lt.rmax.and.x(l1).ge.rmax) iflag = l1
   82    continue
      if (iflag.eq.0) then 
         nn = nn+1
         x(nn) = rmax
         write(*,*) 'Extending model from ',x(nn-1),' to ',rmax
      else
         write(*,*) 'Truncating model from ',x(nn),' to ',rmax
         nn = iflag
         x(nn) = rmax
      endif

   84 continue
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
      write(*,*) 'about to compute gamma derivative'
      if (icase.ge.11.and.icase.le.17) then
c        first compute d(ln gamma)/dr
         call derivk(x,aa(3,1),dq,1,nn,5,1,2)
         write(*,*) 'back from computing derivative'
            do 255 l1=1,nn
            dq(l1) = dq(l1)/aa(3,l1)
  255       continue
c        then compute d(ln c^2)/dr
         write(*,*) 'about to divide by x'
         write(*,*) 'inonz = ',inonz
         if (icase.ne.12) then
            do 253 l1=inonz,nn
            dq(l1) = dq(l1)+(1.-aa(3,l1))*aa(2,l1)/x(l1)
     +                    + aa(4,l1)/x(l1)
  253       continue
         endif
      endif
      write(*,*) 'done all that'
c
c     read in gamma derivatives if computing Y kernels
c     inserted by Earl Bellinger on Nov 16 2016
c     uses dgam.py to generate dgam.txt which is fed in through ofiles
c
      if (icase.eq.5.or.icase.eq.6) then
         call openf(9,'o','f')
c         call skpcom(9)
c        open(9, file='dgam.txt',status='old')
         nnd = 0
            do 270 l1=1,nnmax
            read(9,*,err=269,end=271) xdum(l1), dgammad(1,l1),
     1         dgammad(2,l1), dgammad(3,l1)
  269       nnd = nnd+1
  270       continue
  271    ii = 3
         id = 3
            do 272 l1=1,nn
            z = x(l1)
            call lir(z,xdum,dgamma(1,l1),dgammad,ii,id,nnd,0,inter)
  272       continue
         write(*,*) 'Finished interpolating dgamma'
            do 280 l1=1,nnmax
            write(*,*) x(l1),xdum(l1),dgammad(1,l1),dgamma(1,l1),
     1        dgammad(2,l1),dgamma(2,l1),dgammad(3,l1),dgamma(3,l1)
  280       continue
         close(9)
      endif
c
c     read in gamma derivatives if computing Y kernels
c
c      if (icase.eq.5.or.icase.eq.6) then
c         call openf(9,'o','f')
c         call skpcom(9)
c         nnd = 0
c            do 270 l1=1,nnmax
c            read(9,*,err=269,end=271) idum,xdum(l1),dum,dum,dum,dum,dum,
c     +        dgammad(1,l1),dgammad(2,l1),dgammad(3,l1)
c  269       nnd = nnd+1
c  270       continue
c  271    ii = 3
c         id = 3
c            do 272 l1=1,nn
c            z = x(l1)
c            call lir(z,xdum,dgamma(1,l1),dgammad,ii,id,nnd,0,inter)
c  272       continue
c         write(*,*) 'Finished interpolating dgamma'
c      endif
c
c     write mesh to kernel files
      call write1(ichkr1,nn,x)
      write(*,*) 'Back from write1'
      call write1(ichkr2,nn,x)
      write(*,*) "write1 ichkr2 nn x", ichkr2, nn
      write(*,*) 'Back from write2'
      call write1(ichpsi,nn,x)
c
c     loop to read in modes and (if suitable degree) compute kernels
c
      icompl = 1
      goto 150
  190 lseek = lmin
  200 continue
      if (iprec1.eq.2) then
         read(ichxi,end=900)
     +      (cs(l1),l1=1,50),(xi0(1,l1),xi0(2,l1),l1=1,nn)
      else
         read(ichxi,end=900)
     +      (csr(l1),l1=1,50),(xi0r(1,l1),xi0r(2,l1),l1=1,nn)
         do 812 l1=1,50
         cs(l1)=csr(l1)
  812    continue
         do 813 l1=1,nn
         xi0(1,l1)=xi0r(1,l1)
         xi0(2,l1)=xi0r(2,l1)
  813    continue
      endif
c---temporary filter in frequency
c     if (cs(27).lt.2.5.or.cs(27).gt.3.5) goto 200
c---end of temporary addition
      el     = cs(18)
      elx    = el+0.5
      iel    = elx
c     
  201 if (iel.gt.lseek) then
         lseek = lseek+lstep
         if (lseek.gt.lmax) goto 900
         goto 201
      endif
c     
      if (iel.lt.lseek) goto 200
c
      ordx = abs(cs(19))+0.1
      iord   = ordx
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
      iteeny=0
         do 256 l1=inonz,nn
         facx = sqrt(aax(1,l1)*x(l1))*x(l1)
         xi0(1,l1)=xi0(1,l1)/facx
         if (iel.ne.0) then
            xi0(2,l1)=xi0(2,l1)/facx/sqrt(el*(el+1.0))
         else
            xi0(2,l1)=0.
         endif
         if (xi0(1,l1)**2+xi0(2,l1)**2.gt.tiny**2.and.iteeny.eq.0) 
     +        iteeny=l1
  256    continue
c.. Following line commented out 22/10/94 after advice from MJT
c..      iteeny=0
      write(*,*) 'iteeny, x(iteeny) = ',iteeny,x(iteeny)
      if (zero1) then
         xi0(1,1) = 0.
c  l=1 case might be improved upon
         if (iel.eq.1) then
            xi0(1,1) = xi0(1,2)
         endif
      xi0(2,1) = xi0(1,1)
      endif
c +++ diagnostic
c      if (iord.lt.10) goto 200
c      diag=.true.
c      if (diag) write(6,702)(x(l1),xi0(1,l1),xi0(2,l1),l1=1,nn)
c  702 format(' r,xi,eta '/(3e15.7))
c      diag=.false.
c      call stop1
c +++ end diagnostic
c
c     calculate integral for denominator
      call calci(x,aa,xi0,nn,el,denom)
      ell1 = el*(el+1.)

      if (icase.eq.-1) then
         write(*,*) 'this option withdrawn - need 2 kernels!'
         write(*,*) '(could just write out from main)'
         stop
c           do 260 l1=1,nn
c           uu = xi0(1,l1)**2 + el*(el+1.)*xi0(2,l1)**2
c           uu = uu * aax(1,l1) * x(l1)**2
c           rk(1,l1) = uu/denom
c 260       continue
c        call write2(x,ichkr1,1,iel,iord,rnu,rk,2,nn)
c        goto 250
      endif
c
c     form d(xi)/dr
c
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
  150 continue

      call kercro(x,nn,aa,aax,xi0,xi0x,rk,2,omv2,
     +                            el,ell1,denom)

      write(*,*) "after kercro"
      if (icase.eq.1) then
         call write2(x,ichkr1,ichkr2,iel,iord,rnu,rk,2,nn,1)
      endif

      if (icase.eq.11) then
         call rholag(x,nn,aax(2,1),3,dq,1,aax,rk,el)
         call write2(x,ichkr1,ichkr2,iel,iord,rnu,rk,2,nn,1)
      endif

      if (icase.eq.2) then
         call kergro(x,nn,aa,aax,rk,rka,2)
         call write2(x,ichkr1,ichkr2,iel,iord,rnu,rka,2,nn,1)
      endif

      if (icase.eq.12) then
         call kergro(x,nn,aa,aax,rk,rka,2)
         call rholag(x,nn,aa(3,1),5,dq,1,aax,rka,el)
         call write2(x,ichkr1,ichkr2,iel,iord,rnu,rka,2,nn,1)
      endif

      if(icase.eq.3.or.icase.eq.4.or.icase.eq.7) then
         call kercg(x,nn,aa,aax,rk,rka)
         if (icase.eq.3.or.icase.eq.4)
     +     call copyrk(rka,rkb,nn)
         if (icase.eq.7) then
              do 851 l1=1,nn
              rkb(1,l1)=2.*rka(1,l1)+rka(2,l1)
  851         continue
         endif
         if (icase.eq.3.or.icase.eq.7)
     +     call copyrk(rka(2,1),rkb(2,1),nn)
         if (icase.eq.4) then
              do 853 l1=1,nn
              rkb(2,l1)=rka(1,l1)+rka(2,l1)
  853         continue
         endif
         call write2(x,ichkr1,ichkr2,iel,iord,rnu,rkb,1,nn,0)
      endif

      if(icase.eq.13.or.icase.eq.14.or.icase.eq.17) then
         write(*,*) 'entering rholag'
         call rholag(x,nn,aax(2,1),3,dq,1,aax,rk,el)
         write(*,*) 'back from rholag'
         call lagcg(x,nn,aa,aax,rk,rka)
         if (icase.eq.13.or.icase.eq.14)
     +     call copyrk(rka,rkb,nn)
         if (icase.eq.17) then
              do 855 l1=1,nn
              rkb(1,l1)=2.*rka(1,l1)+rka(2,l1)
  855         continue
         endif
         if (icase.eq.13.or.icase.eq.17)
     +     call copyrk(rka(2,1),rkb(2,1),nn)
         if (icase.eq.14) then
              do 857 l1=1,nn
              rkb(2,l1)=rka(1,l1)+rka(2,l1)
  857         continue
         endif
         call write2(x,ichkr1,ichkr2,iel,iord,rnu,rkb,1,nn,0)
      endif

      if (icase.eq.5) then
         call kergro(x,nn,aa,aax,rk,rka,2)
         call keruy(x,nn,aa,aax,rka,rk,dgamma,3,psi)
         call write2(x,ichkr1,ichkr2,iel,iord,rnu,rk,2,nn,0)
         call writepsi(ichpsi,iel,iord,rnu,psi,nn)
      endif

      if (icase.eq.6) then
         call kergro(x,nn,aa,aax,rk,rka,2)
         call kerroy(x,nn,aa,aax,rka,rk,dgamma,3)
         call write2(x,ichkr1,ichkr2,iel,iord,rnu,rk,2,nn,1)
      endif


c     return to start of loop, to read in another mode
  250 nk = nk+1
      if (icompl.eq.1) then
         icompl = 0
         goto 190
      else
         goto 200
      endif

c
c     finished computing kernels
c
  900 write(6,*) 'Number of modes for which kernels stored = ',nk
      stop

c
 1000 format(i5)
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
      parameter(nnmax=10000)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension xi0(2,nn), xi0x(2,nn)
      dimension rk(ik,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
      common/comple/rkcomp(2,nnmax),rkcomi,psic(2,nnmax),icompl
      common/work/dummy(nnmax,6)
      iopt = 0
      write(*,*) 'In kercro with l = ',el
c
c     icompl=1 indicates only complementary function required
c
      if (icompl.eq.1) then
         do 90 l1=1,nn
         rk(1,l1) = 0.
         rk(2,l1) = aax(1,l1)*x(l1)**2
   90    continue
         return
      endif
      
         do 100 l1=inonz,nn
         lback=nn-l1+1
         if (l1.lt.iteeny) then
            dummy(lback,1) = 0
            dummy(lback,2) = 0
         else
            dummy(lback,1) = xi0(1,l1) * 
     +        ( 2.*aax(1,l1)*xi0x(2,l1)+aax(3,l1)*xi0(1,l1) )
            dummy(lback,2) = aax(1,l1)/x(l1)**el*(xi0(1,l1)
     1                            -el*xi0(2,l1))
         endif
         dummy(lback,3) = x(l1)
  100    continue
      write(*,*) 'Finished 100 loop'
c   linear extrapolation of dummy(.,2) might be improved upon
      if (zero1) then
        if (iteeny.gt.1) then
           dummy(nn,1) = xi0(1,1) * 
     +        ( 2.*aax(1,1)*xi0x(2,1)+aax(3,1)*xi0(1,1) )
           dummy(nn,2) = (dummy(nn-2,2)*dummy(nn-1,3) 
     +      - dummy(nn-1,2)*dummy(nn-2,3))/(dummy(nn-1,3)-dummy(nn-2,3))
        else
           dummy(nn,1) = 0
           dummy(nn,2) = 0
        endif
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
      write(*,*) 'Got here'
         do 115 l1=inonz,nn
         if (l1.lt.iteeny) then
            dummy(l1,5) = 0.
         else
            dummy(l1,5) = x(l1)**(el+1.)*aax(1,l1)*
     +             (xi0(1,l1)+(el+1.)*xi0(2,l1))
         endif
  115    continue
      if (zero1) dummy(1,5) = 0.
      write(*,*) 'And here'
      call zero(dummy(1,4),nn)
c
      call vintk(x,dummy(1,5),dummy(1,4),1,nn,1,1,2)
c
         do 116 l1=inonz,nn
         if (l1.lt.iteeny) then
            dummy(l1,2) = 0.
            dummy(l1,4) = 0.
         else
            dummy(l1,2) = (el+1.)*dummy(l1,2)
     +        -aax(1,l1)*xi0(1,l1)/x(l1)**(el-1.)
     1        +aax(1,nn)*xi0(1,nn)/x(nn)**(el-1.)
            dummy(l1,4) = -el*dummy(l1,4)
     +        +aax(1,l1)*xi0(1,l1)*x(l1)**(el+2.)
         endif
  116    continue
c   linear extrapolation below might be improved upon
      if (zero1) then
         call interp(x(2),dummy(2,2),1,x(1),dummy(1,2))
         call interp(x(4),dummy(2,4),1,x(1),dummy(1,4))
      endif
c
c
c
      write(*,*) 'Even got here'
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
         if ((zero1.and.l1.eq.1).or.(l1.lt.iteeny)) then
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

c                                                           ************
c                                                           *          *
c                                                           *  interp  *
c                                                           *          *
c                                                           ************
      subroutine interp(xi,yi,iorder,x,y)
c
c     given xi(j) and yi(j), j=1,...,iorder+1 and a value x
c        fits a polynomial of order iorder through xi,yi and evaluate it
c        at x  -  returns value in y
c
      implicit double precision (a-h,o-z)
      dimension xi(*), yi(*)
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
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rkcro(ik,nn),rkgro(ik,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny

      call kergr1(x,nn,aa,aax,rkcro,rkgro(2,1),2,2)
c
c     finally assemble the kernels
c
         do 150 l1=1,nn
         rkgro(1,l1) = rkcro(1,l1)
         rkgro(2,l1) = rkcro(2,l1) - rkcro(1,l1) + rkgro(2,l1)
  150    continue
c    
      return
      end

c                                                           ************
c                                                           *          *
c                                                           *  kergr1  *
c                                                           *          *
c                                                           ************
      subroutine kergr1(x,nn,aa,aax,rk1,rk2,ik1,ik2)
c
c     calculates kernels for ln(gamma1) and ln(rho)
c       given equilibrium state and kernels for ln(c^2) and ln(rho)
c
      implicit double precision(a-h,o-z)
      logical zero1
      parameter(nnmax=10000)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rk1(ik1,nn),rk2(ik2,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
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
         dummy(l1,1) = rk1(1,l1)/press(l1)
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
         do 150 l1=1,nn
         rk2(1,l1) = dummy(l1,1) 
     +                + gconst*dummy(l1,2)*aa(1,l1)*aax(1,l1)*x(l1)
  150    continue
c
      return
      end
c                                                           ************
c                                                           *          *
c                                                           *  kercg   *
c                                                           *          *
c                                                           ************
      subroutine kercg(x,nn,aa,aax,rkcro,rkcg)
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
      parameter (nnmax=10000)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rkcro(2,nn),rkcg(2,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
      common/work/dummy(nnmax,6)

      
c     kercg1 returns psi and (rho r^2)^-1 d(psi)/dr
      call kercg1(x,nn,aa,aax,rkcro(2,1),rkcg,2)
      
         do 100 l1=inonz,nn
         dummy(l1,1) = rkcg(1,l1)
         dummy(l1,2) = rkcg(2,l1)*aax(1,l1)*x(l1)**2
  100    continue
         
         do 110 l1=inonz,nn
         rkcg(2,l1) = dummy(l1,2) + aa(2,l1)*aa(3,l1)/x(l1)*dummy(l1,1)
  110    continue
      if (zero1) rkcg(2,l1) = 0.

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
      subroutine kercg1(x,nn,aa,aax,rkrhs,rkcg,ik)
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
      parameter (nnmax=10000)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rkrhs(ik,nn),rkcg(2,nn)
      dimension ea(2,3),v(2)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23
      common/comple/rkcomp(2,nnmax),rkcomi,psic(2,nnmax),icompl
      
      
      
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
         dummy(l1,6) = rkrhs(1,l1)/aax(1,l1)/x(l1)**2
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
         dummy(l1,3) = dummy(l1,2)-rkrhs(1,l1)
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
         dummy(l1,3) = rkrhs(1,l1)**2
  260    continue
      call zero(dummy(1,5),nn)
      call vintk(x,dummy(1,3),dummy(1,5),1,nn,1,1,2) 
      intrhs = dummy(nn,5)
         do 270 l1=1,nn
         dummy(l1,3) = dummy(l1,2)-rkrhs(1,l1)
     +                    - alpha*aax(1,l1)*x(l1)**2
         dummy(l1,3) = dummy(l1,3)**2
  270    continue
      call zero(dummy(1,5),nn)
      call vintk(x,dummy(1,3),dummy(1,5),1,nn,1,1,2) 
      intres = dummy(nn,5)
      write(*,*) 'Integral of squared residue = ',intres
      write(*,*) 'compared with integrals of squared'
      write(*,*) 'LHS, RHS respectively = ',intlhs,intrhs

      if (icompl.eq.1) then
         do 280 l1=1,nn
         psic(1,l1) = rkcg(1,l1)
         psic(2,l1) = rkcg(2,l1)
  280    continue
      write(*,*) 'stored psic'
      write(*,*) 'z1,z2 at surface = ',psic(1,nn),psic(2,nn)
      else
         if (abs(psic(1,nn)).lt.1.e-15) then
            write(*,*) 'Problem - psic tiny at surface!!!'
            write(*,*) 'psic1,2 = ',psic(1,nn),psic(2,nn)
            stop
         endif
         aaa = -rkcg(1,nn)/psic(1,nn)
         do 290 l1=1,nn
         rkcg(1,l1) = rkcg(1,l1) + aaa*psic(1,l1)
         rkcg(2,l1) = rkcg(2,l1) + aaa*psic(2,l1)
  290    continue
      endif


      return
      end

     
      subroutine rhs(x,y,zk,ap,aq,f,fd,h,hd,ifd,n)
      implicit real*8 (a-h,o-z)
      logical zero1
      parameter (nnmax=10000)
      dimension y(*),f(*),fd(ifd,*),h(1),hd(ifd,1),zk(1),ap(1),aq(1)
      common/centre/zero1,inonz,iteeny
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
      parameter (nnmax=10000)
      dimension y1(*),y2(*),g(*),gd(ig,2),zk(1),ap(1),aq(1)
      common/centre/zero1,inonz,iteeny
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23
      common/comple/rkcomp(2,nnmax),rkcomi,psic(2,nnmax),icompl

      gd(1,1) = -3.
      gd(1,2) = x1*a12(inonz)
      g(1) = gd(1,1)*y1(1) + gd(1,2)*y1(2)
c        gd(2,1) = 1.0
c        gd(2,2) = 0.
c        g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) - 0.0
         gd(2,1) = 0.
         gd(2,2) = 1.
         g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) - g23
      return
      end

c                                                           ************
c                                                           *          *
c                                                           *  rholag  *
c                                                           *          *
c                                                           ************
      subroutine rholag(x,nn,q,idq,dq,iddq,aax,rk,el)
c
c     convert Eulerian kernels for (q,rho) into Lagrangian kernels for (q,rho)
c
c     dq contains d(ln q)/dr
c
      implicit double precision (a-h,o-z)
      logical diag,zero1
      parameter (nnmax=10000)
      dimension x(nn),q(idq,nn),dq(iddq,nn),aax(3,nn),rk(2,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
      common/comple/rkcomp(2,nnmax),rkcomi,psic(2,nnmax),icompl
      common/work/dummy(nnmax,6)
  
      diag = .false.

c
c     icompl=1 indicates only complementary function required
c
      if (icompl.eq.1) then
         do 90 l1=1,nn
         rk(1,l1) = 0.
         rk(2,l1) = x(l1)**2
   90    continue
         return
      endif
      
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
c                                                           ************
c                                                           *          *
c                                                           *  lagcg   *
c                                                           *          *
c                                                           ************
      subroutine lagcg(x,nn,aa,aax,rkcro,rkcg)
c
c     calculates Lagrangian kernels for ln(c^2) and ln(gamma1) 
c       given equilibrium state and kernels for ln(c^2) and ln(rho)
c
c     solves p d(psi)/dr  -4 G r^2 int_r^R (m rho /(r^5) psi) dr  
c                = K_{rho,c^2} 
c
c     (for more details see s/r lagcg1)
c
c     and then K_{gamma,c^2} = p d(psi)/dr
c              K_{c^2,gamma} = K_{c^2,rho}-K_{gamma,c^2}
c
      implicit double precision(a-h,o-z)
      logical zero1
      parameter (nnmax=10000)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rkcro(2,nn),rkcg(2,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
      common/work/dummy(nnmax,6)

      
c     lagcg1 returns psi and (p / r^2) d(psi)/dr
      call lagcg1(x,nn,aa,aax,rkcro(2,1),rkcg,2)
      
         do 100 l1=inonz,nn
         rkcg(2,l1) = rkcg(2,l1)*x(l1)**2
  100    continue
      if (zero1) rkcg(2,l1) = 0.

         do 120 l1=1,nn
         rkcg(1,l1) = rkcro(1,l1)-rkcg(2,l1)
  120    continue

      return
      end
c                                                           ************
c                                                           *          *
c                                                           *  lagcg1  *
c                                                           *          *
c                                                           ************
      subroutine lagcg1(x,nn,aa,aax,rkrhs,rkcg,ik)
c
c     solves p d(psi)/dr  -4 G r^2 int_r^R (m rho /(r^5) psi) dr  
c                = K_{rho,c^2} 
c
c     the dependent variables used are 
c          z1 = psi,     z2 = (p / r^2) d(psi)/dr
c
c     these are returned in rkcg(1,.) and rkcg(2,.) respectively
c
      implicit double precision(a-h,o-z)
c.... double precision intlhs,intrhs,intres
      integer*4 v
      logical zero1
      parameter (nnmax=10000)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rkrhs(ik,nn),rkcg(2,nn)
      dimension ea(2,3),v(2)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23
      common/comple/rkcomp(2,nnmax),rkcomi,psic(2,nnmax),icompl
      common/work/dummy(nnmax,6)
      external rhs1,bc1

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
         a12(l1) = x(l1)**2 / press
         a21(l1) = -4.0*gconst*aa(1,l1)*aax(1,l1)/(x(l1)**2)
         a22(l1) = 0.
  100    continue
c     put d(Kroc/(r^2))/dr into b2
         do 105 l1=inonz,nn
         dummy(l1,6) = rkrhs(1,l1)/x(l1)**2
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
      call nrk(x(inonz),rkcg(1,inonz),zk,ap,aq,rhs1,bc1,ii,kk,ka,kb,ki,
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

      if (icompl.eq.1) then
         do 280 l1=1,nn
         psic(1,l1) = rkcg(1,l1)
         psic(2,l1) = rkcg(2,l1)
  280    continue
      write(*,*) 'stored psic'
      write(*,*) 'z1,z2 at surface = ',psic(1,nn),psic(2,nn)
      else
         if (abs(psic(1,nn)).lt.1.e-15) then
            write(*,*) 'Problem - psic tiny at surface!!!'
            write(*,*) 'psic1,2 = ',psic(1,nn),psic(2,nn)
            stop
         endif
         aaa = -rkcg(1,nn)/psic(1,nn)
         do 290 l1=1,nn
         rkcg(1,l1) = rkcg(1,l1) + aaa*psic(1,l1)
         rkcg(2,l1) = rkcg(2,l1) + aaa*psic(2,l1)
  290    continue
      endif

      return
      end

     
      subroutine rhs1(x,y,zk,ap,aq,f,fd,h,hd,ifd,n)
      implicit real*8 (a-h,o-z)
      logical zero1
      parameter (nnmax=10000)
      dimension y(*),f(*),fd(ifd,*),h(1),hd(ifd,1),zk(1),ap(1),aq(1)
      common/centre/zero1,inonz,iteeny
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
 
      subroutine bc1(x1,x2,y1,y2,zk,ap,aq,g,gd,ig,id,n)
      implicit real*8(a-h,o-z)
      logical zero1
      parameter (nnmax=10000)
      dimension y1(*),y2(*),g(*),gd(ig,*),zk(1),ap(1),aq(1)
      common/centre/zero1,inonz,iteeny
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23
      common/comple/rkcomp(2,nnmax),rkcomi,psic(2,nnmax),icompl

      gd(1,1) = -3.
      gd(1,2) = x1*a12(inonz)
      g(1) = gd(1,1)*y1(1) + gd(1,2)*y1(2)
c        gd(2,1) = 0.
c        gd(2,2) = 1.
c        g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) - 0.0
         gd(2,1) = 0.
         gd(2,2) = 1.
         g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) - g23
      return
      end
c                                                           ************
c                                                           *          *
c                                                           *  keruy   *
c                                                           *          *
c                                                           ************
      subroutine keruy(x,nn,aa,aax,rkgro,rkuy,dgamma,idgam,psi)
c
c     calculates kernels for ln(u) and Y
c       given equilibrium state and kernels for ln(gamma1) and ln(rho)
c
c     let gam,p   = d(ln gam)/d(ln p)
c         gam,rho = d(ln gam)/d(ln rho)
c         gam,Y   = d(ln gam)/d(Y)
c     
c     solves d(psi)/dr  -4 pi G rho r^2 int_r^R (rho /(r^2 p) psi) dr
c                = - ( K_{rho,gam} + (gam,p + gam,rho)*K_{gam,rho} )
c
c     warning! 
c     this subroutine overwrites K_{rho,gam}, but this could
c     easily be modified given another array of length nn
c
c     then K_{Y,u} = gam,Y * K_{gam,rho}
c          K_{u,Y} = p d(p^-1 psi)/dr + gam,p*K_{gam,rho}
c
      implicit double precision(a-h,o-z)
      logical zero1
      parameter (nnmax=10000)
      dimension x(nn), psi(nn), aa(5,nn), aax(3,nn), dgamma(idgam,nn)
      dimension rkgro(2,nn),rkuy(2,nn)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
      common/work/dummy(nnmax,6)
 
 
c
c     set up right hand side
         do 100 l1=1,nn
         rkgro(2,l1) = -( rkgro(2,l1) + 
     +          (dgamma(1,l1)+dgamma(2,l1))*rkgro(1,l1) )
  100    continue

c     kercg1 returns psi and (rho r^2)^-1 d(psi)/dr
      call kercg1(x,nn,aa,aax,rkgro(2,1),rkuy,2)
      
         do 101 l1=inonz,nn
         psi(l1) = rkuy(1,l1)
  101     continue 
      
         do 105 l1=inonz,nn
         dummy(l1,1) = rkuy(1,l1)
         dummy(l1,2) = rkuy(2,l1)*aax(1,l1)*x(l1)**2
  105    continue
         
         do 110 l1=inonz,nn
         rkuy(1,l1) = dummy(l1,2) + aa(2,l1)*aa(3,l1)/x(l1)*dummy(l1,1)
     +                  + dgamma(1,l1)*rkgro(1,l1)
  110    continue
      if (zero1) rkuy(1,l1) = 0.
 
         do 120 l1=1,nn
         rkuy(2,l1) = dgamma(3,l1)*rkgro(1,l1)
  120    continue
 
      return
      end
c                                                           ************
c                                                           *          *
c                                                           *  kerroy  *
c                                                           *          *
c                                                           ************
      subroutine kerroy(x,nn,aa,aax,rkgro,rkroy,dgamma,idgam)
c
c     calculates kernels for ln(rho) and Y
c       given equilibrium state and kernels for ln(gamma1) and ln(rho)
c
c     (see header to s/r keruy)
c
      implicit double precision(a-h,o-z)
      logical zero1
      parameter (nnmax=10000)
      dimension x(nn), aa(5,nn), aax(3,nn), dgamma(idgam,nn)
      dimension rkgro(2,nn),rkroy(2,nn)
      dimension wspace(nnmax)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz,iteeny
      common/work/dummy(nnmax,6)

         do 100 l1=1,nn
         wspace(l1) = dgamma(1,l1)*rkgro(1,l1)
  100    continue

      call kergr1(x,nn,aa,aax,wspace,rkroy,1,2)
c
c     finally assemble the kernels
c
         do 150 l1=1,nn
         rkroy(1,l1) = rkgro(2,l1) + dgamma(2,l1)*rkgro(1,l1) 
     +                                               + rkroy(1,l1)
         rkroy(2,l1) = dgamma(3,l1)*rkgro(1,l1)
  150    continue
c
      return
      end
      subroutine write1(ich,nn,x)
      implicit double precision (a-h,o-z)
      parameter(nnmax=10000)
      real rx(nnmax)
      dimension x(nn)
         do 100 l1=1,nn
         rx(l1)=x(l1)
  100    continue
      write(ich) nn,(rx(l1),l1=1,nn)
      return
      end
      
      subroutine writepsi(ichpsi,iel,iord,rnu,psi,nn)
      implicit double precision (a-h,o-z)
      parameter(nnmax=10000)
      real ppsi(nnmax), rrnu
      dimension psi(nn)
      write(*,*) 'Entered writepsi'
      rrnu = rnu
            do 134 l1=1,nn
            ppsi(l1) = psi(l1)
  134       continue
      write(ichpsi) iel,iord,rrnu,(ppsi(l1),l1=1,nn)
      return
      end
      
      subroutine write2(x,ich1,ich2,iel,iord,rnu,rk,ik,nn,itype)
      implicit double precision (a-h,o-z)
      parameter(nnmax=10000)
      common/comple/rkcomp(2,nnmax),rkcomi,psic(2,nnmax),icompl
      common/work/dummy(nnmax,6)
      real rrk(nnmax),rrnu
      dimension rk(2,nn),x(nnmax)
      write(*,*) 'Entered write2'
c
c     write out kernel to disk (but not before taking out a multiple of the 
c     complementary function if itype=1
c
      if (icompl.eq.1.and.itype.eq.0) return

      if (icompl.eq.1.and.itype.eq.1) then
         do 100 l1=1,nn
         dummy(l1,1) = rk(1,l1)**2+rk(2,l1)**2
  100    continue
         dummy(1,2)=0.0
         write(*,*) "before vintk"
         write(*,*) "nn", nn
         call vintk(x,dummy,dummy(1,2),1,nn,1,1,2)
         write(*,*) "after vintk"
         rkcomi = dummy(nn,2)
         do 110 l1=1,nn
         rkcomp(1,l1) = rk(1,l1)
         rkcomp(2,l1) = rk(2,l1)
  110    continue
         rrnu = -1.0
         write(*,*) 'written compl. fn.'
         write(*,*) '   int. sq. = ',rkcomi
         do 115 l1=1,nn
         rrk(l1) = rk(1,l1)
  115    continue
         write(ich1) -1,-1,rrnu,(rrk(l1),l1=1,nn)
         do 117 l1=1,nn
         rrk(l1) = rk(2,l1)
  117    continue
         write(ich2) -1,-1,rrnu,(rrk(l1),l1=1,nn)
      else 
         if (itype.eq.1) then
            do 120 l1=1,nn
            dummy(l1,1) = rkcomp(1,l1)*rk(1,l1)
     +                    + rkcomp(2,l1)*rk(2,l1)
  120       continue
            dummy(1,2) = 0.0
            call vintk(x,dummy,dummy(1,2),1,nn,1,1,2)
            alpha = dummy(nn,2)/rkcomi

            write(*,*) 'alpha = ',alpha
            do 130 l1=1,nn
            rrk(l1) = rk(1,l1)-alpha*rkcomp(1,l1)
  130       continue
         else
            do 131 l1=1,nn
            rrk(l1) = rk(1,l1)
  131       continue
         endif
         rrnu = rnu
         write(ich1) iel,iord,rrnu,(rrk(l1),l1=1,nn)
         if (itype.eq.1) then
            do 132 l1=1,nn
            rrk(l1) = rk(2,l1)-alpha*rkcomp(2,l1)
  132       continue
         else
            do 133 l1=1,nn
            rrk(l1) = rk(2,l1)
  133       continue
         endif
         write(ich2) iel,iord,rrnu,(rrk(l1),l1=1,nn)
      endif
      return
      end

      subroutine copyrk(rka,rkb,nn)
      implicit double precision (a-h,o-z)
      dimension rka(2,nn),rkb(2,nn)
         do 100 l1=1,nn
         rkb(1,l1) = rka(1,l1)
  100    continue
      return
      end
