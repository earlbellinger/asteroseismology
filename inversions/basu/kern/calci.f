      subroutine calci(x,aa,xi0,nn,el,denom)
c
c     calculate integral I of ( xi^2 + L^2 eta^2 ) rho r^2 dr 
c     (everything dimensionless) and return answer in denom
c
      implicit double precision (a-h,o-z)
      logical diag
      parameter (nnmax=1500)
      dimension x(nn), aa(5,nn), xi0(2,nn)
      common/work/t1(nnmax), t2(nnmax)
      diag = .false.
      ell1 = el*(el+1.)
      fourpi = acos(0.d0)*8.
         do 100 l1=1,nn
         rho = aa(1,l1)*aa(5,l1)/fourpi
         t1(l1) = ( xi0(1,l1)**2 + ell1*xi0(2,l1)**2 )*rho*x(l1)**2
  100    continue
      call vinta(x,t1,t2,nn,1,1)
      denom = t2(nn)
      diag = .true.
c +++ diagnostic
      if(diag) write(6,701) denom
  701 format('Integral I = ',e15.7)
c +++ end diagnostic
      diag = .false.
      return
      end
