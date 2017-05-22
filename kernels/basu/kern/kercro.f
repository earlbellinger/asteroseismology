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
         dummy(2000+nn-l1+1) = x(l1)
  100    continue
      call vinta(dummy(2001),dummy,dummy(4001),nn,1,1)
         do 110 l1=1,nn
         dummy(l1) = dummy(4000+nn-l1+1)
  110    continue
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
c        write(6,*) r1,r2,r12,r3,r4,r5
         rk(irow,l1) = (r12+r3+r4+r5)/2.
  140    continue
c     do 150 l1=1,nn
c     write(6,*) x(l1),rk(1,l1),rk(2,l1)
c 150 continue
      return
      end
