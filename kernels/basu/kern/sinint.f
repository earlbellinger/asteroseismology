      subroutine sinint(x,f,x0,anu,y,nn,if,iy)
c     sets integral of the type
c
c     (1)  int(f(x)(x - x0)**nu)dx,    or
c     (2)  int(f(x)(x0 - x)**nu)dx, 
c
c     from x(1) to x(n), where f(x(n))*(x - x0)**nu = f(1,n) or
c     f(x(n))*(x0 - x)**nu = f(1,n), into  real y(1,n),  n=1,nn   
c
c     When nu is a non-negative integer, there are no restrictions
c     on how x0 is located realtive to the integration interval.
c     Otherwise x0 must be outside the interval of integration and
c     the appropriate of types (1) or (2) is assumed.
c
c     If x0 is within the integration interval in this case, an error
c     message is wriotten, and the integral is evaluated using the
c     trapezoidal rule.
c
c     the independent variable x, which need not be uniformly divided
c     or increasing with n (but must be monotonic), must be supplied 
c     by the calling programme.
c
      logical intnu
      dimension x(1),f(if,1),y(iy,1)
c
c  test for integration type, etc.
c
      if(anu.ge.0.and.abs(nint(anu)-anu).lt.1.e-6) then
        nu=nint(anu)
        anu=nu
        intnu=.true.
c..        write(6,*) 'Integer formulation'
      else
        intnu=.false.
      end if
c
      if(.not.intnu) then
c
c  test for location of x0 relative to integration interval
c
        if((x(1)-x0)*(x(nn)-x0).le.0) then
          write(6,110) anu, x0, x(1), x(nn)
          y(1,1)=0
          do 20 n=2,nn
          n1=n-1
   20     y(1,n)=y(1,n1)+0.5*(x(n)-x(n1))*(f(1,n)+f(1,n1))
          return
        end if
      end if
c
      y(1,1)=0
      delx2=x(1)-x0
      adelx2=abs(delx2)
      if(intnu) then
        if(nu.eq.0) then
          xnu2=1
        else
          xnu2=delx2**nu
        end if
        if(xnu2.eq.0) then
          ft2=0
        else
          ft2=f(1,1)/xnu2
        end if
      else
        xnu2=abs(delx2)**anu
        ft2=f(1,1)/xnu2
      end if
c
c  start loop in x
c
      do 50 n=2,nn
c..      write(6,*) 'n, x(n), f(1,n)',n, x(n), f(1,n)
      n1=n-1
      ft1=ft2
      xnu1=xnu2
      delx1=delx2
      adelx1=adelx2
      delx2=x(n)-x0
      adelx2=abs(delx2)
      if(intnu) then
        if(nu.eq.0) then
          xnu2=1
        else
          xnu2=delx2**nu
        end if
        if(xnu2.eq.0) then
          ft2=0
        else
          ft2=f(1,n)/xnu2
        end if
      else
        xnu2=abs(delx2)**anu
        ft2=f(1,n)/xnu2
      end if
c
c  test for expansion
c
      deltax=x(n)-x(n1)
      if(delx1.ne.0) then
        eps=deltax/delx1
      else
        eps=1
      end if
      if(abs(eps).lt.0.01) then
c
        aint0=deltax*xnu1*(1+anu*eps*(0.5+eps*(anu-1)))
        aint1=0.5*deltax*xnu1*(1+0.66666667*anu*eps)
c
c  test for different cases of nu
c
      else if(anu.eq.-1) then
c
        aint0=alog((x(n)-x0)/(x(n1)-x0))
        aint1=1-alog((x(n)-x0)/(x(n1)-x0))/eps
c
      else if(anu.eq.-2) then
c
        aint0=(xnu2*delx2/(xnu1*delx1)-1)*delx1*xnu1/(anu+1)
        aint1=(alog((x(n)-x0)/(x(n1)-x0))/deltax-1./delx2)
c
      else if(intnu) then
c
c  non-negative integral nu
c
        aint0=(xnu2*delx2-xnu1*delx1)/(anu+1)
        aint1=(xnu2*delx2-(xnu2*delx2*delx2-xnu1*delx1*delx1)/
     *    ((anu+2)*deltax))/(anu+1)
c
      else
c
c  general case
c
        rat1=xnu2*adelx2/(xnu1*adelx1)
        rat2=rat1*adelx2/adelx1
        aint0=(rat1-1)*delx1*xnu1/(anu+1)
        aint1=(rat1-(rat2-1)*delx1/((anu+2)*deltax))*delx1*xnu1/(anu+1)
c
      end if
c
c  now partial integrals are set. set final contribution
c
   50 y(1,n)=y(1,n1)+ft1*aint0+(ft2-ft1)*aint1
c
      return
  110 format(//' **** error in sinint. for nu =',1pe13.5,
     *  ' x0 =',e13.5,' is within'/
     *  '      integration interval x(1), x(nn) =',2e13.5/
     *  '      integral evaluated with trapezoidal rule')
      end
