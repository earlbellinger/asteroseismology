      subroutine sqinta(x,a,alpha,y,nn,ia,id)
c
c  integrates function of the form sqrt(1 - (alpha/x)**2)*a(1,x).
c
c  result is put in y(1,x).
c  nn is the number of mesh points.
c  ia and id are first dimensions of a and y.
c
c  original version: 1/11/86
c
      dimension x(1),a(ia,1),y(id,1)
      save pi2
      data pi2 /1.5707963/
c
      y(1,1)=0
      alpha2=alpha*alpha/2
c
      do 40 n=1,nn
      if(x(n).le.alpha) then
        ac=0
        bc=0
      else
        xa=x(n)/alpha
        xa2=xa*xa
        sxa2=sqrt(xa2-1)
        ac=alpha*(sxa2-pi2+asin(1./xa))
        bc=alpha2*(xa*sxa2-alog(xa+sxa2))
      end if
c
      if(n.ne.1) then
c
c  set contribution to integral
c
        n1=n-1
        da=(a(1,n)-a(1,n1))/(x(n)-x(n1))
        y(1,n)=y(1,n1)+(a(1,n1)-x(n1)*da)*(ac-acp)+da*(bc-bcp)
      end if
c
      acp=ac
      bcp=bc
c
   40 continue
      return
      end
