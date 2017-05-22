      subroutine linint(x,y,rhs,ii,iy,ig,isn,nd1,nn)
c
c  Initial value integrator.
c  =========================
c  Integrates a set of first-order ordinary
c  linear homogeneous differential equations, by means of second-order
c  centred difference approximation, from given initial values.
c
c  The equations are on the form
c
c       d y(i; x)/dx = sum ( a(i,j; x) * y(j; x))
c                       j
c
c  The order of the equations is ii.
c  
c  The coefficients are set up in the routine rhs, which must be
c  supplied by the user as external.
c  
c  The routine has the option for determining in a single call the
c  solutions to a given set of equations for several initial conditions;
c  ig specifies the number of such sets.
c  
c  For use with Richardson extrapolation the routine finds the solution
c  at every isn-th point in the input mesh. For ordinary use 
c  isn is set to 1.
c  
c  nd1 specifies the initial value of the mesh index passed into the
c  routine rhs (see below).
c
c  iy is the first dimension of y in the calling programme.
c  
c  On input x(1+isn*(n-1)), n = 1,nn, must contain the mesh in the 
c  independent variable. The initial conditions for set k must be
c  set into y(i + ii*(k-1),1), i = 1,ii, for k = 1,ig.
c  
c  The routine returns the solution for initial value set k in
c  y(i+ii*(k-1),1+isn*(n-1)), i=1,ii, n=1,nn, k=1,ig.
c  
c  The right hand side subroutine rhs:
c  -----------------------------------
c  
c  This is called by linint and should be defined as
c  
c        subroutine rhs(x,aa,iaa,nd)
c        dimension aa(iaa,1)
c           .
c           .
c           .
c  
c  A single call must set the coefficient matrix at a single point.
c  On input x (a scalar) gives the value of the independent variable
c  at the given point; iaa is the first dimension of aa;
c  nd is passed from linint and may be used to address a data
c  array in rhs. 
c  
c  The routine should return
c  
c  aa(i,j) = a(i,j; x)
c  
c  where a(i,j; x) as defined above is the coefficient matrix in
c  the equations.
c  
c  The definition of nd is slightly convoluted. It is assumed that the
c  use might need certain variables to define the right hand side
c  of the equations. Assume, for example, that these are passed into 
c  the routine in
c  
c        common/rhsdat/ data(10,200)
c  
c  When rhs is called at the meshpoint x(1+isn*(n-1)), nd is
c  set to nd1+isn*(n-1). The meshpoint should therefore correspond
c  to the data in data(i,nd). In most cases presumably x and the
c  data would be given on the same mesh, and nd1 would be 1.
c   
      dimension x(nn),y(iy,nn),fd(20,20),w(20,20),y1(20)
      common/modfac/ r21,ntld   
      external rhs  
      ifd=20
      iw=20 
c  right hand sides at first point  
    5 n=1   
      nc=1  
      nd=nd1
      ntst=nn-2 
      iig=ii*ig 
      x2=x(1)   
      call rhs(x2,fd,ifd,nd)
      go to 20  
c  right hand sides at next point   
   10 nd=nd+isn 
      call rhs(x2,fd,ifd,nd)
c  set coefficient matrix for linear equations  
      do 15 i=1,ii  
      do 13 j=1,ii  
   13 w(i,j)=dx*fd(i,j) 
   15 w(i,i)=1+w(i,i)   
c  solve linear equations   
      call leq(w,y(1,n),ii,ig,iw,ii,err)
      if(nc.eq.nn) return   
c  set up right hand side of linear equation at next point  
   20 n1=n  
      n=n+isn   
      nc=nc+1   
      x1=x2 
      x2=x(n)   
      dx=0.5*(x1-x2)
      kg=-ii
      do 25 l=1,ig  
      kg=kg+ii  
      do 25 i=1,ii  
      k=kg+i
      sum=0 
      do 22 j=1,ii  
   22 sum=sum+fd(i,j)*y(kg+j,n1)
   25 y(k,n)=y(k,n1)-dx*sum 
      if(ntld.ne.1) go to 10
      if(nc.ne.ntst.or.ig.ne.2) go to 10
c  test for linear dependence (only implemented for ig = 2) 
      do 30 i=1,iig 
   30 y1(i)=y(i,n1) 
      aym=0 
      im=0  
      aym2=0
      do 32 i=1,ii  
      ay=abs(y1(i)) 
      aym2=amax1(aym2,abs(y1(i+ii)))
      if(ay.le.aym) go to 32
      aym=ay
      im=i  
   32 continue  
c   
      r21=y1(ii+im)/y1(im)  
      s1=0  
      s2=0  
      j=ii  
      do 34 i=1,ii  
      j=j+1 
      yi=y1(i)  
      if(yi.ne.0) go to 33  
      if(abs(y1(j)).gt.1.e-5*aym2) go to 10 
      go to 34  
   33 ri=y1(j)/y1(i)
      s1=s1+abs(r21-ri) 
      s2=s2+abs(ri) 
   34 continue  
c  test 
      if(s1/s2.gt.1.e-7) go to 10   
c  nearly linear dependence. modify initial values and start again  
      j=0   
      aym=0 
      do 36 i=1,ii  
      y1(i)=y(i+ii,1)-r21*y(i,1)
   36 aym=amax1(aym,abs(y1(i))) 
c  normalize initial y's
      aym=1./aym
      do 38 i=1,ii  
   38 y(ii+i,1)=aym*y1(i)   
c   
      write(6,100) (y(i,1),i=1,iig) 
      go to 5   
c   
  100 format(//' initial solution has been modified. new y(n=1):'/  
     *  (1p10e13.5))
c   
      end   
