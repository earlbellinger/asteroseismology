      subroutine surf1 (m,n,x,y,z,iz,zx1,zxm,zy1,zyn,zxy11,zxym1,zxy1n,
     1                  zxymn,zp,temp,sigma)   
c  
c  
c dimension of           x(m),y(n),z(iz,n),zx1(n),zxm(n),zy1(m),zyn(m),
c arguments              zp(3*m*n),temp(2*m+n) 
c  
c latest revision        february 7, 1974  
c  
c purpose                surf1 determines the parameters necessary to  
c                        compute an interpolatory surface passing  
c                        through a rectangular grid of functional  
c                        values.  the surface determined can be
c                        represented as the tensor product of splines  
c                        under tension.  the x- and y-derivatives  
c                        around the boundary and the xy derivatives at 
c                        the four corners may be specified or omitted. 
c                        for actual mapping of points onto the surface 
c                        it is necessary to call the function surf2.   
c  
c access cards           *fortran,s=ulib,n=surf
c                        *cosy 
c  
c usage                  call surf1 (m,n,x,y,z,iz,zx1,zxm,zy1,zyn,zxy11,   
c                                    zxym1,zxy1n,zxymn,zp,temp,sigma)  
c  
c arguments
c  
c on input               m 
c                          is the number of grid lines in the  
c                          x-direction, i.e., lines parallel to the
c                          y-axis (m .ge. 2).  
c  
c                        n 
c                          is the number of grid lines in the  
c                          y-direction, i.e., lines parallel to the
c                          x-axis (n .ge. 2).  
c  
c                        x 
c                          is an array containing the x-coordinates of 
c                          the m grid lines in the x-direction.
c  
c                        y 
c                          is an array containing the y-coordinates of 
c                          the n grid lines in the y-direction.
c  
c                        z 
c                          is an array containing the m*n functional   
c                          values at the grid points, i.e., z(i,j) 
c                          contains the functional value at (x(i),y(j))
c                          i = 1,...,m j = 1,...,n.
c  
c                        iz
c                          is the fortran row dimension of the matrix z
c                          used in the calling program (iz .ge. m).
c  
c                        zx1   
c                          is an array of x-derivatives of the function
c                          along the x(1) grid line, i.e., zx1(j)  
c                          contains the x-partial derivative at the
c                          point (x(1),y(j)) j = 1,...,n.  (this   
c                          parameter is ignored if sigma .lt. 0.). 
c  
c                        zxm   
c                          is an array of x-derivatives of the function
c                          along the x(m) grid line, i.e., zxm(j)  
c                          contains the x-partial derivative at the
c                          point (x(m),y(j)) j = 1,...,n.  (this   
c                          parameter is ignored if sigma .lt. 0.). 
c  
c                        zy1   
c                          is an array of y-derivatives of the function
c                          along the y(1) grid line, i.e., zy1(i)  
c                          contains the y-partial derivative at the
c                          point (x(i),y(1)) i = 1,...,m.  (this   
c                          parameter is ignored if sigma .lt. 0.). 
c  
c                        zyn   
c                          is an array of y-derivatives of the function
c                          along the y(n) grid line, i.e., zyn(i)  
c                          contains the y-partial derivative at the
c                          point (x(i),y(n)) i = 1,...,m.  (this   
c                          parameter is ignored if sigma .lt. 0.). 
c  
c                        zxy11, zxym1, zxy1n and zxymn 
c                          are the xy-derivatives of the function at the   
c                          four corners of the grid, i.e., 
c                          zxy11 (zxym1,zxy1n,zxymn) contains the x-y  
c                          second partial derivative at the point  
c                          (x(1),y(1)) ((x(m),y(1)),(x(1),y(n)),   
c                          (x(m),y(n)), respectively).  (these 
c                          parameters are ignored if sigma .lt. 0.).   
c  
c                        zp
c                          is an array of at least 3*m*n locations.
c  
c                        temp  
c                          is an array of at least 2*m+n locations used
c                          for internal temporary storage. 
c  
c                        sigma 
c                          contains the tension factor.  this is   
c                          non-zero and indicates the curviness desired.   
c                          if abs(sigma) is nearly zero (e.g., .001) the   
c                          resulting surface is a bi-cubic spline.  if 
c                          abs(sigma) is large (e.g., 50.) the resulting   
c                          surface is nearly bi-linear.  if sigma is   
c                          negative, values for the boundary and corner
c                          partial derivatives will be determined  
c                          internally and only place-holding (dummy)   
c                          parameters need be given in the call for
c                          arguments zx1, zxm, zy1, zyn, zxy11, zxym1, 
c                          zxy1n and zxymn.  a standard value for sigma
c                          is approximately 1. in absolute value.  
c  
c on output              zp
c                          contains values proportional to xx-, yy- and
c                          xxyy-partial derivatives at the grid points 
c                          to be used for the actual interpolation.
c  
c                        m, n, x, y, z, iz, zx1, zxm, zy1, zyn, zxy11, 
c                        zxym1, zxy1n, zxymn and sigma are unaltered.  
c  
c entry points           surf1 
c  
c special conditions     none  
c  
c common blocks          none  
c  
c i/o                    none  
c  
c precision              single
c  
c required ulib          none  
c routines 
c  
c specialist             russell k. rew, ncar, boulder, colorado  80302
c  
c language               fortran   
c  
c history                originally written by a. k. cline,
c                        september 1973.   
c  
c  
c  
c  
      integer         n          ,m          ,iz   
      real            x(m)       ,y(n)       ,z(iz,n)    ,zx1(n)     , 
     1                zxm(n)     ,zy1(m)     ,zyn(m)     ,zxy11      , 
     2                zxym1      ,zxy1n      ,zxymn      ,zp(m,n,3)  , 
     3                temp(1)    ,sigma
c  
      fxp(t)=exp(amax1(-87.,amin1(87.,t))) 
c  
      mm = m   
      nn = n   
      mm1 = mm-1   
      mp1 = mm+1   
      nm1 = nn-1   
      np1 = nn+1   
      npm = nn+mm  
c  
c check for given derivatives  
c  
      if (sigma .lt. 0.) go to 103 
c  
c fill given derivatives into zp and temp  
c  
      do 101 i=1,mm
         zp(i,1,1) = zy1(i)
         npi = nn+i
         temp(npi) = zyn(i)
  101 continue 
      do 102 j=1,nn
         zp(1,j,2) = zx1(j)
         npmpj = npm+j 
         temp(npmpj) = zxm(j)  
  102 continue 
      zp(1,1,3) = zxy11
      zp(mm,1,3) = zxym1   
      zxy1ns = zxy1n   
      zxymns = zxymn   
      go to 113
c  
c derivatives not given, use second order interpolation
c polynomial on input values for derivatives   
c  
  103 if (nn .eq. 2) go to 106 
      del1 = y(2)-y(1) 
      del2 = y(3)-y(2) 
      del12 = y(3)-y(1)
      c1 = -(del12+del1)/del12/del1
      c2 = del12/del1/del2 
      c3 = -del1/del12/del2
      do 104 i=1,mm
         zp(i,1,1) = c1*z(i,1)+c2*z(i,2)+c3*z(i,3) 
  104 continue 
      deln = y(nn)-y(nm1)  
      delnn = y(nn)-y(nn-2)
      delnm1 = y(nm1)-y(nn-2)  
      c1 = deln/delnn/delnm1   
      c2 = -delnn/deln/delnm1  
      c3 = (delnn+deln)/delnn/deln 
      do 105 i=1,mm
         npi = nn+i
         temp(npi) = c1*z(i,nn-2)+c2*z(i,nm1)+c3*z(i,nn)   
  105 continue 
      go to 108
  106 c1 = 1./(y(2)-y(1))  
      do 107 i=1,mm
         zp(i,1,1) = c1*(z(i,2)-z(i,1))
         npi = nn+i
         temp(npi) = zp(i,1,1) 
  107 continue 
  108 if (mm .eq. 2) go to 111 
      del1 = x(2)-x(1) 
      del2 = x(3)-x(2) 
      del12 = x(3)-x(1)
      c1 = -(del12+del1)/del12/del1
      c2 = del12/del1/del2 
      c3 = -del1/del12/del2
      do 109 j=1,nn
         zp(1,j,2) = c1*z(1,j)+c2*z(2,j)+c3*z(3,j) 
  109 continue 
      zp(1,1,3) = c1*zp(1,1,1)+c2*zp(2,1,1)+c3*zp(3,1,1)   
      zxy1ns = c1*temp(nn+1)+c2*temp(nn+2)+c3*temp(nn+3)   
      deln = x(mm)-x(mm1)  
      delnm1 = x(mm1)-x(mm-2)  
      delnn = x(mm)-x(mm-2)
      c1 = deln/delnn/delnm1   
      c2 = -delnn/deln/delnm1  
      c3 = (delnn+deln)/delnn/deln 
      do 110 j=1,nn
         npmpj = npm+j 
         temp(npmpj) = c1*z(mm-2,j)+c2*z(mm1,j)+c3*z(mm,j) 
  110 continue 
      zp(mm,1,3) = c1*zp(mm-2,1,1)+c2*zp(mm1,1,1)+c3*zp(mm,1,1)
      zxymns = c1*temp(npm-2)+c2*temp(npm-1)+c3*temp(npm)  
      go to 113
  111 c1 = 1./(x(2)-x(1))  
      do 112 j=1,nn
         zp(1,j,2) = c1*(z(2,j)-z(1,j))
         npmpj = npm+j 
         temp(npmpj) = zp(1,j,2)   
  112 continue 
      zp(1,1,3) = c1*(zp(2,1,1)-zp(1,1,1)) 
      zp(mm,1,3) = zp(1,1,3)   
      zxy1ns = c1*(temp(nn+2)-temp(np1))   
      zxymns = zxy1ns  
c  
c denormalize tension factor in y direction
c  
  113 sigmay = abs(sigma)*float(nn-1)/(y(nn)-y(1)) 
      del1 = 1./(y(2)-y(1))
      do 114 i=1,mm
         zp(i,2,1) = del1*(z(i,2)-z(i,1))  
  114 continue 
      zp(1,2,3) = del1*(zp(1,2,2)-zp(1,1,2))   
      zp(mm,2,3) = del1*(temp(npm+2)-temp(npm+1))  
c  
c set up right hand sides and tridiagonal system for y grid
c perform forward elimination  
c  
      dels = sigmay/del1   
      exps = fxp(dels) 
      sinhs = .5*(exps-1./exps)
      sinhin = del1/sinhs  
      diag1 = sinhin*(dels*.5*(exps+1./exps)-sinhs)
      diagin = 1./diag1
      do 115 i=1,mm
         zp(i,1,1) = diagin*(zp(i,2,1)-zp(i,1,1))  
  115 continue 
      zp(1,1,3) = diagin*(zp(1,2,3)-zp(1,1,3)) 
      zp(mm,1,3) = diagin*(zp(mm,2,3)-zp(mm,1,3))  
      spdiag = sinhin*(sinhs-dels) 
      temp(1) = diagin*spdiag  
      if (nn .eq. 2) go to 119 
      do 118 j=2,nm1   
         jm1 = j-1 
         jp1 = j+1 
         del2 = 1./(y(jp1)-y(j))   
         do 116 i=1,mm 
            zp(i,jp1,1) = del2*(z(i,jp1)-z(i,j))   
  116    continue  
         zp(1,jp1,3) = del2*(zp(1,jp1,2)-zp(1,j,2))
         npmpj = npm+j 
         zp(mm,jp1,3) = del2*(temp(npmpj+1)-temp(npmpj))   
         dels = sigmay/del2
         exps = fxp(dels)  
         sinhs = .5*(exps-1./exps) 
         sinhin = del2/sinhs   
         diag2 = sinhin*(dels*(.5*(exps+1./exps))-sinhs)   
         diagin = 1./(diag1+diag2-spdiag*temp(jm1))
         do 117 i=1,mm 
            zp(i,j,1) = diagin*
     1                  (zp(i,jp1,1)-zp(i,j,1)-spdiag*zp(i,jm1,1)) 
  117    continue  
         zp(1,j,3) = diagin*(zp(1,jp1,3)-zp(1,j,3)-spdiag*zp(1,jm1,3)) 
         zp(mm,j,3) = diagin*  
     1                (zp(mm,jp1,3)-zp(mm,j,3)-spdiag*zp(mm,jm1,3))
         spdiag = sinhin*(sinhs-dels)  
         temp(j) = diagin*spdiag   
         diag1 = diag2 
  118 continue 
  119 diagin = 1./(diag1-spdiag*temp(nm1)) 
      do 120 i=1,mm
         npi = nn+i
         zp(i,nn,1) = diagin*(temp(npi)-zp(i,nn,1)-spdiag*zp(i,nm1,1)) 
  120 continue 
      zp(1,nn,3) = diagin*(zxy1ns-zp(1,nn,3)-spdiag*zp(1,nm1,3))   
      temp(nn) = diagin*(zxymns-zp(mm,nn,3)-spdiag*zp(mm,nm1,3))   
c  
c perform back substitution
c  
      do 122 j=2,nn
         jbak = np1-j  
         jbakp1 = jbak+1   
         t = temp(jbak)
         do 121 i=1,mm 
            zp(i,jbak,1) = zp(i,jbak,1)-t*zp(i,jbakp1,1)   
  121    continue  
         zp(1,jbak,3) = zp(1,jbak,3)-t*zp(1,jbakp1,3)  
         temp(jbak) = zp(mm,jbak,3)-t*temp(jbakp1) 
  122 continue 
      del1 = 1./(x(2)-x(1))
      do 123 j=1,nn
         zp(2,j,2) = del1*(z(2,j)-z(1,j))  
         zp(2,j,3) = del1*(zp(2,j,1)-zp(1,j,1))
  123 continue 
c  
c denormalize tension factor in x direction
c  
      sigmax = abs(sigma)*float(mm-1)/(x(mm)-x(1)) 
c  
c set up right hand sides and tridiagonal system for x grid
c perform forward elimination  
c  
      dels = sigmax/del1   
      exps = fxp(dels) 
      sinhs = .5*(exps-1./exps)
      sinhin = del1/sinhs  
      diag1 = sinhin*(dels*.5*(exps+1./exps)-sinhs)
      diagin = 1./diag1
      do 124 j=1,nn
         zp(1,j,2) = diagin*(zp(2,j,2)-zp(1,j,2))  
         zp(1,j,3) = diagin*(zp(2,j,3)-zp(1,j,3))  
  124 continue 
      spdiag = sinhin*(sinhs-dels) 
      temp(nn+1) = diagin*spdiag   
      if (mm .eq. 2) go to 128 
      do 127 i=2,mm1   
         im1 = i-1 
         ip1 = i+1 
         del2 = 1./(x(ip1)-x(i))   
         do 125 j=1,nn 
            zp(ip1,j,2) = del2*(z(ip1,j)-z(i,j))   
            zp(ip1,j,3) = del2*(zp(ip1,j,1)-zp(i,j,1)) 
  125    continue  
         dels = sigmax/del2
         exps = fxp(dels)  
         sinhs = .5*(exps-1./exps) 
         sinhin = del2/sinhs   
         diag2 = sinhin*(dels*(.5*(exps+1./exps))-sinhs)   
         npi = nn+i
         diagin = 1./(diag1+diag2-spdiag*temp(npi-1))  
         do 126 j=1,nn 
            zp(i,j,2) = diagin*
     1                  (zp(ip1,j,2)-zp(i,j,2)-spdiag*zp(im1,j,2)) 
            zp(i,j,3) = diagin*
     1                  (zp(ip1,j,3)-zp(i,j,3)-spdiag*zp(im1,j,3)) 
  126    continue  
         spdiag = sinhin*(sinhs-dels)  
         temp(npi) = diagin*spdiag 
         diag1 = diag2 
  127 continue 
  128 diagin = 1./(diag1-spdiag*temp(npm-1))   
      do 129 j=1,nn
         npmpj = npm+j 
         zp(mm,j,2) = diagin*(temp(npmpj)-zp(mm,j,2)-spdiag*zp(mm1,j,2))   
         zp(mm,j,3) = diagin*(temp(j)-zp(mm,j,3)-spdiag*zp(mm1,j,3))   
  129 continue 
c  
c perform back substitution
c  
      do 131 i=2,mm
         ibak = mp1-i  
         ibakp1 = ibak+1   
         npibak = nn+ibak  
         t = temp(npibak)  
         do 130 j=1,nn 
            zp(ibak,j,2) = zp(ibak,j,2)-t*zp(ibakp1,j,2)   
            zp(ibak,j,3) = zp(ibak,j,3)-t*zp(ibakp1,j,3)   
  130    continue  
  131 continue 
      return   
      end  
      function surf2 (xx,yy,m,n,x,y,z,iz,zp,sigma,noder,dzx,dzy)   
c  
c  
c dimension of           x(m),y(n),z(iz,n),zp(3*m*n)   
c arguments
c  
c latest revision        october 22, 1973  
c  
c purpose                surf2 interpolates a surface at a given   
c                        coordinate pair using a bi-spline under   
c                        tension.  the subroutine surf1 should be called   
c                        earlier to determine certain necessary
c                        parameters.   
c  
c access cards           *fortran,s=ulib,n=surf
c                        *cosy 
c  
c usage                  zxy = surf2 (xx,yy,m,n,x,y,z,iz,zp,sigma) 
c  
c arguments
c  
c on input               xx and yy 
c                          contain the x- and y-coordinates of the point   
c                          to be mapped onto the interpolating surface.
c  
c                        m and n   
c                          contain the number of grid lines in the x-  
c                          and y-directions respectively, of the   
c                          rectangular grid which determined the   
c                          surface.
c  
c                        x and y   
c                          are arrays containing x- and y-grid values  
c                          respectively, each in increasing order. 
c  
c                        z 
c                          is a matrix containing the m*n functional   
c                          values corresponding to the grid values 
c                          (i.e., z(i,j) is the surface value at the   
c                          point (x(i),y(j)) i = 1,...,m j = 1,...,n). 
c  
c                        iz
c                          contains the fortran row dimension of the   
c                          array z declared in the calling program.
c  
c                        zp
c                          is an array of 3*m*n locations stored with  
c                          the various surface derivative values   
c                          determined by surf1.
c  
c                        sigma 
c                          contains the tension factor (its sign is
c                          ignored).   
c  
c                        noder 
c                          if noder is .true. derivatives of z with
c                          respect to x and y are set in dzx and dzy   
c  
c                        the parameters m, n, x, y, z, iz, zp and sigma
c                        should be input unaltered from the output of  
c                        surf1.
c  
c on output              surf2 
c                          contains the interpolated surface value.
c  
c                        none of the input parameters are altered. 
c  
c entry points           surf2 
c  
c special conditions     none  
c  
c common blocks          none  
c  
c i/o                    none  
c  
c precision              single
c  
c required ulib          none  
c routines 
c  
c specialist             russell k. rew, ncar, boulder, colorado  80302
c  
c language               fortran   
c  
c history                originally written by a. k. cline,
c                        september 1973.   
c  
c  
c  
c  
      logical       noder  
      integer         m          ,n          ,iz   
      real            xx         ,yy         ,x(m)       ,y(n)       , 
     1                z(iz,n)    ,zp(m,n,3)  ,sigma
      data i1/2/,j1/2/ 
c  
c inline one dimensional spline under tension interpolation
c  
      herm(f1,f2,fp1,fp2) = (fp2*sinhd1+fp1*sinhd2)/sinhs+ 
     1                      ((f2-fp2)*del1+(f1-fp1)*del2)/dels 
c  
      fxp(t)=exp(amax1(-87.,amin1(87.,t))) 
c  
c denormalize tension factor in x and y directions 
c  
      sigmax = abs(sigma)*float(m-1)/(x(m)-x(1))   
      sigmay = abs(sigma)*float(n-1)/(y(n)-y(1))   
c  
c determine x interval 
c  
  101 do 102 i=i1,m
         if (x(i)-xx) 102,102,103  
  102 continue 
      i = m
  103 if (xx.ge.x(i-1) .or. xx.lt.x(1)) go to 104  
      i1 = 2   
      go to 101
c  
c determine y interval 
c  
  104 do 105 j=j1,n
         if (y(j)-yy) 105,105,106  
  105 continue 
      j = n
  106 if (yy.ge.y(j-1) .or. yy.lt.y(1)) go to 107  
      j1 = 2   
      go to 104
  107 del1 = yy-y(j-1) 
      del2 = y(j)-yy   
      dels = y(j)-y(j-1)   
      if(dels.lt.87./sigmay) go to 110 
      ex1=fxp(-sigmay*del2)
      ext1=fxp(sigmay*(yy+y(j-1)-2*y(j)))  
      ex2=fxp(-sigmay*del1)
      ext2=fxp(sigmay*(2*y(j-1)-yy-y(j)))  
      sinhs=1  
      go to 115
  110 ex1 = exp(sigmay*del1)   
      ext1=1/ex1   
      ex2 = exp(sigmay*del2)   
      ext2=1/ex2   
      exps = ex1*ex2   
      sinhs = exps-1./exps 
c  
  115 sinhd1=ex1-ext1  
      sinhd2=ex2-ext2  
c  
c interpolate in y direction   
c  
      zim1 = herm(z(i-1,j-1),z(i-1,j),zp(i-1,j-1,1),zp(i-1,j,1))   
      zi = herm(z(i,j-1),z(i,j),zp(i,j-1,1),zp(i,j,1)) 
      zxxim1 = herm(zp(i-1,j-1,2),zp(i-1,j,2),zp(i-1,j-1,3),zp(i-1,j,3))   
      zxxi = herm(zp(i,j-1,2),zp(i,j,2),zp(i,j-1,3),zp(i,j,3)) 
c  
      if(noder) go to 119  
      sinhd1=sigmay*(ex1+ext1) 
      sinhd2=-sigmay*(ex2+ext2)
      del1=1   
      del2=-1  
c  coefficients for y-derivative   
      dzim1 = herm(z(i-1,j-1),z(i-1,j),zp(i-1,j-1,1),zp(i-1,j,1))  
      dzi = herm(z(i,j-1),z(i,j),zp(i,j-1,1),zp(i,j,1))
      dzxim1 = herm(zp(i-1,j-1,2),zp(i-1,j,2),zp(i-1,j-1,3),zp(i-1,j,3))   
      dzxi = herm(zp(i,j-1,2),zp(i,j,2),zp(i,j-1,3),zp(i,j,3)) 
c  
  119 del1 = xx-x(i-1) 
      del2 = x(i)-xx   
      dels = x(i)-x(i-1)   
      if(dels.lt.87./sigmax) go to 120 
      ex1=fxp(-sigmax*del2)
      ext1=fxp(sigmax*(xx+x(i-1)-2*x(i)))  
      ex2=fxp(-sigmax*del1)
      ext2=fxp(sigmax*(2*x(i-1)-xx-x(i)))  
      sinhs=1  
      go to 125
  120 ex1 = exp(sigmax*del1)   
      ext1=1/ex1   
      ex2 = exp(sigmax*del2)   
      ext2=1/ex2   
      exps = ex1*ex2   
      sinhs = exps-1./exps 
c  
  125 sinhd1=ex1-ext1  
      sinhd2=ex2-ext2  
c  
c final interpolation in x direction   
c  
      surf2 = herm(zim1,zi,zxxim1,zxxi)
c  
      if(noder) go to 130  
c  y-derivative
      dzy=herm(dzim1,dzi,dzxim1,dzxi)  
c  
      sinhd1=sigmax*(ex1+ext1) 
      sinhd2=-sigmax*(ex2+ext2)
      del1=1   
      del2=-1  
c  x-derivative
      dzx=herm(zim1,zi,zxxim1,zxxi)
  130 i1 = i   
      j1 = j   
      return   
      end  
