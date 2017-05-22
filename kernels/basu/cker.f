	program nonasy
	implicit real*8(a-h,o-z)
	parameter(npt=8500,nmd=10500,ndim=200,pi=3.1415926535d0)
	parameter (npt2=4500)
	real*4  rd(npt), yy(4,npt), oma, enera,a1,a2,omm,rll
	real*4 akc(npt2,nmd), akr(npt2,nmd),stds,ooo, rrp(npt2)
	character*100 filam,fileig, filck, fildk
	parameter(GG=6.67232d-8)
	dimension fb(10)
	dimension csq(npt),rad(npt),rho(npt),g(npt),hrho(npt)
	dimension xir(npt),p1(npt),psi(npt)
	dimensionxih(npt),rho1(npt),divxi(npt),dxir(npt),f1(npt),f2(npt)
	dimension n(nmd),rl(nmd),om(nmd),energy(nmd),omf(nmd)
c	new arrays defined 
	dimension aan(5,npt)
	dimension datan(8)
	dimension  amass(npt)
	dimension rin1(npt),rin2(npt), rin3(npt)
	dimension  ddr(npt),ddx(npt), stds(nmd)


51	format(' spline coef:'/(1p5e12.3))
52	format(1p9e13.5)
53	format(1p30e13.5)
54	format(i4,f8.1,1p4e13.5)
55	format(1x, 2i3,1x,1p9e14.6)
59	format(1x, 5e14.6,1x,4i3)

	read(*,*)filam
	read(*,*)fileig
	read(*,*)filck
	read(*,*)fildk
c	###################################
        open(unit=12,file=filam ,form='unformatted')
	open(unit=14,file=fileig,   form='unformatted')

	open(21,file=filck ,  form='unformatted')
	open(22,file=fildk,  form='unformatted')


c       Constants to normalise various quantities to dimensionless units
	anorml=4.*pi/3.

c	READING IN REFERNCE MODEL QUANTITIES
	read(14)
	read(14)a1,a2
	rs=a1
	rms=a2
	print*,rs,rms
	write(21)a1,a2
	write(22)a1,a2


	np=0
	do 20 i=1,20000
	read(12,end=21)rad(i),(aan(j,i),j=1,5)
	np=np+1
20	continue
21	print*, 'np ', np, rad(1),rad(np)
	np=np-1
	print*,' printed np1'

	csu=3.*GG*rms/(4.*pi*rs)
	rhou=3.*rms/(4.*pi*rs**3)
	gu=3.*GG*rms/(4.*pi*rs**2)
	dgru=3.*GG*rms/(4.*pi*rs**3)
	omu=sqrt(dgru)
	anorml=4*pi/3.
	omu1=sqrt(anorml)
	pu=GG*(3.*rms/(4.*pi*rs**2))**2
	omu=omu*omu1

	print *,csu,rhou,gu,dgru,omu,pu,enu
	print*,np,rad(1),rad(np)

	do 800 i=1,np
	if(rad(i).lt.1.d-80)rad(i)=1.d-80
	amass(i)=aan(1,i)*rad(i)**3
	csq(i)=aan(1,i)*rad(i)**2/aan(2,i)
	hrho(i)=rad(i)/(aan(4,i)+aan(2,i))
	rho(i)=aan(1,i)*aan(5,i)
	g(i)=aan(1,i)*rad(i)
cc	if(hrho(i).gt.1.d10)then
cc	hrho(i)=hrho(i-1)
cc	endif
cc	write(1,59)rad(i),amass(i),csq(i), rho(i), g(i), hrho(i),aan(2,i)
800	continue	
	print*,np, 'rad(1)=',rad(1),' rad(np)=',rad(np)

	print*,'CSU ', csu, anorml, csq(1),rad(1)
 	
c  ############################
	np0=1


	
	nmode=0
	do 4000 i=1,200000
	read(14,end=4100)ni,rll, oma, omm, enera, enera
	print*,ni,rll,oma*(5.d5*omu)/pi,omm
	write(77,*)i,ni,rll
	read(14) ((yy(j,jj),j=1,4),jj=1,np)
c	print*,yy(1,1),yy(2,1),yy(3,1),yy(4,1)
cc	if(rll.eq.3)go to 4100
	nmode=nmode+1
	rl(nmode)=rll
	rhl=rll*(rll+1.)

	         do 4001 j=1,np
	xir(j)=yy(1,j)
	p1(j)=yy(2,j)
	psi(j)=yy(3,j)
cc	dpsi(j)=yy(3,j)
	rr=rad(j)
	rin1(j)=0.
	rin2(j)=0.
	rin3(j)=0.
4001	continue

	li=rll+0.1


	om(nmode)=oma
	omf(nmode)=omm
	n(nmode)=ni
	l=rll+0.1

	do 3300 j=np0,np
	rr=rad(j)
cc	print*, '3300 ', j, rr
c################################
cc		 if(rli.ne.0)then
cc                psi(j)=yy(3,j)*g(j)
cc                psij=psi(j)
cc                p1(j)=rho(j)*(rr*omi**2*xih(j)+psij)
cc                else
cc                p1(j)=yy(2,j)*rho(j)*omi**2
cc                endif
	psij=psi(j)
	xih(j)=(p1(j)/rho(j)-psi(j))/(rr*om(nmode)**2)
	rho1(j)=(p1(j)-xir(j)*g(j)*rho(j))/csq(j)+xir(j)*rho(j)/hrho(j)
	divxi(j)=-rho1(j)/rho(j)+xir(j)/hrho(j)
	dxir(j)=divxi(j)-2.*xir(j)/rr+rhl*xih(j)/rr

	f1(j)=rr**(l+2)*rho(j)*(divxi(j)-(l+2)*xir(j)/rr-dxir(j))

	if(l*log(rr).lt.-200.)then
	f2(j)=0.
	else
	f2(j)=rr**(1-l)*rho(j)*(divxi(j)-(1-l)*xir(j)/rr-dxir(j))
	endif

	ddr(j)=rr**2*csq(j)*rho(j)*xir(j)*divxi(j)
	ddx(j)= xir(j)*(rho1(j)-rho(j)*divxi(j))
3300	continue
c	print*, 'out of 3300'
cc	enorm=rho(81)*(xir(81)**2+rhl*xih(81)**2)*rad(81)**2
cc	print*,'ENORM= ',enorm
	erri=1.0


c	---kernels only---
	s2k=0.0
	s1k=0.0
	s3k=0.0
c	---kernels only---
        if(l*log(rad(np)).lt.-200)then
	rin2(np)=0.0
	else
	rin2(np)=s2k/rad(np)**(l+1)
	endif
	rin2(np)=0.0
	rin1(np0)=0.0
	rin3(np0)=0.0

c	print*,' before 3320'

	do 3320 j=np0,np-1
c		---
c	print*,'3320 ', j, l

        if(l*log(rad(j+1)).lt.-200.)then
        rin1(j+1)=0.
	else
	s1k=s1k+0.5*(rad(j)-rad(j+1))*(rho1(j)*rad(j)**(1-l)+
     1        rho1(j+1)*rad(j+1)**(1-l))
        sgn=s1k/abs(s1k)
	rin1(j+1)=s1k
	
        endif

c	if(s1k.gt.1.d150.or.s1k.lt.-1.d150)s1k=sgn*1.d150
cc	if(rin1(j+1).gt.1.d10)print*,'rin1 ',rin1(j+1)
c	print*,'rin1 ', rin1(j+1), s1k

	s3k=s3k+0.5*(rad(j)-rad(j+1))*(ddx(j)+ddx(j+1))
	rin3(j+1)=rho(j+1)*s3k
c	print*, ' 3320 end ', j, s3k, rin3(j+1)
c		---
3320	continue
c	print*, 'after 3320'

	do 3330 j=np,np0+1,-1
c		---

	s2k=s2k+0.5*(rad(j-1)-rad(j))*(rho1(j-1)*rad(j-1)**(l+2)+
     1     rho1(j)*rad(j)**(l+2))
	rin2(j-1)=s2k
3330	continue
c	print*,'rin ', rint1(np0+1,1),rint2(np0+1,1)
c	print*, 'after 3330'


c		---- 
	std=0.
c		---- 
	sint=0.0
	ss=0.5*(rad(np0)-rad(np0+1))

	np2=0.
	do 3600 j=np0,np
	rr=rad(j)
	ssr=ss*rr**2
	rrr=rr**2
	r1=-rrr*csq(j)*divxi(j)*rho(j)*divxi(j)
	r2=-ssr*rin3(j)
	nuse=4
	reps=0.
	x0=rad(j)
      call DIVDIF(X0,rad,ddr,NUSE,np,FB,REPS,IER1,DFB,DDFB)
	r3=-dfb*ss
	r5=ss*f1(j)*rin1(j)
	r52=ss*f2(j)*rin2(j)
	r6=ssr*rho(j)*(xir(j)**2+rhl*xih(j)**2)
c	rem=(j/2.)-(j/2)
c	rem=(j/4.)-(j/4)
	rem=0.

	if(rem.eq.0)then
	np2=np2+1
	rrp(np2)=rad(j)
	akc(np2,nmode)=-(r1)/(2.*om(nmode)**2)
	akr(np2,nmode)=-((r2+r3-(r5+r52)/(2.*rll+1.))/ss)
     1  /(2.*om(nmode)**2)
	endif

	lll=rl(nmode)+0.5
	omega=om(nmode)

	std=std+r6

c		---  

	ss=0.5*(rad(j)-rad(j+2))
	if(j.ge.np-1)ss=0.5*(rad(j)-rad(j-1))
3600	continue



	energy(nmode)=enera
	print*,enera
	stds(nmode)=std
	if(std.eq.0) then
	print*,'problem'
	print*, ni,l
	stop
	endif

c	endif
4000	continue

4100	continue
	write(21)np2,(rrp(ji),ji=1,np2)
	write(22)np2,(rrp(ji),ji=1,np2)
	print*,rrp(1),rrp(np2)
	do 1210 ii=1,nmode
	lll=rl(ii)+.1
	ooo=om(ii)*(5.d5*omu)/pi
cc	ooo=omf(ii)
cc	if(lll.le.2)then
cc	write(21)lll, n(ii)-1, ooo,((akc(i,ii)/stds(ii)),i=1,np2)
cc	write(22)lll, n(ii)-1, ooo,((akr(i,ii)/stds(ii)),i=1,np2)
cc	else
	write(21)lll, n(ii), ooo,((akc(i,ii)/stds(ii)),i=1,np2)
	write(22)lll, n(ii), ooo,((akr(i,ii)/stds(ii)),i=1,np2)
cc	endif
	write(42,*)lll,n(ii),ooo, rl(ii)
1210	continue
	close(19)
	print*, 'nmode= ', nmode, ' np2 ',np2


c        STOP HERE IF ONLY FINDING KERNELS
c		STOP

8000	continue
	end

C	-------------------------------------------------


      SUBROUTINE DIVDIF(XB,X,F,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	implicit real*8(a-h,o-z)
      PARAMETER(NMAX=10)
      DIMENSION X(NTAB),F(NTAB),FB(*),XN(NMAX),XD(NMAX)

      NEXT=NEARST(XB,X,NTAB)
      FB(1)=F(NEXT)
      XD(1)=F(NEXT)
      XN(1)=X(NEXT)
      IER=0
      PX=1.0

      DFB=0.0
      DDFB=0.0
      DPX=0.0
      DDPX=0.0

      IP=NEXT
      IN=NEXT

      NIT=MIN0(NMAX,NUSE,NTAB)
      IF(NUSE.GT.NMAX.OR.NUSE.GT.NTAB) IER=12
      IF(NUSE.LT.1) THEN
        IER=11
        NIT=MIN0(6,NTAB,NMAX)
      ENDIF
      NUSE=1

      DO 5000 J=2,NIT

        IF(IN.LE.1) GO TO 2200
        IF(IP.GE.NTAB) GO TO 2000
        IF(ABS(XB-X(IP+1)).LT.ABS(XB-X(IN-1))) GO TO 2200
2000    IN=IN-1
        NEXT=IN
        GO TO 2800
2200    IP=IP+1
        NEXT=IP

2800    XD(J)=F(NEXT)
        XN(J)=X(NEXT)
        DO 3000 K=J-1,1,-1
3000    XD(K)=(XD(K+1)-XD(K))/(XN(J)-XN(K))

        DDPX=DDPX*(XB-XN(J-1))+2.*DPX
        DPX=DPX*(XB-XN(J-1))+PX
        DFB=DFB+DPX*XD(1)
        DDFB=DDFB+DDPX*XD(1)

        PX=PX*(XB-XN(J-1))
        ERR=XD(1)*PX
        FB(J)=FB(J-1)+ERR
        NUSE=J

        IF(ABS(ERR).LT.REPS) RETURN
5000  CONTINUE

      IER=24
      END

c	-------------------------------------------

      FUNCTION NEARST(XB,X,NTAB)
	implicit real*8(a-h,o-z)
      DIMENSION X(NTAB)

      LOW=1
      IGH=NTAB
      IF(.NOT.(XB.LT.X(LOW).EQV.XB.LT.X(IGH))) THEN


1500    IF(IGH-LOW.GT.1) THEN
          MID=(LOW+IGH)/2
          IF(XB.LT.X(MID).EQV.XB.LT.X(LOW)) THEN
            LOW=MID
          ELSE
            IGH=MID
          ENDIF
          GO TO 1500
        ENDIF
      ENDIF

      IF(ABS(XB-X(LOW)).LT.ABS(XB-X(IGH))) THEN
        NEARST=LOW
      ELSE
        NEARST=IGH
      ENDIF
      END

c	--------------------------------------------------

      real*8 FUNCTION GASDEV(IDUM)
	implicit real*8(a-h,o-z)
	save
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
1       V1=2.*RAN1(IDUM)-1.
        V2=2.*RAN1(IDUM)-1.
        R=V1**2+V2**2
        IF(R.GE.1.)GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END

c	---------------------------


      real*8 FUNCTION RAN1(IDUM)
	implicit real*8(a-h,o-z)
	save
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END
c	--------------------

      SUBROUTINE SPLINE(X,F,N,C,IER)
	implicit real*8(a-h,o-z)
      DIMENSION X(N),F(N),C(3,N)

      IF(N.LE.1) THEN
        IER=111
        RETURN

      ELSE IF(N.EQ.2) THEN
        C(1,1)=(F(2)-F(1))/(X(2)-X(1))
        C(2,1)=0.0
        C(3,1)=0.0
        RETURN

      ELSE IF(N.EQ.3) THEN
        DIV12=(F(2)-F(1))/(X(2)-X(1))
        DIV23=(F(3)-F(2))/(X(3)-X(2))
        C(3,1)=0.0
        C(3,2)=0.0
        C(2,1)=(DIV23-DIV12)/(X(3)-X(1))
        C(2,2)=C(2,1)
        C(1,1)=DIV12+C(2,1)*(X(1)-X(2))
        C(1,2)=DIV23+C(2,1)*(X(2)-X(3))
        RETURN

      ELSE

        C(3,N)=(F(N)-F(N-1))/(X(N)-X(N-1))
        DO 1000 I=N-1,2,-1
          C(3,I)=(F(I)-F(I-1))/(X(I)-X(I-1))
          C(2,I)=2.*(X(I+1)-X(I-1))
1000    C(1,I)=3.*(C(3,I)*(X(I+1)-X(I))+C(3,I+1)*(X(I)-X(I-1)))

        C1=X(3)-X(1)
        C(2,1)=X(3)-X(2)
        C(1,1)=C(3,2)*C(2,1)*(2.*C1+X(2)-X(1))+C(3,3)*(X(2)-X(1))**2
        C(1,1)=C(1,1)/C1
        CN=X(N)-X(N-2)
        C(2,N)=X(N-1)-X(N-2)
        C(1,N)=C(3,N)*C(2,N)*(2.*CN+X(N)-X(N-1))
        C(1,N)=(C(1,N)+C(3,N-1)*(X(N)-X(N-1))**2)/CN

        G=(X(3)-X(2))/C(2,1)
        C(2,2)=C(2,2)-G*C1
        C(1,2)=C(1,2)-G*C(1,1)
        DO 2000 J=2,N-2
          G=(X(J+2)-X(J+1))/C(2,J)
          C(2,J+1)=C(2,J+1)-G*(X(J)-X(J-1))
2000    C(1,J+1)=C(1,J+1)-G*C(1,J)
        G=CN/C(2,N-1)
        C(2,N)=C(2,N)-G*(X(N-1)-X(N-2))
        C(1,N)=C(1,N)-G*C(1,N-1)

        C(1,N)=C(1,N)/C(2,N)
        DO 3000 I=N-1,2,-1
3000    C(1,I)=(C(1,I)-C(1,I+1)*(X(I)-X(I-1)))/C(2,I)
        C(1,1)=(C(1,1)-C(1,2)*C1)/C(2,1)

        DO 4000 I=1,N-1
          C(2,I)=(3.*C(3,I+1)-2.*C(1,I)-C(1,I+1))/(X(I+1)-X(I))
          C(3,I)=(C(1,I)+C(1,I+1)-2.*C(3,I+1))/(X(I+1)-X(I))**2
4000    CONTINUE
        C(2,N)=0.0
        C(3,N)=0.0
      ENDIF
      END

C	-------------------------------------------

      FUNCTION SPLEVL(XB,N,X,F,C,DFB,DDFB,IER)
	implicit real*8(a-h,o,p,r-z)
      IMPLICIT LOGICAL(Q)
      DIMENSION X(N),F(N),C(3,N)
      SAVE
      DATA LOW/0/

      IF(N.LE.1) THEN
        IER=111
        RETURN
      ENDIF

      QASCND=X(N).GT.X(1)
      IER=0

      IF(LOW.LT.1.OR.LOW.GE.N) THEN
        LOW=1
        IGH=N
      ENDIF

1000  IF((XB.LT.X(LOW).AND.XB.LT.X(IGH)).OR.
     1    (XB.GT.X(LOW).AND.XB.GT.X(IGH))) THEN
        IF(XB.GT.X(LOW).EQV.QASCND) THEN
          IF(IGH.GE.N) THEN
            IER=11
            LOW=N-1
          ELSE
            NIGH=MIN0(N,IGH+2*(IGH-LOW))
            LOW=IGH
            IGH=NIGH
            GO TO 1000
          ENDIF
        ELSE
          IF(LOW.LE.1) THEN
            IER=12
          ELSE
            NIGH=LOW
            LOW=MAX0(1,LOW-2*(IGH-LOW))
            IGH=NIGH
            GO TO 1000
          ENDIF
        ENDIF
      ELSE

1500    IF(IGH-LOW.GT.1.AND.XB.NE.X(LOW)) THEN
          MID=(LOW+IGH)/2
          IF(XB.LE.X(MID).EQV.XB.LE.X(LOW)) THEN
            LOW=MID
          ELSE
            IGH=MID
          ENDIF
          GO TO 1500
        ENDIF
      ENDIF

      DX=XB-X(LOW)
      SPLEVL=((C(3,LOW)*DX+C(2,LOW))*DX+C(1,LOW))*DX+F(LOW)
      DFB=(3.*C(3,LOW)*DX+2.*C(2,LOW))*DX+C(1,LOW)
      DDFB=6.*C(3,LOW)*DX+2.*C(2,LOW)
      END
