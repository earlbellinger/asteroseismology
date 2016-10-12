	program sola
        implicit real*8(a-h,o-z)
        character*100 filenm,filek1,filek2,filmod, fileng,filobs
	character*100 ftest
        parameter(npt=5500,nmd=5500,ndim=300,pi=3.1415926535d0)
        parameter(nmmd=5500,nmd1=5500)
        parameter(npt2=5500,nx0=100)
        parameter(GG=6.67232d-8)
	real*4 arr(npt2), ak(npt2), rd(npt2), ooo,rs,rms

cc	dimension arr(npt2), ak(npt2), rd(npt2)
cc	dimension
cc	dimension
        dimension rad(npt), datan(8), aan(5,npt2),rdd(npt2)
        dimension  fdif(50,2000), rho(npt2), csq(npt2)
        dimension std(nmd), enorm(nmd), xc(nx0)
        dimension  aa1(nmd,nmd)
        dimension   fakc(npt,nmd), fakr(npt,nmd), rhl(npt2)
        dimension cov(nmd1,nmd1), w(nmd), error(nmd)
        dimension n(nmd),rl(nmd),om(nmd),energy(nmd),dom(nmd)
     1     ,err(50,2000),ener(50,2000)
        dimension erorg(nmd),rhol(npt2), cs(npt2), freq(50,2000)
	dimension weight(npt), ov1(nmd,nmd), ov2(nmd,nmd), gmg(npt2),gam(npt2)
	dimension surcon(nmd,nx0), vv1(nmd,nx0), target(npt,nx0), x0(nx0)
	dimension delta(nx0), sumrc(nx0)
	dimension ctil1(nmd,nx0),  cenc(nx0), cenr(nx0)
	dimension pp(0:nx0,nmd),ipiv(nmd), widthc(3,nx0)
	dimension sumcc(nx0), sumr(nx0), avc(npt,nx0)
	dimension ckint(nx0), rkint(nx0), avrc(npt,nx0), acros(nx0)
	dimension  fb(10),uuu(npt2),uu(npt2)
	dimension sol(nx0)
c
	dimension enl(50),ome(50)
c

cc        data ob,oe,pb,pe/1.0,5.25,5.00,17.3/
        data ob,oe,pb,pe/1.0,3.01,5.00,17.3/
	data num/76/

51      format(' spline coef:'/(1p5e12.3))
52      format(1p100e14.6)
53      format(1x,'theta = ',f13.9,3x,'beta = ',f11.8,3x,
     1  'Lambda = ',i3, ' deldel =', e13.5)
54      format(f7.5,6(1x,f6.4),1p9e13.5)
55      format(1x, 2i6,1x,1p9e14.6)
59      format(1x, 5e14.6,1x,4i3)

	read(*,*)ftest
        read(*,*)filmod
        read(*,*)fileng
        read(*,*)filek1
        read(*,*)filek2
        read(*,*)filobs
	read(*,*)ob,oe
cc	open(32,file=ftest, form='unformatted')
	open(32,file=ftest)
	open(unit=12,file=filmod,form='unformatted')
        open(unit=21, file=filek1,form='unformatted')
        open(unit=22, file=filek2,form='unformatted')
cc	open(unit=24,file=fileng,form='unformatted')
	open(unit=24,file=fileng)
        print*,'file= ',filenm,filek1,filek2,filmod,fileng
	open(11,file=filobs)

	open(35,file='../tmp/surcon',
     1 form='unformatted')

	
	open(19,file='../tmp/fort19')
	open(20,file='../tmp/fort20')

c	-------------
	open(25,file=
     1  '../tmp/ov.u-Y.small'
     1  ,form='unformatted')
	open(26,file=
     1  '../tmp/ov.Y-u.small'
     1  ,form='unformatted')
	open(27,file='../tmp/coef.u-Y',
     1  form='unformatted')
	open(28,file='../tmp/coef.Y-u',
     1  form='unformatted')
	open(37,file='../tmp/coef')

c	%np%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c	Constants to normalise various quantities to dimensionless units
        anorml=4.d0*pi/3.d0
 
c       READING IN REFERNCE MODEL QUANTITIES
	np0=1
	np=0
	do 800 i=1,npt2
	read(12,end=801)rdd(i),(aan(j,i),j=1,5)
	np=np+1
        term=aan(1,i)*aan(5,i)/3.d0
	rhl(i)=log10(term)
        cs(i)=aan(1,i)*rdd(i)**2d0*anorml/aan(2,i)
	uu(i)=cs(i)/aan(3,i)
	gmg(i)=aan(3,i)
	write(33,52)rdd(i),cs(i), aan(3,i)
800     continue
801	print*,np
        print*,'rad', rdd(np), np
	write(33,*)
	ntabb=np
 
	read(21)rs,rms
	read(22)
	print*,rs,rms
cc	rs=6.9578e10
cc	rms=1.98900E+33
        csu=3.d0*GG*rms/(4.d0*pi*rs)
        rhou=3.d0*rms/(4.d0*pi*rs**3)
        gu=3.d0*GG*rms/(4.d0*pi*rs**2)
        dgru=3.d0*GG*rms/(4.d0*pi*rs**3)
        omu=sqrt(dgru)
        omu1=sqrt(anorml)
        pu=GG*(3.d0*rms/(4.d0*pi*rs**2))**2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do 1400 i=1,50
        do 1400 j=1,2000
        fdif(i,j)=0.0
        freq(i,j)=0.0
	ener(i,j)=0.0
1400    err(i,j)=0.0

        do 9092 i=1,nmd1
                do 9091 j=1,nmd1
        cov(j,i)=0.
	ov1(j,i)=0.
	ov2(j,i)=0.
9091    continue
	do 9090 j=1,num
	vv1(i,j)=0.
9090	continue
9092	continue
 
c       READ IN OBSERVED errors
cc	read(11,*)
cc	read(11,*)
cc	read(11,*)
	llimit=0
        do 1600 i=1,9000
cc        read(11,*,end=1700) ni,l,fre,er
        read(11,*,end=1700) l,ni,fre,er
	if(l.gt.llimit)llimit=l
        err(ni+1,l+1)=er
1600    continue
1700	continue

	
 
c       READING IN Test FREQUENCY 
        do 900 i=1,20000
        read(32,*,end=910)l,ni,fre
cc        read(32,*,end=910)l,ni,fre
	freq(ni+1,l+1)=fre
900     continue
910     continue
	print*,'i frequwncies', i
        print*, 'MESSAGE: READ IN FREQUENCY DIFFERENCES'

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
c	READ ENERGY
	nsp=0
	do 911 i=1,20000
	read(24,*,end=919)li,ni, fre, ene
cc	print*, li, ni, fre, ene
	ener(ni+1,li+1)=ene
	fdif(ni+1,li+1)=fre

	if(li.eq.0)then
	nsp=nsp+1
	enl(nsp)=log(ene)
	ome(nsp)=fre*pi/(5.d5*omu)
	endif

911	continue
919	continue


c	READ IN THE KERNELS
        read(21)np,(rd(ji),ji=1,np)
	print*,'read radius', rd(np), np
        read(22)
	nmode=0

	omax=0.
	omin=20000
	do 8100 i=1,20000
	read(21,end=999)lll,nnn, ooo,(ak(ii),ii=1,np)	
	read(22)lll,nnn, ooo,(arr(ii),ii=1,np)	
cc	if(nnn.eq.0)go to 8110
cc	if(nnn.ne.0.and.lll.gt.4)go to 8110
	

	if(abs(fdif(nnn+1,lll+1)).lt.1e-13)then
	print*, 'no model frequency', lll, nnn
	go to 8110
	endif

	if(ooo/1000.lt.ob)go to 8110

cc	print *, lll, nnn, fdif(nnn+1,lll+1)
	dif=(freq(nnn+1,lll+1)-fdif(nnn+1,lll+1))/fdif(nnn+1,lll+1)
        erri=err(nnn+1,lll+1)
cc	dif=dif+8.618676657e-05  695.82 1
cc	dif=dif-0.0002372438757 
cc695.67 2

	

	if(lll.eq.int(rl(nmode)+.1).and.nnn.eq.n(nmode))then
	print *,'lll', lll,nnn,ooo
	print*,int(rl(nmode)+.1),n(nmode)
cc	    nmode=nmode-1
	go to 8110
	endif

	if(lll.lt.0)go to 8110



	if(abs(erri).le.1e-11)then
	print*,'no obs. ', lll,nnn,ooo
	go to 8110
	endif

	if(ener(nnn+1,lll+1).le.0)then
	print*, 'no energy', lll,nnn
	go to 8110
	endif

	ot=freq(nnn+1,lll+1)
	if(ot/1000. gt.oe.or.ot/1000.lt.ob)go to 8110
81	nmode=nmode+1
	if(ooo.gt.omax)omax=ooo
	if(ooo.lt.omin)omin=ooo
cc	print*, lll,ot, ooo

		inum=0
		do 8111 jj=1,np,1
		inum=inum+1
		fakc(inum,nmode)=ak(jj)
		fakr(inum,nmode)=arr(jj)
		rad(inum)=rd(jj)
cc		if(lll.eq.100.and.nnn.eq.3)then
cc		write(33,52)rd(jj),ak(jj),arr(jj)
cc		endif
8111		continue

	n(nmode)=nnn
	rl(nmode)=lll
	om(nmode)=fdif(nnn+1,lll+1)*pi/(5.d5*omu)
	energy(nmode)=ener(nnn+1,lll+1)
	error(nmode)=erri*pi/(5.d5*omu)
	erorg(nmode)=erri/ooo
	dom(nmode)=dif
cc	write(23,*)lll,nnn,ooo/1000.,energy(nmode)
	write(23,55)lll,nnn,ooo, dif, erorg(nmode)
	print*,lll
cc	write(35)lll,nnn, fdif(nnn+1,lll+1)/1000., dif, erri
	go to 8100

8110	continue
8100	continue
999	continue

	im=-1
	ami=-10.
cc	write(35)im, im, ami, ami, ami
	np=inum
	print*,'np ',np

	sume=0.d0
	do 8101 i=1,nmode
	error(i)=error(i)/om(i)
cc	error(i)=1.d0
	sume=sume+1./(error(i)**2)
8101	continue
	
	print*,sume,sume/nmode,' ave. error'

	print*,'read kernels'
 
	print*,'omax,omin', omax, omin
	omax=omax*pi/(5.d5*omu)
	omin=omin*pi/(5.d5*omu)
	print*,'omax,omin', omax, omin

        print*,'read kernels',nmode
c %%%%%%%%%%%%%%%%%%%%%%
c	INTERPOLATING
        reps=1.e-7
        ntab=ntabb
	do 993 i=1,np
	xx=rad(i)
        nuse=4
      call DIVDIF(xx,rdd,cs,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	csq(i)=fb(nuse)
	nuse=4
      call DIVDIF(xx,rdd,rhl,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	rho(i)=10.d0**fb(nuse)
	rhol(i)=fb(nuse)
	nuse=4
      call DIVDIF(xx,rdd,uu,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	uuu(i)=fb(nuse)
	nuse=4
      call DIVDIF(xx,rdd,gmg,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	gam(i)=fb(nuse)
cc	write(33,52)xx,rh(i), csq(i) 
cc	write(33,52)xx,rh(i), csq(i) 
993	continue



	read*, deldel
        read *, idum, lam, icase, invcas,irep, ifil
c	read covariance matrix


	sumc=0.
	do 840 j=1,nmode
	cov(j,j)=error(j)*error(j)
	sumc=sumc+cov(j,j)
840	continue

	if(ifil.eq.1)then

	read(3)n1,n1,((cov(i,j),i=1,n1),j=1,n1)
		if(n1.ne.nmode)then
	print*,'there is something inconsistent in the cov. matrix'
		stop
		endif
	sumc=0.
	do 991 i=1,nmode
	sumc=sumc+cov(i,i)
991	continue
	print*,'read cov filterd matric'
	endif

	if(ifil.eq.2)then
	sumc=0.
	read(4,*)
	do 841 j=1,nmode
	read(4,*)l1,n1,er
	cov(j,j)=er*er
	sumc=sumc+cov(j,j)
841	continue
	endif

	sumc=sumc/nmode
	print*, 'sumc', sumc

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	print*,'type idum'
c	adding errors to data
        if(idum.ne.0)then
        do 404 i=1,nmode
        dom(i)=dom(i)+erorg(i)*gasdev(idum)
404     continue
        endif

        do 9010 i=1,nmode
 
 
 
	xxx=dom(i)
        w(i)=xxx

9010    continue

c	conservation of mass condition in case of rho inversion
	if(icase.eq.0)then
	do 312 j=1,np
	fakr(j,nmode+1)=0.
	fakc(j,nmode+1)=rho(j)*rad(j)**2
312	continue
	nmode=nmode+1
	error(nmode)=0.d0
	w(nmode)=0
	n(nmode)=-1
	om(nmode)=-1*pi/(5.d5*omu)
	rl(nmode)=-1d0
	endif


c %%%%%%%%%%%%%%%%%%%%
	print*,'calc. integration weights'
c	get integration weights

	weight(1)=0.5d0*(rad(1)-rad(2))	
	do 100 j=1,np-2
	weight(j+1)=0.5d0*(rad(j)-rad(j+2))
100	continue
	weight(np)=0.5d0*(rad(np-1)-rad(np))


c%%%%%%%%%%%%%%%%%%%%%%
c NOW THE SOLA RELATED STUFF 
	print*,'calc overlap matrix'

	if(invcas.eq.0)then
c	set the overlap matrix for both variables
c	----------------------------------------

	do 200 i=1,nmode
	   do 200 j=i,nmode

		sum1=0
		sum2=0
		
		do 202 k=1,np
		sum1=sum1+weight(k)*fakc(k,i)*fakc(k,j)
	sum2=sum2+((rad(k)+1.)**4)*weight(k)*fakr(k,i)*fakr(k,j)
202	continue
	if(sum1.eq.0)then
	print*,'000000 outout in overlap'
	print*,sum1, sum2, i, j
	pause
	endif


	ov1(i,j)=sum1
	ov1(j,i)=sum1

	ov2(i,j)=sum2
	ov2(j,i)=sum2
200	continue

	do 203 i=1,nmode
		sum1=0.
		sum2=0.
		do 204 k=1,np
		sum1=sum1+weight(k)*fakc(k,i)
		sum2=sum2+weight(k)*fakr(k,i)
204	continue
	ov1(i,nmode+1)=0.5*sum1
	ov2(i,nmode+1)=0.5*sum2
	ov1(nmode+1,i)=ov1(i,nmode+1)
	ov2(nmode+1,i)=ov2(i,nmode+1)
203	continue

	numker=nmode+1

	do 351 i=1,numker
	write(25)(ov1(i,j),j=1,numker)
	write(26)(ov2(i,j),j=1,numker)
351	continue

	ELSE

	numker=nmode+1
	do 352 i=1,numker
	read(25)(ov1(i,j),j=1,numker)
	read(26)(ov2(i,j),j=1,numker)
352	continue
	print*,'read', ov1(1,1), ov2(3,3)
	
	ENDIF

c	addtional constraints (surface)
c	------------------------------

	if(lam.ge.0)then
	print*,'calc. constraints'
	call constr(nmd, nmode,lam,om, energy, surcon,omax,omin,pp, 
     1  ome, enl, nsp, omu)
cc	call constr1(nmd, nmode,lam,om, energy, surcon,omax,omin,pp, 
cc     1  ome, enl, nsp, omu)
	do 350 j=1,lam+1
	do 300 i=1,nmode
	ov1(i,j+nmode+1)=surcon(i,j)
	ov1(j+nmode+1,i)=ov1(i,j+nmode+1)
	ov2(i,j+nmode+1)=surcon(i,j)
	ov2(j+nmode+1,i)=ov2(i,j+nmode+1)
300	continue
350	continue
	numker=numker+lam+1


	write(35)nmode,lam+1
	do 301 ij=1,nmode
	write(35)(surcon(ij,j),j=1,lam+1)
301	continue
	endif
	print*,'numker ',numker



c	set target kernels
c	------------------
	print*,'calc. targets'
	ntab=np
	rsamp=0.2d0
	iii=nearst(rsamp,rad,ntab)
	do 360 i=1,num
	x0(i)=0.05+(i-1)*.0126666666
cc	x0(i)=0.045+(i-1)*0.005
        ii= NEARST(x0(i),rad,NTAB)
	delta(i)=deldel*sqrt(csq(ii))/sqrt(csq(iii))
360	continue
	
	call targetf(npt,np,num,x0, delta, rad, weight, target)
cc	call tarff(npt,np,num,x0, delta, rad, weight, target)
	do 1111 i=1,np
	write(2,52)rad(i), (target(i,j),j=1,num)
1111	continue
	

c	set the rhs
c	-----------
	print*,'calc. rhs'
	
	do 400 j=1,num

		do 405 i=1, nmode
		sum1=0.
		sum2=0.
		do 410 k=1,np
		sum1=sum1+weight(k)*target(k,j)*fakc(k,i)
410	continue
	vv1(i,j)=sum1
405	continue
	vv1(nmode+1,j)=0.5

400	continue



c	NOW FOR THE SOLUTIONS
c	=====================

	do 700 KK=1,40
	read(*,*)beta, amu11, beta2, amu22
	amu1=tan(amu11)
	if(beta.lt.0)stop

c	set matrix for lhs
c	------------------
	print*,'calc. lhs'
	do 710 j=1,nmode
		do 720 i=1,nmode
		cc=cov(i,j)/sumc
		aa1(i,j)=ov1(i,j)+beta*ov2(i,j)+amu1*cc
720	continue
710	continue

	do 711 j=1,numker
	do 712 i=nmode+1, numker
	aa1(i,j)=ov1(i,j)
	aa1(j,i)=ov1(j,i)
712	continue
		do 730 k=1,num
		ctil1(j,k)=vv1(j,k)
730	continue
711	continue

c	do 7101 ijk=1,numker
c	write(55)(aa1(ijk,i),i=1,numker)
c7101	continue
cc	write(37)numker,numker,((aa1(i,j),i=1,numker),j=1,numker)



	if(irep.eq.0)then
c	solving the equations
c	--------------------

	print*,'solving eqn'
	call solve(nmd,nmd,numker,num,aa1,ctil1,ierr,ipiv)
	print*,'done c'
	do 1029 j=1,nmode+1
	write(37,52)(ctil1(j,kjj),kjj=1,num)
1029	continue
	write(27)numker,num,((ctil1(i,j),i=1,numker),j=1,num)
	write(27)nmode,(cov(icov,icov), icov=1,nmode)
	write(27) num, (x0(in), in=1,num)
	ELSE
	read(27)numker,num,((ctil1(i,j),i=1,numker),j=1,num)
	read(27)
	read(27)
	ENDIF


c	Calculating averaging kernel
c	-------------------------------
cc	write(19,*)beta,amu11
cc	write(20,*)beta2,amu22
	do 1112 j=1,np

		do 1113 jj=1,num
		sumcc(jj)=0.
		sumrc(jj)=0.
		do 1114 i=1,nmode
		sumcc(jj)=sumcc(jj)+ctil1(i,jj)*fakc(j,i)
		sumrc(jj)=sumrc(jj)+ctil1(i,jj)*fakr(j,i)
1114	continue
	avc(j,jj)=sumcc(jj)
	avrc(j,jj)=sumrc(jj)
1113	continue


	write(19,52)rad(j),(sumcc(jj),jj=1,num)
	write(20,52)rad(j),(sumrc(jj),jj=1,num)
1112	continue

c	Finding cg of averaging kernel
c	------------------------------

	do 1115 jj=1,num
		sumc1=0.
		sumc2=0.
		akint=0.0
		do 1116 j=1,np
		sumc1=sumc1+weight(j)*rad(j)*avc(j,jj)
		sumc2=sumc2+weight(j)*avc(j,jj)
		akint=akint+weight(j)*(avrc(j,jj)**2)
1116		continue
	print*, 'int. ', sumc2, sumr2
	cenc(jj)=sumc1/sumc2
	ckint(jj)=sumc2
	acros(jj)=sqrt(akint)
1115	continue
cc	print 52, cenc
	print 52,cenr
	
	print*,'calling width'
	call fwidth(npt,np,num,rad,avc,widthc)

c	the solution
c	------------
	write(10,53)amu11, beta, lam, deldel
	write(40,53)amu11, beta, lam, deldel
	write(10,*)
	write(30,*)
cc	if(ifil.eq.1)then
cc	read(3)n1,n1,((cov(i,jj),i=1,n1),jj=1,n1)
cc	print*,'reading full cov matrix'
cc	endif
	do 415 j=1,num
	sum1=0.
	er1=0.

	do 420 i=1,nmode
	sum1=sum1+ctil1(i,j)*w(i)
	if(j.eq.63)then
	lll=rl(i)+.3
	nnn=n(i)
	write(43,55)lll,nnn, ctil1(i,j)*w(i)
	endif
	er1=er1+ctil1(i,j)**2*cov(i,i)
	oll=ctil1(i,j)**2*cov(i,i)
cc		if(j.eq.65)write(37,*)i,oll,er1
	
420	continue

	er1=sqrt(er1)
	s1=abs(widthc(1,j)-widthc(3,j))

	if(ifil.eq.1)then
	er1=0.d0
	er2=0.d0
	print*,' calculating error '
	do 9111 k=1,nmode
		si1=0.
		si2=0.
		do 912 i=1,nmode
		si1=si1+ctil1(i,j)*cov(k,i)
912		continue
		er1=er1+ctil1(k,j)*si1
9111		continue
	er1=sqrt(er1)
	endif

	diff1=0.
	do 421 i=1,np
	diff1=diff1+weight(i)*(avc(i,j)-target(i,j))**2
421	continue

	print*, x0(j), sum1, sum2
	write(10,54)x0(j), delta(j), cenc(j), (widthc(i,j),i=3,1,-1),
     1  s1, sum1, er1, acros(j),diff1
	xc(j)=widthc(2,j)
	sol(j)=sum1

	reps=1.e-7
	ntab=np
	nuse=4
	if(icase.eq.0)then
      call DIVDIF(Xc(j),rad,csq,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	csaa=sqrt(fb(nuse)*(1.+sum1)*csu)
	fbb=fb(nuse)*csu
	nuse=4
      call DIVDIF(Xc(j),rad,rhol,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
cc	summ=sum2*sqrt(1.+sum2*sum2/4.)-sum2*sum2/2
cc	csaa=fb(nuse)*(1.+summ)*rhou
	fb(nuse)=10.d0**fb(nuse)
	fbb=fb(nuse)*rhou*(1.+sum1)

	elseif (icase.eq.1)then
      call DIVDIF(Xc(j),rad,csq,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	csaa=(fb(nuse)*(1.+sum1)*csu)
	csa=fb(nuse)*sum1*csu
	nuse=4
      call DIVDIF(Xc(j),rad,uuu,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	fbb=(fb(nuse)*(1.+sum1)*csu)

	elseif (icase.eq.2)then
      call DIVDIF(Xc(j),rad,gam,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	fbb=(fb(nuse)*(1.+sum1))
	nuse=4
      call DIVDIF(Xc(j),rad,csq,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	csaa=(fb(nuse)*(1.+sum1)*csu)

	endif

	write(40,52)xc(j),fbb,sqrt(csaa)

415	continue

	write(40,*)
	write(40,*)
	do 713 i=1,10
	xx=0.005*(i-1)
	if(icase.ne.0)then
	nuse=3
	call DIVDIF(Xx,x0,sol,NUSE,num,FB,REPS,IER,DFB,DDFB)
	sum=fb(nuse)
	nuse=3
	call DIVDIF(Xx ,rad,csq,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	 csaa=(fb(nuse)*(1.+sum)*csu)
	write(40,*)xx,sum,sqrt(csaa)
	else
	nuse=3
        call DIVDIF(Xx,x0,sol,NUSE,num,FB,REPS,IER,DFB,DDFB)
        sum=fb(nuse)
        nuse=3
	call DIVDIF(Xx,rdd,rhl,NUSE,NTABb,FB,REPS,IER,DFB,DDFB)
        fb(nuse)=10.d0**fb(nuse)
        fbb=fb(nuse)*rhou*(1.+sum)
	write(40,*)xx,sum,fbb
	endif
	
713	continue

	
700	continue
	end
		
c	=======================================================
c	       SUBROUTINES AND FUNCTIONS
c	=======================================================

	subroutine constr(nmd, nmode,lam,om,energy,surcon,omax,omin,
     1  pp,ome, enl, nsp,omu)
	implicit real*8(a-h,o-z)
	parameter(nx0=100)
	parameter(pi=3.1415926535d0)
	dimension surcon(nmd,nx0), om(nmode),energy(nmode)
	dimension pp(nx0,nmode), fb(10)
	dimension enl(50),ome(50),cnl(3,50)
	dimension co(4,200),xo(-3:200)

	call spline(ome, enl, nsp, cnl, ier)
        print*,'spline ier=', ier, nsp

	print*, omax,omin,' omax, omin'
	den=omax-omin
cc	omax=omax*(5.d2*omu)/pi
cc	omin=omin*(5.d2*omu)/pi
	print*, omax,omin,' omax, omin'
	write(33,*) omax,omin,' omax, omin'
	doo=(omax-omin)/(lam-1)

	do 150 i=1,lam
	xo(i)=omin+(i-1)*doo
150	continue

	print*,'lam ',lam
	ic=4
	call bspcof(lam,xo,co,ic)
	print*, 'called bspcof'

	do 100 i=1,nmode-1
	x=om(i)*(5.d2*omu)/pi
	x=om(i)
cc	print*,'i in loop', i
	  do 120 j=1,lam+2
	  pp(j,i)=bspev(lam,xo,co,ic,x,j)
120	continue
100	continue
	print*,'called bspev'

	do 200 i=1,nmode-1
cc	oo=om(i)*(5.d5*omu)/pi
	oo=om(i)
	en0=splevl(oo,nsp,ome,enl,cnl,dfb,ddfb,ier)
	en0=exp(en0)
        nuse=4
	reps=1.e-5
      call DIVDIF(oo,ome,enl,NUSE,Nsp,FB,REPS,IER,DFB,DDFB)
	en0=fb(nuse)
	en0=exp(en0)
	print *,en0
cC	en0=1.
	eee=energy(i)/en0
	print*,oo, eee
	write(33,*)oo, eee, energy(i)
	do 200 j=1,lam+2
	surcon(i,j)=pp(j,i)/eee
200	continue
	
	do 300 j=1,lam+2
	surcon(nmode,j)=0.
300	continue

	return
	end
c	-------------------------
	
	subroutine targetf(npt,np,num,x0,delta, rad,weight, target)
	implicit real*8(a-h,o-z)
	dimension target(npt,num), x0(num), delta(num)
	dimension weight(np), rad(np),sum(300)

	print*,(rad(i),i=1,np)
	print*,num,(x0(i),i=1,num),(delta(i),i=1,num)

	do 150 j=1,num
	sum(j)=0.
	tt=delta(j)/sqrt(2.)

	print*, x0(j), tt

	del=delta(j)

		if(x0(j).lt.0.94)then
	term= (2.*(x0(j)**2)-del*del)/(2.*x0(j))

	do 111 i=1,np
	xxx=rad(i)
	arg=(xxx-term)/del
	arg=arg*arg
	if(arg.gt.120d0) then
	target(i,j)=0
	else
	target(i,j)=xxx*exp(-arg)
	endif

cc	if(xxx.gt..95.and.x0(j).gt..25)then
cc	target(i,j)=target(i,j)+(xxx-.95)
cc	elseif(xxx.gt..95.and.x0(j).le..25)then
cc	target(i,j)=target(i,j)+(xxx-.95)/4.
cc	endif
	sum(j)=sum(j)+weight(i)*target(i,j)
111	continue

		else
        if(del.lt.0.012)del=0.012
	if(x0(j).eq.1.)x0(j)=0.99d0
        term=x0(j)-(del*del)/(2*(1.-x0(j)))


        do 101 i=1,np
        xxx=rad(i)
        arg=(xxx-term)/del
        arg=arg*arg
        if(arg.gt.120) then
        target(i,j)=0
        else
        target(i,j)=(1-xxx)*exp(-arg)
        endif
        sum(j)=sum(j)+weight(i)*target(i,j)
101     continue

        endif


150	continue

cc	do 150 j=1,num
cc	sum(j)=0.
cc	do 100 i=1,np
cc	arg=(rad(i)-x0(j))/delta(j)
cc	arg=arg**2
cc	if(arg.gt.120)then
cc	target(i,j)=0.d0
cc	else
cc	target(i,j)=exp(-arg)
cc	endif
cc	sum(j)=sum(j)+weight(i)*target(i,j)
cc100	continue
cc150	continue

	print*,(sum(j),j=1,num)
	do 200 j=1,num
		do 250 i=1,np
	target(i,j)=target(i,j)/sum(j)
250	continue
200	continue
	
	return
	end
c	----
c	-ve side lobe targets

	subroutine tarff(npt,np,num,x0,delta, rad,weight, target)
	implicit real*8(a-h,o-z)
	dimension target(npt,num), x0(num), delta(num)
	dimension weight(np), rad(np),sum(300)

	print*,(rad(i),i=1,np)
	print*,num,(x0(i),i=1,num),(delta(i),i=1,num)

	do 150 j=1,num
	sum(j)=0.
	tt=delta(j)/sqrt(2.)

	print*, x0(j), tt

	del=delta(j)

		if(x0(j).lt.0.94)then
	term= x0(j)

	do 111 i=1,np
	xxx=rad(i)
	arg=(xxx-term)/del
	arg=arg*arg
	if(arg.gt.120) then
	target(i,j)=0
	else
	target(i,j)=(1.5-arg)*exp(-arg)
	endif

	sum(j)=sum(j)+weight(i)*target(i,j)
111	continue

		else
        if(del.lt.0.012)del=0.012
	if(x0(j).eq.1.)x0(j)=0.99d0

	term=x0(j)

        do 101 i=1,np
        xxx=rad(i)
        arg=(xxx-term)/del
        arg=arg*arg
        if(arg.gt.120) then
        target(i,j)=0
        else
        target(i,j)=(1.5-arg)*exp(-arg)
        endif
        sum(j)=sum(j)+weight(i)*target(i,j)
101     continue

        endif


150	continue

cc	do 150 j=1,num
cc	sum(j)=0.
cc	do 100 i=1,np
cc	arg=(rad(i)-x0(j))/delta(j)
cc	arg=arg**2
cc	if(arg.gt.120)then
cc	target(i,j)=0.d0
cc	else
cc	target(i,j)=dexp(-arg)
cc	endif
cc	sum(j)=sum(j)+weight(i)*target(i,j)
cc100	continue
cc150	continue

	print*,(sum(j),j=1,num)
	do 200 j=1,num
		do 250 i=1,np
	target(i,j)=target(i,j)/sum(j)
250	continue
200	continue
	
	return
	end
c	-----------------------------
	subroutine fwidth(npt,np,num,rad,avc,widthc)
	implicit real*8(a-h,o-z)
	parameter(nn=3500)
	parameter(nx0=100)
	dimension rad(npt), avc(npt,nx0), widthc(3,nx0)
	dimension dx(nn), sumker(nn,nx0), qq(3)
	dimension fb(10),y(nn)
	
	q1=0.25d0
	q2=0.5d0
	q3=0.75d0

	do 100 j=1,np-1
		dx(j)=rad(j)-rad(j+1)
100	continue

	do 200 k=1,num
	sum=0.d0
	sumker(1,k)=0.d0
		do 300 j=1,np-1
		sker=0.5d0*(avc(j,k)+avc(j+1,k))
		sum=sum+dx(j)*sker
		sumker(j+1,k)=sum
300	continue
200	continue

	do 410 j=1,num
	anor=sumker(np,j)
	do 400 i=1,np
	sumker(i,j)=sumker(i,j)/anor
400	continue
	print*, sumker(i,j),sumker(np,j)
410	continue

	ntab=np
	reps=1.d-5

	do 1000 j=1,num
		do 1010 i=1,np
		y(i)=sumker(i,j)
1010	continue
	nuse=4
      call DIVDIF(q1,y,rad,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	widthc(1,j)=fb(nuse)
	nuse=4
      call DIVDIF(q2,y,rad,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	widthc(2,j)=fb(nuse)
	nuse=4
      call DIVDIF(q3,y,rad,NUSE,NTAB,FB,REPS,IER,DFB,DDFB)
	widthc(3,j)=fb(nuse)
1000	continue

	return
	end
c	-------------------------
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
c	SUBROUTINE WHICH CALLS LAPACK ROUTINES
c	=====================================

	subroutine solve(la,lb,norder,nrhs,a,b,ierr,ipiv)
	implicit real*8(a-h,o-z)
	character*1 trans, uplo

	dimension a(la,norder), b(lb,nrhs), ipiv(norder)
	dimension work (10000)
	
	lwork=10000
	uplo='u'
	call dsytrf(uplo, norder, a ,la, ipiv, work, lwork, ier)
	if(ier.ne.0)then
	print*,'error in transformation', ier
	endif

	trans='N'
	uplo='u'
	call dsytrs(uplo,norder,nrhs,a,la,ipiv,b,lb,ierr)
	print*, 'solution ier', ierr
	return
	end


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

	subroutine bspcof(n,x,c,ic)
	implicit real*8(a-h,o-z)
	dimension c(ic,n+2),a(4,4),int(4),x(-3:n+3)

	x(n+1)=2.*x(n)-x(n-1)
	x(n+2)=3.*x(n)-2.*x(n-1)
	x(n+3)=4.*x(n)-3.*x(n-1)
	x(0)=2*x(1)-x(2)
	x(-1)=3.*x(1)-2.*x(2)
	x(-2)=4.*x(1)-3.*x(2)
	x(-3)=5.*x(1)-4.*x(2)
	do 1000 i=1,n+2
	a10=x(i-2)-x(i-3)
	a34=x(i)-x(i+1)
	a21=x(i-1)-x(i-2)
	a23=x(i-1)-x(i)
	a(1,1)=a10*(a10**2+3.*a10*a21+3.*a21**2)
	a(1,2)=0.0
	a(1,3)=a21**3
	a(1,4)=0.0
	c(1,i)=1.0
	a(2,1)=0.0
	a(2,2)=a34*(a34**2+3.*a34*a23+3.*a23**2)
	a(2,3)=0.0
	a(2,4)=a23**3
	c(2,i)=1.0
	a(3,1)=a10*(a10+2.*a21)
	a(3,2)=-a34*(a34+2.*a23)
	a(3,3)=a21**2
	a(3,4)=-a23**2
	c(3,i)=0.0
	a(4,1)=a10
	a(4,2)=-a34
	a(4,3)=a21
	a(4,4)=-a23
	c(4,i)=0.0

	m=4
	num=1
	lj=4
	iflg=0
      call GAUELM(m,NUM,A,c(1,i),DET,INT,LJ,IER,IFLG)
	if(ier.gt.0) print *,'error in gauelm, ier=',ier,x(i)
1000	continue
	end

C	-------------------------------------------------

	function bspevl(n,x,c,ic,wt,x0)
	implicit real*8(a-h,o-z)
	dimension x(-3:n),c(ic,n),wt(n)

	bspevl=0.0
	do 1000 i=1,n+2
	a0=x(i-3)
	a1=x(i-2)
	a2=x(i-1)
	a3=x(i)
	a4=x(i+1)
	if(x0.le.a0.or.x0.ge.a4) then
	bsp=0.0
	else if(x0.le.a1) then
	bsp=c(1,i)*(x0-a0)**3
	else if(x0.le.a2) then
	bsp=c(3,i)*(x0-a1)**3+c(1,i)*(a1-a0)*
     1     (3.*(x0-a1)**2+(a1-a0)*(3.*(x0-a1)+(a1-a0)))
	else if(x0.le.a3) then
	bsp=c(4,i)*(x0-a3)**3+c(2,i)*(a3-a4)*
     1      (3.*(x0-a3)**2+(a3-a4)*(3.*(x0-a3)+(a3-a4)))
	else
	bsp=c(2,i)*(x0-a4)**3
	endif
	bspevl=bspevl+wt(i)*bsp
1000	continue

	end

C	-------------------------------------------------

	function bspd1(n,x,c,ic,wt,x0)
	implicit real*8(a-h,o-z)
	dimension x(-3:n),c(ic,n),wt(n)

	bspd1=0.0
	do 1000 i=1,n+2
	a0=x(i-3)
	a1=x(i-2)
	a2=x(i-1)
	a3=x(i)
	a4=x(i+1)
	if(x0.le.a0.or.x0.ge.a4) then
	bsp=0.0
	else if(x0.le.a1) then
	bsp=3.*c(1,i)*(x0-a0)**2
	else if(x0.le.a2) then
	bsp=3.*c(3,i)*(x0-a1)**2+c(1,i)*(a1-a0)*
     1     (6.*(x0-a1)+3.*(a1-a0))
	else if(x0.le.a3) then
	bsp=3.*c(4,i)*(x0-a3)**2+c(2,i)*(a3-a4)*
     1      (6.*(x0-a3)+3.*(a3-a4))
	else
	bsp=3.*c(2,i)*(x0-a4)**2
	endif
	bspd1=bspd1+wt(i)*bsp
1000	continue

	end
C	------------------------------------------------------------
	function bspd2a(n,x,c,ic,wt,x0)
	implicit real*8(a-h,o-z)
	dimension x(-3:n),c(ic,n),wt(n+2)

	bspd2a=0.0
	do 1000 i=1,n+2
	a0=x(i-3)
	a1=x(i-2)
	a2=x(i-1)
	a3=x(i)
	a4=x(i+1)
	if(x0.le.a0.or.x0.ge.a4) then
	bsp=0.0
	else if(x0.le.a1) then
	bsp=6.*c(1,i)*(x0-a0)
	else if(x0.le.a2) then
	bsp=6.*c(3,i)*(x0-a1)+c(1,i)*(a1-a0)*6
	else if(x0.le.a3) then
	bsp=6.*c(4,i)*(x0-a3)+c(2,i)*(a3-a4)*6
	else
	bsp=6.*c(2,i)*(x0-a4)
	endif
	bspd2a=bspd2a+wt(i)*bsp
1000	continue

	end

C	-------------------------------------------------

      SUBROUTINE GAUELM(N,NUM,A,X,DET,INT,LJ,IER,IFLG)
	implicit real*8(a-h,o-z)
      DIMENSION A(LJ,N),INT(N),X(LJ,NUM)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=111
        RETURN
      ENDIF

      IER=122
      IF(IFLG.LE.1) THEN
        DET=1.0
        DO 2600 K=1,N-1
          R1=0.0
          DO 2200 L=K,N
            IF(ABS(A(L,K)).GT.R1) THEN
              R1=ABS(A(L,K))
              KM=L
            ENDIF
2200      CONTINUE

          INT(K)=KM
          IF(KM.NE.K) THEN
            DO 2300 L=K,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
2300        A(KM,L)=T1
            DET=-DET
          ENDIF

          DET=DET*A(K,K)
          IF(A(K,K).EQ.0.0) RETURN
C         IF(ABS(A(K,K)).LT.REPS) RETURN
          DO 2500 L=K+1,N
            A(L,K)=A(L,K)/A(K,K)
            DO 2500 L1=K+1,N
2500      A(L,L1)=A(L,L1)-A(L,K)*A(K,L1)
2600    CONTINUE
        DET=DET*A(N,N)
        INT(N)=N
        IF(A(N,N).EQ.0.0) RETURN
C         IF(ABS(A(N,N)).LT.REPS) RETURN

        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
          IF(K.NE.INT(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INT(K),J)
            X(INT(K),J)=T1
          ENDIF
          DO 3000 L=K+1,N
3000    X(L,J)=X(L,J)-A(L,K)*X(K,J)

        X(N,J)=X(N,J)/A(N,N)
        DO 3300 K=N-1,1,-1
          DO 3200 L=N,K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L)
3300    X(K,J)=X(K,J)/A(K,K)
5000  CONTINUE
      END

	function bspev(n,x,c,ic,x0,i)
	implicit real*8(a-h,o-z)
	dimension x(-3:n+3),c(ic,n+2)

	a0=x(i-3)
	a1=x(i-2)
	a2=x(i-1)
	a3=x(i)
	a4=x(i+1)
	if(x0.le.a0.or.x0.ge.a4) then
	bspev=0.0
	return
	else if(x0.le.a1) then
	bspev=c(1,i)*(x0-a0)**3
	else if(x0.le.a2) then
	bspev=c(3,i)*(x0-a1)**3+c(1,i)*(a1-a0)*(3.*(x0-a1)**2
     1   +(a1-a0)*(3.*(x0-a1)+(a1-a0)))
	else if(x0.le.a3) then
	bspev=c(4,i)*(x0-a3)**3+c(2,i)*(a3-a4)*(3.*(x0-a3)**2
     1   +(a3-a4)*(3.*(x0-a3)+(a3-a4)))
	else
	bspev=c(2,i)*(x0-a4)**3
	endif
	end

	subroutine constr1(nmd, nmode,lam,om,energy,surcon,omax,omin,
     1  pp,ome, enl, nsp,omu)
	implicit real*8(a-h,o-z)
	parameter(nx0=100)
	parameter(pi=3.1415926535d0)
	dimension surcon(nmd,nx0), om(nmode),energy(nmode)
	dimension pp(0:nx0,nmode)
	dimension enl(50),ome(50),cnl(3,50)

	call spline(ome, enl, nsp, cnl, ier)
        print*,'spline ier=', ier

	den=omax-omin
	do 100 i=1,nmode-1
	x=(2.d0*om(i)-omin-omax)/den
	pp(0,i)=1.d0

	if(lam.ge.0)then
	pp(1,i)=x
		do 150 j=1,lam-1
c	Legendre pol
	pp(j+1,i)=((2.d0*j+1)*x*pp(j,i)-j*pp(j-1,i))/(j+1.d0)
c	chebyshev
cc	pp(j+1,i)=2.d0*x*pp(j,i)-pp(j-1,i)
150	continue
	endif
100	continue

	do 200 i=1,nmode-1
	oo=om(i)*(5.d5*omu)/pi
	oo=om(i)
	en0=splevl(oo,nsp,ome,enl,cnl,dfb,ddfb,ier)
	en0=exp(en0)
cc	print *,en0
	eee=energy(i)/en0
	write(33,*)oo, eee
	do 200 j=1,lam+1
	surcon(i,j)=1.d4*pp(j-1,i)/eee
200	continue
	
	do 300 j=1,lam+1
	surcon(nmode,j)=0.
300	continue

	return
	end
c	-------------------------
