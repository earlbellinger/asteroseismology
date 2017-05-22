C     PROGRAM TO FIND ALL COMPLEX ZEROS OF AN ANALYTIC FUNCTION INSIDE A
C     GIVEN CIRCLE USING QUADRATURE BASED METHOD

      PROGRAM ROOT
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      EXTERNAL CF
      DIMENSION CZERO(10)

C     EXAMPLE 7.8

51    FORMAT('   IER =',I4,5X,'NO. OF ZEROS =',I2/5X,'CENTRE =',
     1   1P2D14.6,5X,'RADIUS =',D14.6/'    ZEROS:',(1P2D14.6,2X,2D14.6))

      REPS=1.D-7
      AEPS=1.D-8

100   PRINT *,'TYPE CZ=CENTRE(COMPLEX VARIABLE), RAD=RADIUS OF CIRCLE'
      PRINT *,'             (QUITS WHEN RAD.LE.0)'
      READ *,CZ,RAD
      IF(RAD.LE.0) STOP
      CALL DELVES(CZ,RAD,CF,NZ,CZERO,IER,REPS,AEPS)

C	If NZ turns out to be very large no zeros are determined but
C	the write statement may give error trying to write nonexistent
C	numbers. This can happen if one of the zero falls very close
C	to the circle.

      WRITE(6,51) IER,NZ,CZ,RAD,(CZERO(I),I=1,NZ)
      GO TO 100
      END

C     -------------------------------------------

C	To find contour integrals over a circle as required by subroutine DELVES
C
C	CF : (input) Name of the function routine to calculate the function value
C	CZ : (input) Complex variable specifying the centre of the required circle
C	RAD : (input) Real variable specifying the radius of the circle.
C	NZ : (output) Number of zeros located by the subroutine inside the circle
C	CS : (output) Complex array of length NMAX containing the values
C		of contour integrals
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=44 implies that the first integration to determine NZ
C			did not converge to satisfactory accuracy, but its
C			magnitude is less than 0.5. NZ is set to zero in this case.
C		IER=45 implies that the first integration to determine NZ
C			did not converge to satisfactory accuracy, but a
C			weaker criterion is satisfied.
C		IER=46 implies that one of the higher integrals to determine
C			CS(I) did not converge to satisfactory accuracy, but a
C			weaker criterion is satisfied.
C		IER=434 implies that NZ>NMAX and no further calculations are done
C		IER=437 implies that the first integration to determine NZ
C			did not converge satisfactorily, no further calculations 
C			are done.
C	NMAX : (input) Maximum number of zeros to be located. If NZ>NMAX
C		no further calculations are done.
C	
C	FUNCTION CF(CX,CDF) must be supplied by the user. Here CF, CX and CDF are
C		complex variables and CDF is the first derivative of CF.
C
C	Required routines : CF

      SUBROUTINE CONTUR(CF,CZ,RAD,NZ,CS,IER,NMAX)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      PARAMETER(REPS0=1D-2,REPS=1D-4,NPT=1024,PI=3.14159265358979324D0)
      PARAMETER(CI=(0.D0,1.D0))
      DIMENSION CFUN(NPT),CZP(NPT),CS(NMAX),CSM(NPT)

      IER=0
      NP=16
      CXP=EXP(CI*2.*PI/NP)
      CZ0=1.
      CSUM=0.0
      CNT=0.0
      II=0

C	Trapezoidal rule for the first integral
C	I=1,7 will try NP from 16 to 1024=NPT
C	To increase number of points increase NPT and the loop count
      DO 2000 I=1,7
        DO 1000 J=1,NP
          II=II+1
          CZP(II)=CZ0
          CFN=CF(RAD*CZ0+CZ,CDF)
          CFUN(II)=CZ0*CDF/CFN
          CZ0=CZ0*CXP
          CSUM=CSUM+CFUN(II)
1000    CONTINUE

        IF(I.EQ.1) THEN
          CZ0=EXP(CI*PI/NP)
        ELSE
          NP=NP*2
          CZ0=EXP(CI*PI/NP)
          CXP=CZ0*CZ0
        ENDIF

C	Latest approximation to the integral
        CINT=CSUM*RAD/NP
        NZ=CINT+0.4D0
        EPS=ABS(CINT-CNT)
        EPS1=ABS(CINT-NZ)
        IF(EPS.LT.REPS0.AND.EPS1.LT.0.1D0) GO TO 2500
        CNT=CINT
2000  CONTINUE

      IF(ABS(CINT).LT.0.5) THEN
C	The integral is approximately zero
        IER=44
        NZ=0
        RETURN
      ENDIF

      IF(EPS1.GT.0.1D0) THEN
C	Integral fails to converge
        IER=437
        RETURN
      ENDIF
C	Accept the value of NZ, although integral hasn't converged properly
      IER=45

2500  NZ=CINT+0.4D0
      IF(NZ.LE.0) RETURN
      IF(NZ.GT.NMAX) THEN
C	Too many zeros inside the circle
        IER=434
        RETURN
      ENDIF

      RMAX=MAX(ABS(CZ),RAD)
      DO 3000 I=1,NZ
C	Contour integrals using trapezoidal rule
        CSUM=0.0
        DO 2800 J=1,NP
          CSUM=CSUM+CFUN(J)*CZP(J)**I
2800    CONTINUE
        CSM(I)=CSUM
        CS(I)=CSUM*RAD**(I+1)/NP
3000  CONTINUE
      II=NP

C	Double the number of points until convergence
3100  CONTINUE
      DO 3400 J=1,NP
        CFN=CF(RAD*CZ0+CZ,CDF)
        CFUNI=CZ0*CDF/CFN
        DO 3200 K=1,NZ
          CSM(K)=CSM(K)+CFUNI*CZ0**K
3200    CONTINUE
        CZ0=CZ0*CXP
3400  CONTINUE

      NP=NP*2
      CZ0=EXP(CI*PI/NP)
      CXP=CZ0*CZ0

C	The convergence test
      QC=.TRUE.
      DO 3600 K=1,NZ
        CINT=CSM(K)*RAD**(K+1)/NP
        EPS=ABS(CINT-CS(K))
        IF(EPS.GT.REPS*MAX(ABS(CINT),RMAX**K)) QC=.FALSE.
        CS(K)=CINT
3600  CONTINUE

      IF(QC) RETURN
      IF(NP.LE.16*NPT) GO TO 3100

C	Integrals fail to converge
      IER=46
      END

C     ------------------------------------------------

C	Complex zeros of analytic function inside a given circle
C
C	CZ : (input) Complex variable specifying the centre of the required circle
C	RAD : (input) Real variable specifying the radius of the circle.
C		All zeros inside the circle with centre CZ and radius RAD
C		need to be found.
C	CF : (input) Name of the function routine to calculate the function value
C	NZ : (output) Number of zeros located by the subroutine
C	CZERO : (output) Complex array of length NMAX (=5) containing the
C		zeros determined by the subroutine
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=405 implies RAD.LE.0, in which case no calculations are done
C		IER=434 implies that NZ>NMAX and zeros are not determined
C		IER=435 implies that iteration failed to converge to specified
C			accuracy for at least one zero
C		IER=436 implies that POLYC failed to find all roots of the
C			required polynomial
C		Other values of IER may be set by CONTUR
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than max(AEPS, REPS*ABS(CZERO))
C	
C	FUNCTION CF(CX,CDF) must be supplied by the user. Here CF, CX and CDF are
C		all complex variables, and CDF is the first derivative of CF.
C	Required routines : CONTUR, NEWRAC, POLYC, LAGITC, CF
C
      SUBROUTINE DELVES(CZ,RAD,CF,NZ,CZERO,IER,REPS,AEPS)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      PARAMETER(NMAX=5)
      EXTERNAL CF
      DIMENSION CZERO(NMAX),CS(NMAX),COF(NMAX+1),CWK(NMAX+1)

      NZ=0
      IER=405
      IF(RAD.LE.0.0) RETURN

C	To calculate the contour integrals
      CALL CONTUR(CF,CZ,RAD,NZ,CS,IER,NMAX)
      IF(IER.GT.100) RETURN
      IF(NZ.LE.0) RETURN
      IF(NZ.GT.NMAX) THEN
C	If the number of zeros is too large, then quit
        IER=434
        RETURN
      ENDIF

C	Coefficients of the equivalent polynomial
      COF(NZ+1)=1.0
      DO 2000 N=NZ,1,-1
        J=NZ+1-N
        C1=CS(J)
        DO 1500 I=1,J-1
1500    C1=C1+CS(J-I)*COF(NZ+1-I)
        COF(N)=-C1/J
2000  CONTINUE

C	Find zeros of the polynomial
      QREFIN=.FALSE.
      CALL POLYC(NZ,COF,CZERO,IER1,QREFIN,CWK)
      IF(IER1.GT.100) IER=436
      DO 3000 I=1,NZ
3000  CZERO(I)=CZERO(I)+CZ

C	Refine the zeros by iterating on the original function
      DO 4000 I=1,NZ
        CALL NEWRAC(CZERO(I),REPS,AEPS,IER1,CF)
        IF(IER1.GT.100) IER=435
4000  CONTINUE
      END

C  ------------------------------------------------------------

C	Root of a polynomial with complex coefficients using Laguerre iteration
C
C	N : (input) The degree of the polynomial
C	COF : (input) Complex array of length N+1 containing the coefficients of
C		the polynomial. COF(I) is the coefficient of X**(I-1)
C	CXI : (input/output) Complex variable containing the initial guess,
C		 after execution it will contain the computed root
C	IER : (output) Error parameter, IER=0 for successful execution
C		IER=438 implies that denominator is zero and iteration cannot
C			be continued further
C		IER=439 implies that iteration has failed to converge
C
C	Required routines : None

      SUBROUTINE LAGITC(N,COF,CXI,IER)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      PARAMETER(ITMAX=50,REPS=1.D-4,AEPS=1.D-6)
      DIMENSION COF(N+1)

      IER=0
      CDX1=ABS(CXI)+1
      QC=.FALSE.
      IC=ITMAX

C	The Laguerre's iteration
      DO 2000 I=1,ITMAX
        CF=COF(N+1)
        CFP=0.0
        CFPP=0.0
        DO 1000 J=N,1,-1
          CFPP=CXI*CFPP+2.*CFP
          CFP=CXI*CFP+CF
          CF=CXI*CF+COF(J)
1000    CONTINUE

        CH=(N-1)*((N-1)*CFP*CFP-N*CF*CFPP)
        CH=SQRT(CH)
        CDEN=CFP+CH
        IF(ABS(CFP-CH).GT.ABS(CDEN)) CDEN=CFP-CH

        IF(CDEN.NE.0.0) THEN
          CDX=-N*CF/CDEN
          IF(ABS(CDX).LT.MAX(REPS*ABS(CXI),AEPS).AND.I.GT.1.AND.
     1        (.NOT.QC)) THEN
            QC=.TRUE.
            IC=I
          ENDIF
          IF(QC.AND.ABS(CDX/CDX1).GT.1.0) RETURN
          IF(I-IC.GT.5.AND.ABS(CDX/CDX1).GT.0.99D0) RETURN
          CDX1=CDX
          IF(CDX.EQ.0.0) RETURN
          CXI=CXI+CDX
        ELSE
          IF(CF.EQ.0.0) RETURN
C	If the denominator vanishes, then quit
          IER=438
          RETURN
        ENDIF
2000  CONTINUE

C	Iteration fails to converge
      IF(.NOT.QC) IER=439
      END

C     ----------------------------------------

C	Complex zero of a given function using Newton-Raphson iteration
C
C	CX : (input/output) Initial guess for the zero, after execution
C		it will contain the computed zero.
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(CX))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=425 implies that iteration failed to converge to specified accuracy
C		IER=426 implies that derivative is zero and the iteration
C			has to be abandoned at this stage.
C	CFUN : (input) Name of the function routine to calculate the function
C		FUNCTION CFUN(CX,CDX) must be supplied by the user.
C		Here CDX is the first derivative of CFUN.
C
C	Required routines : CFUN

      SUBROUTINE NEWRAC(CX,REPS,AEPS,IER,CFUN)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      PARAMETER(NIT=75)

      IER=0
      DO 1000 L=1,NIT
        CF=CFUN(CX,CDF)
        IF(CDF.EQ.0.0) THEN
          IF(CF.EQ.0.0) RETURN
          IER=426
          RETURN
        ENDIF

C	Newton-Raphson iteration
        CDX=-CF/CDF
        CX=CX+CDX
        IF(ABS(CDX).LT.MAX(REPS*ABS(CX),AEPS).AND.L.GT.2) RETURN
1000  CONTINUE

      IER=425
      END

C     -------------------------------------------

C	Roots of a polynomial with complex coefficients using Laguerre iteration
C
C	N : (input) The degree of the polynomial
C	COF : (input) Complex array of length N+1 containing the coefficients of
C		the polynomial. COF(I) is the coefficient of X**(I-1)
C	CX : (output) Complex array of length N, containing the computed roots
C	IER : (output) Error parameter, IER=0 for successful execution
C		IER=k*11 implies that iteration for refining the roots failed
C			to converge for k of the roots
C		IER=406 implies that N.LE.0 and no calculations are done
C		IER=408 implies that A(N+1)=0 and no calculations are done
C		IER=430 implies that iteration failed to converge for some root
C	QREFIN : (input) Logical parameter to decide if roots need to be refined
C		If QREFIN=.TRUE. the roots are refined using original polynomial
C		otherwise no refinement is done.
C	CWK : Complex array of length N+1 used as a scratch space
C
C	Required routines : LAGITC
C
      SUBROUTINE POLYC(N,COF,CX,IER,QREFIN,CWK)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      DIMENSION COF(N+1),CX(N),CWK(N+1)

      IF(N.LE.0) THEN
        IER=406
        RETURN
      ENDIF
      IF(COF(N+1).EQ.0.0) THEN
        IER=408
        RETURN
      ENDIF

      DO 1000 I=1,N+1
1000  CWK(I)=COF(I)
      NP=N
      CXR=0.0
      IER=0

C	Find the next root
2000  CALL LAGITC(NP,CWK,CXR,IER1)
      IF(IER1.NE.0) THEN
C	If iteration fails to converge, try once more
        CXR=1.123456
        CALL LAGITC(NP,CWK,CXR,IER1)
        IF(IER1.NE.0) THEN
C	If iteration fails again, then quit
          IER=430
          RETURN
        ENDIF
      ENDIF

      CXRT=CXR
      IF(NP.LT.N.AND.QREFIN) THEN
C	Refine the roots with original polynomial
        CALL LAGITC(N,COF,CXRT,IER1)
        IF(IER1.NE.0) THEN
C	If iteration fails to converge then retain the old value
          IER=IER+11
          CXRT=CXR
        ENDIF
      ENDIF

      IF(NP.GT.1) THEN
C	Perform deflation using unrefined root
        CN0=CWK(NP)
        CWK(NP)=CWK(NP+1)
        DO 2400 I=NP-1,1,-1
          CN=CXR*CWK(I+1)+CN0
          CN0=CWK(I)
          CWK(I)=CN
2400    CONTINUE
      ENDIF

      NP=NP-1
      CX(N-NP)=CXRT
C	If any more roots are left, find them
      IF(NP.GT.0) GO TO 2000
      END

C     -----------------------

      FUNCTION CF(CZ,CDF)
      IMPLICIT COMPLEX*16(C)

C     THE ANALYTIC FUNCTION AND ITS DERIVATIVE

      CF=CZ+SIN(CZ)
      CDF=1.+COS(CZ)
      END


