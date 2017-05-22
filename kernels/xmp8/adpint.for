C     PROGRAM FOR ADAPTIVE INTEGRATION OVER A FINITE INTERVAL
C	This program uses GAUS16 instead of KRONRD and hence
C	is less efficient as compared to the version used in quad.f

      PROGRAM QUAD
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUN

C     EXAMPLE 6.1: INTEGRATE SQRT(X) OVER [0,1]

53    FORMAT('   A =',1PD14.6,5X,'B =',D14.6,5X,'IER =',I4/
     1       '    NO. OF FUNCTION EVALUATIONS =',I7,5X,'INTEGRAL =',
     2       D14.6/5X,'ESTIMATED ERROR =',D14.6)

      REPS=1.D-13
      AEPS=1.D-18

100   PRINT *,'TYPE   A=LOWER LIMIT,  B=UPPER LIMIT,',
     1        '    (QUITS WHEN A.EQ.B)'
      READ *,A,B
      IF(A.EQ.B) STOP
      NMAX=16385
      CALL ADPINT(RI,A,B,REPS,AEPS,DIF,FUN,IER,NPT,NMAX)
      WRITE(6,53) A,B,IER,NPT,RI,DIF
      GO TO 100
      END
 

C     ----------------------------------------------------------

C	To integrate a function over finite interval using adaptive control
C	of step size
C
C	RINT : (output) Calculated value of the integral
C	XL : (input) The lower limit
C	XU : (input) The upper limit
C	REPS : (input) The required relative accuracy
C	AEPS : (input) The required absolute accuracy
C		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
C	DIF : (output) estimated (absolute) error achieved by the subroutine
C	F : (input) Name of the function routine to calculate the integrand
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=31 implies specified accuracy was not achieved on
C			at least one subinterval
C		IER=32 implies that this failure occurred more than IFMAX (=5) times
C		IER=325 implies that subroutine failed to attain required
C			accuracy using NMAX function evaluations
C		In all cases DIF will contain the estimated accuracy
C	NPT : (output) Number of function evaluations used by subroutine
C	NMAX : (input/output) Maximum number of function evaluations to be tried
C		If NMAX.LE.0 it is set to MAXPT (=100000)
C
C		FUNCTION F(X) must be supplied by the user.
C
C	Required routines : KRONRD (or GAUS16), F

      SUBROUTINE ADPINT(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(IPMAX=100,IFMAX=5,MAXPT=100000)
      EXTERNAL F
      DIMENSION XU1(IPMAX)

      IER=0
      IFAIL=0
      RINT=0.0
      DIF=0.0
      IF(XL.EQ.XU) RETURN
      IF(NMAX.LE.0) NMAX=MAXPT
      AEPSL=AEPS
      NPT=0
      RL=XL
      RU=XU
      IU=0

C	To evaluate the integral over [RL,RU]
C1000  CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
1000  CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RM=0.5*(RL+RU)
C	Q=.TRUE. if the interval cannot be divided further
      Q=IU.GE.IPMAX.OR.RM.EQ.RL.OR.RM.EQ.RU

      IF(DIF0.LT.MAX(ABS(FINT)*REPS,AEPSL).OR.Q) THEN
C	Accept the value of FINT if adequate convergence or if the interval
C	cannot be subdivided further
        RINT=RINT+FINT
        DIF=DIF+DIF0
        IF(Q.AND.DIF0.GT.MAX(ABS(RINT)*REPS,AEPSL)) THEN
C	Integration fails to converge on this subinterval. Go to the next subinterval
          IER=31
          IFAIL=IFAIL+1
          IF(IFAIL.GT.IFMAX) THEN
C	If failure is frequent then adjust the convergence criterion.
            IER=32
            AEPSL=DIF*0.5
          ENDIF
        ENDIF

C	If all subintervals are exhausted then return
        IF(IU.LE.0) RETURN

C	otherwise try next subinterval
        RL=RU
        RU=XU1(IU)
        IU=IU-1
      ELSE

C	Subdivide the current interval and try again
        IU=IU+1
        XU1(IU)=RU
        RU=RM
      ENDIF

      IF(NPT.LT.NMAX) GO TO 1000
C	If the number of function evaluations has exceeded the limit then
C	try a last call to estimate the integral over the remaining interval
      IER=325
      RU=XU
C      CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
      CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END

C     ----------------------------------------------------------

C	To integrate a function over a finite interval using 16 point
C	Gauss-Legendre formula, for use with ADPINT
C
C	RI : (output) Calculated value of the integral
C	A : (input) The lower limit
C	B : (input) The upper limit
C	DIF : (output) estimated (absolute) error achieved by the subroutine
C	N : (output) Number of function evaluations used by subroutine
C	F : (input) Name of the function routine to calculate the integrand
C
C	FUNCTION F(X) must be supplied by the user
C
C	Required routines : F

      SUBROUTINE GAUS16(RI,A,B,DIF,N,F)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  W8(4),A8(4),W16(8),A16(8)

C	W8 and A8 are the weights and abscissas for the 8-point Gauss formula
C	W16 and A16 are the weights and abscissas for the 16-point Gauss formula
C	Because of symmetry only half the points are given.

      DATA W8 /0.10122853629037625915D0, 0.22238103445337447054D0,
     *         0.31370664587788728734D0, 0.36268378337836198297D0/
      DATA A8 /0.96028985649753623168D0, 0.79666647741362673959D0,
     *         0.52553240991632898582D0, 0.18343464249564980494D0/

      DATA W16/0.02715245941175409485D0, 0.06225352393864789286D0,
     *         0.09515851168249278481D0, 0.12462897125553387205D0,
     *         0.14959598881657673208D0, 0.16915651939500253819D0,
     *         0.18260341504492358887D0, 0.18945061045506849629D0/
      DATA A16/0.98940093499164993260D0, 0.94457502307323257608D0,
     *         0.86563120238783174388D0, 0.75540440835500303390D0,
     *         0.61787624440264374845D0, 0.45801677765722738634D0,
     *         0.28160355077925891323D0, 0.09501250983763744019D0/

      AT=(B-A)/2.
      BT=(B+A)/2.
      R1=0.0
C	8-point Gauss-Legendre formula
      DO 2000 K=1,4
        R1=R1+W8(K)*(F(AT*A8(K)+BT)+F(BT-AT*A8(K)))
2000  CONTINUE

      RI=0.0
C	16-point Gauss-Legendre formula
      DO 2500 K=1,8
        RI=RI+W16(K)*(F(AT*A16(K)+BT)+F(BT-AT*A16(K)))
2500  CONTINUE

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=24
      END

C     ----------------------------------------------------------
 
      FUNCTION FUN(X)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE INTEGRAND

      FUN=SQRT(X)
      END
