C     PROGRAM TO INTEGRATE OSCILLATORY FUNCTIONS USING FILON'S FORMULA

      PROGRAM INTEG
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      EXTERNAL FUN

C     EXERCISE 6.21 : INTEGRAL I1
 
51    FORMAT('   IER =',I4,5X,'A =',1PD14.6,5X,'B =',D14.6,5X,
     1       'RK =',D14.5/'  NO. OF FUNCTION EVALUATIONS =',I7,5X,
     2       'INTEGRAL =',D14.6/9X,'ESTIMATED ERROR =',D14.6,5X,
     3       'EXACT VALUE =',D14.6)

      REPS=1.D-13
      AEPS=1.D-18
      QSIN=.TRUE.
C      QSIN=.FALSE.
      B=4.*ACOS(0.D0)

100   PRINT *, 'TYPE  A=LOWER LIMIT,  B=UPPER LIMIT,'
      READ *,A,B
      PRINT *, 'RK=COEFFICIENT OF X IN OSCILLATORY FACTOR'
      PRINT *, '             (QUITS WHEN RK.EQ.0)'
      READ *,RK
      IF(RK.EQ.0) STOP
      CALL FILON(RI,A,B,RK,QSIN,REPS,AEPS,DIF,IER,N,FUN)

C	The exact value for integral values of RK
      REX=-B*RK/(RK*RK-1)
      WRITE(6,51) IER,A,B,RK,N,RI,DIF,REX
      GO TO 100
      END
 
C     ---------------------------------------------------
 
C	To calculate integrals with oscillatory integrand of the form
C	FUN(x)*SIN(RK*x)  or  FUN(x)*COS(RK*x)
C
C	RI : (output) Computed value of the integral
C	XL : (input) The lower limit
C	XU : (input) The upper limit
C	RK : (input) Coefficient of x in the oscillatory term in the integrand
C	QSIN : (input) Logical variable to specify the form of integrand
C		If QSIN=.TRUE. the oscillatory factor is SIN(RK*x)
C		If QSIN=.FALSE. the oscillatory factor is COS(RK*x)
C	REPS : (input) The required relative accuracy
C	AEPS : (input) The required absolute accuracy
C		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
C	DIF : (output) estimated (absolute) error achieved by the subroutine
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=30 implies specified accuracy was not achieved
C			DIF will contain the estimated accuracy
C	N : (output) Number of function evaluations used by subroutine
C	FUN : (input) Name of the function routine to calculate the integrand
C		excluding the oscillatory factor
C		FUNCTION FUN(X) must be supplied by the user.
C	
C	Required routines : FUN

      SUBROUTINE FILON(RI,XL,XU,RK,QSIN,REPS,AEPS,DIF,IER,N,FUN)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
C	THC should be about (100*hcross)**(1/6)
      PARAMETER(NMIN=5,NMAX=15,THC=0.005D0)
C	For REAL*4 version use this value
C      PARAMETER(NMIN=5,NMAX=15,THC=0.15)

      IF(QSIN) THEN
C	For Sin(kx)
        FEND=FUN(XL)*COS(RK*XL)-FUN(XU)*COS(RK*XU)
        EVEN=0.5*(FUN(XL)*SIN(RK*XL)+FUN(XU)*SIN(RK*XU))
      ELSE
C	For Cos(kx)
        FEND=FUN(XU)*SIN(RK*XU)-FUN(XL)*SIN(RK*XL)
        EVEN=0.5*(FUN(XL)*COS(RK*XL)+FUN(XU)*COS(RK*XU))
      ENDIF

      ODD=0.0
      IER=0
      RI=0.0
      DIF=0.0
      N=2
      H=(XU-XL)
      IF(H.EQ.0.0) RETURN

      DO 3000 I=1,NMAX
        H=H/2.
        EVEN=EVEN+ODD
        ODD=0.0
        X1=XL+H
        H2=2.*H

C	Starting with 3 points subdivide the intervals into 2 until convergence
        DO 1000 J=1,N/2
          X=X1+H2*(J-1)
          IF(QSIN) THEN
            ODD=ODD+FUN(X)*SIN(RK*X)
          ELSE
            ODD=ODD+FUN(X)*COS(RK*X)
          ENDIF
1000    CONTINUE

        T=RK*H
        IF(ABS(T).GT.THC) THEN
C	Use normal functions
          ALPHA=(T*T+SIN(T)*(T*COS(T)-2.*SIN(T)))/T**3
          BETA=2.*(T+COS(T)*(T*COS(T)-2.*SIN(T)))/T**3
          GAMMA=4.*(SIN(T)-T*COS(T))/T**3
        ELSE
C	Use Taylor series expansion
          ALPHA=2.*T**3*(1.+T*T*(-1.+T*T/15.)/7.)/45.
          BETA=2.D0/3.+2.*T*T*(1.D0/15.+T*T*(-2.D0/105.+T*T/567.))
          GAMMA=4.D0/3.+T*T*(-2.D0/15.+T*T*(1.D0/210.-T*T/11340.))
        ENDIF
        R1=H*(ALPHA*FEND+BETA*EVEN+GAMMA*ODD)

        DIF=R1-RI
        RI=R1
        IF(I.LE.NMIN) GO TO 3000
        IF(ABS(DIF).LT.MAX(REPS*ABS(R1),AEPS)) RETURN
3000  N=N*2

      N=N/2
      IER=30
      END
 
C     ---------------------------------------------
 
      FUNCTION FUN(X)
      IMPLICIT REAL*8(A-H,O,P,R-Z)

C     SPECIFY THE INTEGRAND (OMIT THE OSCILLATORY FACTOR COS(KX) OR SIN(KX))

      FUN=X*COS(X)
      END
