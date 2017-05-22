C     PROGRAM FOR INTEGRATING F(X)*SQRT((X-A)/(B-X)) OVER (A,B)

      PROGRAM INTEG
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F

C     EXERCISE 6.18 

51    FORMAT('   IER =',I4,5X,'A =',1PD14.6,5X,'B=',D14.6/5X,
     1  'NO. OF FUNCTION EVALUATIONS =',I7/5X,'INTEGRAL =',D14.6,5X,
     2       'ESTIMATED ERROR =',3D14.6)

      REPS=1.D-14
      AEPS=1.D-19

100   PRINT *,'TYPE  A=LOWER LIMIT,  B=UPPER LIMIT   QUITS WHEN(A.EQ.B)'
      READ *,A,B
      IF(A.EQ.B) STOP
      CALL GAUCB1(RI,A,B,REPS,AEPS,DIF,IER,N,F)
      WRITE(6,51) IER,A,B,N,RI,DIF
      GO TO 100
      END
 
C     ---------------------------------------------

C	To integrate a function over finite interval using Gauss-Chebyshev formulas
C	Calculates the integral of FUN(X)*SQRT(X-A)/SQRT(B-X)
C
C	RINT : (output) Calculated value of the integral
C	A : (input) The lower limit
C	B : (input) The upper limit
C	REPS : (input) The required relative accuracy
C	AEPS : (input) The required absolute accuracy
C		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
C	DIF : (output) estimated (absolute) error achieved by the subroutine
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=30 implies specified accuracy was not achieved
C			DIF will contain the estimated accuracy
C	NPT : (output) Number of function evaluations used by the subroutine
C	FUN : (input) Name of the function routine to calculate the
C		integrand multiplied by SQRT((B-X)/(X-A))
C
C		FUNCTION FUN(X) must be supplied by the user.
C
C	Required routines : FUN

      SUBROUTINE GAUCB1(RINT,A,B,REPS,AEPS,DIF,IER,NPT,FUN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=13,PI=3.14159265358979324D0)
 
      N=1
      IER=0
      RINT=0.0
      NPT=0
      DA=(B-A)
 
      DO 5000 I=1,NMAX
        N=N*2
        DX=PI/(4*N+2)
 
C	Apply N-point formula after transforming the range to (0,1)
        R1=0.0
        DO 3000 J=1,N
          A1=COS((2*J-1)*DX)**2
          R1=R1+A1*FUN(A+DA*A1)
3000    CONTINUE
        R1=R1*DX*4.*DA
 
        DIF=R1-RINT
        RINT=R1
        NPT=NPT+N
        IF(I.GT.3.AND.ABS(DIF).LT.MAX(AEPS,ABS(RINT)*REPS)) RETURN
5000  CONTINUE
 
C	Integral fails to converge
      IER=30
      END

C     ---------------------------------------------

      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE INTEGRAND (WITHOUT THE WEIGHT FUNCTION) FOR SUBROUTINE GAUCB1

      F=EXP(X)
      END
