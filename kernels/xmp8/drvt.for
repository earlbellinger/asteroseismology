C     PROGRAM FOR DIFFERENTIATION USING H --> 0 EXTRAPOLATION

      PROGRAM DIFF
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F

C     EXAMPLE 5.2: DERIVATIVES OF EXP(X) 

51    FORMAT('   IER =',I4,I6,'TH DERIVATIVE AT X =',1PD14.6,
     1  3X,1PD14.6/ '   EXACT VALUE =',D14.6,5X,'INITIAL SPACING ='
     2  ,D14.6) 
 
      REPS=1.D-5
      AEPS=1.D-6

100   PRINT *,'TYPE  ID=ORDER OF DERIVATIVE,    (QUITS WHEN ID.LE.0)'
      READ *,ID
      IF(ID.LE.0) STOP
      PRINT *,'A=REQUIRED POINT,   HH=INITIAL SPACING'
      READ *,A,HH
      DF=DRVT(A,ID,HH,REPS,AEPS,F,IER)
      WRITE(6,51) IER,ID,A,DF,EXP(A),HH
      GO TO 100
      END
 
C     --------------------------------------------------
 
C	To calculate derivative of a function
C
C	A : (input) Value of x at which the derivative needs to be calculated
C	ID : (input) Order of derivative required, ID may be 1,2,3 or 4
C	HH0 : (input) Initial value of step length to be used.
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy.
C		Calculations will continue till estimated error is less
C		than MAX(AEPS, REPS*ABS(DRVT))
C	F : (input) The name of function routine to calculate the given function
C		FUNCTION F(X) must be supplied by the user.
C	IER : (output) Error parameter, IER=0 for successful execution
C		IER=28 implies required accuracy is not achieved
C		IER=29 implies roundoff error appears to be dominating
C		IER=208 implies ID<1 or ID>4, no calculations are done.
C
C	DRVT is the calculated derivative at x=A
C		
C	FUNCTION F(X) must be supplied by the user
C
C	Required routines : F

      FUNCTION DRVT(A,ID,HH0,REPS,AEPS,F,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=15,NCMAX=6,HFAC=1.5)
      DIMENSION H(NMAX),T(NMAX,NCMAX)

C	Set spacing to the specified initial value
      HH=HH0
      ERR=0.0
      IER=0

      DO 2000 I=1,NMAX
        H(I)=HH

        IF(ID.EQ.1) THEN
          T(I,1)=(F(A+HH)-F(A-HH))/(2.*HH)
        ELSE IF(ID.EQ.2) THEN
          T(I,1)=(F(A-HH)-2.*F(A)+F(A+HH))/(HH*HH)
        ELSE IF(ID.EQ.3) THEN
          T(I,1)=(F(A+2.*HH)-F(A-2.*HH)-2.*(F(A+HH)-F(A-HH)))/(2.*HH**3)
        ELSE IF(ID.EQ.4) THEN
          T(I,1)=(F(A+2*HH)+F(A-2*HH)-4*(F(A+HH)+F(A-HH))+6*F(A))/HH**4
        ELSE
          IER=208
          RETURN
        ENDIF

C	Reduce spacing by a factor of HFAC
        HH=HH/HFAC
        IF(I.GT.1) THEN
C	Calculate at most NCMAX columns of the T-table
          JU=MIN(I,NCMAX)
          DO 1200 J=2,JU
           T(I,J)=T(I,J-1)+(T(I,J-1)-T(I-1,J-1))/((H(I-J+1)/H(I))**2-1.)
1200      CONTINUE

          DRVT=T(I,JU)
          ERR1=ERR
          ERR=MIN(ABS(DRVT-T(I,JU-1)),ABS(DRVT-T(I-1,MIN(I-1,NCMAX))))
          IF(ERR.LT.MAX(ABS(DRVT)*REPS,AEPS)) RETURN

          IF(I.GT.5.AND.ERR.GT.2.*ERR1) THEN
C	Roundoff error dominates
            IER=29
            DRVT=T(I-1,MIN(I-1,NCMAX))
            RETURN
          ENDIF
        ENDIF

2000  CONTINUE
      IER=28
      END
 
C     ------------------------------------------------------
 
      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE GIVEN FUNCTION

      F=EXP(X)
      END
