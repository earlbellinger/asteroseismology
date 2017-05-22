C     PROGRAM FOR MINIMISATION IN ONE DIMENSION
C     USING CUBIC HERMITE INTERPOLATION

      PROGRAM MINMIS
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F

C     EXAMPLE 8.3

51    FORMAT('   INITIAL APPROXIMATIONS =',1P2D14.6)
52    FORMAT('  IER =',I4,3X,'MINIMISER =',1PD14.6,3X,
     1       'F =',D14.6,3X,4HF''=,D14.6)

      REPS=1.D-7
      AEPS=1.D-8

100   PRINT *,'TYPE X1,X2=INITIAL APPROXIMATIONS  (QUITS WHEN X1.EQ.X2)'
      READ *,X1,X2
      IF(X1.EQ.X2) STOP
      WRITE(6,51) X1,X2
      CALL DAVIDM(X1,X2,FX,D2F,REPS,AEPS,IER,F)
      WRITE(6,52) IER,X2,FX,D2F
      GO TO 100
      END

C     -------------------------------------------

C	To minimise a function in one dimension using Hermite cubic interpolation
C
C	X1,X2 : (input/output) Two starting values for iteration
C		These will be updated during execution
C		with X2 giving the best estimate for minimiser.
C	F2 : (output) The function value at X2, which should be the minimum
C	D2F : (output) Estimated value of second derivative at X2, to detect
C		if the extremum is a minimum or maximum. If D2F>0 then
C		X2 should be minimum, otherwise it will be maximum
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The Minimiser will be calculated with accuracy MAX(AEPS, REPS*ABS(X2))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=52 implies that iteration has converged to a maximum
C		IER=502 implies that initial values X1=X2 and no calculations are done
C		IER=524 implies that iteration cannot be continued 
C		IER=525 implies that iteration failed to converge to specified accuracy
C	F : (input) Name of the function routine to calculate the function
C		which is to be minimised
C
C	FUNCTION F(X,DX) to calculate the required function, must be supplied
C		by the user. Here DX is the first derivative of F(X) at X
C
C	Required routines : F

      SUBROUTINE DAVIDM(X1,X2,F2,D2F,REPS,AEPS,IER,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=75)

      IER=0
      IF(X1.EQ.X2) THEN
        IER=502
        RETURN
      ENDIF

      F1=F(X1,DF1)
      F2=F(X2,DF2)
      DX1=X2-X1

C	Loop for iteration
      DO 2000 I=1,NIT
        DF12=(F2-F1)/DX1
        R=2.*DF2+DF1-3.*DF12
        R1=(3.*DF12-DF2-DF1)**2-DF1*DF2
        IF(R1.GT.0.0) THEN
C	Perform Hermite cubic interpolation
          R1=SIGN(SQRT(R1),R)
          R=R+R1
          DX=-DF2*DX1/R
        ELSE
C	Perform parabolic interpolation
          R2=2.*(DF12-DF2)
          IF(R2.NE.0.0) THEN
            DX=DX1*DF2/R2
          ELSE
C	Both cubic and parabolic interpolation fail, hence quit
            IER=524
            RETURN
          ENDIF
        ENDIF

        X1=X2
        X2=X2+DX
        D2F=(DF2-DF1)/DX1
        DX1=DX
        F1=F2
        DF1=DF2
        F2=F(X2,DF2)

C	Convergence check
        IF(X1.EQ.X2.OR.(ABS(DX).LT.MAX(REPS*ABS(X2),AEPS).AND.I.GT.2))
     1     THEN
C	If it is maximum set error flag
          IF(D2F.LE.0) IER=52
          RETURN
        ENDIF
2000  CONTINUE

C	Iteration fails to converge
      IER=525
      END

C     ----------------------------------------------------

      FUNCTION F(X,DF)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE REQUIRED FUNCTION AND ITS DERIVATIVE

      F=(X*X-0.01D0)*EXP(-10.0*X)
      DF=(0.1D0+2.*X-10.0*X*X)*EXP(-10.0*X)
      END

