C     PROGRAM TO FIND REAL ROOT OF A NONLINEAR EQUATION USING
C     NEWTON-RAPHSON METHOD WITH MODIFICATIONS FOR TREATING MULTIPLE ROOTS

      PROGRAM ROOT
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F

C     EXAMPLE 7.4

51    FORMAT('   IER =',I4,5X,'STARTING VALUE =',1PD14.6,5X,
     1       'ROOT =',D14.6)

      REPS=1.D-14
      AEPS=1.D-19

100   PRINT *,'TYPE X0=STARTING VALUE, XL=LOWER LIMIT, XU=UPPER LIMIT'
      PRINT *,'           (QUITS WHEN XL.EQ.XU)'
      READ *,X0,XL,XU
      IF(XL.EQ.XU) STOP
      CALL NEWRAP(X0,XL,XU,X,REPS,AEPS,IER,F)
      WRITE(6,51) IER,X0,X
      GO TO 100
      END

C     -------------------------------------------

C	Real zero of a given function using Newton-Raphson iteration
C
C	X0 : (input) Initial guess for the zero
C	XL : (input) Lower limit of interval where zero is expected
C	XU : (input) Upper limit of interval where zero is expected
C	X : (output) Computed value of the zero
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=-k implies that the zero is multiple and iteration
C			is modified to account for multiplicity k.
C		IER=126-100*k implies that the derivative is zero and the zero is
C			detected to be multiple. Iteration is terminated.
C		IER=403 implies XL>X0 or XU<X0, in which case no calculations are done
C		IER=424 implies that iteration goes outside the specified limits
C		IER=425 implies that iteration failed to converge to specified accuracy
C		IER=426 implies that the derivative is zero and iteration
C			is terminated
C	FUN : (input) Name of the function routine to calculate the function
C		FUNCTION FUN(X,DX) must be supplied by the user.
C		Here DX is the first derivative of FUN at X.
C
C	Required routines : FUN

      SUBROUTINE NEWRAP(X0,XL,XU,X,REPS,AEPS,IER,FUN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=75)

      IER=0
      IF(XL.GT.X0.OR.XU.LT.X0) THEN
C	If X0 is outside the specified limits then quit
        IER=403
        RETURN
      ENDIF

      MR=1
      X=X0
      DX=1
      DXR=1
      L0=1

      DO 1000 L=1,NIT
        F=FUN(X,DF)
        DX1=DX
        DXR1=DXR

        IF(DF.EQ.0.0) THEN
C	If the derivative is zero, then quit
          IF(F.EQ.0.0) RETURN
C	If the function is nonzero, then set the error flag
          IER=426
          IF(MR.GT.1) IER=126-100*IER
          RETURN
        ENDIF

C	Newton-Raphson iteration for root with multiplicity MR
        DX=-MR*F/DF
        DXR=DX/DX1
        X=X+DX

        IF(ABS(DX).LT.MAX(REPS*ABS(X),AEPS).AND.L.GT.2) RETURN
        IF(X.LT.XL.OR.X.GT.XU) THEN
C	If the iteration goes outside the specified limits, then quit
          IER=424
          RETURN
        ENDIF

        IF(L-L0.GT.3) THEN
C	If 3 iterations have been done with same multiplicity, then
C	get a new estimate for multiplicity MR
          MR1=MR
          MR2=MR
          IF(DXR.LT.0.99D0) MR1=MR/(1.-DXR)+0.5
          IF(DXR1.LT.0.99D0) MR2=MR/(1.-DXR1)+0.5
C	Accept the new value of MR only if both estimates match
          IF(MR1.EQ.MR2.AND.MR1.NE.MR) THEN
            IER=-MR1
            L0=L
            MR=MR1
            IF(MR.EQ.1) IER=0
          ENDIF
        ENDIF
1000  CONTINUE

C	Iteration fails to converge
      IER=425
      END

C     ----------------------------------------------------

      FUNCTION F(X,DF)
      IMPLICIT REAL*8(A-H,O-Z)

C	The required function and its derivative

      F=X-TAN(X)
      DF=1.-1./COS(X)**2
      END
