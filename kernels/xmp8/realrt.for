C     PROGRAM TO FIND A REAL ROOT OF A NONLINEAR EQUATION
C     USING BISECTION OR SECANT ITERATION OR BRENT'S METHOD

C	BISECT and bRENT may converge to a singularity if sign change
C        is across a singularity, while SECANT will not do so


      PROGRAM ROOT
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F

C     EXAMPLE 7.2

51    FORMAT('   IER =',I4,5X,'IT =',I2,5X,'X0 =',1PD14.6,5X,
     1       'ROOT =',D14.6)
52    FORMAT(' IER =',I4,2X,'IT =',I2,2X,'ROOT =',1PD14.6,2X,'XL =',
     1       D14.6,2X,'XU =',D14.6)
53    FORMAT('       NB =',I4,5X,'XL =',1PD14.6,5X,'XU =',D14.6)

      REPS=1.D-7
      AEPS=1.D-8

100   PRINT *,'XL=LOWER LIMIT,  XU=UPPER LIMIT    (QUITS WHEN XL.EQ.XU)'
      READ *,XL,XU
      IF(XL.EQ.XU) STOP
      PRINT *,'TYPE IT=1/2/3 FOR  BISECT/SECANT/BRENT'
      READ *,IT

      IF(IT.EQ.1) THEN
        PRINT *,'TYPE NB=NO. OF BISECTIONS REQUIRED'
        READ *,NB
        WRITE(6,53) NB,XL,XU
        CALL BISECT(X,XL,XU,NB,IER,F)
        WRITE(6,52) IER,IT,X,XL,XU
      ELSE IF(IT.EQ.2) THEN
        PRINT *,'TYPE X0=INITIAL APPROXIMATION'
        READ *,X0
        CALL SECANT(X0,XL,XU,X,REPS,AEPS,IER,F)
        WRITE(6,51) IER,IT,X0,X
      ELSE
        CALL BRENT(XL,XU,X,REPS,AEPS,IER,F)
        WRITE(6,52) IER,IT,X,XL,XU
      ENDIF
      GO TO 100
      END

C     ---------------------------------------------

C	Real zero of a continuous function using bisection
C
C	X : (output) Computed value of the zero using interpolation
C	XL : (input/output) Lower limit of interval containing the zero
C	XU : (input/output) Upper limit of interval containing the zero
C		These limits will be replaced by refined limits after execution
C	NB : (input) Number of bisections to be performed
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=-1 implies that function value is zero at some trial point
C		IER=401 implies NB.LE.0 in which case no calculations are done
C		IER=421 implies that function has same sign at both limits
C		and it is not possible to apply bisection
C	F : (input) Name of the function routine to calculate the function
C		FUNCTION F(X) must be supplied by the user.
C
C	Required routines : F

      SUBROUTINE BISECT(X,XL,XU,NB,IER,F)
      IMPLICIT REAL*8(A-H,O-Z)

      IF(NB.LE.0) THEN
        IER=401
        RETURN
      ENDIF

      FL=F(XL)
      FU=F(XU)
      IER=-1
C	If the function value is zero at either point then quit
      IF(FL.EQ.0.0) THEN
        X=XL
        RETURN
      ELSE IF(FU.EQ.0.0) THEN
        X=XU
        RETURN
      ENDIF

      IF(FL.GT.0.0.EQV.FU.GT.0.0) THEN
C	If the function has the same sign at both end points then quit
        IER=421
        RETURN
      ENDIF

C	Loop for bisection
      DO 1000 K=1,NB
        X=(XL+XU)/2.
        FX=F(X)
C	If function is zero then quit
        IF(FX.EQ.0.0) RETURN
        IF(FX.GT.0.0.EQV.FU.GT.0.0) THEN
          XU=X
          FU=FX
        ELSE
          XL=X
          FL=FX
        ENDIF
1000  CONTINUE

C	linear interpolation between XL and XU
      X=(XL*FU-XU*FL)/(FU-FL)
      IER=0
      END

C     -------------------------------------------

C	Real zero of a given function using Brent's method
C
C	A : (input/output) Lower limit of interval where zero is located
C	B : (input/output) Upper limit of interval where zero is located
C		The limits A, B are updated by the subroutine
C	X : (output) Computed value of the zero
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=427 implies that function has same sign at x=A and B
C			in which case no calculations are done
C		IER=428 implies that iteration failed to converge to specified
C			accuracy
C	F : (input) Name of the function routine to calculate the function
C		FUNCTION F(X) must be supplied by the user.
C
C	Required routines : F

      SUBROUTINE BRENT(A,B,X,REPS,AEPS,IER,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=75)

      IER=0
      FA=F(A)
      FB=F(B)
      IF(FA.GT.0.0.EQV.FB.GT.0.0) THEN
C	If the function has the same sign at both the end points, then quit
        IER=427
        RETURN
      ENDIF

      C=A
      FC=FA
      D=B-A
      E=D

C	The iteration loop
      DO 2000 I=1,NIT
        IF(ABS(FC).LT.ABS(FB)) THEN
C	If F(C) < F(B), then interchange B and C
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF

        EPS=MAX(REPS*ABS(B),AEPS)
        DX=0.5*(C-B)
        IF(ABS(DX).LT.EPS.OR.FB.EQ.0) THEN
C	Iteration has converged
          IER=0
          X=B
          RETURN
        ENDIF

        IF(ABS(E).LT.EPS.OR.ABS(FA).LE.ABS(FB)) THEN
C	If the previous step is too small or if F(A).LE.F(B), perform bisection
          D=DX
          E=DX
        ELSE
          R3=FB/FA
          IF(A.EQ.C) THEN
C	Try linear interpolation
            P=2.*DX*R3
            P1=1.-R3
          ELSE
C	Try inverse quadratic interpolation
            R1=FA/FC
            R2=FB/FC
            P=R3*(2.*DX*R1*(R1-R2)-(B-A)*(R2-1))
            P1=(R1-1.)*(R2-1.)*(R3-1.)
          ENDIF
          IF(P.GT.0) THEN
            P1=-P1
          ELSE
            P=-P
          ENDIF

          E=D
          IF(2.*P.LT.3.*DX*P1-ABS(EPS*P1).AND.P.LT.ABS(0.5*E*P1)) THEN
C	Accept the interpolated value
            D=P/P1
          ELSE
C	otherwise, perform bisection
            D=DX
            E=DX
          ENDIF
        ENDIF

        A=B
        FA=FB
        IF(ABS(D).GT.EPS) THEN
          B=B+D
        ELSE
C	If the change is too small, shift by EPS
          B=B+SIGN(EPS,DX)
        ENDIF
        FB=F(B)
        IF(FB.GT.0.0.EQV.FC.GT.0.0) THEN
C	If F(B) and F(C) have the same sign, then replace C by A
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
2000  CONTINUE

C	Iteration fails to converge
      IER=428
      X=B
      END

C     ----------------------------------------------------

C	Real zero of a given function using secant iteration
C
C	X0 : (input) Initial guess for the zero
C	XL : (input) Lower limit of interval where zero is expected
C	XU : (input) Upper limit of interval where zero is expected
C	X : (output) Computed value of the zero
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=40 implies that function value is equal at two points
C			and it is not possible to continue the iteration
C		IER=402 implies XL>X0 or XU<X0, in which case no calculations are done
C		IER=422 implies that iteration goes outside the specified limits
C		IER=423 implies that iteration failed to converge to specified accuracy
C	FUN : (input) Name of the function routine to calculate the function
C		FUNCTION FUN(X) must be supplied by the user.
C
C	Required routines : FUN

      SUBROUTINE SECANT(X0,XL,XU,X,REPS,AEPS,IER,FUN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIS=75)

      IER=0
      IF(XL.GT.X0.OR.XU.LT.X0) THEN
C	If X0 is outside the specified interval (XL,XU) then quit
        IER=402
        RETURN
      ENDIF

      X=X0
C	Select the increment for the next point X+DX
      DX=(XU-X0)/100.
      IF(ABS(DX).LT.5.*MAX(REPS*ABS(X),AEPS)) DX=(XL-X0)/5.
      IF(ABS(DX).GT.0.1D0*MAX(ABS(X),100.*AEPS))
     1    DX=SIGN(0.1D0*MAX(ABS(X),100.*AEPS),DX)
      F1=0.0

      DO 1000 L=1,NIS
        F=FUN(X)
        DX1=DX

        IF(F1-F.EQ.0.0) THEN
          IF(F.EQ.0.0) RETURN
C	If F1=F and F.NE.0, then quit
          IER=40
          RETURN
        ENDIF

C	The secant iteration
        IF(L.GT.1) DX=DX1*F/(F1-F)
        X=X+DX
        F1=F

        IF(ABS(DX).LT.MAX(REPS*ABS(X),AEPS).AND.L.GT.2) RETURN
        IF(X.LT.XL.OR.X.GT.XU) THEN
C	If iteration goes outside the specified limits (XL,XU), then quit
          IER=422
          RETURN
        ENDIF
1000  CONTINUE

C	The iteration fails to converge
      IER=423
      END

C     ----------------------------------------------------

      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)

C	The required function

      F=X-TAN(X)
      END
