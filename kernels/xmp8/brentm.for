C     PROGRAM FOR MINIMISATION IN ONE DIMENSION USING BRENT'S METHOD

      PROGRAM MINMIZ
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F

C     EXAMPLE 8.2

51    FORMAT('  INITIAL BRACKETING TRIPLET :',1P3D14.6)
52    FORMAT('    FINAL BRACKETING TRIPLET :',1P3D14.6/
     1     5X,'IER =',I4,5X,'MINIMISER =',D14.6,5X,'MINIMUM =',D14.6)

      REPS=1.D-7
      AEPS=1.D-9

100   PRINT *,'TYPE XL,XU,X=BRACKETING TRIPLET WITH X BETWEEN XL AND XU'
      PRINT *,'           (QUITS WHEN XL.EQ.XU)'
      READ *,XL,XU,X
      IF(XL.EQ.XU) STOP
      WRITE(6,51) XL,X,XU
      CALL BRENTM(XL,XU,X,FX,REPS,AEPS,IER,F)
      WRITE(6,52) XL,X,XU,IER,X,FX
      GO TO 100
      END

C     -------------------------------------------

C	To minimise a function in one dimension using Brent's method
C
C	A,B,X : (input/output) Triplet which brackets the minimum.
C		After execution these will contain the final triplet
C		with X giving the best estimate for minimiser.
C		X must be between A and B; and F(X)<MIN(F(A),F(B))
C	FX : (output) The function value at X.
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The bracketing interval is subdivided until
C		ABS(B-A) < MAX(AEPS, REPS*ABS(X))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=51 implies that subroutine failed to reduce the
C			bracketing interval to required level
C		IER=523 implies that initial values of A,X,B do not bracket
C			the minimum
C	F : (input) Name of the function routine to calculate the function
C		which is to be minimised
C
C	FUNCTION F(X) to calculate the required function, must be supplied
C		by the user.
C
C	Required routines : F

      SUBROUTINE BRENTM(A,B,X,FX,REPS,AEPS,IER,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=75,GR=0.381966D0)

      IER=0
      IF(A.GT.B) THEN
C	Interchange A and B
        T=A
        A=B
        B=T
      ENDIF
      FA=F(A)
      FB=F(B)
      FX=F(X)
      IF(FA.LT.FX.OR.FB.LT.FX.OR.X.LE.A.OR.X.GE.B) THEN
        IER=523
        RETURN
      ENDIF

      E=0
      V=X
      W=X
      FV=FX
      FW=FX

C	Loop for iteration
      DO 2000 I=1,NIT
        XM=0.5*(A+B)
        EPS2=MAX(REPS*ABS(X),AEPS)
        EPS=0.5*EPS2

C	The convergence test
        IF(ABS(X-XM).LT.EPS2-0.5*(B-A)) THEN
          IER=0
          RETURN
        ENDIF

        P=0
        T=0
        R=0
        IF(ABS(E).GT.EPS) THEN
C	Parabolic interpolation
          R=(X-W)*(FX-FV)
          T=(X-V)*(FX-FW)
          P=(X-V)*T-(X-W)*R
          T=2*(T-R)
          IF(T.GT.0) P=-P
          T=ABS(T)
          R=E
          E=D
        ENDIF

        IF(ABS(P).LT.ABS(.5*T*R).AND.P.GT.T*(A-X).AND.P.LT.T*(B-X)) THEN
C	accept the interpolated point
          D=P/T
          U=X+D
          IF(U-A.LT.EPS2.OR.B-U.LT.EPS2) THEN
C	If it is too close to end points shift it by EPS at least
            D=EPS
            IF(X.GE.XM) D=-EPS
          ENDIF
        ELSE
C	Perform golden section search
          E=B-X
          IF(X.GE.XM) E=A-X
          D=GR*E
        ENDIF
        IF(ABS(D).GE.EPS) THEN
          U=X+D
        ELSE
C	Shift the point by at least EPS
          U=X+SIGN(EPS,D)
        ENDIF
        FU=F(U)

C	Updating the bracketing triplet
        IF(FU.LE.FX) THEN
          IF(U.LT.X) THEN
C	(A, U, X) is the triplet
            B=X
          ELSE
C	(X, U, B) is the triplet
            A=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
C	(U, X, B) is the triplet
            A=U
          ELSE
C	(A, X, U) is the triplet
            B=U
          ENDIF
          IF(FU.LE.FW.OR.W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF

2000  CONTINUE

C	Iteration fails to converge
      IER=51
      END

C     ----------------------------------------------------

      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE FUNCTION TO BE MINIMISED

      F=(X*X-0.01D0)*EXP(-10.0*X)
      END
