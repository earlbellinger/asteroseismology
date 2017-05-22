C     PROGRAM FOR MINIMISATION IN ONE DIMENSION 
C     IT FIRST FINDS A BRACKETING TRIPLET TO LOCATE A MINIMUM AND
C     THEN USES GOLDEN SECTION SEARCH TO FIND THE ACCURATE VALUE

      PROGRAM MINMIZ
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F

C     EXAMPLE 8.1

51    FORMAT('   INITIAL CHOICE FOR BRACKETING :',1P2D14.6)
52    FORMAT('   IER =',I4,5X,'BRACKETING TRIPLET:',1P3D14.6)
53    FORMAT('   MINIMISER =',1PD14.6,5X,'MINIMUM =',D14.6)

      REPS=1.D-6
      AEPS=1.D-7

100   PRINT *,'TYPE A,B=INITIAL PTS FOR BRACKETING  (QUITS WHEN A.EQ.B)'
      READ *,A,B
      IF(A.EQ.B) STOP
      WRITE(6,51) A,B
      CALL BRACKM(A,B,X,F,IER)
      WRITE(6,52) IER,A,X,B

C     USE THE BRACKETING TRIPLET A,X,B FOR GOLDEN SECTION SEARCH

      CALL GOLDEN(A,B,X,FX,REPS,AEPS,IER,F)
      WRITE(6,52) IER,A,X,B
      WRITE(6,53)  X,FX
      GO TO 100
      END

C     ------------------------------------------------------------

C	To bracket a minimum in one dimension
C
C	A,B : (input/output) Starting values for the search. After
C		execution they will contain the limits of bracketing triplet
C	X : (output) Real variable which will contain the point X between
C		A and B, such that F(X) < MIN(F(A),F(B))
C	F : (input) Name of the function routine to calculate the function
C		which is to be minimised
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=501 implies that A=B, in which case no calculations are done
C		IER=521 implies that subroutine failed to find a bracketing triplet
C
C	FUNCTION F(X) to calculate the required function, must be supplied
C		by the user.
C
C	Required routines : F

      SUBROUTINE BRACKM(A,B,X,F,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(GR=1.618034D0,RMAX=10.,NTRY=200)

      IF(A.EQ.B) THEN
        IER=501
        RETURN
      ENDIF
      IER=0

      FA=F(A)
      FB=F(B)
      IF(FB.GT.FA) THEN
C	Exchange A and B
        T=A
        A=B
        B=T
        T=FA
        FA=FB
        FB=T
      ENDIF
C	First guess for X
      X=B+GR*(B-A)
      FX=F(X)

      DO 2000 I=1,NTRY
        IF(FB.LT.FX) THEN
C	The minimum is bracketed by the triplet (A,B,X), exchange B and X
          T=B
          B=X
          X=T
          RETURN
        ENDIF

C	Next guess by parabolic interpolation
        R=(B-A)*(FB-FX)
        T=(B-X)*(FB-FA)
        P=(B-X)*T-(B-A)*R
        T=2.*(T-R)
        IF(T.GT.0) P=-P
        T=ABS(T)

        IF(T.GT.0.AND.(P.GT.0.EQV.(X-B).GT.0)) THEN
          IF(ABS(P).LT.ABS(T*(X-B))) THEN
C	Interpolated point is between B and X
            XP=B+P/T
            FP=F(XP)
            IF(FP.LT.FX) THEN
C	Minimum is bracketed by (B,XP,X)
              A=B
              B=X
              X=XP
              RETURN
            ELSE IF(FP.GT.FB) THEN
C	Minimum is bracketed by (A,B,XP)
              X=B
              B=XP
              RETURN
            ENDIF
C	If the minimum is not bracketed, then reject the interpolated point
            XP=X+GR*(X-B)
            FP=F(XP)

          ELSE IF(ABS(P).LT.ABS(T*(X-B)*RMAX)) THEN
C	Interpolated point is between X and the allowed limit
            XP=B+P/T
            FP=F(XP)
            IF(FP.LT.FX) THEN
C	If the minimum is not bracketed, then reject the interpolated point
              XP=X+GR*(X-B)
              FP=F(XP)
            ENDIF
          ELSE
C	Reject the interpolated point if it is beyond the allowed limit
            XP=B+RMAX*(X-B)
            FP=F(XP)
          ENDIF
        ELSE
C	Reject the interpolated point if it is on the other side of B
          XP=X+GR*(X-B)
          FP=F(XP)
        ENDIF

        A=B
        B=X
        X=XP
        FA=FB
        FB=FX
        FX=FP
2000  CONTINUE

C	Fails to bracket the minimum
      IER=521
      END

C     ----------------------------------------------------

C	To minimise a function in one dimension using Golden section search
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
C		IER=50 implies that function values at A, X, B are not
C			distinct, most likely due to roundoff errors
C		IER=51 implies that subroutine failed to reduce the
C			bracketing interval to required level
C		IER=522 implies that initial values of A,X,B do not bracket the minimum
C	F : (input) Name of the function routine to calculate the function
C		which is to be minimised
C
C	FUNCTION F(X) to calculate the required function, must be supplied
C		by the user.
C
C	Required routines : F

      SUBROUTINE GOLDEN(A,B,X,FX,REPS,AEPS,IER,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(GR=0.381966D0,NIT=100)

      FA=F(A)
      FB=F(B)
      FX=F(X)
      IF(FX.GT.MIN(FA,FB).OR.(X-A.GT.0.EQV.X-B.GT.0)) THEN
C	(A, X, B) does not bracket a minimum
        IER=522
        RETURN
      ENDIF
      IER=0

      DO 2000 I=1,NIT
C	Choose the next point in the larger section
        IF(ABS(X-A).GT.ABS(X-B)) THEN
          XP=X+GR*(A-X)
          FP=F(XP)
          IF(FP.LT.FX) THEN
C	(A, XP, X) bracket the minimum
            B=X
            X=XP
            FB=FX
            FX=FP
          ELSE
C	(XP, X, B) bracket the minimum
            A=XP
            FA=FP
          ENDIF
        ELSE
          XP=X+GR*(B-X)
          FP=F(XP)
          IF(FP.LT.FX) THEN
C	(X, XP, B) bracket the minimum
            A=X
            X=XP
            FA=FX
            FX=FP
          ELSE
C	(A, X, XP) bracket the minimum
            B=XP
            FB=FP
          ENDIF
        ENDIF

        IF(ABS(B-A).LT.MAX(ABS(X)*REPS,AEPS)) RETURN
        IF(FX.GE.FA.AND.FX.GE.FB) THEN
C	Roundoff errors may be dominating
          IER=50
          RETURN
        ENDIF
2000  CONTINUE

C	Fails to converge to required accuracy
      IER=51
      END

C     --------------------------------------------------

      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE FUNCTION TO BE MINIMISED

      F=(X*X-0.01D0)*EXP(-10.0*X)
      END
