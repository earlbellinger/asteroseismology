C     PROGRAM TO OBTAIN RATIONAL FUNCTION MINIMAX APPROXIMATION FOR
C     MATHEMATICAL FUNCTIONS USING REMES ALGORITHM

      PROGRAM MINMAX
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(20),WK(400),IWK(20),F(200),E(25),X(200)

C     EXAMPLE 9.10: RELATIVE MINIMAX APPROXIMATION TO ARC TAN(X)

51    FORMAT('   IER =',I4,5X,'MAX. ERROR =',1PD14.6/
     1       '   COEF. IN NUMERATOR',4D16.8/(2X,5D16.8))
52    FORMAT('   COEF. IN DENOMINATOR',1P4D16.8/(2X,5D16.8))
53    FORMAT('  EXTREMA IN ERROR CURVE',1P4D14.6/(2X,5D14.6))
54    FORMAT('   M =',I3,5X,'K =',I3,5X,'N =',I3,5X,'IFLG =',I3)
55    FORMAT('  INITIAL GUESS FOR COEFFICIENTS :',1P2D14.6/
     1       (2X,1P5D14.6))
56    FORMAT('  INITIAL GUESS FOR EXTREMA OF ERROR CURVE :',1P2D14.6/
     1       (2X,1P5D14.6))

100   PRINT *,'TYPE M=DEGREE OF NUMERATOR,  K=DEGREE OF DENOMINATOR'
      PRINT *,'     (QUITS WHEN (M+K).LE.0)'
      READ *,M,K
      IF(M+K.LE.0) STOP
      PRINT *,'TYPE N=NO. OF PTS FOR INITIAL SEARCH OF EXTREMA'
      READ *,N
      PRINT *,'TYPE XL=LOWER LIMIT,  XU=UPPER LIMIT, IFLG=0/1/2'
      READ *,XL,XU,IFLG
      MM=M
      KK=K
      IF(XL.EQ.XU) STOP

      EPS=1.D-8
      EPSM=1.D-5
      IE=M+K+2
      WRITE(6,54) M,K,N,IFLG
      IF(IFLG.EQ.0) THEN
        PRINT *,'TYPE INITIAL GUESS FOR COEF'
        READ *,(A(I),I=1,M+K+1)
        WRITE(6,55)(A(I),I=1,M+K+1)
      ELSE IF(IFLG.EQ.2) THEN
        PRINT *,'TYPE INITIAL GUESS FOR EXTREMA'
        READ *,(E(I),I=1,M+K+2)
        WRITE(6,56) (E(I),I=1,M+K+2)
      ENDIF

      CALL REMES(M,K,N,XL,XU,A,X,F,E,IE,EMAX,EPS,EPSM,IFLG,IER,WK,IWK)
      WRITE(6,51) IER,EMAX,(A(I),I=K+1,M+K+1)
      WRITE(6,52) (A(I),I=1,K)
      WRITE(6,53) (E(I),I=1,IE)
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

C     --------------------------------------------------

C	To calculate the error in rational function approximation
C	For use with subroutine REMES. It is called by BRENTM to find
C	extrema of error curve. It may not be used for any other purpose.
C
C	X : (input) the value at which error is to be calculated
C	FM = (FUN(X) - FUND(X)*R_mk(X))*SI  is the calculated error
C	The parameters for rational function are passed through common block
C	M,K are degree of numerator and denominator
C	A is a real array containing the coefficient of rational function
C		approximation
C	SI is -1 for maximum and +1 for minimum. The error is multiplied
C		by SI to make it a minimum.
C		For initial scan SI>10 and function is not evaluated.
C
C	Functions FUN(x) and FUND(x) must be provided by the user to
C		seek approximation of form FUN(X) = FUND(X)*RMK(X)
C
C	Required routines : FUN, FUND
C
      FUNCTION FM(X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=50)
      COMMON/ZZFN/A(NMAX),SI,M,K

      NK=M+K+2
      FM=0.0
      IF(SI.LT.10) FM=FUN(X)

C	Calculate the numerator using nested multiplication
      FN=A(K+M+1)
      DO 2200 J=1,M
2200  FN=FN*X+A(NK-J-1)

      IF(K.GT.0) THEN
C	Calculate the denominator using nested multiplication
        FD=A(K)
        DO 2300 J=1,K-1
2300    FD=FD*X+A(K-J)
        FD=FD*X+1
      ELSE
        FD=1.
      ENDIF

      FM=(FM-FUND(X)*FN/FD)*SIGN(1.D0,SI)
      END

C     -------------------------------------------

C	Solution of a system of linear equations using Gaussian elimination
C	with partial pivoting
C
C	N : (input) Number of equations to be solved
C	NUM : (input) Number of different sets (each with N equations) of
C	        equations to be solved
C	A : (input/output) The matrix of coefficient of size LJ*N
C	        A(I,J) is the coefficient of x_J in Ith equation
C	     	at output it will contain the triangular decomposition
C	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
C	        X(I,J) is the Ith element of Jth right hand side
C	     	at output it will contain the solutions
C	DET : (output) The determinant of the matrix
C	INC : (output) Integer array of length N containing information about
C		interchanges performed during elimination
C	LJ : (input) First dimension of arrays A and X in calling program
C	IER : (output) Error flag, IER=0 signifies successful execution
C		IER=101 implies (N.LE.0 or N.GT.LJ) 
C		IER=121 implies some pivot turned out to be zero and hence
C			matrix must be nearly singular
C	IFLG : (input) Integer parameter to specify the type of computation required
C		If IFLG.LE.0, both elimination and solution are
C			done and IFLG is set to 2
C		If IFLG=1, only elimination is done and IFLG is set to 2
C		If IFLG.GE.2 only solution is calculated, the triangular
C			decomposition should have been calculated earlier
C
C	Required routines : None

      SUBROUTINE GAUELM(N,NUM,A,X,DET,INC,LJ,IER,IFLG)
      IMPLICIT REAL*8(A-H,O-Z)
C	For complex matrices use the following statements instead
C      IMPLICIT REAL*8(R)
C      IMPLICIT COMPLEX*16(A-H,S-Z)

      DIMENSION A(LJ,N),INC(N),X(LJ,NUM)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=101
        RETURN
      ENDIF

      IER=121
      IF(IFLG.LE.1) THEN
C	Perform elimination

        DET=1.0
        DO 2600 K=1,N-1
C	Find the maximum element in the Kth column
          R1=0.0
          KM=K
          DO 2200 L=K,N
            IF(ABS(A(L,K)).GT.R1) THEN
              R1=ABS(A(L,K))
              KM=L
            ENDIF
2200      CONTINUE

          INC(K)=KM
          IF(KM.NE.K) THEN
C	Interchange the rows if needed
            DO 2300 L=K,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
2300        A(KM,L)=T1
            DET=-DET
          ENDIF

          DET=DET*A(K,K)
          IF(A(K,K).EQ.0.0) RETURN
C	To check for singular or nearly singular matrices replace this
C	statement by, where REPS is approximately \hcross*Max(A(I,J))
C         IF(ABS(A(K,K)).LT.REPS) RETURN
          DO 2500 L=K+1,N
            A(L,K)=A(L,K)/A(K,K)
            DO 2500 L1=K+1,N
2500      A(L,L1)=A(L,L1)-A(L,K)*A(K,L1)
2600    CONTINUE
        DET=DET*A(N,N)
        INC(N)=N
C	If pivot is zero then return, IER has been set to 121
        IF(A(N,N).EQ.0.0) RETURN
C	To check for singular or nearly singular matrices replace this
C	statement by, where REPS is approximately \hcross*Max(A(I,J))
C         IF(ABS(A(N,N)).LT.REPS) RETURN

        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
C	Solution for the NUM different right-hand sides
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
C	Forward substitution
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF
          DO 3000 L=K+1,N
3000    X(L,J)=X(L,J)-A(L,K)*X(K,J)

C	back-substitution
        X(N,J)=X(N,J)/A(N,N)
        DO 3300 K=N-1,1,-1
          DO 3200 L=N,K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L)
3300    X(K,J)=X(K,J)/A(K,K)
5000  CONTINUE
      END

C     ----------------------------------

C	To calculate coefficients of Rational function minimax approximation
C	for a given function over a finite interval using Remes algorithm
C
C	M : (input) Required degree of polynomial in the numerator
C	K : (input) Required degree of polynomial in the denominator
C	N : (input/output) Number of points to be used for scanning the
C		extrema of error curves. If N<3(M+K+2) it will be set to
C		this value.
C	XL : (input) Lower limit of the interval over which approximation is required
C	XU : (input) Upper limit of the interval over which approximation is required
C	A : (input/output) Real array of length M+K+1 containing the coefficients
C		of approximation. A(I) is the coefficient of x**I in
C		the denominator, the constant term being 1. A(K+J+1) is
C		the coefficient of x**J in the numerator. If IFLG=0 the
C		initial guess to coefficients must be supplied.
C	X : (output) Real array of length N containing the points at which
C		function value is calculated for scanning the extrema in
C		error curve
C	F : (output) Real array of length N containing the function values at X(I)
C	EX : (input/output) Real array of length M+K+5 containing the
C		extrema of error curve. If IFLG=2, then the initial guess
C		for these extrema must be supplied
C	IE : (input/output) Number of extrema in error curve. This value
C		must be supplied if IFLG=2
C	EMAX : (output) Maximum error in the calculated approximation
C	EPS : (input) Required accuracy, the iteration for calculating the
C		coefficients of approximations is continued until
C		difference in maximum error is less than EPS
C	EPSM : (input) Required accuracy to which extrema in error curve
C		are determined.
C	IFLG : (input) Integer parameter specifying the nature of initial
C		approximation that is supplied. 
C		If IFLG=0 then iteration is started from initial guess for
C			coefficients supplied in array A
C		If IFLG=1 then no initial guess is required. Iteration is
C			started by assuming that the extrema of error curve
C			coincides with those of T_{M+K+1}(x).
C		If IFLG=2 then iteration is started from initial guess for
C			position of extrema in error curve. These must be
C			supplied in array EX and IE should be set to the number
C			of extrema
C	IER : (output) error parameter, IER=0 implies successful execution
C		IER=614 implies that M<0, K<0, XU.LE.XL or M+K+2>NMAX (=50)
C			in this case no calculations are done
C		IER=632 implies that the Remes iteration did not converge
C			to the specified accuracy
C		IER=633 implies that at some stage the error curve does not
C			have the required number of extrema.
C		Other values of IER may be set by GAUELM or BRENTM
C	WK : Real array of length (K+M+2)**2 used as scratch space
C	IWK : Integer array of length K+M+2 used as scratch space
C
C	Functions FUN(X) and FUND(X) must be supplied by the user.
C	Here X, FUN and FUND are real variables and we seek approximation
C	of the form FUN(x) ~ FUND(X)*R_mk(x)
C	To obtain minimax approximation with respect to absolute error
C	set FUND(X)=1 and FUN(X) to required function
C	To obtain minimax approximation with respect to relative error
C	set FUN(X)=1 and FUND(X) to reciprocal of required function.
C
C	Required routines : GAUELM, BRENTM, FM, FUN, FUND
C
      SUBROUTINE REMES(M,K,N,XL,XU,A,X,F,EX,IE,EMAX,EPS,EPSM,IFLG,IER,
     1                 WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FM
      PARAMETER(NIT=30,NJT=10,NMAX=50,PI=3.14159265358979324D0)
      COMMON/ZZFN/AA(NMAX),SI,MM,KK
      DIMENSION A(M+K+2),F(N),IWK(K+M+2),WK(K+M+2,*),EX(*),X(N)

      NK=M+K+2
      IF(M.LT.0.OR.K.LT.0.OR.NK.GT.NMAX.OR.XU.LE.XL) THEN
        IER=614
        RETURN
      ENDIF

C	Copy the value in the common block for use by function FM
      MM=M
      KK=K
      IF(N.LT.3*NK) N=3*NK

C	Generating the mesh points for crude scan of extrema
      H=(XU-XL)/(N-1.)
      DO 2000 I=1,N
        X(I)=XL+(I-1)*H
        F(I)=FUN(X(I))
2000  CONTINUE

      IER=0
      NUM=1
      LJ=NK
      REPS=0.0

C	Loop for Remes iteration
      DO 5000 IT=1,NIT
        E1=0.0
C	Copy the coefficients to common block for function FM
        DO 2200 I=1,M+K+1
2200    AA(I)=A(I)

        IF(IT.EQ.1.AND.IFLG.GE.1) THEN
          IF(IFLG.EQ.1) THEN
C	Use the extrema of Chebyshev polynomial T_{M+K+1} as initial guess
            EX(1)=X(1)
            IE=1
            DO 2300 I=2,NK-1
2300        EX(I)=0.5*(XL+XU)+0.5*(XU-XL)*COS((NK-I)*PI/(NK-1.))
            EX(NK)=XU
          ENDIF
          EI=0.0
          EMAX=0.0
          EMIN=0.0
        ELSE

C	Locate the extrema of the error curve
          DO 2400 I=1,N
C	Flag to avoid calculating the function repeatedly
            SI=22.0
            EI=F(I)+FM(X(I))
            IF(I.GT.2) THEN
              IF(E1.GE.E2.EQV.E1.GE.EI) THEN
C	An extrema is bracketed
                SI=1.
C	To convert maxima to minima
                IF(E1.GT.E2.OR.E1.GT.EI) SI=-1.
                AI=X(I-2)
                B=X(I)
                XI=X(I-1)
                CALL BRENTM(AI,B,XI,FX,REPS,EPSM,IER,FM)
                IF(IER.GT.0) RETURN

                IE=IE+1
                EX(IE)=XI
                WK(IE,1)=FX*SI
                IF(ABS(FX).GT.EMAX) EMAX=ABS(FX)
                IF(ABS(FX).LT.EMIN) EMIN=ABS(FX)
C	To ensure that dimensions of EX are not exceeded
                IF(IE.GT.M+K+4) GO TO 2500
              ENDIF
            ELSE IF(I.EQ.1) THEN
C	The end point is always included in the list of extrema
              IE=1
              EX(1)=X(1)
              EMAX=ABS(EI)
              EMIN=EMAX
              WK(I,1)=EI
              E0=EI
            ENDIF
            E2=E1
            E1=EI
2400      CONTINUE

          IE=IE+1
C	The end point is always included in the list of extrema
          EX(IE)=X(N)
          WK(IE,1)=EI
          EMAX=MAX(EMAX,ABS(EI))
          EMIN=MIN(EMIN,ABS(EI))

          IF(IE.LT.NK) THEN
C	If the number of extrema is less than required then quit
            IER=633
            RETURN
          ENDIF
2500      IF(IE.GT.NK) THEN

C	Remove one extrema from the list
            IE=IE-1
            IF(ABS(WK(IE+1,1)).GT.ABS(WK(1,1))) THEN
C	remove the first extrema
              EMIN=ABS(WK(2,1))
              E0=WK(2,1)
              DO 2600 I=1,IE
                EX(I)=EX(I+1)
                IF(ABS(WK(I+1,1)).LT.EMIN) EMIN=ABS(WK(I+1,1))
2600          WK(I,1)=WK(I+1,1)
            ELSE
C	remove the last extrema
              EMIN=ABS(WK(1,1))
              DO 2700 I=2,IE
                IF(ABS(WK(I,1)).LT.EMIN) EMIN=ABS(WK(I,1))
2700          CONTINUE
            ENDIF
C	Repeat until the number of extrema = NK
            GO TO 2500
          ENDIF

C	Convergence check, quit if difference between various extrema is
C	less than 1%
          IF(EMAX-EMIN.LT.1.D-2*EMAX) RETURN
          EI=MIN(0.5*(EMAX+EMIN),2.*EMIN)
          IF(E0.LT.0.0) EI=-EI
        ENDIF

C	Iteration to determine the coefficients
        DO 4000 JT=1,NJT
          AI=1
C	Setting up the equation matrix
          DO 3600 I=1,NK
            FI=FUN(EX(I))
            FD=FUND(EX(I))
            DO 3000 J=1,K
3000        WK(I,J)=-(FI-AI*EI)*EX(I)**J
            WK(I,K+1)=FD
            DO 3200 J=1,M
3200        WK(I,J+K+1)=FD*EX(I)**J
            WK(I,NK)=AI
            A(I)=FI
            AI=-AI
3600      CONTINUE

          IFL=0
          CALL GAUELM(NK,NUM,WK,A,DET,IWK,LJ,IER,IFL)
          IF(IER.GT.0) RETURN
          DIF=ABS(A(NK)-EI)
          EI=A(NK)
C	convergence check
          IF(DIF.LT.EPS.OR.DIF.LT.0.3D0*(EMAX-EMIN)) GO TO 5000
          IF(IT.EQ.1.AND.IFLG.GE.1) EMAX=0.3D0*ABS(A(NK))
4000    CONTINUE

C	Even if iteration for calculating coefficients fails continue further

5000  CONTINUE

C	Remes iteration fails to converge
      IER=632
      END

C     ------------------------------------

      FUNCTION FUND(X)
      IMPLICIT REAL*8(A-H,O-Z)
C     SET FUND=1/F(X) FOR RELATIVE MINIMAX APPROXIMATION

      IF(X.EQ.0.0) THEN
C     TAKE THE LIMITING VALUE
      FUND=1.
      ELSE
      X1=SQRT(X)
      FUND=X1/ATAN(X1)
      ENDIF
      END

C     ------------------------------------

      FUNCTION FUN(X)
      IMPLICIT REAL*8(A-H,O-Z)

C     SET FUN=1 FOR RELATIVE MINIMAX APPROXIMATION 
      FUN=1.0
      END
