C     PROGRAM TO CALCULATE COEFFICIENTS OF NONLINEAR LEAST SQUARES
C     APPROXIMATION USING QUASI-NEWTON METHOD WITH BFGS FORMULA
C     OR THE DIRECTION SET METHOD

      PROGRAM NLFIT
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL NLLSQ,NLLSQ_F
      DIMENSION H(400),G(20),WK(6000),B(20)
      COMMON/ZZFUN/FX(100),X(100),EF(100),FX1(100),N

C     EXAMPLE 9.4
C	With most initial guess the iteration does not converge to
C	specified accuracy, but the final value may still be acceptable
C	The sequence of iterates depends on roundoff error and it is
C	unlikely that the output given in the file will be reproduced.

C     FUNCTION WITH RANDOM ERROR TO GENERATE INPUT DATA

      F(Y)=EXP(-Y)+2.*EXP(-2.*Y)+3.0*EXP(-3.*Y)+RANGAU(IS)*1.D-6

51    FORMAT('   IER =',I4,5X,'NO. OF DATA PTS =',I4,5X,
     1     'NO. OF PARAMETERS =',I3/5X,'NO. OF FUNC. EVALUATIONS =',I5,
     2     5X,'SQUARED DIFF =',1PD14.6,5X,'IT =',I2/'   COEF =',
     3     5D14.6/(4X,5D14.6))
52    FORMAT('  INITIAL GUESS FOR COEF =',1P4D13.5/(2X,5D13.5))

      REPS=1.D-7
      AEPS=1.D-8

C     THE NO. OF PARAMETERS M SHOULD BE EVEN

100   PRINT *, 'TYPE N=NO. OF DATA PTS, M=NO. OF PARAMETERS'
      PRINT *,'           (QUITS WHEN N.LE.1)'
      READ *,N,M
      IF(N.LE.1) STOP
      IS=2

C     GENERATING INPUT DATA WITH RANDOM ERROR

      HH=1./(N-1.)
      DO 1000 I=1,N
        X(I)=(I-1)*HH
        FX(I)=F(X(I))
        EF(I)=1.D-2
1000  CONTINUE

      PRINT *,'TYPE INITIAL GUESS FOR PARAMETERS'
      READ *,(B(I),I=1,M)
      WRITE(6,52) (B(I),I=1,M)
      PRINT *,'TYPE IT=1/2  FOR BFGS/NMINF'
      READ *,IT

      IF(IT.EQ.1) THEN
        CALL BFGS(M,B,F0,G,H,NUM,REPS,AEPS,IER,NLLSQ,WK)
      ELSE
        CALL NMINF(M,B,F0,NUM,REPS,AEPS,IER,NLLSQ_F,WK)
      ENDIF
      WRITE(6,51) IER,N,M,NUM,F0,IT,(B(I),I=1,M)
      GO TO 100
      END

C     -------------------------------------------

C	To minimise a function of several variables using quasi-Newton method
C
C	N : (input) Number of variables
C	X : (input/output) Real array of length N containing the initial
C		guess for the minimum.
C		After execution it should contain the coordinates of minimiser
C	F : (output) The function value at X
C	G : (output) Real array of length N containing the gradient vector at X
C	H : (output) Real array of length N*N containing the estimated
C		Hessian matrix at X. The first dimension of H should be N.
C	NUM : (output) Number of function evaluations used by the subroutine
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy, all components of the
C		Minimiser will be calculated with accuracy MAX(AEPS, REPS*ABS(X(I)))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=53 implies that Hessian matrix at final point is probably singular
C			The iteration may have converged to a saddle point
C		IER=503 implies that N < 1, in which case no calculations are done
C		IER=526 implies that iteration failed to converge to specified accuracy
C		Other values may be set by LINMIN
C	FCN : (input) Name of the subroutine to calculate the function value
C		and its derivatives
C	WK : Real array of length 3N used as scratch space
C
C	SUBROUTINE FCN(N,X,F,G) to calculate the required function, must be supplied
C		by the user. Here N is the number of variables, F is the
C		function value at X and G is the gradient vector. X and G
C		are real arrays of length N. F and G must be calculated by FCN.
C
C	Required routines : LINMIN, FLNM, FCN
C
      SUBROUTINE BFGS(N,X,F,G,H,NUM,REPS,AEPS,IER,FCN,WK)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL FCN
      PARAMETER(NIT=200)
      DIMENSION X(N),WK(3*N),H(N,N),G(N)

      IER=0
      IF(N.LT.1) THEN
        IER=503
        RETURN
      ENDIF

      DO 2000 I=1,N
C	Initialise the Hessian matrix to unit matrix
        DO 1800 J=1,N
1800    H(J,I)=0.0
        H(I,I)=1.0
2000  CONTINUE
C	If some variable needs to be kept fixed the corresponding
C	diagonal elements, H(I,I) should be set to zero.

      CALL FCN(N,X,F,G)
      DF=ABS(F)
      IF(DF.EQ.0.0) DF=1
      N2=2*N
      NUM=1
      H2=1.

C	The iteration loop
      DO 4000 IT=1,NIT
        DF1=0.0
        H1=H2
        H2=0.0
C	Calculating the search direction WK =S^(k)
        DO 2400 I=1,N
          WK(I)=0.0
          DO 2200 J=1,N
            H2=H2+ABS(H(I,J))
2200      WK(I)=WK(I)-H(I,J)*G(J)
          DF1=DF1+WK(I)*G(I)
2400    CONTINUE

        IF(DF1.EQ.0.0) THEN
C	If gradient vanishes, then quit
C	If Hessian matrix appears to be singular, set the error flag
          IF(ABS(H2/H1).GT.1.3D0) IER=53
          RETURN
        ENDIF

C	Initial guess for line search
        X1=0
        X2=MIN(1.D0,-2.*DF/DF1)
        F1=F
        IF(X2.LE.0.0) X2=1
        CALL LINMIN(X1,X2,F1,DF1,REPS,AEPS,IER1,FCN,WK,X,N,NUM)
        IF(IER1.GT.0) IER=IER1
C	If line search fails, then quit
        IF(IER1.GT.100) RETURN

C	The convergence test
        QC=.TRUE.
        DO 2800 I=1,N
          X(I)=WK(N+I)
          WK(N+I)=X1*WK(I)
          IF(ABS(WK(N+I)).GT.MAX(REPS*ABS(X(I)),AEPS)) QC=.FALSE.
          GI=WK(N2+I)
          WK(N2+I)=WK(N2+I)-G(I)
          G(I)=GI
2800    CONTINUE
C	It is possible to apply convergence check on Function value using DF
C	instead of X(I)
        DF=F-F1
        F=F1
        IF(QC) THEN
          IF(ABS(H2/H1).GT.1.3D0) IER=53
          RETURN
        ENDIF

C	Update the matrix using BFGS formula
        DO 3200 I=1,N
          WK(I)=0.0
          DO 3200 J=1,N
3200    WK(I)=WK(I)+H(I,J)*WK(N2+J)
        GHG=0.0
        DG=0.0
        DO 3400 I=1,N
          DG=DG+WK(N2+I)*WK(N+I)
3400    GHG=GHG+WK(I)*WK(N2+I)
        R1=(1.+GHG/DG)/DG
        DO 3600 J=1,N
          DO 3600 I=1,N
3600   H(I,J)=H(I,J)+R1*WK(N+I)*WK(N+J)-(WK(N+I)*WK(J)+WK(I)*WK(N+J))/DG

4000  CONTINUE

C	Iteration fails to converge
      IER=526
      END

C     --------------------------------------------

C	Function routine to calculate the function value as required for
C	line search without derivatives
C
C	FCN : (input) Name of subroutine to calculate the required function
C	X : (input) Parameter along the line to specifying the point where
C		function evaluation is required
C	V : (input/output) Real array of length 2N, first N elements specify the
C		direction of line search. After execution next N elements will contain
C		the coordinates of the point at which function is evaluated.
C	X0 : (input) Real array of length N, containing the coordinates
C		of the starting point for line search
C	N : (input) Number of variables in the function to be minimised
C	NUM : (input/output) Integer variable to keep count of function evaluations
C
C	SUBROUTINE FCN(N,X,FX) to calculate the required function, must be supplied
C		by the user. Here N is the number of variables, FX is the
C		function value at X. X is a real array of length N.
C
C	Required routines : FCN

      FUNCTION FLN(FCN,X,V,X0,N,NUM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(N),X0(2*N)

      NUM=NUM+1
C	coordinates of the required points
      DO 1000 I=1,N
1000  X0(N+I)=X0(I)+V(I)*X
      CALL FCN(N,X0(N+1),FLN)
      END

C     --------------------------------------------------

C	Function routine to calculate the function value and its derivative
C	as required for line search
C
C	FCN : (input) Name of subroutine to calculate the required function
C	X : (input) Parameter along the line to specify the point where
C		function evaluation is required
C	DF : (output) First derivative of function along the line at X
C	V : (input/output) Real array of length 3N, first N elements specify the
C		direction of line search. Next N elements will contain the
C		coordinates of the point at which function is evaluated,
C		while the last N elements contain the gradient vector at the point
C	X0 : (input) Real array of length N, containing the coordinates
C		of starting point for line search
C	N : (input) Number of variables in the function to be minimised
C	NUM : (input/output) Integer variable to keep count of function evaluations
C
C	SUBROUTINE FCN(N,X,FX,G) to calculate the required function, must be supplied
C		by the user. Here N is the number of variables, FX is the
C		function value at X and G is the gradient vector. X and G
C		are real arrays of length N.
C
C	Required routines : FCN

      FUNCTION FLNM(FCN,X,DF,V,X0,N,NUM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(3*N),X0(N)

      NUM=NUM+1
      N2=2*N
C	The coordinates of the required point
      DO 1000 I=1,N
1000  V(N+I)=X0(I)+V(I)*X

      CALL FCN(N,V(N+1),FLNM,V(N2+1))
C	The first derivative along the search direction
      DF=0.0
      DO 2000 I=1,N
2000  DF=DF+V(I)*V(N2+I)
      END

C     ----------------------------------------------------

C	To perform a line search for minimum of several variables as required
C	by Quasi-Newton methods
C	This routine should not be used for any other purpose
C
C	X1 : (input/output) Starting value for the line search. After execution
C		it should contain the distance to minimiser along the line
C	X2 : (input/output) Initial estimate for the minimum along the line.
C		 This value will be modified by the subroutine.
C	F1 : (input/output) The function value at X1, this value must be supplied
C	DF1 : (output) The first derivative along the search direction at X1
C	REPS : (input) Required relative accuracy 
C	AEPS : (input) Required absolute accuracy,
C		These criterion is only used to terminate line search under
C		certain conditions and not generally applicable to the line
C		search. These values should be same as what is used in subroutine BFGS
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=54 implies that subroutine failed to find acceptable point
C			but the function value at last point is less than
C			that at beginning.
C		IER=55 implies that subroutine failed to find acceptable point
C			even though the interval has been reduced to required accuracy
C			but the function value at last point is less than
C			that at beginning.
C		IER=527 implies that iteration failed to find any point where
C			the function value is less than the starting value
C		IER=528 implies that iteration failed to find any point where
C			the function value is less than the starting value
C			even though the interval has been reduced to required accuracy
C	F : (input) Name of the subroutine to calculate the function value
C		and its derivatives
C	V : (input/output) Real array of length 3N. First N element specify
C		the direction in which minimisation is required. The next
C		N elements will contain the coordinates of the minimiser
C		found by LINMIN. The last 3N elements will contain the gradient
C		vector at the minimiser.
C	XI : (input) Real array of length N containing the coordinates of
C		starting point for line search
C	N : (input) Number of variables in the function to minimised
C	NUM : (output) Integer variable to keep count of the number of
C		function evaluations used so far
C
C	SUBROUTINE F(N,X,FX,G) to calculate the required function, must be supplied
C		by the user. Here N is the number of variables, FX is the
C		function value at X and G is the gradient vector. X and G
C		are real arrays of length N.
C
C	Required routines :  FLNM, F
C
      SUBROUTINE LINMIN(X1,X2,F1,DF1,REPS,AEPS,IER,F,V,XI,N,NUM)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL F
      PARAMETER(NIT=15,RHO=.01D0,SIGMA=.1D0,T1=9.,T2=.9D0,T3=.5D0)
      DIMENSION V(3*N),XI(N)

      IER=0
C	Select the bracketing phase
      QB=.FALSE.
      F2=FLNM(F,X2,DF2,V,XI,N,NUM)
      DX1=X2-X1
      F0=F1
      DF0=DF1
      X0=X1

      DO 2000 I=1,NIT
        FC=F0+DF0*RHO*(X2-X0)
        IF(ABS(DF2).LE.-SIGMA*DF0.AND.F2.LE.FC) THEN
C	Found an acceptable point
          X1=X2
          F1=F2
          DF1=DF2
          RETURN
        ENDIF
C	Test for bracketing
        IF(.NOT.QB) THEN
          IF(F2.GT.FC.OR.F2.GT.F1.OR.DF2.GE.0) QB=.TRUE.
        ENDIF

C	Hermite cubic interpolation
        DF12=(F2-F1)/DX1
        R=2.*DF2+DF1-3.*DF12
        R1=(3.*DF12-DF2-DF1)**2-DF1*DF2
        DX=0.0
        IF(R1.GT.0.0) THEN
          R1=SIGN(SQRT(R1),R)
          R=R+R1
          DX=-DF2*DX1/R
        ELSE
C	try parabolic interpolation
          R2=2.*(DF12-DF2)
          IF(R2.NE.0.0) DX=DX1*DF2/R2
        ENDIF

        IF(QB) THEN
C	Minimum is bracketed and hence improve on the bracket
          IF(DX.LT.-T2*DX1) DX=-T2*DX1
          IF(DX.GT.-T3*DX1) DX=-T3*DX1
          XA=X2+DX
          FA=FLNM(F,XA,DFA,V,XI,N,NUM)
          FC=F0+DF0*RHO*(XA-X0)

          IF(ABS(DFA).LE.-SIGMA*DF0.AND.FA.LE.FC) THEN
C	The new point is acceptable
            X1=XA
            F1=FA
            DF1=DFA
            RETURN
          ELSE IF(FA.GT.FC.OR.FA.GE.F1.OR.DFA.GT.0.0) THEN
            X2=XA
            F2=FA
            DF2=DFA
          ELSE
            X1=XA
            F1=FA
            DF1=DFA
          ENDIF

          DX1=X2-X1
          IF(ABS(DX1).LT.MAX(REPS*ABS(X2),AEPS)) THEN
C	If the interval is too small, then quit
            IER=528
            IF(F2.LE.F0) THEN
C	Accept the last point in any case
              X1=X2
              F1=F2
              DF1=DF2
              IER=55
            ENDIF
            RETURN
          ENDIF
        ELSE
C	Minimum hasn't been bracketed, choose the point further down.
          IF(DX.LT.X2-X1.AND.DX.GT.0) DX=X2-X1
          IF(DX.GT.T1*(X2-X1).OR.DX.LE.0.0) DX=T1*(X2-X1)
          X1=X2
          X2=X2+DX
          DX1=DX
          F1=F2
          DF1=DF2
          F2=FLNM(F,X2,DF2,V,XI,N,NUM)
        ENDIF
2000  CONTINUE

C	Iteration has failed to find an acceptable point
      IF(F2.LE.F0) THEN
C	accept this point if function value is smaller
        F1=F2
        X1=X2
        DF1=DF2
        IER=54
        RETURN
      ELSE IF(F1.LE.F0.AND.X1.NE.X0) THEN
        IER=54
        RETURN
      ENDIF

C	No acceptable point found
      IER=527
      END

C     ----------------------------------------------------

C	To perform a line search for minimum of several variables as required
C	by direction set method
C	This routine should not be used for any other purpose
C
C	X0 : (input/output) Starting value for the line search. After execution
C		it should contain the distance to minimiser along the line
C	X1 : (input/output) Initial estimate for the minimum along the line.
C		 This value will be modified by the subroutine.
C	F0 : (input/output) The function value at X1, this value must be supplied
C	REPS : (input) Required relative accuracy 
C	AEPS : (input) Required absolute accuracy,
C		This criterion is only used to terminate line search under
C		certain conditions and not generally applicable to the line
C		search. These values should be same as what is used in subroutine NMINF
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=56 implies that subroutine failed to find acceptable point
C			in this case the starting point itself is accepted
C	F : (input) Name of the subroutine to calculate the function value
C	V : (input/output) Real array of length 2N. First N element specify
C		the direction in which minimisation is required. The next
C		N elements will contain the coordinates of the minimiser
C		found by LINMNF.
C	XI : (input) Real array of length N containing the coordinates of
C		starting point for line search
C	N : (input) Number of variables in the function to be minimised
C	NUM : (output) Integer variable to keep count of the number of
C		function evaluations used so far
C
C	SUBROUTINE F(N,X,FX) to calculate the required function, must be supplied
C		by the user. Here N is the number of variables, FX is the
C		function value at X. X is a real array of length N.
C
C	Required routines :  FLN, F
C
      SUBROUTINE LINMNF(X0,X1,F0,REPS,AEPS,IER,F,V,XI,N,NUM)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      PARAMETER(NIT=25,T1=9.D0,GR=1.618034D0,GC=0.381966D0)
      DIMENSION V(N),XI(2*N)

      IER=0
      F1=FLN(F,X1,V,XI,N,NUM)
      IF(F1.LT.F0) THEN
C	Choose the next point further down on the same side
        X2=X0+GR*(X1-X0)
      ELSE
C	choose the next point on opposite side
        X2=X0-GR*(X1-X0)
      ENDIF
      F2=FLN(F,X2,V,XI,N,NUM)

      DO 2000 I=1,NIT
C	Parabolic interpolation
        R=(X2-X1)*(F2-F0)
        T=(X2-X0)*(F2-F1)
        P=(X2-X0)*T-(X2-X1)*R
        T=2.*(T-R)
        IF(T.GT.0) P=-P
        T=ABS(T)
        XP=X0
        IF(ABS(P).LT.ABS(T1*T*(X2-X1))) THEN
C	Try the interpolated value
          XP=X2+P/T
          FP=FLN(F,XP,V,XI,N,NUM)
          IF(FP.LT.F0) THEN
C	accept the point
            X0=XP
            F0=FP
            RETURN
          ENDIF
        ENDIF

        IF(F1.LT.F0) THEN
          X0=X1
          F0=F1
          RETURN
        ELSE IF(F2.LT.F0) THEN
          X0=X2
          F0=F2
          RETURN
        ELSE IF(XP.NE.X0) THEN
C	subdivide the interval
          IF(XP-X0.GT.0.0.EQV.X1-X0.GT.0.0) THEN
            IF(ABS(XP-X0).LT.0.5*ABS(X1-X0)) THEN
              X1=XP
              F1=FP
            ELSE
C	use golden section
              X1=X0+GC*(X1-X0)
              F1=FLN(F,X1,V,XI,N,NUM)
            ENDIF
          ELSE
            IF(ABS(XP-X0).LT.0.5*ABS(X2-X0)) THEN
              X2=XP
              F2=FP
            ELSE
C	use golden section
              X2=X0+GC*(X2-X0)
              F2=FLN(F,X2,V,XI,N,NUM)
            ENDIF
          ENDIF

C	If interpolated point is not acceptable use golden section
        ELSE IF(ABS(X2-X0).GT.ABS(X1-X0)) THEN
          X2=X0+GC*(X2-X0)
          F2=FLN(F,X2,V,XI,N,NUM)
        ELSE
          X1=X0+GC*(X1-X0)
          F1=FLN(F,X1,V,XI,N,NUM)
        ENDIF

C	If the change in function value is too small, then quit
        IF(MIN(F1-F0,F2-F0).LT.MAX(REPS*ABS(F0),AEPS)) RETURN
2000  CONTINUE

C	fails to find an acceptable point
      IER=56
      END

C     ----------------------------------------------------

C	Subroutine to calculate the Chi square function for a nonlinear
C	least squares fit, for use with subroutine BFGS.
C
C	N : (input) Number of parameters to be fitted
C	A : (input) Real array of length N containing the parameter values
C		at which the function needs to be evaluated
C	F : (output) The function value at specified parameter values
C	G : (output) Real array of length N containing the gradient vector
C		G(I) will contain dF/dA(I)
C
C	The data points are passed through common block ZZFUN, which
C	must be initialised in the calling program before calling BFGS
C
C	FX : (input) Real array of length NP containing the function values
C	X : (input) Real array of length NP containing the abscissas at which
C		function values are available
C	EF : (input) Real array of length NP containing the estimated errors in 
C		function values, for use in deciding the weights for each point
C		Although EF should contain the error estimate, but in many
C		cases it is found that multiplying all values by a suitable
C		constant can improve the convergence of BFGS dramatically
C	FX1 : (output) Real array of length NP containing the fitted value
C		of the function at each tabular point
C	NN : (input) The number of data points in the table to be fitted
C
C	The parameter NP must be equal to the dimension of arrays as
C	declared in the calling program.
C
C	This routine requires subroutine FCN to calculate the required
C	function which has to be fitted. There is no provision to pass on
C	the name of this subroutine and hence it must be changed explicitly
C	to the required name. SUBROUTINE FCN(N,A,X,F,DF) must be supplied
C	by the user. Here N is the number of parameters, A is a real array
C	of length N containing the values of parameters and X is the value
C	of independent variable where the fitting function needs to be evaluated.
C	F is the calculated function value and DF is a real array of length N
C	containing the calculated derivatives.
C
C	Required routines : FCN

      SUBROUTINE NLLSQ(N,A,F,G)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NP=100)
      COMMON/ZZFUN/FX(NP),X(NP),EF(NP),FX1(NP),NN
      DIMENSION A(N),G(N),DF(NP)
 
C     GENERATING THE FUNCTION FOR MINIMISATION
 
      IF(N.GT.NP) STOP 129
      F=0.0
      DO 1000 I=1,N
1000  G(I)=0.0
      DO 3000 I=1,NN
        CALL FCN(N,A,X(I),FA,DF)
        FX1(I)=FA
 
C     SUM OF SQUARED DIFFERENCE AND ITS GRADIENT
 
        F=F+((FX(I)-FA)/EF(I))**2
        DO 2400 J=1,N
          G(J)=G(J)+2*(FA-FX(I))*DF(J)/EF(I)**2
2400    CONTINUE
3000  CONTINUE
      END

C     ----------------------------------------------------

C	Subroutine to calculate the Chi square function for a nonlinear
C	least squares fit, for use with subroutine NMINF.
C	Version of NLLSQ for use with NMINF, when derivatives are not available
C
C	N : (input) Number of parameters to be fitted
C	A : (input) Real array of length N containing the parameter values
C		at which the function needs to be evaluated
C	F : (output) The function value at specified parameter values
C
C	The data points are passed through common block ZZFUN, which
C	must be initialised in the calling program before calling NMINF
C
C	FX : (input) Real array of length NP containing the function values
C	X : (input) Real array of length NP containing the abscissas at which
C		function values are available
C	EF : (input) Real array of length NP containing the estimated errors in 
C		function values, for use in deciding the weights for each point
C	FX1 : (output) Real array of length NP containing the fitted value
C		of the function at each tabular point
C	NN : (input) The number of data points in the table to be fitted
C
C	The parameter NP must be equal to the dimension of arrays as
C	declared in the calling program.
C
C	This routine requires subroutine FCN to calculate the required
C	function which has to be fitted. There is no provision to pass on
C	the name of this subroutine and hence it must be changed explicitly
C	to the required name. SUBROUTINE FCN(N,A,X,F) must be supplied
C	by the user. Here N is the number of parameters, A is a real array
C	of length N containing the values of parameters and X is the value
C	of independent variable where the fitting function needs to be evaluated.
C	F is the calculated function value.
C
C	Required routines : FCN

      SUBROUTINE NLLSQ_F(N,A,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NP=100)
      COMMON/ZZFUN/FX(NP),X(NP),EF(NP),FX1(NP),NN
      DIMENSION A(N)
 
C     GENERATING THE FUNCTION FOR MINIMISATION
 
      F=0.0
      DO 3000 I=1,NN
C     Change the subroutine name explicitly to FCN1
        CALL FCN1(N,A,X(I),FA)
        FX1(I)=FA
 
C     SUM OF SQUARED DIFFERENCE
 
        F=F+((FX(I)-FA)/EF(I))**2
3000  CONTINUE
      END

C     ----------------------------------------------------

C	To minimise a function of several variables using direction set method
C
C	N : (input) Number of variables
C	X : (input/output) Real array of length N containing the initial
C		guess for the minimum.
C		After execution it should contain the coordinates of minimiser
C	F : (output) The function value at X
C	NUM : (output) Number of function evaluations used by the subroutine
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy, iteration will stop when
C		change in function value is less than MAX(AEPS, REPS*ABS(F))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=504 implies that N < 2, in which case no calculations are done
C		IER=529 implies that iteration failed to converge to specified accuracy
C		Other values may be set by LINMNF
C	FCN : (input) Name of the subroutine to calculate the function value
C	WK : Real array of length N*(2N+2) used as scratch space
C
C	SUBROUTINE FCN(N,X,F) to calculate the required function, must be supplied
C		by the user. Here N is the number of variables, F is the
C		function value at X. X is a real array of length N.
C
C	Required routines : LINMNF, FLN, SVD, FCN
C
      SUBROUTINE NMINF(N,X,F,NUM,REPS,AEPS,IER,FCN,WK)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL FCN
      PARAMETER(NIT=200)
      DIMENSION X(N),WK(N,2*(N+1))

      IER=0
      IF(N.LE.1) THEN
        IER=504
        RETURN
      ENDIF

C	Initialise the direction set matrix to identity matrix
      DO 2000 I=1,N
        WK(I,N+3)=0.0
        WK(I,N+4)=0.0
        DO 1800 J=1,N
1800    WK(J,I)=0.0
        WK(I,I)=1.0
2000  CONTINUE

      CALL FCN(N,X,F)
      NUM=1

C	The main iteration loop
      DO 6000 IT=1,NIT
        FI=F
        DO 4000 K=1,N

C	The starting point for line search
          DO 2400 I=1,N
2400      WK(I,N+1)=X(I)
          F0=F
          DFMAX=0.0
          KMAX=1

          DO 2800 I=1,N
            X1=0.0
C	Use previous value as initial approximation to minimum
            X2=WK(I,N+3)
            FUN=F0
            IF(X2.EQ.0.0) X2=1
            CALL LINMNF(X1,X2,F0,REPS,AEPS,IER1,FCN,WK(1,I),WK(1,N+1),
     1                  N,NUM)
            IF(IER1.GT.0) IER=IER1

            WK(I,N+3)=X1
C	Estimate of second derivative along the line
            IF(X1.NE.0.0) WK(I,N+4)=ABS((FUN-F0)/X1**2)
C	The new starting point
            DO 2600 J=1,N
2600        WK(J,N+1)=WK(J,N+1)+X1*WK(J,I)
            IF(FUN-F0.GE.DFMAX.AND.I.LE.N-K+1) THEN
              DFMAX=FUN-F0
              KMAX=I
            ENDIF
2800      CONTINUE

C	Remove the KMAX th direction
          DO 3000 I=KMAX,N-1
            WK(I,N+3)=WK(I+1,N+3)
            WK(I,N+4)=WK(I+1,N+4)
            DO 3000 J=1,N
3000      WK(J,I)=WK(J,I+1)
C	Add a new direction
          DO 3200 I=1,N
3200      WK(I,N)=WK(I,N+1)-X(I)

          X1=0.0
          X2=1
          FUN=F
          WK(N,N+4)=0.0
C	Starting point for the final line search in the loop
          DO 3400 I=1,N
3400      WK(I,N+1)=X(I)
          CALL LINMNF(X1,X2,F,REPS,AEPS,IER1,FCN,WK(1,N),WK(1,N+1),
     1                N,NUM)
          IF(IER1.GT.0) IER=IER1

          WK(N,N+3)=X1
          IF(X1.NE.0.0) WK(N,N+4)=ABS((FUN-F)/X1**2)
          DO 3600 J=1,N
3600      X(J)=X(J)+X1*WK(J,N)
4000    CONTINUE

        IF(ABS(F-FI).LT.MAX(REPS*ABS(F),AEPS)) RETURN

C	The matrix V for SVD
        DO 4500 J=1,N
          IF(WK(J,N+4).GT.0.0) THEN
            DO 4400 I=1,N
4400        WK(I,J)=WK(I,J)/SQRT(WK(J,N+4))
          ENDIF
4500    CONTINUE
        M=N
        LA=N
        LV=N
        CALL SVD(N,M,WK,WK(1,N+3),WK(1,N+1),LA,LV,WK(1,N+2),IER1)
        IF(IER1.GT.0) THEN
          IER=IER1
          RETURN
        ENDIF

        DO 4600 I=1,N
          WK(I,N+3)=SQRT(FUN-F)*WK(I,N+1)
          WK(I,N+4)=0.0
          IF(WK(I,N+1).NE.0.0) WK(I,N+4)=1./WK(I,N+1)**2
4600    CONTINUE

6000  CONTINUE
      IER=529
      END

C     ----------------------------------------------------

C	To calculate the Singular Value Decomposition of a matrix A=U D V-transpose
C
C	N : (input) Number of variables
C	M : (input) Number of equations
C	A : (input/output) Matrix of coefficients of size LA*N
C		After execution it will contain the matrix U
C	V : (output) The matrix V of size LV*N
C	SIGMA : (output) Array of length N, containing the singular values
C	LA : (input) Actual value of first dimension of A in the calling program
C	LV : (input) Actual value of first dimension of V in the calling program
C	E : Scratch array of length N
C	IER : (output) Error parameter, IER=0 if execution is successful
C		IER=12 QR iteration failed to converge to required accuracy
C		IER=105 implies N.LE.0, N.GT.LV, M.LE.0, M.GT.LA, N.GT.M
C
C	Required routines : None

      SUBROUTINE SVD(N,M,A,V,SIGMA,LA,LV,E,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ITMAX=30,REPS=1.D-16)
C	For REAL*4 use REPS=6.E-8
C      PARAMETER(ITMAX=30,REPS=6.E-8)
      DIMENSION A(LA,N),V(LV,N),SIGMA(N),E(N)

      IF(N.GT.M.OR.N.LE.0.OR.M.LE.0.OR.M.GT.LA.OR.N.GT.LV) THEN
        IER=105
        RETURN
      ENDIF

      IER=0
C	Reduction to Bidiagonal form using Householder transformations
      G=0
      RMAX=0

      DO 3000 I=1,N
C	Off-diagonal elements of bidiagonal form
        E(I)=G
        S=0
        DO 1200 J=I,M
1200    S=S+A(J,I)**2

        IF(S.LE.0.0) THEN
C	transformation not required
          G=0
        ELSE
          F=A(I,I)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I)=F-G

          DO 1800 J=I+1,N
            S=0
            DO 1400 K=I,M
1400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 1600 K=I,M
1600        A(K,J)=A(K,J)+F*A(K,I)
1800      CONTINUE
        ENDIF

C	Diagonal elements of bidiagonal form
        SIGMA(I)=G
        S=0
        DO 2000 J=I+1,N
2000    S=S+A(I,J)**2

        IF(S.LE.0.0) THEN
C	Transformation not required
          G=0
        ELSE
          F=A(I,I+1)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I+1)=F-G
          DO 2200 J=I+1,N
C	Temporary storage of intermediate results
2200      E(J)=A(I,J)/H

          DO 2800 J=I+1,M
            S=0
            DO 2400 K=I+1,N
2400        S=S+A(J,K)*A(I,K)
            DO 2600 K=I+1,N
2600        A(J,K)=A(J,K)+S*E(K)
2800      CONTINUE
        ENDIF
        R1=ABS(SIGMA(I))+ABS(E(I))
        IF(R1.GT.RMAX) RMAX=R1
3000  CONTINUE

C	Accumulation of right hand transformation in array V
      DO 4000 I=N,1,-1
        IF(G.NE.0.0) THEN
          H=A(I,I+1)*G
          DO 3200 J=I+1,N
3200      V(J,I)=A(I,J)/H

          DO 3800 J=I+1,N
            S=0
            DO 3400 K=I+1,N
3400        S=S+A(I,K)*V(K,J)
            DO 3600 K=I+1,N
3600        V(K,J)=V(K,J)+S*V(K,I)
3800      CONTINUE
        ENDIF

        DO 3900 J=I+1,N
          V(I,J)=0.0
          V(J,I)=0.0
3900    CONTINUE
        V(I,I)=1
        G=E(I)
4000  CONTINUE

C	Accumulation of left hand transformation overwritten on matrix A
      DO 5000 I=N,1,-1
        G=SIGMA(I)
        DO 4200 J=I+1,N
4200    A(I,J)=0
        IF(G.NE.0.0) THEN
          H=A(I,I)*G

          DO 4700 J=I+1,N
            S=0
            DO 4400 K=I+1,M
4400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 4600 K=I,M
4600        A(K,J)=A(K,J)+F*A(K,I)
4700      CONTINUE

          DO 4800 J=I,M
4800      A(J,I)=A(J,I)/G
        ELSE
          DO 4900 J=I,M
4900      A(J,I)=0.0
        ENDIF
        A(I,I)=A(I,I)+1
5000  CONTINUE

C	Diagonalisation of the bidiagonal form
      AEPS=REPS*RMAX
C	Loop over the singular values
      DO 8000 K=N,1,-1
C	The QR transformation
        DO 7500 ITR=1,ITMAX

C	Test for splitting
          DO 5200 L=K,1,-1
            IF(ABS(E(L)).LT.AEPS) GO TO 6000
            IF(ABS(SIGMA(L-1)).LT.AEPS) GO TO 5400
5200      CONTINUE

C	cancellation of E(L) if L>1
5400      C=0.0
          S=1.0
          DO 5800 I=L,K
            F=S*E(I)
            E(I)=C*E(I)
            IF(ABS(F).LT.AEPS) GO TO 6000
            G=SIGMA(I)
            SIGMA(I)=SQRT(F*F+G*G)
            C=G/SIGMA(I)
            S=-F/SIGMA(I)

            DO 5600 J=1,M
              R1=A(J,L-1)
              R2=A(J,I)
              A(J,L-1)=R1*C+R2*S
              A(J,I)=C*R2-S*R1
5600        CONTINUE
5800      CONTINUE

6000      Z=SIGMA(K)
          IF(L.EQ.K) THEN
C	QR iteration has converged
            IF(Z.LT.0.0) THEN
              SIGMA(K)=-Z
              DO 6200 J=1,N
6200          V(J,K)=-V(J,K)
            ENDIF
            GO TO 8000
          ENDIF

          IF(ITR.EQ.ITMAX) THEN
            IER=12
            GO TO 7500
          ENDIF

C	calculating shift from bottom 2x2 minor
          X=SIGMA(L)
          Y=SIGMA(K-1)
          G=E(K-1)
          H=E(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.*H*Y)
          G=SQRT(1.+F*F)
          IF(F.LT.0.0) G=-G
          F=((X-Z)*(X+Z)+H*(Y/(F+G)-H))/X

C	next QR transformation
          C=1.0
          S=1.0
C	Given's rotation
          DO 7000 I=L+1,K
            G=E(I)
            Y=SIGMA(I)
            H=S*G
            G=C*G
            E(I-1)=SQRT(F*F+H*H)
            C=F/E(I-1)
            S=H/E(I-1)
            F=C*X+S*G
            G=C*G-S*X
            H=S*Y
            Y=C*Y

            DO 6400 J=1,N
              X=V(J,I-1)
              Z=V(J,I)
              V(J,I-1)=C*X+S*Z
              V(J,I)=C*Z-S*X
6400        CONTINUE

            SIGMA(I-1)=SQRT(F*F+H*H)
            IF(SIGMA(I-1).NE.0.0) THEN
              C=F/SIGMA(I-1)
              S=H/SIGMA(I-1)
            ENDIF
            F=C*G+S*Y
            X=C*Y-S*G
            DO 6600 J=1,M
              Y=A(J,I-1)
              Z=A(J,I)
              A(J,I-1)=C*Y+S*Z
              A(J,I)=C*Z-S*Y
6600        CONTINUE
7000      CONTINUE

          E(L)=0
          E(K)=F
          SIGMA(K)=X
7500    CONTINUE
8000  CONTINUE
      END

C     ----------------------------------------------------

C	To generate random numbers with Gaussian probability distribution
C	It generates random numbers with zero mean and variance of 1.
C	
C	ISEED : (input/output) integer seed, it should be positive and
C		less than AN. It is updated by the routine and should
C		not be modified between two calls, unless a fresh
C		sequence is required
C
C	Required routines : None

      FUNCTION RANGAU(ISEED)
      IMPLICIT REAL*8(A-H,O-Z)
C	Retain the following declaration even for REAL*4 version
C	otherwise AN and A will be rounded
      REAL*8 AN,A,B
      PARAMETER(AN=199017.0,A=24298.0,B=99991.,PI=3.14159265358979324D0)
 
      N1=1+MOD(A*ISEED+B,AN)
      RANGAU=SQRT(2.D0*LOG(AN/ISEED))*COS(2.0*PI*N1/AN)
      ISEED=N1
      END

C     ----------------------------------------------------

      SUBROUTINE FCN1(N,A,X,F)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(*)

C     GENERATING THE FITTING FUNCTION FOR USE WITH NMINF

        F=0.0
        DO 2200 J=1,N,2
2200    F=F+A(J)*EXP(-X*A(J+1))

      END

C     ----------------------------------------------------

      SUBROUTINE FCN(N,A,X,F,G)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(*),G(9)

C     GENERATING THE FITTING FUNCTION FOR USE WITH BFGS

        F=0.0
        DO 2200 J=1,N,2
2200    F=F+A(J)*EXP(-X*A(J+1))

C     THE GRADIENT

        DO 2400 J=1,N,2
          G(J)=EXP(-X*A(J+1))
          G(J+1)=-X*A(J)*EXP(-X*A(J+1))
2400    CONTINUE
      END
