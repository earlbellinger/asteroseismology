C     PROGRAM TO SOLVE NONLINEAR PARABOLIC EQUATIONS USING THE
C     METHOD OF LINES

      PROGRAM PARBLC
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL LINES,F
      PARAMETER(NP=101,NQ=5)
      COMMON/ZZLINE/DX0,XX0,XN,XI(NP),U0(NQ),UN(NQ),UX(NQ,2),NX,NE
      COMMON/FN/TT,XX
      DIMENSION X(350),DX(350),WK(16000),YF(90),IWK(90)

C     EXAMPLE 13.2      

51    FORMAT('   INITIAL VALUES =',1P4D14.6/(2X,5D14.6))
52    FORMAT('   IER =',I4,5X,'TIME =',1PD14.6,5X,'TIME STEP =',D14.6/
     1       5X,'NO. OF PTS =',I5,5X,'NO. OF FUNCTION EVALUATIONS =',I8/
     2       5X,'SOLUTION =',4D14.6/(3X,5D14.6))

      NE=1
C	KEEP STEP-SIZE FIXED
      IFLG=1
      IST=1
      NMAX=100000
      PRINT *,'TYPE T0=INITIAL TIME, H=TIME STEP, NX=NO. OF PTS'
      READ *,T0,H,NX
      N=NX-2
      DX0=1.D0/(NX-1)
      XX0=0.0
      XN=1.0
      TT=T0
      REPS=1.D-6
      AEPS=1.D-7

C     FINDING THE INITIAL VALUES FOR THE NONLINEAR WAVE BY SOLVING
C     IMPLICIT EQUATION

      DO 500 I=1,NX
        XI(I)=DX0*I
        XX=XI(I)
        IF(100.*XX.GT.5.0) THEN
          XS=1.+EXP(-300.*XX)
        ELSE
          XL=1.0
          XU=1.D5
          CALL BRENT(XL,XU,XS,REPS,AEPS,IER,F)
          IF(IER.GT.0) STOP
        ENDIF
        X(I)=XS
500   CONTINUE
      WRITE(6,51) (X(I),I=1,NX-1)
      U0(1)=X(1)
      UN(1)=X(NX-1)

C     GO ON TYPING INTERMEDIATE VALUES OF T AT WHICH SOLUTION IS REQUIRED

600   PRINT *,'TYPE TN=FINAL TIME        (QUITS WHEN TN.LT.-10)'
      READ *,TN
      IF(TN.LT.-10.0) STOP

      REPS=1.D-4
      CALL MSTEP(N,X,DX,LINES,H,T0,TN,YF,REPS,NSTEP,NMAX,IER,IFLG,IST,WK
     1           ,IWK)
      WRITE(6,52) IER,TN,H,NX,NSTEP,U0(1),(YF(I),I=1,N),UN(1)
      IF(IER.EQ.0) GO TO 600
      END

C     --------------------------------------------------------

C	To perform one step of solution of initial value problems in ordinary
C	differential equations using fourth order Adams-Bashforth-Moulton
C	predictor corrector method. It is called by subroutine MSTEP.
C
C	N : (input) Number of first order differential equations to be solved
C	Y : (input/output) Real array of length 7N containing the solution
C		at last four points. After execution new point will be added.
C	DY : (input/output) Real array of length 7N containing the derivatives
C		of Y at the points in Y
C	DIF : (input) Name of subroutine to calculate the right hand side
C		of differential equation y'=f(t,y)
C	H : (input) Step length to be used.
C	T : (input) Value of independent variable t, at which the solution
C		is required.
C	REPS : (input) Required accuracy in each component of the solution.
C		The subroutine only controls error in iteration on corrector.
C	NSTP : (input/output) Number of calls to DIF required so far.
C		The count is updated by the subroutine.
C	IJ : (input) The index j+1 in corrector formula
C	IJM1 : (input) The index j in corrector formula
C	IJM2 : (input) The index j-1 in corrector formula
C	IJM3 : (input) The index j-2 in corrector formula
C	IJM4 : (input) The index j-3 in corrector formula
C	IER : (output) Error parameter; IER=0 implies successful execution
C		IER=729 implies that iteration on corrector failed to converge
C	WK : Real array of length 2N used to pass the predicted values
C
C	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
C	the differential equation. T is the value of independent variable,
C	N is the number of variables, Y is a real array of length N containing
C	the values of variables. DY is a real array of length N which should
C	contain the calculated values of derivatives at (T,Y).
C
C	Required routines : DIF
C	
C	
      SUBROUTINE ADAMS(N,Y,DY,DIF,H,T,REPS,NSTP,IJ,IJM1,IJM2,
     1                 IJM3,IJM4,IER,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(EPS=1.D-30,NIT=10,CFAC=0.2D0)
      DIMENSION Y(N,7),DY(N,7),WK(N,2)

      IER=0
      DO 2000 J=1,N
C	The predictor
        T1=Y(J,IJM1)+H*(55.*DY(J,IJM1)-59.*DY(J,IJM2)+37.*DY(J,IJM3)
     1     -9.*DY(J,IJM4))/24.
C	Modifier
        Y(J,IJ)=T1+251.*(Y(J,IJM1)-WK(J,2))/270.
2000  WK(J,1)=T1

C	Iteration on corrector
      DO 3000 J=1,NIT
        NSTP=NSTP+1
        CALL DIF(T,N,Y(1,IJ),DY(1,IJ))
        ERR=0.0
        DO 2600 K=1,N
C	The corrector
          T2=Y(K,IJM1)+H*(9.*DY(K,IJ)+19.*DY(K,IJM1)-5.*DY(K,IJM2)
     1       +DY(K,IJM3))/24.
          R1=ABS(T2)+ABS(T2-Y(K,IJM1))+EPS
          T1=ABS((Y(K,IJ)-T2)/R1)
          IF(T1.GT.ERR) ERR=T1
2600    Y(K,IJ)=T2
C	The convergence test
        IF(ERR.LT.CFAC*REPS) RETURN
3000  CONTINUE

C	iteration on corrector fails to converge
      IER=729

      END

C   ----------------------------------------------------------------

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

C     -------------------------------------------------

C	To perform one step of solution of initial value problems in ordinary
C	differential equations using fourth order stiffly stable method.
C	It is called by subroutine MSTEP.
C
C	N : (input) Number of first order differential equations to be solved
C	Y : (input/output) Real array of length 7N containing the solution
C		at last four points. After execution new point will be added.
C	DY : (input/output) Real array of length 7N containing the derivatives
C		of Y at the points in Y
C	DIF : (input) Name of subroutine to calculate the right hand side
C		of differential equation y'=f(t,y)
C	H : (input) Step length to be used.
C	T : (input) Value of independent variable t, at which the solution
C		is required.
C	REPS : (input) Required accuracy in each component of the solution.
C		The subroutine only controls error in iteration on corrector.
C	NSTP : (input/output) Number of calls to DIF required so far.
C		The count is updated by the subroutine.
C	IJ : (input) The index j+1 in corrector formula
C	IJM1 : (input) The index j in corrector formula
C	IJM2 : (input) The index j-1 in corrector formula
C	IJM3 : (input) The index j-2 in corrector formula
C	IJM4 : (input) The index j-3 in corrector formula
C	IFLAG : (input/output) Integer variable used as a flag
C		If IFLAG=0 initial approximation to Jacobian is generated
C			and IFLAG is set to 1
C		Otherwise old approximation to J^{-1} is used.
C	IER : (output) Error parameter; IER=0 implies successful execution
C		IER=729 implies that iteration on corrector failed to converge
C			or cannot be continued because the estimated Jacobian is singular.
C	WK1 : (input/output) Real array of length 3N used to pass the
C		predicted values
C	WK : Real array of length N*MAX(2N, N+4) used as scratch space
C	IWK : Integer array of length N used as scratch space
C
C	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
C	the differential equation. T is the value of independent variable,
C	N is the number of variables, Y is a real array of length N containing
C	the values of variables. DY is a real array of length N which should
C	contain the calculated values of derivatives at (T,Y).
C
C	Required routines : GAUELM, DIF
C	
C	
      SUBROUTINE GEAR(N,Y,DY,DIF,H,T,REPS,NSTP,IJ,IJM1,IJM2,
     1                IJM3,IJM4,IFLAG,IER,WK1,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(EPS=1.D-30,NIT=20,CFAC=1.D-2)
      DIMENSION Y(N,7),DY(N,7),WK1(N,3),WK(N,*),IWK(N)

      IER=0
      DO 1000 J=1,N
C	predictor
        T1=Y(J,IJM1)+H*(55.*DY(J,IJM1)-59.*DY(J,IJM2)+37.*DY(J,IJM3)
     1     -9.*DY(J,IJM4))/24.
C	Modifier
        Y(J,IJ)=T1+6275.*(Y(J,IJM1)-WK1(J,2))/8003.
1000  WK1(J,1)=T1

      LJ=N
      IPAS=0
      IER=0
      NSTP=NSTP+1
      CALL DIF(T,N,Y(1,IJ),DY(1,IJ))
      DO 1200 J=1,N
C	Residual in the corrector
        WK1(J,3)=Y(J,IJ)-(48.*Y(J,IJM1)-36.*Y(J,IJM2)+16.*Y(J,IJM3)-
     1           3.*Y(J,IJM4)+12.*H*DY(J,IJ))/25.
1200  CONTINUE

      IF(IFLAG.EQ.0) THEN
C	Generate initial approximation to Jacobian J 
        DO 2000 IP=1,N
          X1=Y(IP,IJ)
          DK=100.*REPS*ABS(X1)
          IF(DK.EQ.0.0) DK=100.*REPS
          Y(IP,IJ)=X1+DK
          NSTP=NSTP+1
          CALL DIF(T,N,Y(1,IJ),DY(1,IJ))
          DO 1800 J=1,N
            WK(J,1)=Y(J,IJ)-(48.*Y(J,IJM1)-36.*Y(J,IJM2)+16.*Y(J,IJM3)-
     1              3.*Y(J,IJM4)+12.*H*DY(J,IJ))/25.
1800      WK(J,N+IP)=(WK(J,1)-WK1(J,3))/DK
          Y(IP,IJ)=X1
2000    CONTINUE

        NUM=N
        IFLG=0
        DO 2500 I=1,N
          DO 2400 J=1,N
2400      WK(J,I)=0.0
2500    WK(I,I)=1.0

C	Calculate inverse of Jacobian
        CALL GAUELM(N,NUM,WK(1,N+1),WK,DET,IWK,LJ,IER1,IFLG)
        IF(IER1.GT.0) THEN
          IER=729
          RETURN
        ENDIF
C	Reset the flag so that Jacobian is not calculated again
        IFLAG=1
      ENDIF

      DO 2800 I=1,N
        WK(I,N+1)=WK1(I,3)
        WK(I,N+2)=Y(I,IJ)
2800  CONTINUE

C	Iteration on parameters using Broyden's method
      DO 5000 IPAS=1,NIT

C	The convergence check
        ERR=0.0
        DO 3500 J=1,N
          S1=0.0
C	Correction to the Jth component
          DO 3200 K=1,N
3200      S1=S1+WK(J,K)*WK(K,N+1)
          R2=ABS(Y(J,IJ))+ABS(Y(J,IJ)-Y(J,IJM1))+EPS
          ERR=MAX(ERR,ABS(S1/R2))
          Y(J,IJ)=WK(J,N+2)-S1
          WK(J,N+3)=-S1
3500    CONTINUE

        NSTP=NSTP+1
        CALL DIF(T,N,Y(1,IJ),DY(1,IJ))
        IF(ERR.LT.CFAC*REPS) RETURN

C	calculate the residuals
        DO 3800 J=1,N
          WK1(J,3)=Y(J,IJ)-(48.*Y(J,IJM1)-36.*Y(J,IJM2)+16.*Y(J,IJM3)-
     1             3.*Y(J,IJM4)+12.*H*DY(J,IJ))/25.
3800    CONTINUE

C	Update J inverse using Broyden's formula
        DO 4000 J=1,N
          WK(J,N+4)=WK1(J,3)-WK(J,N+1)
          WK(J,N+1)=WK1(J,3)
          WK(J,N+2)=Y(J,IJ)
4000    CONTINUE
        SS=0.0
        DO 4300 J=1,N
          S1=0.0
          S2=0.0
          DO 4100 K=1,N
            S1=S1+WK(K,N+3)*WK(K,J)
            S2=S2+WK(J,K)*WK(K,N+4)
4100      CONTINUE
          WK1(J,3)=S1
          Y(J,IJ)=S2-WK(J,N+3)
          SS=SS+S1*WK(J,N+4)
4300    CONTINUE

        IF(SS.EQ.0.0) THEN
          IER=729
          RETURN
        ENDIF
        DO 4400 J=1,N
          DO 4400 K=1,N
            WK(K,J)=WK(K,J)-Y(K,IJ)*WK1(J,3)/SS
4400    CONTINUE
5000  CONTINUE

C	Iteration fails to converge
      IER=729
      END

C     -------------------------------------------------

C	To solve a system of nonlinear parabolic differential equation using
C	the method of lines
C	The differential equation is assumed to be of the form
C
C	du/dt = f_i( d^2u/dx^2, du/dx, u, x, t),    i=1,2,...,NE
C
C	with Dirichlet boundary conditions
C
C	u(x0,t) = g1(t);  	u(xN,t) = g2(t)
C
C	This subroutine is used through subroutine MSTEP or RKM for solving
C	a system of ordinary differential equations. This routine specifies
C	the equivalent system of ordinary differential equations.
C	The parameters for the partial differential equations are passed
C	through the common block ZZLINE. Elements of this common block must
C	be initialised before calling MSTEP or RKM.
C
C	T : (input) Value of "time" where the derivatives need to be calculated.
C	N : (input) Number of ordinary differential equations to be solved
C		after applying the method of lines, N=NE*(NX-2)
C	U : (input) Real array of length N containing the estimated
C		solution at t.
C	DU : (output) Real array of length N which will contain the derivatives
C		of U at T.
C
C	The arguments in common block are as follows:
C
C	DX : (input) The step length in X to be used for computations.
C		It is assumed to be uniform.
C	X0 : (input) Lower limit on X where the solution is required
C	XN : (input) Upper limit on X where the solution is required.
C		Solution is computed in the interval (X0,XN)
C	X : (input) Real array of length NX-2 containing the mesh points used
C		in X direction. The end points X0 and XN are not included.
C		This array must be set before calling MSTEP or RKM.
C	U0 : (output) Real array of length NE containing the solution
C		at the first boundary X0. It is calculated using the
C		boundary conditions calculated by subroutine BC.
C	UN : (output) Real array of length NE containing the solution
C		at the second boundary XN. It is calculated using the
C		boundary conditions calculated by subroutine BC.
C	UX : Real array of length at least 2*NE used as scratch space
C	NE : (input) Number of differential equations in the parabolic system
C	NX : (input) Number of mesh points in the X direction.
C
C	The dimensions of all arrays in the common block must match those
C	in the calling program and should be sufficient to define the
C	equations. No checks are made and there is no error exit.
C	It is the responsibility of user to ensure consistency.
C
C	This subroutine requires subroutines FCN to calculate the derivatives
C	defining the differential equations and subroutine BC to calculate
C	the boundary conditions. The names of these routines are fixed and
C	cannot be passed on as arguments.
C
C	SUBROUTINE FCN(NE,X,T,U,UX,UXX,DU) and SUBROUTINE BC(NE,T,X0,XN,U0,UN)
C	must be supplied by the user 
C	Subroutine FCN should calculate the derivatives DU = dU/dt at
C	specified T, X, U, UX (=dU/dX), UXX (=d^2U/dX^2)
C	Here NE is the number of parabolic equations in the system, T and X
C	are the values of t and x, U, UX, UXX and DU are real arrays of
C	length NE containing the solution and its derivative as described
C	above.
C
C	Subroutine BC should calculate the solution at boundaries at x=X0 and XN
C	Here NE is the number of parabolic equations in the system,
C	T is the value of t at which boundary conditions are required.
C	X0 and XN are the two boundaries in X where boundary conditions are
C	applied. U0 and UN are real arrays of length NE containing the
C	solution at X0 and XN. U0 and UN must be calculated by subroutine BC.
C
C	Required routines : FCN, BC

      SUBROUTINE LINES(T,N,U,DU)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NP=101,NQ=5)
      COMMON/ZZLINE/DX,X0,XN,X(NP),U0(NQ),UN(NQ),UX(NQ,2),NX,NE
      DIMENSION U(N),DU(N)

C	Calculate the boundary values
      CALL BC(NE,T,X0,XN,U0,UN)
      DO 2000 I=1,N,NE
        J=(I-1)/NE+1

        DO 1500 K=1,NE
          IF(I.EQ.1) THEN
            UM=U0(K)
          ELSE
            UM=U((J-2)*NE+K)
          ENDIF

          IF(I+NE.GT.N) THEN
            UP=UN(K)
          ELSE
            UP=U(J*NE+K)
          ENDIF

          UX(K,1)=(UP-UM)/(2.*DX)
          UX(K,2)=(UP-2.*U(I-1+K)+UM)/DX**2
1500    CONTINUE

        CALL FCN(NE,X(J),T,U(I),UX,UX(1,2),DU(I))
2000  CONTINUE
      END

C     --------------------------------------------------

C	To solve initial value problems in ordinary differential equations
C	using a fourth-order multistep method with adaptive step size control
C	It can be used with ADAMS (Adams-Bashforth-Moulton predictor corrector method)
C	or with GEAR (Gear's stiffly stable method)
C
C	N : (input) Number of first order differential equations to be solved
C	Y : (input/output) Real array of length 7N, the first N elements
C		should contain the initial values of variables. During
C		execution, the solution at various intermediate points will
C		be stored in this array. The contents of this array must
C		be preserved between two calls to MSTEP if the integration
C		needs to be continued further.
C	DY : (output) Real array of length 7N, containing the derivatives
C		of Y the points stored in Y. This array should also be
C		preserved between two calls to MSTEP.
C	DIF : (input) Name of subroutine to calculate the right hand side
C		of differential equation y'=f(t,y)
C	H : (input/output) Initial guess for the step size. After execution
C		it will contain the step size used by the program
C	T0 : (input/output) Initial value of independent variable t, at
C		which initial values are specified. After execution it will
C		be set to the point up to which integration has been successful
C	TN : (input) The final value of t at which the solution is required.
C		Intermediate values will not be preserved so if solution
C		is required at intermediate points, TN must be set to first
C		such value and multiple calls will be needed to calculate
C		all required values. For each subsequent call only TN needs
C		to be updated. Other variables including the scratch arrays
C		in the call statement must be preserved to their old values.
C	YF : (output) Real array of length N, containing the solution at
C		the last point up to which integration is successful.
C		If execution is successful it will contain the required
C		solution at TN.
C	REPS : (input) Required accuracy in each component of the solution.
C		The subroutine only controls local truncation error and hence
C		actual error could be larger
C	NSTP : (input/output) Number of calls to DIF required since the starting
C		of integration. The count is updated by the subroutine.
C	NMAX : (input/output) Maximum number of calls to DIF to be tried.
C		If NMAX.LE.0 it will be set to a default value of NMX=100000.
C	IER : (output) Error parameter; IER=0 implies successful execution
C		IER=702 implies N.LE.0, no calculations are done
C		IER=724 or 725 implies STRT4 failed to generate starting values
C		IER=726 implies that step-size has become smaller than
C			REPS*|TN-T0|
C		IER=727 implies that step size is too small for arithmetic used
C		IER=728 implies that integration could not be completed in
C			the specified number of steps.
C	IFLG : (input/output) Integer variable used as flag to decide the
C		type of computation.
C		If IFLG=0 the computation is started afresh with starting
C			values being generated and the step size is adjusted
C			to meet the specified accuracy. NSTP is also initialised
C			to zero in the beginning. After execution IFLG is set to 2
C		If IFLG=1 the computation is started afresh with starting
C			values being generated but the step size is kept fixed.
C			NSTP is also initialised to zero in the beginning.
C			After execution IFLG is set to 3
C		If IFLG=2 the integration is continued using already available
C			solution at previous steps. The step size is adjusted
C			to meet the specified accuracy.
C		If IFLG=3 the integration is continued using already available
C			solution at previous steps. The step size is kept fixed.
C		IFLG should not be changed after the first call, unless
C		fresh starting values need to be generated.
C	IST : (input/output) Integer variable used as flag to decide the
C		multistep method to be used:
C		If IST=0 Adams method is used, which may be good for
C			general equations
C		If IST=1 Gear's method is used, which may be used for
C			stiff equations
C	WK : Real array used as scratch space. For use with ADAMS the size
C		of WK should be 2N, while for use with GEAR it should be
C		3N+N*MAX(N+4, 2N). This scratch space must be preserved
C		between two calls to the routine
C	IWK : Integer array of length N used as scratch space.
C
C	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
C	the differential equation. T is the value of independent variable,
C	N is the number of variables, Y is a real array of length N containing
C	the values of variables. DY is a real array of length N which should
C	contain the calculated values of derivatives at (T,Y).
C
C	Required routines : ADAMS, GEAR, STRT4, RK4, GAUELM, DIF
C	
C	
      SUBROUTINE MSTEP(N,Y,DY,DIF,H,T0,TN,YF,REPS,NSTP,NMAX,IER,IFLG,
     1                 IST,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL DIF
      PARAMETER(NMX=100000,EPS=1.D-30)
      DIMENSION Y(N,7),DY(N,7),WK(N,*),YF(N),DEL(4),IWK(*)
      SAVE

      IF(N.LE.0) THEN
        IER=702
        RETURN
      ENDIF
      IF(NMAX.LE.0) NMAX=NMX
      IER=0
      TSTEP=TN-T0

      IF(IFLG.LE.1) THEN
        IF(TSTEP.EQ.0.0) RETURN
C	Generate the starting values
        NSTP=0
        IF(H.EQ.0.0) H=TSTEP/8
C	Change the sign of H if needed
        IF(H.LT.0.0.EQV.TN.GT.T0) H=-H

        T=T0
C	Generate the starting values
        CALL STRT4(N,Y,DY,DIF,H,T,REPS,IFLG,TSTEP,NSTP,IER,WK)
        IF(IER.GT.0) RETURN

C	Initialising the array indices
        IJ=4
        IJM1=3
        IJM2=2
        IJM3=1
        T0=T
C	Set the flag for subroutine GEAR
        IFLAG=0
C	Update IFLG
        IFLG=IFLG+2
        NS=4
      ENDIF

      IF(TN.LE.T0.EQV.H.GT.0) GO TO 6000
      IF(TSTEP.EQ.0.0) GO TO 6000
C	Updating the array indices
3400  IJM4=IJM3
      IJM3=IJM2
      IJM2=IJM1
      IJM1=IJ
      IJ=IJ+1
      IF(IJ.GT.7) IJ=1
      T=T0+H
      REPS1=ABS(REPS*H/TSTEP)

C	To perform one step of integration using Adams method
      IF(IST.EQ.0) THEN
        CALL ADAMS(N,Y,DY,DIF,H,T,REPS1,NSTP,IJ,IJM1,IJM2,
     1             IJM3,IJM4,IER,WK)
      ELSE

C	To perform one step of integration using stiffly stable method
        CALL GEAR(N,Y,DY,DIF,H,T,REPS1,NSTP,IJ,IJM1,IJM2,
     1            IJM3,IJM4,IFLAG,IER,WK,WK(1,4),IWK)
      ENDIF

C	Estimating the truncation error
      ERR=0
      DO 4500 I=1,N
        R2=ABS(Y(I,IJ))+ABS(Y(I,IJ)-Y(I,IJM1))+EPS
        R1=ABS(WK(I,1)-Y(I,IJ))/R2
        IF(R1.GT.ERR) ERR=R1
4500  CONTINUE
      ERR=ABS(0.25*ERR*TSTEP/H)

      IF(T.EQ.T0) THEN
C	Step size is too small for the arithmetic used
        IER=727
        DO 4600 I=1,N
4600    YF(I)=Y(I,IJ)
        RETURN
      ENDIF

      IF(IFLG.LE.2) THEN
        IF(ERR.LT.REPS.AND.IER.EQ.0) THEN
C	Integration is successful
          T0=T0+H
          DO 4700 I=1,N
4700      WK(I,2)=WK(I,1)
          NS=NS+1
          IF(TN.LE.T0.EQV.H.GT.0) GO TO 6000
          IF(ABS((TN-T0)/TSTEP).LE.REPS) GO TO 6000

          IF(ERR.LT.REPS/32.AND.NS.GT.7) THEN
C	Double the step size
            H=2.*H
            I2=IJM4-2
            IF(I2.LE.0) I2=I2+7
            DO 4800 I=1,N
              Y(I,IJM1)=Y(I,IJM2)
              Y(I,IJM2)=Y(I,IJM4)
              Y(I,IJM3)=Y(I,I2)
              DY(I,IJM1)=DY(I,IJM2)
              DY(I,IJM2)=DY(I,IJM4)
              DY(I,IJM3)=DY(I,I2)
4800        CONTINUE
            NS=4
          ENDIF

        ELSE
C	If integration has failed, then reduce the step size
          IF(NS.GT.-4.AND.ERR.LT.16.*REPS) THEN
C	H is halved
            DO 5000 I=1,N
              Y1=(45.*Y(I,IJM1)+72.*Y(I,IJM2)+11.*Y(I,IJM3))/128.+
     1            H*(-9.*DY(I,IJM1)+36.*DY(I,IJM2)+3.*DY(I,IJM3))/128.
              Y2=(11.*Y(I,IJM1)+72.*Y(I,IJM2)+45.*Y(I,IJM3))/128.+
     1            H*(-3.*DY(I,IJM1)-36.*DY(I,IJM2)+9.*DY(I,IJM3))/128.
              Y(I,IJ)=Y(I,IJM1)
              Y(I,IJM1)=Y1
              Y(I,IJM3)=Y2
              DY(I,IJ)=DY(I,IJM1)
5000        CONTINUE
            T=T0-H/2
            CALL DIF(T,N,Y(1,IJM1),DY(1,IJM1))
            T=T-H
            CALL DIF(T,N,Y(1,IJM3),DY(1,IJM3))
            NSTP=NSTP+2
            H=H/2.
          ELSE

C	If error is too large or the halving has failed once, then
C	generate fresh starting values with smaller h
            IFLG=IFLG-2
            T=T0
            H=H/8.
            DO 5200 I=1,N
5200        Y(I,1)=Y(I,IJM1)
            CALL STRT4(N,Y,DY,DIF,H,T,REPS,IFLG,TSTEP,NSTP,IER,WK)
            IF(IER.GT.0) THEN
C	If STRT4 fails, then quit
              DO 5300 I=1,N
5300          YF(I)=Y(I,IJM1)
              RETURN
            ENDIF
            IJ=4
            IJM1=3
            IJM2=2
            IJM3=1
            T0=T
            IFLG=IFLG+2
            IFLAG=0
          ENDIF

          NS=-4
          IF(ABS(H/TSTEP).LT.REPS) THEN
C	The step size is too small
            IER=726
            DO 5400 I=1,N
5400        YF(I)=Y(I,IJM1)
            RETURN
          ENDIF

        ENDIF
      ELSE
        IF(IER.GT.0) THEN
C	For integration with fixed H, quit if corrector fails
          DO 5600 I=1,N
5600      YF(I)=Y(I,IJM1)
          RETURN
        ENDIF
        T0=T0+H
        DO 5700 I=1,N
5700    WK(I,2)=WK(I,1)
        IF(TN.LE.T0.EQV.H.GT.0) GO TO 6000
        IF(ABS((TN-T0)/TSTEP).LE.REPS) GO TO 6000
      ENDIF

      IF(NSTP.LT.NMAX) GO TO 3400
C	Quit if the specified number of function evaluations are exceeded
      IER=728
      DO 5800 I=1,N
5800  YF(I)=Y(I,IJ)
      RETURN

6000  TX=(TN-T0)/H
C	Use interpolation to calculate solution at TN
      DO 7000 I=1,N
        DEL(1)=DY(I,IJ)*H
        DEL(2)=Y(I,IJ)-Y(I,IJM1)
        DEL(3)=DY(I,IJM1)*H
        DEL(4)=Y(I,IJM1)-Y(I,IJM2)
        DO 6500 J1=2,3
          DO 6500 J2=4,J1,-1
6500    DEL(J2)=DEL(J2-1)-DEL(J2)
        DEL(4)=DEL(3)-0.5*DEL(4)
        YF(I)=Y(I,IJ)+TX*(DEL(1)+TX*(DEL(2)+(TX+1)*(DEL(3)+
     1        (TX+1)*DEL(4)/2.)))
7000  CONTINUE
      END

C     ---------------------------------------------------

C	To perform one step of integration of ordinary differential equations
C	using a fourth-order Runge-Kutta method 
C
C	N : (input) Number of first order differential equations to be solved
C	T : (input) Initial value of independent variable t, at
C		which initial values are specified. This value is not updated.
C	Y0 : (input) Real array of length N containing the initial
C		values of variables.
C	DY0 : (input) Real array of length N containing the derivatives
C		of Y at the initial point Y0
C	H : (input) The step size to be used for integration
C	Y1 : (output) Real array of length N containing the solution at t=T+H
C	DIF : (input) Name of subroutine to calculate the right hand side
C		of differential equation y'=f(t,y)
C	WK : Real array of length 2N used as scratch space
C
C	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
C	the differential equation. T is the value of independent variable,
C	N is the number of variables, Y is a real array of length N containing
C	the values of variables. DY is a real array of length N which should
C	contain the calculated values of derivatives at (T,Y).
C
C	Required routines : DIF
C	
C	
      SUBROUTINE RK4(N,T,Y0,DY0,H,Y1,DIF,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y0(N),DY0(N),Y1(N),WK(N,2)

      H2=H/2.
      DO 1000 I=1,N
        Y1(I)=H2*DY0(I)
C	y_0+0.5k_1
1000  WK(I,1)=Y0(I)+Y1(I)

      T1=T+H2
      CALL DIF(T1,N,WK,WK(1,2))
      DO 1200 I=1,N
        Y1(I)=Y1(I)+H*WK(I,2)
C	y_0+0.5k_2
1200  WK(I,1)=Y0(I)+H2*WK(I,2)

      CALL DIF(T1,N,WK,WK(1,2))
      DO 1400 I=1,N
        Y1(I)=Y1(I)+H*WK(I,2)
C	y_0+k_3
1400  WK(I,1)=Y0(I)+H*WK(I,2)

      T1=T+H
      CALL DIF(T1,N,WK,WK(1,2))
      DO 1600 I=1,N
1600  Y1(I)=Y0(I)+(Y1(I)+H2*WK(I,2))/3.
      END

C     --------------------------------------------------------------

C	To generate the starting values for fourth order multistep methods
C	in solution of initial value problem in ordinary differential equations.
C	It is used by MSTEP.
C
C	N : (input) Number of first order differential equations to be solved
C	Y : (input/output) Real array of length 4N, the first N elements
C		should contain the initial values of variables. During
C		execution, the solution at the next three points will
C		be stored in this array. 
C	DY : (output) Real array of length 4N, containing the derivatives
C		of Y the points stored in Y.
C	DIF : (input) Name of subroutine to calculate the right hand side
C		of differential equation y'=f(t,y)
C	H : (input/output) Initial guess for the step size. After execution
C		it will contain the step size used by the program
C	T : (input/output) Initial value of independent variable t, at
C		which initial values are specified. After execution it will
C		be set to T+3H if execution is successful.
C	REPS : (input) Required accuracy in each component of the solution.
C		The subroutine only controls local truncation error and hence
C		actual error could be larger
C	IFLG : (input) Integer variable used as flag to decide the
C		type of computation.
C		If IFLG=0 the step size is adjusted to meet the specified accuracy.
C		If IFLG=1 the step size is kept fixed.
C	TSTEP : (input) The size of interval over which integration is requested
C			It is used only for convergence check.
C	NSTP : (input/output) Number of calls to DIF required since the starting
C		of integration. The count is updated by the subroutine.
C	IER : (output) Error parameter; IER=0 implies successful execution
C		IER=724 implies that step size becomes too small
C		IER=725 implies that iteration to generate starting values
C			failed to converge
C	WK : Real array of length 2N used to pass on the value for modifier
C
C	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
C	the differential equation. T is the value of independent variable,
C	N is the number of variables, Y is a real array of length N containing
C	the values of variables. DY is a real array of length N which should
C	contain the calculated values of derivatives at (T,Y).
C
C	Required routines : RK4, DIF
C	
C	
      SUBROUTINE STRT4(N,Y,DY,DIF,H,T,REPS,IFLG,TSTEP,NSTP,IER,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(EPS=1.D-30,SFAC=0.9D0,NIT=10)
      EXTERNAL DIF
      DIMENSION Y(N,4),DY(N,4),WK(N,2)

      T0=T
      IER=0
      DO 4000 IT=1,NIT
        CALL DIF(T0,N,Y(1,1),DY(1,1))
C	Generate y_1
        CALL RK4(N,T0,Y(1,1),DY(1,1),H,Y(1,2),DIF,WK)
        T=T0+H
        CALL DIF(T,N,Y(1,2),DY(1,2))
C	Generate y_2
        CALL RK4(N,T,Y(1,2),DY(1,2),H,Y(1,3),DIF,WK)
        H2=2.*H
C	Calculate y_2 using double step
        CALL RK4(N,T0,Y(1,1),DY(1,1),H2,Y(1,4),DIF,WK)
        NSTP=NSTP+11

C	Estimate the truncation error in y_2
        ERR=0
        DO 2000 I=1,N
          R2=ABS(Y(I,1))+ABS(Y(I,2))+ABS(Y(I,3))+EPS
          R1=ABS(Y(I,3)-Y(I,4))/R2
          IF(R1.GT.ERR) ERR=R1
2000    CONTINUE
        ERR=ERR*TSTEP/H

        IF(ERR.LE.REPS.OR.IFLG.GT.0) THEN
C	Accept the computed values
          T=T+H
          CALL DIF(T,N,Y(1,3),DY(1,3))
C	Generate y_3
          CALL RK4(N,T,Y(1,3),DY(1,3),H,Y(1,4),DIF,WK)
          T=T+H
          CALL DIF(T,N,Y(1,4),DY(1,4))
          NSTP=NSTP+5
C	Store Y_3 for use by modifier
          DO 3000 I=1,N
3000      WK(I,2)=Y(I,4)
          RETURN
        ELSE

C	Reduce the step size
          H=SFAC*H*(REPS/ERR)**0.25D0
          IF(ABS(H/TSTEP).LT.REPS.OR.T.EQ.T0) THEN
C	Step size is too small
            IER=724
            RETURN
          ENDIF
        ENDIF
4000  CONTINUE

C	Iteration to generate starting values fails to converge
      IER=725
      END

C     --------------------------------------------------------------

      SUBROUTINE FCN(NE,X,T,U,UX,UXX,DU)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(*),UX(*),UXX(*),DU(*)

C     THE DIFFERENTIAL EQUATION  DU/DT=5U**4D2U/DX2+20U**3(DU/DX)**2

      DU(1)=5.*U(1)**4*UXX(1)+20.*U(1)**3*UX(1)**2
      END

C     --------------------------------------------------------------

      SUBROUTINE BC(NE,T,X0,XN,U0,UN)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      COMMON/FN/TT,X
      DIMENSION U0(*),UN(*)

C     THE BOUNDARY CONDITIONS : SOLVE THE EQUATION FOR NONLINEAR
C     WAVE TO OBTAIN THE BOUNDARY VALUES

      TT=T
      XL=1.0
      XU=1.D5
      REPS=1.D-6
      AEPS=1.D-6
      X=X0
      CALL BRENT(XL,XU,U0(1),REPS,AEPS,IER,F)
      IF(IER.GT.0) STOP  11

      X=XN
      IF(300.*(XN-1500.*T).GT.8) THEN
        UN(1)=1.+EXP(-300.*(XN-1500.*T))
      ELSE
        XL=1.0
        XU=1.D5
        CALL BRENT(XL,XU,UN(1),REPS,AEPS,IER,F)
        IF(IER.GT.0) PRINT *,IER,XL,XU,UN(1)
        IF(IER.GT.0) STOP  12
      ENDIF
      END

C     ------------------------------

      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/FN/T,XX
      DATA X0,V/0.0,1500/

C     THE IMPLICIT EQUATION DEFINING THE SOLUTION

      IF(X.GT.1) THEN
        F=1.25*(X-1.)**4+20.*(X-1.)**3/3.+15.*(X-1.)**2
     1     +20.*(X-1.)+5.*LOG(X-1.)-V*(V*T-XX+X0)
      ELSE
        F=-5.D8
      ENDIF
      END
