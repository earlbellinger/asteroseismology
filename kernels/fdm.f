C     PROGRAM TO SOLVE BOUNDARY VALUE PROBLEMS IN 
C     ORDINARY DIFFERENTIAL EQUATIONS USING FINITE DIFFERENCE METHOD

      PROGRAM BVP
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL EQN,BCS
      DIMENSION PAR(7446),X(2,2482),WK(125000),IWK(74460),T(2482),
     1          XC(2,2482)

51    FORMAT('   IER =',I4,5X,'NO. OF PTS =',I5,5X,'LAMDA =',1PD14.6/
     1       11X,'T',20X,'SOLUTION',18X,'CORRECTED SOLUTION')
52    FORMAT(I4,2X,1PD14.6,2X,2D14.6,2X,2D14.6)

      N=2482
      M=2
      ML=1
      IFLAG=3
      REPS=1.D-8

      OPEN(13,FILE='PAR',STATUS='OLD',BLANK='ZERO',ERR=21)
      DO 20 I=1,7446
        READ(13,130) PAR(I)
130     FORMAT(E25.10)
20    CONTINUE
21    CLOSE(13)

      OPEN(14,FILE='T',STATUS='OLD',BLANK='ZERO',ERR=31)
      DO 30 I=1,N
        READ(14,140) T(I)
140     FORMAT(E25.10)
30    CONTINUE
31    CLOSE(14)

      DO 500 I=1,N
        X(1,I)=0.
        X(2,I)=0.
500   CONTINUE

      CALL FDM(N,M,ML,PAR,X,XC,T,EQN,BCS,IWK,WK,IFLAG,REPS,IER)
C      WRITE(6,51)IER,N
      DO 1000 I=1,N
        WRITE(6,52) I,T(I),X(1,I),X(2,I),XC(1,I),XC(2,I)
1000  CONTINUE
      END
 
C     ------------------------------------------------------

C	To solve two-point boundary value problem in ordinary differential
C	equations using finite difference method
C
C	N : (input) Number of mesh points to be used.
C	M : (input) Number of first order differential equations in the system
C	ML : (input) Number of boundary conditions at the first boundary t=T(1)
C	PAR : (input) Real array to be passed on to EQN and BCS for calculating
C		the equations and boundary conditions. This array is not used
C		by FDM, but is merely used to pass on any extra parameters
C		that may be required to specify the equations
C	X : (input/output) Real array of length M*N containing the solution.
C		It should contain the initial guess to solution at the time
C		of calling. After execution it will contain the calculated
C		solution. X(i,j) is the ith component of solution at jth mesh point.
C		First dimension of X in calling program must be M
C	XC : (output) Real array of length M*N containing the solution after
C		applying the deferred correction. It is stored in the same
C		format as X.
C	T : (input) Real array of length N containing the mesh points.
C		These points must be in ascending or descending order.
C		For calculating the deferred correction the mesh spacing
C		must be uniform.
C	EQN : (input) Name of subroutine to specify the differential equation
C		to be solved.
C	BCS : (input) Name of subroutine to calculate the boundary conditions
C		at t=T(1) and T(N)
C		Outline of a sample routine for EQN and BCS can be found
C		at the end of this file
C	IWK : Integer array of length M*N used as scratch space
C	WK : Real array of length (M+ML)*2M*(N+1) used as scratch space
C	IFLAG : (input) Integer variable used as a flag to decide the type of
C		computation required.
C		IFLAG=0 implies that equations are nonlinear and deferred
C			correction is to be calculated. In this case X must
C			contain the initial guess and mesh spacing must be uniform
C		IFLAG=1 implies that equations are nonlinear and deferred
C			correction is not required. In this case X must contain
C			the initial guess and mesh spacing could be arbitrary
C		IFLAG=2 implies that equations are linear and deferred correction
C			is required. In this case initial guess is not required,
C			but mesh spacing must be uniform
C		IFLAG=3 implies that equations are linear and deferred correction
C			is not required. In this case initial guess is not
C			required and mesh spacing can be arbitrary
C	REPS : (input) Required accuracy. This is only used to check convergence
C		of Newton's method for solving finite difference equations.
C		The truncation error depends on mesh spacing and is not
C		controlled in this routine.
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=704 implies that N<3, M.LE.ML, or ML.LE.0, in which case
C			no calculations are done
C		IER=734 implies that N<5 and deferred correction is requested.
C			In this case deferred correction is not calculated
C		IER=735 implies that the finite difference matrix is singular
C		IER=736 implies that mesh spacing is not uniform and
C			deferred correction is not calculated
C		IER=737 implies that Newton's iteration for solving the
C			finite difference equations failed to converge.
C
C	SUBROUTINE EQN and BCS must be supplied by the user.
C	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) calculates the right hand
C		sides for differential equations By'_i=f_i(t,y,par)
C		J is the serial number of mesh point at which calculation
C		is required. M is the number of first order differential
C		equations in the system, ML the number of boundary conditions
C		at the first point, PAR is a real array which can be used
C		to pass on any required parameters to define the equations.
C		A and B are real arrays of length (M+ML)*M defining the
C		differential equation B Y'=f(T,PAR,Y) and
C		A(I,K)=dF_I/dY(K) is the Jacobian matrix
C		Y is real array of length M specifying the solution at t=T
C		F is a real array of length M containing the right hand
C		sides of differential equations f_I as defined above.
C		F, A and B must be calculated by the subroutine.
C
C	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) calculates the boundary
C		conditions at both boundaries. M is the number of differential
C		equations, ML is the number of boundary condition at the
C		first mesh point t=T1. PAR is a real array which can be used
C		to pass any required parameters. BC is a real array of length
C		(M+ML)*M which should contain the coefficients of boundary
C		conditions. First ML rows will specify the boundary conditions
C		at t=T1, while remaining rows will specify those at t=TN.
C		G is a real array of length M specifying the boundary conditions
C		G(I)=g_i(T1,PAR,Y1) (I.LE.ML) are the boundary conditions
C		at T1 (g_i=0), while G(I)=g_i(TN,PAR,YN) (I.GT.ML) are the
C		boundary conditions at TN. BC is the Jacobian matrix dg_i/dY(K)
C		Y1 and YN are real arrays of length M specifying the
C		solution at t=T1 and TN.
C
C	Required routines : SETMAT, GAUBLK, EQN, BCS
C
      SUBROUTINE FDM(N,M,ML,PAR,X,XC,T,EQN,BCS,IWK,WK,IFLAG,REPS,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL EQN,BCS
      PARAMETER(NIT=20,EPS=1.D-30)
      DIMENSION PAR(*),X(M,N),XC(M,N),IWK(M,N),WK(M+ML,2*M,N+1),T(N)

      IF(N.LT.3.OR.M.LE.ML.OR.ML.LE.0) THEN
        IER=704
        RETURN
      ENDIF

      DO 1100 I=1,NIT
C	Set up the matrix of finite difference equations
        CALL SETMAT(N,M,ML,WK,WK(1,1,N+1),X,XC,T,PAR,EQN,BCS)
        IFLG=0
C	Solve the system of linear equations
        CALL GAUBLK(N,M,ML,WK,IFLG,DET,IDET,IWK,XC,IER)
        IF(IER.GT.0) RETURN

C	Checking for convergence
        RERR=0.0
        DO 1000 J=1,N
          J1=J+1
          IF(J.EQ.N) J1=J-1
          DO 1000 K=1,M
            XJ=X(K,J)+XC(K,J)
            R2=ABS(XJ)+ABS(X(K,J)-X(K,J1))+EPS
            RE=ABS(XC(K,J)/R2)
            IF(RE.GT.RERR) RERR=RE
            X(K,J)=XJ
1000    CONTINUE
        IF(RERR.LT.REPS.OR.IFLAG.GT.1) GO TO 1150
1100  CONTINUE

C	Newton's iteration fails to converge
      IER=737
      RETURN

1150  IF(IFLAG.EQ.1.OR.IFLAG.GT.2) RETURN
      IF(N.LE.4) THEN
        IER=734
        RETURN
      ENDIF

C	Calculate the deferred correction
      HH=T(2)-T(1)
      DO 2000 J=1,N-1
        TJ=0.5*(T(J)+T(J+1))
        H=T(J+1)-T(J)
        IF(ABS(H-HH).GT.1.D-4*ABS(HH)) THEN
C	Mesh spacing is not uniform, hence quit
          IER=736
          RETURN
        ENDIF

        JJ=J
        CALL EQN(JJ,M,ML,PAR,WK(1,1,N+1),WK(1,M+1,N+1),X(1,J),XC(1,J+1),
     1           TJ)
        DO 1200 I=1,M
          IF(J.EQ.1) THEN
            WK(M+1,I,N+1)=(-2*X(I,1)+7*X(I,2)-9*X(I,3)+5*X(I,4)-X(I,5))
     1                     /24.
            WK(M+1,I+M,N+1)=(43*X(I,1)-112*X(I,2)+102*X(I,3)-40*X(I,4)+
     1                       7*X(I,5))*H/192.
          ELSE IF(J.EQ.N-1) THEN
            WK(M+1,I,N+1)=(2*X(I,N)-7*X(I,N-1)+9*X(I,N-2)-5*X(I,N-3)
     1                     +X(I,N-4))/24.
            WK(M+1,I+M,N+1)=(43*X(I,N)-112*X(I,N-1)+102*X(I,N-2)-
     1                       40*X(I,N-3)+7*X(I,N-4))*H/192.
          ELSE
            WK(M+1,I,N+1)=(-X(I,J-1)+3.*X(I,J)-3.*X(I,J+1)+X(I,J+2))/24.
            WK(M+1,I+M,N+1)=(X(I,J-1)-X(I,J)-X(I,J+1)+X(I,J+2))*H/16.
          ENDIF
1200    CONTINUE

C	Set up the RHS for deferred correction
        DO 1500 I=1,M
          XD=0.0
          DO 1400 K=1,M
           XD=XD+WK(I,K+M,N+1)*WK(M+1,K,N+1)-WK(I,K,N+1)*WK(M+1,K+M,N+1)
1400      CONTINUE
          XC(ML+I,J)=XD
1500    CONTINUE
2000  CONTINUE

      DO 2200 I=1,ML
2200  XC(I,1)=0.0
      DO 2300 I=ML+1,M
2300  XC(I,N)=0.0

C	Calculate deferred correction
      CALL GAUBLK(N,M,ML,WK,IFLG,DET,IDET,IWK,XC,IER)
      DO 2500 J=1,N
        DO 2500 I=1,M
2500  XC(I,J)=XC(I,J)+X(I,J)
      END

C	--------------------------------
C	
C	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T)
C       IMPLICIT REAL*8(A,B,D-H,O-Z)
C	DIMENSION A(M+ML,M),B(M+ML,M),Y(M),PAR(*),F(M)
C
C	DO 1000 I=1,M
C	F(I)=f_i(T,PAR,Y)
C	  DO 1000 K=1,M
C	    A(K,I)=df_k/dY(i)
C	    B(K,I)=b_{ki}(T,PAR)
C 1000  CONTINUE
C       END
C
C	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN)
C       IMPLICIT REAL*8(A,B,D-H,O-Z)
C	DIMENSION PAR(*),BC(M+ML,M),G(M),Y1(M),YN(M)
C
C	DO 1000 I=1,M
C	  IF(I.LE.ML) THEN
C	    G(I)=g_i(T1,PAR,Y1)
C	  ELSE
C	    G(I)=g_i(TN,PAR,YN)
C         ENDIF
C	  DO 1000 K=1,M
C	    BC(I,K)=dg_I/dY(K)
C 1000  CONTINUE
C       END

C     ----------------------------------------------------------

C	To solve a system of linear equations arising from finite difference
C	approximation of ordinary differential equations
C
C	N : (input) Number of mesh points used.
C	M : (input) Number of first order differential equations in the system
C	ML : (input) Number of boundary conditions at the first boundary t=T(1)
C	A : (input/output) Real array of length (M+ML)*2M*N containing
C		the matrix of equations. After execution it will contain
C		the triangular decomposition
C	IFLG : (input/output) Integer variable used as a flag to decide
C		the nature of computation.
C		If IFLG=0 the triangular decomposition, determinant and
C			solution of equations is computed. IFLG is set to 2
C		If IFLG=1 only the triangular decomposition and determinant
C			are calculated and IFLG is set to 2
C		If IFLG=2 then it is assumed that triangular decomposition
C			is already done and available in A and only the solution
C			of system of equations is solved
C	DET : (output) Scaled value of determinant of the matrix
C	IDET : (output) exponent of determinant, the value of determinant
C		is DET*2**IDET
C	INC : (input/output) Integer array of length M*N containing the
C		information about interchanges used by Gaussian elimination.
C		It is calculated if IFLG=0 or 1, while for IFLG=2 it must
C		be supplied from previous calculation.
C	X : (input/output) Real array of length M*N containing the right hand
C		side of equation at input. After execution it will be overwritten
C		by the solution if IFLG=0 or 2. For IFLG=1 it is not used.
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=735 implies that the matrix is singular
C
C	Required routines : None

      SUBROUTINE GAUBLK(N,M,ML,A,IFLG,DET,IDET,INC,X,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+ML,2*M,N),INC(M,N),X(M,N)

      IF(IFLG.LE.1) THEN
        IDET=0
        DET=1.
C	The number of rows in each block of the matrix
        MR=M+ML
C	The number of columns in each block of the matrix
        MC=2*M
        IER=735

        DO 3000 J=1,N
          IF(J.EQ.N) THEN
            MR=M
            MC=M
          ENDIF

          DO 2000 K=1,MIN(M,MR-1)
            RMAX=ABS(A(K,K,J))
            KMAX=K
C	Find the pivot
            DO 1200 KI=K+1,MR
              R1=ABS(A(KI,K,J))
              IF(R1.GT.RMAX) THEN
                RMAX=R1
                KMAX=KI
              ENDIF
1200        CONTINUE
            INC(K,J)=KMAX

            IF(KMAX.NE.K) THEN
C	exchange rows K and KMAX
              DET=-DET
              DO 1300 KI=K,MC
                AT=A(K,KI,J)
                A(K,KI,J)=A(KMAX,KI,J)
1300          A(KMAX,KI,J)=AT
            ENDIF

            DET=DET*A(K,K,J)
C	If the pivot is zero, then quit
            IF(A(K,K,J).EQ.0.0) RETURN

C	Gaussian elimination
            DO 1500 KI=K+1,MR
              A(KI,K,J)=A(KI,K,J)/A(K,K,J)
              DO 1500 KJ=K+1,MC
                A(KI,KJ,J)=A(KI,KJ,J)-A(KI,K,J)*A(K,KJ,J)
1500        CONTINUE
2000      CONTINUE

          IF(DET.NE.0.0) THEN
C	Scale the determinant if necessary
2350        IF(ABS(DET).GT.32.) THEN
              DET=DET*0.03125D0
              IDET=IDET+5
              GO TO 2350
            ENDIF

2370        IF(ABS(DET).LT.0.03125D0) THEN
              DET=DET*32.
              IDET=IDET-5
              GO TO 2370
            ENDIF
          ENDIF

C	Copy the overlapping elements into the next block
          IF(J.LT.N) THEN
            DO 2600 K=1,ML
              DO 2600 KI=1,M
2600        A(K,KI,J+1)=A(K+M,KI+M,J)
          ENDIF
3000    CONTINUE
        INC(M,N)=M
        DET=DET*A(M,M,N)
        IF(A(M,M,N).EQ.0.0) RETURN
        IER=0

        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

C	Solve the system of linear equations
      IER=0
      MR=M+ML
      DO 3100 J=1,N
        IF(J.EQ.N) MR=M
        DO 3100 K=1,MIN(M,MR-1)
          KK=INC(K,J)
          IF(K.NE.KK) THEN
C	exchange the corresponding elements of RHS
            XT=X(K,J)
            X(K,J)=X(KK,J)
            X(KK,J)=XT
          ENDIF

C	Gaussian elimination
          DO 2800 L=K+1,MR
2800      X(L,J)=X(L,J)-A(L,K,J)*X(K,J)
3100  CONTINUE

C	back-substitution
      MC=M
      DO 3500 J=N,1,-1
        DO 3300 K=M,1,-1
          D1=X(K,J)
            DO 3200 L=MC,K+1,-1
3200        D1=D1-X(L,J)*A(K,L,J)
3300    X(K,J)=D1/A(K,K,J)
        MC=2*M
3500  CONTINUE
      END

C     -------------------------------------------------------

C	To setup the matrix of system of linear equations arising from
C	finite difference approximation of ordinary differential equations
C	This routine is called by FDM or GEVP
C
C	N : (input) Number of mesh points used.
C	M : (input) Number of first order differential equations in the system
C	ML : (input) Number of boundary conditions at the first boundary t=T(1)
C	A : (output) Real array of length (M+ML)*2M*N containing
C		the matrix of equations.
C	BC : (output) Real array of length (M+ML)*(M+1) containing the
C		coefficients of boundary conditions
C	X : (input) Real array of length M*N containing the current approximation
C		to the solution.
C	XC : (output) Real array of length M*N containing the right hand
C		side of finite difference equations calculated by the routine
C	T : (input) Real array of length N containing the mesh points.
C	PAR : (input) Real array containing the parameters to be passed
C		on to subroutines EQN and BCS. This array is not used by
C		the subroutine.
C	EQN : (input) Name of subroutine to calculate the equation matrix
C	BCS : (input) Name of subroutine to calculate the boundary conditions
C
C	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) and
C	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) must be supplied by the user.
C	The form of these routines is described in documentation for FDM or GEVP
C
C	Required routines : EQN, BCS
C
      SUBROUTINE SETMAT(N,M,ML,A,BC,X,XC,T,PAR,EQN,BCS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+ML,2*M,N),BC(M+ML,M+1),X(M,N),XC(M,N),PAR(*),T(N)

C	Loop over the mesh points
      DO 1500 K=1,N-1
C	t_{k+1/2}
        TK=0.5*(T(K)+T(K+1))
        H=T(K+1)-T(K)
        DO 800 I=1,M
C	y_{k+1/2}
800     BC(I,1)=0.5*(X(I,K)+X(I,K+1))
        KK=K
C	Calculate the equation matrix at t_{k+1/2}
        CALL EQN(KK,M,ML,PAR,A(1,1,K),A(1,M+1,K),BC,XC(ML+1,K),TK)

C	Setup the RHS of finite difference equations
        DO 1100 J=1,M
          XK=XC(ML+J,K)*H
          DO 1000 I=1,M
1000      XK=XK-A(J,M+I,K)*(X(I,K+1)-X(I,K))
          XC(ML+J,K)=XK
1100    CONTINUE

C	Setup the finite difference matrix
        DO 1500 J=1,M
          DO 1200 I=M,1,-1
            A(ML+I,J,K)=-A(I,J+M,K)-0.5*H*A(I,J,K)
1200      A(ML+I,J+M,K)=A(I,J+M,K)-0.5*H*A(I,J,K)
          DO 1500 I=1,ML
1500  A(I,J+M,K)=0.0

C	The boundary conditions
      CALL BCS(M,ML,PAR,BC,BC(1,M+1),T(1),T(N),X(1,1),X(1,N))

C	Boundary conditions at the first boundary
      DO 3000 I=1,ML
        XC(I,1)=-BC(I,M+1)
        DO 3000 J=1,M
3000  A(I,J,1)=BC(I,J)

C	Boundary conditions at the second boundary
      DO 4000 I=ML+1,M
        XC(I,N)=-BC(I,M+1)
        DO 4000 J=1,M
4000  A(I,J,N)=BC(I,J)
      END

C     -------------------------------------------

      SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+ML,M),B(M+ML,M),Y(M),PAR(*),F(*)

C     SOLVES By' = vec f(t,y)
C     THE ARRAY A STORES THE JACOBIAN

C     WE HAVE THE DIFFERENTIAL EQUATION
C     f(1) = -(K_rho_Gamma1 + (Gamma_1,P + Gamma_1,rho) * K_Gamma1_rho) - U y_2
C     f(1) = U y_2 -(K_rho_Gamma1 + (Gamma_1,rho + Gamma_1,P) * K_rho_Gamma1)
C     f(2) = - V y_1 - U y_2 

C     WITH THE JACOBIAN
C     df_1/dy_1 = 0
C     df_1/dy_2 = -U
C     df_2/dy_1 = -V
C     df_2/dy_2 = U

C     LET
C     PAR(1..N) = (K_rho_Gamma1 - (Gamma_1,rho + Gamma_1,P) * K_rho_Gamma1)(T)
C     PAR(N+1..2*N) = U(T)
C     PAR(2*N+1..3*N) = V(T)

      F(1)=PAR(N+J)*Y(2)+PAR(J)
      F(2)=PAR(2*N+J)*Y(1)+PAR(N+J)*Y(2)
      
      F(1)=-F(1)/T
      F(2)=F(2)/T

      A(1,1)=0.
      A(1,2)=-PAR(N+J)/T
      A(2,1)=PAR(2*N+J)/T
      A(2,2)=PAR(N+J)/T

      B(1,1)=1.
      B(1,2)=0.
      B(2,1)=0.
      B(2,2)=1.

C     THE DIFFERENTIAL EQ.  Y1'=Y2,  Y2'=LAMDA*SINH(LAMDA*Y1)
C     WITH LAMDA=PAR(1)

C      DO 1000 I=1,M
C        F(I)=0.0
C        DO 1000 K=1,M
C          A(K,I)=0.0
C1000  B(K,I)=0.0

C      B(1,1)=1.
C      A(1,2)=1.
C      B(2,2)=1.
C      A(2,1)=PAR(1)**2*COSH(PAR(1)*Y(1))
      
C      F(1)=Y(2)
C      F(2)=PAR(1)*SINH(PAR(1)*Y(1))
      END

C     ----------------------------------------------

      SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PAR(*),BC(M+ML,M),G(M),Y1(*),YN(*)

C     JACOBIAN OF BOUNDARY CONDITIONS
      BC(1,1)=1.
      BC(1,2)=0.
      BC(2,1)=0.
      BC(2,2)=1.

C     BOUNDARY CONDITIONS    Y1(T=T1)=0,  YN(T=T2)=0
      G(1)=Y1(1)-1
      G(2)=YN(2)-1
      
      
C     THE BOUNDARY CONDITIONS  X1(T=T1)=0,  X1(T=T2)=1

C      BC(1,1)=1.0
C      BC(1,2)=0.0
C      BC(2,1)=1.0
C      BC(2,2)=0.0

C      G(1)=X1(1)
C      G(2)=X2(1)-1.
      END

