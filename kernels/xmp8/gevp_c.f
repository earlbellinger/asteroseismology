C     PROGRAM FOR SOLVING GENERALISED EIGENVALUE PROBLEM IN
C     ORDINARY DIFFERENTIAL EQUATIONS USING FINITE DIFFERENCE METHOD
C     THIS IS A COMPLEX VERSION OF GEVP, AND CAN BE USED FOR
C     FINDING COMPLEX EIGENVALUES AND CORRESPONDING EIGENFUNCTIONS

      PROGRAM EIGEN
      IMPLICIT REAL*8(H,O,Q-T)
      IMPLICIT COMPLEX*16(A-G,P,U-Z)
      EXTERNAL EQN,BCS,EQND,BCSD
      PARAMETER(M0=2,N0=1001,ML0=1)
      DIMENSION PAR(10),T(N0),X(M0,N0),XC(M0,N0),IWK(M0,N0),
     1          WK(M0+ML0,2*M0,N0+1)

C     EXAMPLE 11.14 : SPHEROIDAL HARMONICS

51    FORMAT(I5,2X,1PD14.6,2X,2D14.6,2X,2D14.6)
52    FORMAT('   IER =',I4,5X,'NO. OF PTS =',I5/
     1    5X,'EIGENVALUE =',1P2D14.6,5X,'CORRECTED EIGENVALUE =',2D14.6/
     2    13X,1HT,18X,'EIGENFUNCTION')
53    FORMAT('   CSQ =',1P2D14.6,5X,'M =',2D14.6/5X,'E0 =',2D14.6)

      M=2
      ML=1
      IFLAG=2
      REPS=1.D-8
      RMAX=1.D5

100   PRINT *,'TYPE E0=INITIAL APP. TO EIGENVALUE, N=NO. OF PTS'
      PRINT *,'            (QUITS WHEN N.LE.0)'
      READ *,E0,N
      IF(N.LE.0) STOP
      PRINT *,'TYPE C SQUARE, M'
      READ *,(PAR(I),I=2,3)
      WRITE(6,53) PAR(2),PAR(3),E0

C     SET UP MESH WITH UNIFORM SPACING

      H=1.D0/(N-1)
      DO 500 I=1,N
500   T(I)=(I-1)*H
      CALL GEVP_C(N,M,ML,PAR,X,XC,T,E0,EQN,BCS,EQND,BCSD,IWK,WK,
     1          IFLAG,REPS,RMAX,IER)

      WRITE(6,52) IER,N,PAR(1),E0
      DO 1000 I=1,N,10
        WRITE(6,51) I,T(I),X(1,I),X(2,I)
1000  CONTINUE
      GO TO 100
      END
 
C     ------------------------------------------------------

C	To solve a system of linear equations arising from finite difference
C	approximation of ordinary differential equations
C	Version for complex matrix
C
C	N : (input) Number of mesh points used.
C	M : (input) Number of first order differential equations in the system
C	ML : (input) Number of boundary conditions at the first boundary t=T(1)
C	A : (input/output) Complex array of length (M+ML)*2M*N containing
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
C			of system of equations is calculated.
C	DET : (output) Value of determinant of the matrix
C	IDET : (output) exponent of determinant, the value of determinant
C		is DET*2**IDET
C	INC : (input/output) Integer array of length M*N containing the
C		information about interchanges used by Gaussian elimination.
C		It is calculated if IFLG=0 or 1, while for IFLG=2 it must
C		be supplied from previous calculation.
C	X : (input/output) Complex array of length M*N containing the right hand
C		side of equation at input. After execution it will be overwritten
C		by the solution if IFLG=0 or 2. For IFLG=1 it is not used.
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=735 implies that the matrix is singular
C
C	Required routines : None

      SUBROUTINE GAUBLK_C(N,M,ML,A,IFLG,DET,IDET,INC,X,IER)
C      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT REAL*8(H,O,P,R,T)
      IMPLICIT COMPLEX*16(A-G,U-Z)
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

C     ----------------------------------------------------------

C	To solve a generalised eigenvalue problem in ordinary differential
C	equations using finite difference method
C	Version for complex eigenvalue
C
C	N : (input) Number of mesh points to be used.
C	M : (input) Number of first order differential equations in the system
C	ML : (input) Number of boundary conditions at the first boundary t=T(1)
C	PAR : (input/output) Complex array to be passed on to EQN, EQND, BCS and BCSD
C		for specifying the equations and boundary conditions. PAR(1)
C		is used for passing the eigenvalue and hence should not be
C		used for any other variable. After execution the eigenvalue
C		will be available in PAR(1). Other elements of this array are
C		not used by GEVP_C, but are merely used to pass on any extra parameters
C		that may be required to specify the equations
C	X : (output) Complex array of length M*N containing the eigenvector.
C		X(i,j) is the ith component of solution at jth mesh point.
C		First dimension of X in calling program must be M
C	XC : (output) Complex array of length M*N containing the left eigenvector.
C		It is stored in the same format as X.
C	T : (input) Real array of length N containing the mesh points.
C		These points must be in ascending or descending order.
C		For calculating the deferred correction the mesh spacing
C		must be uniform.
C	E0 : (output) The calculated eigenvalue including deferred correction
C		PAR(1) contains the uncorrected eigenvalue. If deferred
C		correction is not applied E0 is not calculated.
C	EQN : (input) Name of subroutine to specify the differential equation
C		to be solved.
C	BCS : (input) Name of subroutine to calculate the boundary conditions
C		at t=T(1) and T(N)
C	EQND : (input) Name of subroutine to calculate derivative of equation
C		matrix with respect to the eigenvalue.
C	BCS : (input) Name of subroutine to calculate derivative of the boundary
C		conditions with respect to the eigenvalue
C	IWK : Integer array of length M*N used as scratch space
C	WK : Complex array of length (M+ML)*2M*(N+1) used as scratch space
C	IFLAG : (input) Integer variable used as a flag to decide the type of
C		computation required.
C		IFLAG=0 implies that only the eigenvalue is calculated.
C		IFLAG=1 implies that both eigenvalue and eigenvector are
C			calculated.
C		IFLAG=2 implies that first order correction to eigenvalue
C			is also calculated in addition to eigenvector
C			In this case the mesh spacing must be uniform
C	REPS : (input) Required accuracy. This is only used to check convergence
C		of Muller's method for finding zeros of the determinant
C		The truncation error depends on mesh spacing and is not
C		controlled in this routine.
C	RMX : (input) Expected maximum limit on the eigenvalue. This parameter
C		is passed on to MULER2 and iteration is terminated if the
C		magnitude of estimated eigenvalue exceeds RMX.
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=704 implies that N<3, M.LE.ML, or ML.LE.0, in which case
C			no calculations are done
C		IER=734 implies that N<5 and deferred correction is requested.
C			In this case deferred correction is not calculated
C		IER=735 implies that the finite difference matrix is singular
C			while calculating the eigenvector. In this case
C			eigenvector is not calculated, but eigenvalue is
C			available in PAR(1).
C		IER=736 implies that mesh spacing is not uniform and
C			deferred correction is not calculated
C		IER=738 implies that eigenvector vanishes.
C		IER=739 implies that inverse iteration for calculating the
C			eigenvector failed to converge.
C		IER=740 implies that inverse iteration for calculating the
C			left eigenvector failed to converge.
C
C	SUBROUTINE EQN, EQND, BCS and BCSD must be supplied by the user.
C	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) calculates the right hand
C		sides for differential equations B y'=A y
C		J is the serial number of mesh point at which calculation
C		is required. M is the number of first order differential
C		equations in the system, ML the number of boundary conditions
C		at the first point, PAR is a Complex array which can be used
C		to pass on any required parameters to define the equations.
C		PAR(1) is the eigenvalue.
C		A and B are Complex arrays of length (M+ML)*M defining the
C		differential equation B Y'= A Y.
C		Y is Complex array of length M specifying the solution at t=T
C		F is a Complex array of length M which should be set to zero.
C		F, A and B must be calculated by the subroutine.
C
C	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) calculates the boundary
C		conditions at both boundaries. M is the number of differential
C		equations, ML is the number of boundary condition at the
C		first mesh point t=T1. PAR is a Complex array which can be used
C		to pass any required parameters. PAR(1) is the eigenvalue.
C		BC is a Complex array of length (M+ML)*M which should contain
C		the coefficients of boundary conditions, BC Y =0.
C		First ML rows will specify the boundary conditions
C		at t=T1, while remaining rows will specify those at t=TN.
C		G is a Complex array of length M which should be set to zero.
C		Y1 and YN are Complex arrays of length M specifying the
C		solution at t=T1 and TN.
C
C	SUBROUTINE EQND(J,M,ML,PAR,A,B,T) calculates the derivatives of
C		matrices A and B (calculated by EQN) with respect to the
C		eigenvalue (PAR(1)). 
C		J is the serial number of mesh point at which calculation
C		is required. M is the number of first order differential
C		equations in the system, ML the number of boundary conditions
C		at the first point, PAR is a Complex array which can be used
C		to pass on any required parameters to define the equations.
C		PAR(1) is the eigenvalue.
C		A and B are Complex arrays of length (M+ML)*M which should
C		give the derivatives of arrays A and B calculated by EQN.
C		T is the value of t at which the derivatives are required.
C
C	SUBROUTINE BCSD(M,ML,PAR,BC,T1,TN) calculates the derivative of
C		matrix BC (calculated by BCS) with respect to the eigenvalue PAR(1).
C		M is the number of differential equations,
C		ML is the number of boundary condition at the
C		first mesh point t=T1. PAR is a Complex array which can be used
C		to pass any required parameters. PAR(1) is the eigenvalue.
C		BC is a Complex array of length (M+ML)*M which should contain
C		the derivative of the matrix BC calculated by BCS.
C		First ML rows will specify the boundary conditions
C		at t=T1, while remaining rows will specify those at t=TN.
C
C	Required routines : SETMAT_C, GAUBLK_C, MULER2, EQN, EQND, BCS, BCSD
C
      SUBROUTINE GEVP_C(N,M,ML,PAR,X,XC,T,E0,EQN,BCS,EQND,BCSD,
     1                IWK,WK,IFLAG,REPS,RMX,IER)
C      IMPLICIT REAL*8(A,B,D-H,O-Z)
C      IMPLICIT COMPLEX*16(C)
C	For complex eigenvalues use the following statements
      IMPLICIT REAL*8(H,O,R,S,T)
      IMPLICIT COMPLEX*16(A-G,P,U-Z)
      EXTERNAL EQN,BCS
      PARAMETER(NIT=20,EPS=1.D-30)
      DIMENSION PAR(*),T(N),X(M,N),XC(M,N),IWK(M,N),WK(M+ML,2*M,N+1),
     1          CZERO(1)

      IF(N.LT.3.OR.M.LE.ML.OR.ML.LE.0) THEN
        IER=704
        RETURN
      ENDIF

      RAPS=0.1D0*REPS
      DE=MAX(1.D-3*ABS(E0),100.*REPS)
      CX1=E0+DE
      CX2=E0-DE
      CX3=E0
      NZ=0
      IER=0

C	First call to MULER2
1000  CALL MULER2(CX1,CX2,CX3,REPS,RAPS,IER,CF,CX,IDET,NZ,CZERO,RMX)
      IF(IER.LT.0) THEN
C	calculate the determinant
        PAR(1)=CX
        CALL SETMAT_C(N,M,ML,WK,WK(1,1,N+1),X,XC,T,PAR,EQN,BCS)
        IFLG=1
        CALL GAUBLK_C(N,M,ML,WK,IFLG,DET,IDET,IWK,XC,IER1)
        CF=DET
C	Call MULER2 again
        GO TO 1000
      ENDIF

C	the eigenvalue
      PAR(1)=CX3
      IF(IER.GT.100.OR.IFLAG.EQ.0) RETURN
      IF(IER1.GT.0) THEN
        IER=735
        RETURN
      ENDIF

C	Initialise the vector for inverse iteration
      DO 1200 I=1,N
        DO 1200 J=1,M
          X(J,I)=1.
1200  XC(J,I)=1.
      IFLG=2

C	Loop for inverse iteration to calculate the eigenvector
      DO 2500 I=1,NIT
        CALL GAUBLK_C(N,M,ML,WK,IFLG,DET,IDET,IWK,X,IER)

C	Normalising the eigenvector
        RMAX=0.0
        DO 1500 J=1,N
          DO 1500 K=1,M
            IF(ABS(X(K,J)).GT.RMAX) THEN
              RMAX=ABS(X(K,J))
              XMAX=X(K,J)
            ENDIF
1500    CONTINUE
        IF(RMAX.EQ.0.0) THEN
C	eigenvector vanishes, hence quit
          IER=738
          RETURN
        ENDIF

        DO 1800 J=1,N
          DO 1800 K=1,M
1800    X(K,J)=X(K,J)/XMAX

C	convergence check
        IF(RMAX*REPS.GT.1.0) GO TO 3000
2500  CONTINUE
      IER=739
      RETURN

3000  IF(IFLAG.LE.1) RETURN
      IF(N.LE.4) THEN
        IER=734
        RETURN
      ENDIF

C	Loop for inverse iteration for the left eigenvector
      DO 5000 IT=1,NIT

C	Forward substitution for U^T
        DO 3500 I=1,N
          DO 3500 J=1,M
            IF(I.GT.1) THEN
              DO 3200 K=1,M
3200          XC(J,I)=XC(J,I)-XC(K,I-1)*WK(K,J+M,I-1)
            ENDIF
              DO 3400 K=1,J-1
3400          XC(J,I)=XC(J,I)-XC(K,I)*WK(K,J,I)
3500    XC(J,I)=XC(J,I)/WK(J,J,I)

C	Back substitution for L^T
        MR=M
        DO 4000 I=N,1,-1
          IF(I.EQ.N-1) MR=M+ML
          DO 4000 J=MIN(M,MR-1),1,-1
            DO 3900 K=J+1,MR
              K0=K
              I1=I
C	Checking for backward interchange
              DO 3600 K1=J+1,K
                IF(K1.GT.M) K0=K-M
                IF(IWK(K1,I).EQ.K0) GO TO 3900
3600          CONTINUE

              K1=K
3700          IF(K1.GT.M) THEN
                K1=K1-M
                I1=I1+1
              ENDIF
              IF(IWK(K1,I1).GT.K1) THEN
                K0=IWK(K1,I1)
                K2=K0
                DO 3800 KK=K1+1,K2
                  IF(KK.GT.M) K0=K2-M
                  IF(IWK(KK,I1).EQ.K0) THEN
                    K1=KK
                    GO TO 3900
                  ENDIF
3800            CONTINUE
                K1=K2
                IF(I1.LE.N) GO TO 3700
C	Hopefully, this statement will not be executed
                STOP
              ENDIF
3900        XC(J,I)=XC(J,I)-XC(K1,I1)*WK(K,J,I)
4000    CONTINUE

C	Exchanges due to interchange matrix N
        DO 4200 I=N,1,-1
          DO 4200 J=M,1,-1
            JJ=IWK(J,I)
            IF(JJ.NE.J) THEN
              XT=XC(J,I)
              XC(J,I)=XC(JJ,I)
              XC(JJ,I)=XT
            ENDIF
4200    CONTINUE

        RMAX=0.0
        DO 4400 I=1,N
          DO 4400 J=1,M
            IF(RMAX.LT.ABS(XC(J,I))) RMAX=ABS(XC(J,I))
4400    CONTINUE
        DO 4600 I=1,N
          DO 4600 J=1,M
4600    XC(J,I)=XC(J,I)/RMAX

        IF(RMAX*REPS.GT.1.0) GO TO 6000
5000  CONTINUE
      IER=740
      RETURN

C	calculating the deferred correction
6000  HH=T(2)-T(1)
      ENUM=0.0
      EDEN=0.0
      DO 7000 J=1,N-1
        TJ=0.5*(T(J)+T(J+1))
        H=T(J+1)-T(J)
        IF(ABS(H-HH).GT.1.D-4*ABS(HH)) THEN
          IER=736
          RETURN
        ENDIF

        JJ=J
        CALL EQN(JJ,M,ML,PAR,WK(1,1,2),WK(1,M+1,2),X(1,J),WK(1,1,3),TJ)
        CALL EQND(JJ,M,ML,PAR,WK(1,1,3),WK(1,M+1,3),TJ)

        DO 6200 I=1,M
          IF(J.EQ.1) THEN
            WK(I,1,1)=(-2*X(I,1)+7*X(I,2)-9*X(I,3)+5*X(I,4)-X(I,5))/24.
            WK(I,2,1)=(43*X(I,1)-112*X(I,2)+102*X(I,3)-40*X(I,4)+
     1                 7*X(I,5))*H/192.
          ELSE IF(J.EQ.N-1) THEN
            WK(I,1,1)=(2*X(I,N)-7*X(I,N-1)+9*X(I,N-2)-5*X(I,N-3)+
     1                  X(I,N-4))/24.
            WK(I,2,1)=(43*X(I,N)-112*X(I,N-1)+102*X(I,N-2)-
     1                 40*X(I,N-3)+7*X(I,N-4))*H/192.
          ELSE
            WK(I,1,1)=(-X(I,J-1)+3.*X(I,J)-3.*X(I,J+1)+X(I,J+2))/24.
            WK(I,2,1)=(X(I,J-1)-X(I,J)-X(I,J+1)+X(I,J+2))*H/16.
          ENDIF
6200    CONTINUE

        DO 6500 I=1,M
          XN=0.0
          XD=0.0
          DO 6400 K=1,M
            XD=XD+WK(I,K+M,3)*(X(K,J+1)-X(K,J))-0.5*H*WK(I,K,3)*
     1         (X(K,J)+X(K,J+1))
            XN=XN+WK(I,K+M,2)*WK(K,1,1)-WK(I,K,2)*WK(K,2,1)
6400      CONTINUE
          ENUM=ENUM+XC(ML+I,J)*XN
          EDEN=EDEN+XC(ML+I,J)*XD
6500    CONTINUE
7000  CONTINUE

      CALL BCSD(M,ML,PAR,WK(1,1,1),T(1),T(N))
      DO 7200 I=1,ML
        XD=0.0
        DO 7100 K=1,M
7100    XD=XD+WK(I,K,1)*X(K,1)
7200  EDEN=EDEN+XD*XC(I,1)
      DO 7500 I=ML+1,M
        XD=0.0
        DO 7400 K=1,M
7400    XD=XD+WK(I,K,1)*X(K,N)
7500  EDEN=EDEN+XD*XC(I,N)

      DEF=ENUM/EDEN
      E0=PAR(1)+DEF
      END
C	--------------------------------
C	
C	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T)
C       IMPLICIT REAL*8(H,O,R,S,T)
C       IMPLICIT COMPLEX*16(A-G,P,U-Z)
C	DIMENSION A(M+ML,M),B(M+ML,M),Y(M),PAR(*),F(M)
C
C	DO 1000 I=1,M
C	F(I)=0.0
C	  DO 1000 K=1,M
C	    A(K,I)=a_{ki}(T,PAR)
C	    B(K,I)=b_{ki}(T,PAR)
C 1000  CONTINUE
C       END
C
C	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN)
C       IMPLICIT REAL*8(H,O,R,S,T)
C       IMPLICIT COMPLEX*16(A-G,P,U-Z)
C	DIMENSION PAR(*),BC(M+ML,M),G(M),Y1(M),YN(M)
C
C	DO 1000 I=1,M
C	    G(I)=0.0
C	  DO 1000 K=1,M
C	  IF(I.LE.ML) THEN
C	    BC(I,K)=bc_{ik}(T1,PAR)
C	  ELSE
C	    BC(I,K)=bc_{ik}(TN,PAR)
C         ENDIF
C 1000  CONTINUE
C       END
C	
C	SUBROUTINE EQND(J,M,ML,PAR,A,B,T)
C       IMPLICIT REAL*8(H,O,R,S,T)
C       IMPLICIT COMPLEX*16(A-G,P,U-Z)
C	DIMENSION A(M+ML,M),B(M+ML,M),PAR(*)
C
C	DO 1000 I=1,M
C	  DO 1000 K=1,M
C	    A(K,I)=da_{ki}/dPAR(1) (T,PAR)
C	    B(K,I)=db_{ki}/dPAR(1) (T,PAR)
C 1000  CONTINUE
C       END
C
C	SUBROUTINE BCSD(M,ML,PAR,BC,T1,TN)
C       IMPLICIT REAL*8(H,O,R,S,T)
C       IMPLICIT COMPLEX*16(A-G,P,U-Z)
C	DIMENSION PAR(*),BC(M+ML,M)
C
C	DO 1000 I=1,M
C	  DO 1000 K=1,M
C	  IF(I.LE.ML) THEN
C	    BC(I,K)=d bc_{ik}/dPAR(1) (T1,PAR)
C	  ELSE
C	    BC(I,K)=d bc_{ik}/dPAR(1) (TN,PAR)
C         ENDIF
C 1000  CONTINUE
C       END

C     -------------------------------------------------------

C	Complex zero of a given function using Muller iteration with deflation
C	Function is assumed to be calculated as CF*2**IX
C	This subroutine uses reverse communication to calculate function
C	values. If IER<0 the function should be evaluated and MULER2
C	should be called back with new function value. Calculation of
C	function value should not change any other variables in the
C	call statement. First call should be with IER=0
C
C	CX1,CX2,CX3 : (input/output) Complex starting values for iteration
C		These will be updated during execution and CX3 should be
C		the best estimate for the zero.
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(CX3))
C	IER : (input/output) Error parameter, IER=0 implies successful execution
C		For the first call IER should be set to 0.
C		IER<0 implies that execution is not over and the subroutine
C			needs a new function evaluation at z=CX. After calculating
C			the function value MULER2 should be called back.
C		IER=42 implies that roundoff errors appear to be dominating
C			and calculations are terminated.
C		IER=43 implies that iteration failed to converge to specified accuracy
C			but a weaker convergence criterion is satisfied
C		IER=404 implies 2 of the 3 starting values are identical and
C			no calculations are done
C		IER=431 implies that iteration goes outside the specified limits
C		IER=432 implies that iteration failed to converge to specified accuracy
C		IER=433 implies that denominator for calculating the iteration
C			vanishes and it is not possible to continue further
C	CF : (input) Calculated value of the function at z=CX
C		If subroutine exits with IER<0, then the calling routine
C		should calculate the function value at CX and call MULER2
C		with this value stored in CF and IX. Other variables should
C		not be changed.
C	CX : (output) Value of z at which the function evaluation is needed
C	IX : (input) The exponent of function value, The function value 
C		should be CF*2**IX
C	NZ : (input) Number of zeros already known (for deflation)
C	CZERO : (input) Complex array of length NZ containing the known
C		zeros for deflation.
C	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
C		when ABS(CX) > RMAX
C	
C	Required routines : None (Function has to be calculated by calling program)
C
      SUBROUTINE MULER2(CX1,CX2,CX3,REPS,AEPS,IER,CF,CX,IX,NZ,
     1                  CZERO,RMAX)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      PARAMETER(REPS0=1.D-3,NIT=50)
      DIMENSION CZERO(NZ)
      SAVE

      IFL=-IER
      IF(IFL.GT.0.AND.IFL.LT.4) THEN
C	Jump to proper location and continue execution
        GO TO (1000,2000,3000), IFL
      ENDIF
C	Initial call to subroutine, start from beginning
      IF(CX2.EQ.CX1.OR.CX2.EQ.CX3.OR.CX3.EQ.CX1) THEN
C	If 2 of the starting values are equal, then quit
        IER=404
        RETURN
      ENDIF

      IER=-1
      CX=CX1
C	To evaluate the function at CX1
      RETURN
1000  CF1=CF
      IF1=IX
C	perform deflation
      DO 1500 J=1,NZ
1500  CF1=CF1/(CX1-CZERO(J))

      IER=-2
      CX=CX2
C	To calculate the function value at CX2
      RETURN
2000  CF2=CF
      IF2=IX
C	perform deflation
      DO 2500 J=1,NZ
2500  CF2=CF2/(CX2-CZERO(J))

      CF1I=CF1*2.D0**(IF1-IF2)
      CH1=(CF2-CF1I)/(CX2-CX1)
      I=0
      DX1=ABS(CX3-CX2)
      CX=CX3
      IER=-3
C	To calculate the function value at CX3
      RETURN

C	Loop for the Muller iteration
3000  I=I+1
      CF3=CF
      IF3=IX
      IF(NZ.GT.0) THEN
C	perform deflation
        DO 3500 J=1,NZ
3500    CF3=CF3/(CX3-CZERO(J))
      ENDIF

      IER=0
      IF(CX3.EQ.CX1.OR.CX3.EQ.CX2) RETURN
      IF(CF3.EQ.0.0) RETURN
      CF2A=CF2*2.D0**(IF2-IF3)
      CH2=(CF3-CF2A)/(CX3-CX2)
      CH1A=CH1*2.D0**(IF2-IF3)
      C2=(CH2-CH1A)/(CX3-CX1)
      CG=CH2+(CX3-CX2)*C2
      CI=SQRT(CG*CG-4.*CF3*C2)
      CD=CG+CI
      CD1=CG-CI
      IF(ABS(CD1).GT.ABS(CD)) CD=CD1

      IF(CD.EQ.0.0) THEN
C	If the denominator is zero, then quit
        IER=433
        RETURN
      ENDIF

      CDX=-2.*CF3/CD
      CX1=CX2
      CX2=CX3
C	The new approximation to zero
      CX3=CX3+CDX
      CF1=CF2
      CF2=CF3
      IF1=IF2
      IF2=IF3
      CH1=CH2
      DX=ABS(CDX)
      IF(I.GT.2.AND.DX.LT.MAX(REPS*ABS(CX3),AEPS)) RETURN
      IF(I.GT.5.AND.DX.LT.MAX(REPS0*ABS(CX3),AEPS)
     1   .AND.DX.GT.DX1) THEN
C	Roundoff errors appear to be dominating, hence quit
        IER=42
        RETURN
      ENDIF

      DX1=DX
      IF(ABS(CX3).GT.RMAX) THEN
C	Iteration goes outside the specified range
        IER=431
        RETURN
      ENDIF

      IF(I.LT.NIT) THEN
        IER=-3
        CX=CX3
C	To calculate the function value at CX3 and continue the loop
        RETURN
      ENDIF

C	Iteration fails to converge
      IER=43
      IF(DX.LT.MAX(REPS0*ABS(CX3),AEPS)) RETURN
      IER=432
      END

C     --------------------------------------------------------

C	To setup the matrix of system of linear equations arising from
C	finite difference approximation of ordinary differential equations
C	This routine is called by GEVP_C, the matrix is complex
C
C	N : (input) Number of mesh points used.
C	M : (input) Number of first order differential equations in the system
C	ML : (input) Number of boundary conditions at the first boundary t=T(1)
C	A : (output) Complex array of length (M+ML)*2M*N containing
C		the matrix of equations.
C	BC : (output) Complex array of length (M+ML)*(M+1) containing the
C		coefficients of boundary conditions
C	X : (input) Complex array of length M*N containing the current approximation
C		to the solution.
C	XC : (output) Complex array of length M*N containing the right hand
C		side of finite difference equations calculated by the routine
C	T : (input) Real array of length N containing the mesh points.
C	PAR : (input) Complex array containing the parameters to be passed
C		on to subroutines EQN and BCS. This array is not used by
C		the subroutine.
C	EQN : (input) Name of subroutine to calculate the equation matrix
C	BCS : (input) Name of subroutine to calculate the boundary conditions
C
C	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) and
C	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) must be supplied by the user.
C	The form of these routines is described in documentation for GEVP_C
C
C	Required routines : EQN, BCS
C
      SUBROUTINE SETMAT_C(N,M,ML,A,BC,X,XC,T,PAR,EQN,BCS)
C      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT REAL*8(H,O,P,R,T)
      IMPLICIT COMPLEX*16(A-G,U-Z)
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

C     ----------------------------------------------------

      SUBROUTINE EQN(J,M,ML,PAR,A,B,X,F,T)
      IMPLICIT REAL*8(H,O,Q-T)
      IMPLICIT COMPLEX*16(A-G,P,U-Z)
      DIMENSION A(M+ML,M),B(M+ML,M),X(M),PAR(*),F(M)

C     EQ. FOR SPHEROIDAL HARMONICS   X1'=X2
C     (1-T*T)X2'=-(LAMDA-C**2*T**2)X1+2(M+1)T*X2
C     LAMDA=PAR(1)   C**2=PAR(2)   M=PAR(3)

      DO 1000 I=1,M
        F(I)=0.0
        DO 1000 K=1,M
          A(K,I)=0.0
1000  B(K,I)=0.0

      B(1,1)=1.
      A(1,2)=1.
      B(2,2)=1.-T**2
      A(2,1)=-(PAR(1)-PAR(2)*T**2)
      A(2,2)=2.*(PAR(3)+1)*T
     
      END

C     ----------------------------------------------

      SUBROUTINE BCS(M,ML,PAR,BC,G,T1,T2,X1,X2)
      IMPLICIT REAL*8(H,O,Q-T)
      IMPLICIT COMPLEX*16(A-G,P,U-Z)
      DIMENSION PAR(*),BC(M+ML,M),G(M),X1(*),X2(*)

C     BOUNDARY CONDITIONS : X2=0  AT T=T1
C     X2-(LAMDA-C**2)/(2(M+1))X1=0  AT T=T2

      G(1)=0.0
      G(2)=0.0
      BC(1,1)=0.0
      BC(1,2)=1.0
      BC(2,1)=-0.5*(PAR(1)-PAR(2))/(PAR(3)+1.)
      BC(2,2)=1.0
      END

C     ----------------------------------------------------

      SUBROUTINE EQND(J,M,ML,PAR,A,B,T)
      IMPLICIT REAL*8(H,O,Q-T)
      IMPLICIT COMPLEX*16(A-G,P,U-Z)
      DIMENSION A(M+ML,M),B(M+ML,M),PAR(*)

C     DERIVATIVE OF EQ. MATRIX W.R.T. EIGENVALUE

      DO 1000 I=1,M
        DO 1000 K=1,M
          A(K,I)=0.0
1000  B(K,I)=0.0
      A(2,1)=-1.
      END

C     ----------------------------------------------

      SUBROUTINE BCSD(M,ML,PAR,BC,T1,T2)
      IMPLICIT REAL*8(H,O,Q-T)
      IMPLICIT COMPLEX*16(A-G,P,U-Z)
      DIMENSION PAR(*),BC(M+ML,M)

C     DERIVATIVE OF BOUNDARY COND. MATRIX W.R.T. EIGENVALUE

      DO 1000 I=1,M
        DO 1000 J=1,M
1000  BC(I,J)=0.0
      BC(2,1)=-0.5/(PAR(3)+1.)
      END
