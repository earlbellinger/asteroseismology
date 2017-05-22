C     PROGRAM TO FIND ONE EIGENVALUE AND THE CORRESPONDING EIGENVECTOR
C     OF A GENERAL REAL MATRIX USING INVERSE ITERATION
C     IT CALCULATES BOTH LEFT AND RIGHT EIGENVECTORS

      PROGRAM INVERS
      IMPLICIT REAL*8(A-H,O-Z)
C     COMPLEX*16 P,U,W,EI,RC,P0
      DIMENSION A(21,21),U(21),V(21),INT(21),W(21,22),B(21,21)
      DATA ((A(I,J),J=1,4),I=1,4)/-2,2,2,2, -3,3,2,2, -2,0,4,2,-1,0,0,5/
C      DATA ((A(I,J),J=1,4),I=1,4)/4,-5,0,3, 0,4,-3,-5, 5,-3,4,0,3,0,5,4/

51    FORMAT(10X,13HTHE MATRIX IS)
52    FORMAT(10F8.3)
53    FORMAT('   IER =',I4,5X,'IFLG =',I2,5X,'EIGENVALUE =',1PD14.6/
     1       5X,'RAYLEIGH QUOTIENT =',D14.6,5X,'FINAL SHIFT =',D14.6)
54    FORMAT('  EIGENVECTOR =',1P4D14.6/(2X,4D14.6))
55    FORMAT('   INITIAL SHIFT =',1P2D14.6)
56    FORMAT('  LEFT EIGENVECTOR =',1P4D14.6/(2X,4D14.6))

      NUM=100
      IA=21
      REPS=1.D-14

100   PRINT *,'TYPE M=ORDER, P=SHIFT,  IFLG=0/1/2'
      PRINT *,'            (QUITS WHEN M.LE.0)'
      READ *,M,P,IFLG
      IF(M.LE.0) STOP

      WRITE(6,51)
      PRINT *,'TYPE IN THE MATRIX ROW-WISE'
      DO 1000 I=1,M
        U(I)=1.0
        V(I)=1.0
        PRINT *,I,'TH ROW'
        READ *,(A(I,J),J=1,M)
        WRITE(6,52) (A(I,J),J=1,M)
1000  CONTINUE

      PRINT *,'TYPE INITIAL VECTOR'
      READ *,(U(I),I=1,M)
      WRITE(6,55) P
      P0=P
      CALL INVIT(A,M,IA,P,U,IFLG,EI,RC,REPS,W,INT,NUM,IER)
      WRITE(6,53) IER,IFLG,EI,RC,P
      WRITE(6,54) (U(I),I=1,M)

C	Calculate the left eigenvector for the same matrix 
      P=P0
      DO 2000 I=1,M
        U(I)=1.0
        V(I)=1.0
2000  CONTINUE
      CALL INVIT_L(A,M,IA,P,U,IFLG,EI,RC,REPS,W,INT,NUM,IER)
      WRITE(6,53) IER,IFLG,EI,RC,P
      WRITE(6,56) (U(I),I=1,M)
      GO TO 100
      END

C     -----------------------------------------------------

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

C     ------------------------------------------------------

C	Real eigenvalue and eigenvector of a real matrix using inverse iteration
C
C	A : (input) Real array of length IA*M containing the matrix elements
C	M : (input) Order of the matrix
C	IA : (input) The first dimension of A as declared in the calling program
C		IA.GE.M
C	P : (input/output) Initial value of the shift. This will be modified
C		by the program if IFLG>0
C	U : (input/output) Real array of length M, which should specify the
C		initial approximation to eigenvector. After execution it
C		will contain the calculated eigenvector.
C	IFLG : (input) Integer variable to specify the type of iteration required
C		If IFLG=0 the shift P is kept fixed
C		If IFLG=1 the shift P is varied using Rayleigh quotient
C		If IFLG=2 the shift P is varied using max(V_s+1)
C	EI : (output) Estimated eigenvalue using simple inverse iteration
C	ERC : (output) Estimated eigenvalue using Rayleigh quotient
C	REPS : (input) Required absolute accuracy. Iteration is terminated
C		when all components of eigenvector and the eigenvalue have
C		converged to REPS.
C	WK : Real array of length M*(M+1) used as scratch space
C	IWK : Integer array of length M used as scratch space
C	NIT : (input/output) Number of iterations required. If it is
C		zero or negative NIT is set to NIT0=100
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=106 implies that M.LE.1 or M>IA, in which case no
C			calculations are done
C		IER=141 implies that vector is zero at some stage and
C			calculations are aborted
C		IER=142 implies that inverse iteration has failed to converge
C
C	Required routines : GAUELM
C
      SUBROUTINE INVIT(A,M,IA,P,U,IFLG,EI,ERC,REPS,WK,IWK,NIT,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(IA,M),U(M),WK(M,M+1),IWK(M)
      PARAMETER(NIT0=100)

      IF(M.LE.1.OR.M.GT.IA) THEN
        IER=106
        RETURN
      ENDIF

C	Copy the matrix to WK and apply the shift
      DO 1100 I=1,M
        WK(I,M+1)=U(I)
        DO 1000 J=1,M
1000    WK(J,I)=A(J,I)
1100  WK(I,I)=A(I,I)-P

      NUM=1
      LJ=M
      IFL=1
      IF(NIT.LE.0) NIT=NIT0
C	Perform Gaussian elimination on A-pI
      CALL GAUELM(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
      IF(IER.GT.0) RETURN

      EPI=0.0
C	Loop for inverse iteration
      DO 5000 J=1,NIT
        CALL GAUELM(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
        IF(IER.GT.0) RETURN

C	Normalising the vector U
        R1=0.0
        KM=0
        DO 3500 K=1,M
          IF(R1.LT.ABS(U(K))) THEN
            R1=ABS(U(K))
            KM=K
          ENDIF
3500    CONTINUE
        UKM=U(KM)
        IF(UKM.EQ.0.0) THEN
C	If the vector is zero, then quit
          IER=141
          RETURN
        ENDIF

C	The eigenvalue
        EI=WK(KM,M+1)/UKM+P
        S1=0.0
        S2=0.0
C	Calculating the Rayleigh quotient
        DO 4000 K=1,M
          S1=S1+U(K)*WK(K,M+1)
C	For complex eigenvalues use the following statement instead of the
C	preceding one
C         S1=S1+CONJG(U(K))*WK(K,M+1)
          S2=S2+ABS(U(K))**2
          U(K)=U(K)/UKM
4000    CONTINUE
        ERC=P+S1/S2

C	Convergence check
        R1=ABS(EI-EPI)
        DO 4500 I=1,M
          R1=MAX(R1,ABS(WK(I,M+1)-U(I)))
          WK(I,M+1)=U(I)
4500    CONTINUE
        IF(ABS(R1).LT.REPS) RETURN
        EPI=EI

        IF(IFLG.GE.1) THEN
C	Update the shift
          P=ERC
          IF(IFLG.EQ.2) P=EI
C	Setting up the new matrix A-pI
          DO 4700 I=1,M
            DO 4600 K=1,M
4600        WK(K,I)=A(K,I)
4700      WK(I,I)=A(I,I)-P
          IFL=1
          CALL GAUELM(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
          IF(IER.GT.0) RETURN
        ENDIF

5000  CONTINUE
C	Iteration fails to converge
      IER=142
      END

C	--------------------------------------------------------------

C	Real eigenvalue and left eigenvector of a real matrix using inverse iteration
C
C	A : (input) Real array of length IA*M containing the matrix elements
C	M : (input) Order of the matrix
C	IA : (input) The first dimension of A as declared in the calling program
C	P : (input/output) Initial value of the shift. This will be modified
C		by the program if IFLG>0
C	U : (input/output) Real array of length M, which should specify the
C		initial approximation to eigenvector. After execution it
C		will contain the calculated left eigenvector.
C	IFLG : (input) Integer variable to specify the type of iteration required
C		If IFLG=0 the shift P is kept fixed
C		If IFLG=1 the shift P is varied using Rayleigh quotient
C		If IFLG=2 the shift P is varied using max(V_s+1)
C	EI : (output) Estimated eigenvalue using simple inverse iteration
C	ERC : (output) Estimated eigenvalue using Rayleigh quotient
C	REPS : (input) Required absolute accuracy. Iteration is terminated
C		when all components of eigenvector and the eigenvalue have
C		converged to REPS.
C	WK : Real array of length M*(M+1) used as scratch space
C	IWK : Integer array of length M used as scratch space
C	NIT : (input/output) Number of iterations required. If it is
C		zero or negative NIT is set to NIT0=100
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=106 implies that M.LE.1 or M>IA, in which case no
C			calculations are done
C		IER=141 implies that vector is zero at some stage and
C			calculations are aborted
C		IER=142 implies that inverse iteration has failed to converge
C
C	Required routines : GAUELM
C
      SUBROUTINE INVIT_L(A,M,IA,P,U,IFLG,EI,ERC,REPS,WK,IWK,NIT,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(IA,M),U(M),WK(M,M+1),IWK(M)
      PARAMETER(NIT0=100)

      IF(M.LE.1.OR.M.GT.IA) THEN
        IER=106
        RETURN
      ENDIF

C	Copy the transposed matrix to WK and apply the shift
      DO 1100 I=1,M
        WK(I,M+1)=U(I)
        DO 1000 J=1,M
1000    WK(J,I)=A(I,J)
1100  WK(I,I)=A(I,I)-P

      NUM=1
      LJ=M
      IFL=1
      IF(NIT.LE.0) NIT=NIT0
C	Perform Gaussian elimination on A-pI
      CALL GAUELM(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
      IF(IER.GT.0) RETURN

      EPI=0.0
C	Loop for inverse iteration
      DO 5000 J=1,NIT
        CALL GAUELM(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
        IF(IER.GT.0) RETURN

C	Normalising the vector U
        R1=0.0
        KM=0
        DO 3500 K=1,M
          IF(R1.LT.ABS(U(K))) THEN
            R1=ABS(U(K))
            KM=K
          ENDIF
3500    CONTINUE
        UKM=U(KM)
        IF(UKM.EQ.0.0) THEN
C	If the vector is zero, then quit
          IER=141
          RETURN
        ENDIF

C	The eigenvalue
        EI=WK(KM,M+1)/UKM+P
        S1=0.0
        S2=0.0
C	Calculating the Rayleigh quotient
        DO 4000 K=1,M
          S1=S1+U(K)*WK(K,M+1)
C	For complex eigenvalues use the following statement instead of the
C	preceding one
C         S1=S1+CONJG(U(K))*WK(K,M+1)
          S2=S2+ABS(U(K))**2
          U(K)=U(K)/UKM
4000    CONTINUE
        ERC=P+S1/S2

C	Convergence check
        R1=ABS(EI-EPI)
        DO 4500 I=1,M
          R1=MAX(R1,ABS(WK(I,M+1)-U(I)))
          WK(I,M+1)=U(I)
4500    CONTINUE
        IF(ABS(R1).LT.REPS) RETURN
        EPI=EI

        IF(IFLG.GE.1) THEN
C	Update the shift
          P=ERC
          IF(IFLG.EQ.2) P=EI
C	Setting up the new matrix A-pI
          DO 4700 I=1,M
            DO 4600 K=1,M
4600        WK(K,I)=A(I,K)
4700      WK(I,I)=A(I,I)-P
          IFL=1
          CALL GAUELM(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
          IF(IER.GT.0) RETURN
        ENDIF

5000  CONTINUE
C	Iteration fails to converge
      IER=142
      END
