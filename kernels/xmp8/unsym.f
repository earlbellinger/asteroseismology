C     PROGRAM FOR FINDING EIGENVALUES AND EIGENVECTORS OF UNSYMMETRIC MATRIX

      PROGRAM UNSYM

      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      DIMENSION A(30,30),D(30),E(30),CV(30,30),EI(30),AB(30,30)
     1          ,INT(30),CWK(300)
      DATA ((A(I,J),J=1,7),I=1,7)/6,0,0,0,0,1,0, 0,4,0,3D-8,1D-4,2D-4,
     1   1D-3, 1,10000,7,0,0,-2,20, 0,2D8,0,1,-40000,30000,-4D5,
     2   -2,-30000,0,1D-4,2,2,40, 0,0,0,0,0,0,0,
     3   0,1000,0,4D-5,0.1D0,-0.2D0,3/

51    FORMAT(10X,13HTHE MATRIX IS)
52    FORMAT(1P7D11.3)
53    FORMAT('   IER =',I4,5X,'LOW =',I4,5X,'IGH=',I4,5X,
     1       'THE BALANCED MATRIX IS')
54    FORMAT('   DIAGONAL TRANSFORMATION MATRIX FOR BALANCING ='
     1       ,1P2D14.6/(2X,5D14.6))
55    FORMAT('   IER =',I4,5X,'THE HESSENBERG MATRIX')
56    FORMAT(1P7D13.5)
57    FORMAT(9X,'EIGENVECTORS')
58    FORMAT(I5,(5X,1P2D14.6,2X,2D14.6))
59    FORMAT('   IER =',I4,5X,'EIGENVALUES =',1P2D14.6/
     1       (5X,2D14.6,2X,2D14.6))

100   PRINT *,'TYPE N=ORDER,      (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0) STOP

      WRITE(6,51)
      PRINT *,'TYPE IN THE MATRIX ROW-WISE'
        DO 500 I=1,N
        PRINT *,I,'TH ROW'
        READ *,(A(I,J),J=1,N)
        WRITE(6,52)(A(I,J),J=1,N)
500   CONTINUE

      REPS=1.D-14
      IA=30
C	USE BASE 10 FOR ILLUSTRATION. IT SHOULD BE 2 NORMALLY
      B=10
      LOW=1
      IGH=N

C     BALANCE THE MATRIX 

      CALL BALANC(A,N,IA,B,LOW,IGH,D,IER)
      WRITE(6,53) IER,LOW,IGH

C     PRESERVE THE MATRIX FOR FINDING EIGENVECTORS

      DO 1000 I=1,N
        DO 800 J=1,N
800     AB(J,I)=A(J,I)
1000  WRITE(6,52) (A(I,J),J=1,N)
      WRITE(6,54) (D(I),I=1,N)

C     REDUCE THE MATRIX TO HESSENBERG FORM

      CALL ELMHES(A,N,IA,LOW,IGH,INT,IER)
      WRITE(6,55) IER
      DO 1500 I=1,N
1500  WRITE(6,56) (A(I,J),J=1,N)

C     FIND EIGENVALUES OF HESSENBERG MATRIX

      CALL HQR(A,N,IA,E,EI,REPS,IER)
      WRITE(6,59) IER,(E(I),EI(I),I=1,N)

C     FINDING THE EIGENVECTORS USING INVERSE ITERATION

      NIT=100
      IFLG=0
      DO 3000 I=1,N
C       PERTURB THE EIGENVALUE TO AVOID ZERO PIVOT
        CP=CMPLX(E(I),EI(I))*(1.+20.*REPS)+REPS**2
        IFLG=0
        DO 2500 J=1,N
2500    CV(J,I)=1.0
        CALL INVIT_C(AB,N,IA,CP,CV(1,I),IFLG,CEI,CR,REPS,CWK,INT,
     1               NIT,IER)
C        WRITE(6,58) IER,(CV(J,I),J=1,N)
3000  CONTINUE

C     BACK-TRANSFORM THE EIGENVECTORS

      M=N
      IZ=30
      CALL BALBAK(N,LOW,IGH,CV,M,IZ,D)

C     NORMALISE THE EIGENVECTOR BEFORE PRINTING

      WRITE(6,57)
      DO 4000 I=1,N
        RMAX=0
        DO 3300 J=1,N
          IF(RMAX.GE.ABS(CV(J,I))) GO TO 3300
          IM=J
          RMAX=ABS(CV(J,I))
3300    CONTINUE
        CJ=CV(IM,I)
4000  WRITE(6,58) I,(CV(J,I)/CJ,J=1,N)
      GO TO 100
      END

C     -----------------------------------------------------

C	To balance a general real matrix by exact diagonal similarity transformation
C
C	A : (input/output) Real array of length IA*N containing the matrix
C		elements. After execution it will contain the balanced matrix
C	N : (input) Order of matrix
C	IA : (input) First dimension of array A as declared in the calling
C		program
C	B : (input) Base of the arithmetic to be used for balancing.
C		Should generally be 2.
C	LOW : (output) Index of column such that in balanced matrix A(i,j)=0
C		if i>j and j<LOW, the first LOW-1 eigenvalues are isolated
C	IGH : (output) Index of row such that in balanced matrix A(i,j)=0
C		if i>j and i>IGH, the last N-IGH eigenvalues are isolated
C		After balancing only the sub-matrix between rows and columns
C		LOW to IGH need to be considered
C	D : (output) Real array of length N containing information about
C		transformations used for balancing. D(LOW) to D(IGH)
C		contain elements used for balancing, while other elements
C		contain information about permutations used
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=112 implies that N.LE.1 or N>IA
C
C	Required routines : None

      SUBROUTINE BALANC(A,N,IA,B,LOW,IGH,D,IER)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(GAM=0.95D0)
      DIMENSION A(IA,N),D(N)

      IF(N.LE.1.OR.N.GT.IA) THEN
        IER=112
        RETURN
      ENDIF
      IER=0
      B2=B*B
      LOW=1
      IGH=N

C	Search for rows isolating an eigenvalue
2000  DO 3000 J=IGH,1,-1
        R=0
        DO 2100 I=1,J-1
2100    R=R+ABS(A(J,I))
        DO 2200 I=J+1,IGH
2200    R=R+ABS(A(J,I))

        IF(R.EQ.0.0) THEN
C	Push the row to bottom
          D(IGH)=J
          IF(J.NE.IGH) THEN
            DO 2400 I=1,IGH
              T=A(I,J)
              A(I,J)=A(I,IGH)
              A(I,IGH)=T
2400        CONTINUE
            DO 2600 I=LOW,N
              T=A(J,I)
              A(J,I)=A(IGH,I)
              A(IGH,I)=T
2600        CONTINUE
          ENDIF
          IGH=IGH-1
          GO TO 2000
        ENDIF
3000  CONTINUE

C	Search for columns isolating an eigenvalue
3200  DO 4000 J=LOW,IGH
        C=0
        DO 3300 I=LOW,J-1
3300    C=C+ABS(A(I,J))
        DO 3400 I=J+1,IGH
3400    C=C+ABS(A(I,J))

C	Move the column to the left end
        IF(C.EQ.0) THEN
          D(LOW)=J
          IF(J.NE.LOW) THEN
            DO 3600 I=1,IGH
              T=A(I,J)
              A(I,J)=A(I,LOW)
              A(I,LOW)=T
3600        CONTINUE
            DO 3800 I=LOW,N
              T=A(J,I)
              A(J,I)=A(LOW,I)
              A(LOW,I)=T
3800        CONTINUE
          ENDIF
          LOW=LOW+1
          GO TO 3200
        ENDIF
4000  CONTINUE

C	Balance the sub-matrix in rows LOW to IGH
      DO 4200 I=LOW,IGH
4200  D(I)=1
4300  QC=.FALSE.
      DO 5000 I=LOW,IGH
        C=0
        R=0
        DO 4400 J=LOW,IGH
          IF(J.NE.I) THEN
            C=C+ABS(A(J,I))
            R=R+ABS(A(I,J))
          ENDIF
4400    CONTINUE
        G=R/B
        F=1
        S=C+R
4500    IF(C.LT.G) THEN
          F=F*B
          C=C*B2
          GO TO 4500
        ENDIF
        G=R*B
4600    IF(C.GE.G) THEN
          F=F/B
          C=C/B2
          GO TO 4600
        ENDIF

C	Apply the transformation
        IF((C+R)/F.LT.GAM*S) THEN
          G=1./F
          D(I)=D(I)*F
          QC=.TRUE.
          DO 4800 J=LOW,N
4800      A(I,J)=A(I,J)*G
          DO 4900 J=1,IGH
4900      A(J,I)=A(J,I)*F
        ENDIF
5000  CONTINUE

C	perform one more cycle of balancing
      IF(QC) GO TO 4300
      END

C     ---------------------------------------------------

C	Perform back-transformation of a set of right eigenvectors from
C	those of balanced matrix to that for original matrix
C
C	N : (input) Order of matrix
C	LOW, IGH : (input) After balancing only rows and columns from
C		LOW to IGH need to be considered, since other rows or
C		columns contain isolated eigenvalues
C	CZ : (input/output) Complex array of length IZ*M containing the
C		eigenvectors of balanced matrix. After execution it will
C		be overwritten by eigenvectors of the original matrix.
C	M : (input) Number of eigenvectors to be balanced
C	IZ : (input) First dimension of array CZ as declared in the calling
C		program
C	D : (input) Real array of length N, containing information about
C		transformation used for balancing.
C
C	Required routines : None

      SUBROUTINE BALBAK(N,LOW,IGH,CZ,M,IZ,D)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16(C)
      DIMENSION CZ(IZ,M),D(N)

      DO 3000 I=LOW,IGH
        S=D(I)
C	For back-transforming left eigenvectors use the following
C	statement instead of the preceding one
C       S=1/D(I)

        DO 2200 J=1,M
2200    CZ(I,J)=CZ(I,J)*S
3000  CONTINUE

      DO 3400 I=LOW-1,1,-1
        K=D(I)
        IF(K.NE.I) THEN
C	Exchange the corresponding rows
          DO 3200 J=1,M
            CS=CZ(I,J)
            CZ(I,J)=CZ(K,J)
3200      CZ(K,J)=CS
        ENDIF
3400  CONTINUE

      DO 4000 I=IGH+1,N
        K=D(I)
        IF(K.NE.I) THEN
C	Exchange the corresponding rows
          DO 3800 J=1,M
            CS=CZ(I,J)
            CZ(I,J)=CZ(K,J)
3800      CZ(K,J)=CS
        ENDIF
4000  CONTINUE
      END

C     -------------------------------------------------------

C	To reduce a general real matrix to Hessenberg form using stabilised
C	elementary transformations
C	It is advisable to balance the matrix before applying this transformations
C
C	A : (input/output) Real array of length IA*N containing the matrix
C		After execution the reduced matrix will be overwritten on
C		the same array
C	N : (input) Order of matrix
C	IA : (input) First dimension of array A as declared in the calling
C		program
C	LOW, IGH : (input) After balancing only rows and columns from
C		LOW to IGH need to be considered, since other rows or
C		columns contain isolated eigenvalues. If the matrix is
C		not balanced use LOW=1, IGH=N
C	INC : (output) Integer array of length N, which will contain information
C		about row and column interchanges during reduction
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=113 implies that N.LE.1 or N>IA
C
C	Required routines : None

      SUBROUTINE ELMHES(A,N,IA,LOW,IGH,INC,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(IA,N),INC(N)

      IF(N.LE.1.OR.N.GT.IA) THEN
        IER=113
        RETURN
      ENDIF
      IER=0
      IF(LOW.GT.IGH-2) RETURN

      DO 5000 M=LOW+1,IGH-1
        I=M
C	Find the pivot
        AMAX=0
        DO 2000 J=M,IGH
          IF(ABS(A(J,M-1)).GT.ABS(AMAX)) THEN
            AMAX=A(J,M-1)
            I=J
          ENDIF
2000    CONTINUE
        INC(M)=I

        IF(I.NE.M) THEN
C	Interchange the corresponding rows and columns
          DO 2200 J=M-1,N
            T=A(I,J)
            A(I,J)=A(M,J)
            A(M,J)=T
2200      CONTINUE
          DO 2400 J=1,IGH
            T=A(J,I)
            A(J,I)=A(J,M)
            A(J,M)=T
2400      CONTINUE
        ENDIF

        IF(AMAX.NE.0.0) THEN
C	Perform Gaussian elimination
          DO 4000 I=M+1,IGH
            T=A(I,M-1)
            IF(T.NE.0.0) THEN
              T=T/AMAX
              A(I,M-1)=T
              DO 3000 J=M,N
3000          A(I,J)=A(I,J)-T*A(M,J)
              DO 3200 J=1,IGH
3200          A(J,M)=A(J,M)+T*A(J,I)
            ENDIF
4000      CONTINUE
        ENDIF
5000  CONTINUE
      END

C     -------------------------------------------------------------

C	Solution of a system of linear equations with complex coefficients
C	using Gaussian elimination with partial pivoting
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
C		    	decomposition should have been calculated earlier
C
C	Required routines : None

      SUBROUTINE GAUELM_C(N,NUM,A,X,DET,INT,LJ,IER,IFLG)
C      IMPLICIT REAL*8(A-H,O-Z)
C	For complex matrices use the following statements instead
      IMPLICIT REAL*8(R)
      IMPLICIT COMPLEX*16(A-H,S-Z)

      DIMENSION A(LJ,N),INT(N),X(LJ,NUM)

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

          INT(K)=KM
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
        INT(N)=N
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
          IF(K.NE.INT(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INT(K),J)
            X(INT(K),J)=T1
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

C     ---------------------------------------------------------------

C	To find all eigenvalues of an upper-Hessenberg matrix using QR algorithm
C
C	H : (input/output) Real array of length IH*NN containing the matrix
C		elements. The contents of array are destroyed during execution
C	NN : (input) Order of matrix
C	IH : (input) First dimension of array H as declared in the calling
C		program
C	ER : (output) Real array of length NN containing the real part of
C		the eigenvalues.
C	EI : (output) Real array of length NN containing the imaginary part of
C		the eigenvalues.
C	REPS : (input) Required relative accuracy in eigenvalues. It should
C		be of the order of machine accuracy
C	IER : (output) error parameter; IER=0 implies successful execution
C		IER=114 implies that N.LE.1 or N>IH, in which case no
C			calculations are done
C		IER=145 implies that QR iteration failed to converge at
C			some stage and calculations are abandoned.
C
C	Required routines : None

      SUBROUTINE HQR(H,NN,IH,ER,EI,REPS,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=30)
      DIMENSION H(IH,NN),ER(NN),EI(NN)

      T=0
      N=NN
      IF(N.LE.1.OR.N.GT.IH) THEN
        IER=114
        RETURN
      ENDIF
      IER=0
C	If all eigenvalues are found, then quit
2000  IF(N.EQ.0) RETURN

C	Loop for double QR iteration
      DO 5000 IT=1,NIT
C	Search for a small sub-diagonal element
        DO 2200 L=N,2,-1
          IF(ABS(H(L,L-1)).LE.REPS*(ABS(H(L-1,L-1))+ABS(H(L,L))))
     1       GO TO 2300
2200    CONTINUE
        L=1
2300    X=H(N,N)

        IF(L.EQ.N) THEN
C	One eigenvalue is isolated
          ER(N)=X+T
          EI(N)=0.0
          N=N-1
          GO TO 2000
        ENDIF

        Y=H(N-1,N-1)
        W=H(N,N-1)*H(N-1,N)
        IF(L.EQ.N-1) THEN
C	A pair of eigenvalues are isolated
          P=(Y-X)/2.
          D=P*P+W
          Y=SQRT(ABS(D))
          X=X+T
          IF(D.GT.0) THEN
C	Pair of real eigenvalues
            IF(P.LT.0) Y=-Y
            Y=P+Y
            ER(N-1)=X+Y
            ER(N)=X-W/Y
            EI(N-1)=0
            EI(N)=0
          ELSE
C	Pair of complex conjugate eigenvalues
            ER(N)=X+P
            ER(N-1)=ER(N)
            EI(N-1)=Y
            EI(N)=-Y
          ENDIF
          N=N-2
C	Find next eigenvalue
          GO TO 2000
        ENDIF

C	Apply special shifts
        IF(IT.EQ.10.OR.IT.EQ.20) THEN
          T=T+X
          DO 2500 I=1,N
2500      H(I,I)=H(I,I)-X
          S=ABS(H(N,N-1))+ABS(H(N-1,N-2))
          X=0.75D0*S
          Y=X
          W=-0.4375D0*S*S
        ENDIF

C	Search for two consecutive small sub-diagonal elements
        DO 2800 M=N-2,L,-1
          Z=H(M,M)
          R=X-Z
          S=Y-Z
          P=(R*S-W)/H(M+1,M)+H(M,M+1)
          U=H(M+1,M+1)-Z-R-S
          R=H(M+2,M+1)
          S=ABS(P)+ABS(U)+ABS(R)
          P=P/S
          U=U/S
          R=R/S
          IF(M.EQ.L) GO TO 3000
          IF(ABS(H(M,M-1))*(ABS(U)+ABS(R)).LE.REPS*ABS(P)*
     1       (ABS(H(M-1,M-1))+ABS(Z)+ABS(H(M+1,M+1)))) GO TO 3000
2800    CONTINUE
3000    DO 3200 I=M+3,N
          H(I,I-3)=0.0
3200    H(I,I-2)=0.0
        H(M+2,M)=0.0

C	Double QR transform for rows L to N and columns M to N
        DO 4000 K=M,N-1
          IF(K.NE.M) THEN
            P=H(K,K-1)
            U=H(K+1,K-1)
            R=0
            IF(K.NE.N-1) R=H(K+2,K-1)
            X=ABS(P)+ABS(U)+ABS(R)
            IF(X.EQ.0.0) GO TO 4000
            P=P/X
            U=U/X
            R=R/X
          ENDIF

          S=SQRT(P*P+U*U+R*R)
          IF(P.LT.0) S=-S
          IF(K.NE.M) THEN
            H(K,K-1)=-S*X
          ELSE IF(L.NE.M) THEN
            H(K,K-1)=-H(K,K-1)
          ENDIF
          P=P+S
          X=P/S
          Y=U/S
          Z=R/S
          U=U/P
          R=R/P

C	Row modification
          DO 3400 J=K,N
            P=H(K,J)+U*H(K+1,J)
            IF(K.NE.N-1) THEN
              P=P+R*H(K+2,J)
              H(K+2,J)=H(K+2,J)-P*Z
            ENDIF
            H(K+1,J)=H(K+1,J)-P*Y
            H(K,J)=H(K,J)-P*X
3400      CONTINUE

C	Column modification
          DO 3600 I=L,MIN(N,K+3)
            P=X*H(I,K)+Y*H(I,K+1)
            IF(K.NE.N-1) THEN
              P=P+Z*H(I,K+2)
              H(I,K+2)=H(I,K+2)-P*R
            ENDIF
            H(I,K+1)=H(I,K+1)-P*U
            H(I,K)=H(I,K)-P
3600      CONTINUE
4000    CONTINUE
5000  CONTINUE

C	QR iteration fails to converge
      IER=145
      END

C     ----------------------------------------------------------

C	Complex eigenvalue and eigenvector of a real matrix using inverse iteration
C
C	A : (input) Real array of length IA*M containing the matrix elements
C	M : (input) Order of the matrix
C	IA : (input) The first dimension of A as declared in the calling program
C	P : (input/output) Initial value of the shift. This will be modified
C		by the program if IFLG>0
C	U : (input/output) Complex array of length M, which should specify the
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
C	WK : Complex array of length M*(M+1) used as scratch space
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
C	Required routines : GAUELM_C
C
      SUBROUTINE INVIT_C(A,M,IA,P,U,IFLG,EI,ERC,REPS,WK,IWK,NIT,IER)
      IMPLICIT REAL*8(A,O,R)
      IMPLICIT COMPLEX*16(B-H,P,S-Z)
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
      CALL GAUELM_C(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
      IF(IER.GT.0) RETURN

      EPI=0.0
C	Loop for inverse iteration
      DO 5000 J=1,NIT
        CALL GAUELM_C(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
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
C          S1=S1+U(K)*WK(K,M+1)
C	For complex eigenvalues use the following statement instead of the
C	preceding one
          S1=S1+CONJG(U(K))*WK(K,M+1)
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
          CALL GAUELM_C(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
          IF(IER.GT.0) RETURN
        ENDIF

5000  CONTINUE
C	Iteration fails to converge
      IER=142
      END
