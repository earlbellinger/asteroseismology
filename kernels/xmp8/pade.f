C     PROGRAM TO CALCULATE COEFFICIENTS OF PADE APPROXIMATION FROM THOSE
C     OF THE MACLAURIN SERIES

      PROGRAM PADEAP
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(20),WK(100),C(20),IWK(10)

C     EXAMPLE 9.8: PADE APPROXIMATION FOR ARC TAN(X)

51    FORMAT('    IER =',I4,5X,'M =',I3,5X,'K =',I3/4X,
     1       'COEF. IN NUMERATOR =',1P4D14.6/(2X,5D14.6))
52    FORMAT('  COEF. IN DENOMINATOR =',1P4D14.6/(2X,5D14.6))

100   PRINT *,'TYPE M=DEGREE OF NUMERATOR, K=DEGREE OF DENOMINATOR'
      PRINT *,'             (QUITS WHEN (M+K).LE.0)'
      READ *,M,K
      IF(M+K.LE.0) STOP

C     COEFFICIENTS OF MACLAURIN SERIES FOR (ARC TAN(X))/X IN TERMS
C     OF X**2

      DO 1000 I=1,M+K+1
        C(I)=(-1.)**(I+1)/(2*I-1.D0)
1000  CONTINUE

      CALL PADE(M,K,A,C,IER,WK,IWK)
      WRITE(6,51) IER,M,K,(A(I),I=K+1,M+K+1)
      IF(K.GT.0) WRITE(6,52) (A(I),I=1,K)
      GO TO 100
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

C     --------------------------------------------------

C	To calculate coefficients of Pade approximation using the known
C	coefficients of Maclaurin series
C
C	M : (input) Degree of polynomial in the numerator
C	K : (input) Degree of polynomial in the denominator
C	A : (output) Real array of length M+K+1 containing the coefficients
C		of Pade approximation. A(I) is the coefficient of x**I in
C		the denominator, the constant term being 1. A(K+J+1) is
C		the coefficient of x**J in the numerator
C	C : (input) Real array of length M+K+1 containing the coefficients
C		of Maclaurin series. C(I+1) is the coefficient of x**I
C	IER : (output) error parameter, IER=0 implies successful execution
C		IER=612 implies that M<0 or K<0, no calculations are done
C		Other values of IER may be set by GAUELM
C	WK : Real array of length K*K used as scratch space
C	IWK : Integer array of length K used as scratch space
C
C	Required routines : GAUELM
C
      SUBROUTINE PADE(M,K,A,C,IER,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+K+1),C(M+K+1),IWK(K),WK(K,K)

      IF(M.LT.0.OR.K.LT.0) THEN
        IER=612
        RETURN
      ENDIF

      IF(K.GT.0) THEN
C	Setting up the coefficients matrix of linear equations to
C	calculate the coefficients in the denominator
        DO 2200 I=1,K
          DO 2000 J=1,K
            WK(J,I)=0.0
            IF(M+J-I+1.GT.0) WK(J,I)=C(M+J-I+1)
2000      CONTINUE
C	The right-hand side vector
          A(I)=-C(M+I+1)
2200    CONTINUE

        NUM=1
        LJ=K
        IFLG=0
C	Solve the system of linear equations
        CALL GAUELM(K,NUM,WK,A,DET,IWK,LJ,IER,IFLG)
        IF(IER.GT.0) RETURN
      ENDIF

C	Calculate the coefficients in the numerator
      A(K+1)=C(1)
      IF(M.GT.0) THEN
        DO 3000 I=2,M+1
          AI=C(I)
          IF(K.GT.0) THEN
            DO 2600 J=1,MIN(K,I-1)
2600        AI=AI+C(I-J)*A(J)
          ENDIF
          A(K+I)=AI
3000    CONTINUE
      ENDIF

      END
