C     Solution of a system of linear equations using Gaussian elimination
C     	for a band matrix
      PROGRAM BLK
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(100,40),X(100),INT(100),WK(100)
 
51    FORMAT(' SOLUTION :',(1P5D14.6))
52    FORMAT('   IER =',I4,5X,'N =',I4,5X,'K =',I3,4X,'DET =',1PD14.6,
     1   3X,'IDET =',I5)
 
      IS=-1
100   PRINT *,'TYPE N=NO. OF EQS, K=BAND WIDTH    (QUITS WHEN N.LE.0)'
      READ *,N,K
      IF(N.LE.0) STOP
 
C     Generate a matrix using random numbers and RHS assuming solution=1
 
      DO 1200 I=1,N
        S1=0.0
        DO 1000 J=1,2*K+1
          A(I,J)=0.0
          IF(J.GT.K-I+1.AND.J.LE.K+1+N-I) THEN
            A(I,J)=RANF(IS)
            S1=S1+A(I,J)
          ENDIF
1000    CONTINUE
        X(I)=S1
1200  CONTINUE
 
      NUM=1
      IFLG=0
      LJ=100
      CALL GAUBND(N,K,NUM,A,X,DET,IDET,INT,LJ,IER,IFLG,WK)
      WRITE(6,52) IER,N,K,DET,IDET
      WRITE(6,51) (X(I),I=1,N)
 
      GO TO 100
 
      END
 
C     ------------------------------------------------------
 
C     Solution of a system of linear equations using Gaussian elimination
C     	for a band matrix
C
C     N : (input) Number of equations to be solved
C     KB : (input) Bandwidth of matrix A(I,J)=0 if ABS(I-J)>KB
C     NUM : (input) Number of different sets (each with N equations) of
C		equations to be solved
C     A : (input/output) The matrix of coefficient of size LJ*(3*KB+1)
C		A(I,J-I+KB+1) is the coefficient of x_j in Ith equation
C		at output it will contain the triangular decomposition
C     X : (input/output) The matrix containing right hand sides (size LJ*NUM)
C		X(I,J) is the Ith element of Jth right hand side
C		at output it will contain the solutions
C     DET, IDET : (output) The determinant of the matrix = DET*2**IDET
C     INC : (output) Integer array of length N containing information about
C		interchanges performed during elimination
C     LJ : (input) First dimension of arrays A and X in calling program
C     IER : (output) Error flag, IER=0 signifies successful execution
C		IER=104 implies (N.LE.0 or N.GT.LJ or KB.GT.N)
C		IER=124 implies some pivot turned out to be zero and hence
C	     		matrix must be nearly singular
C     IFLG : (input) Integer variable used as a flag to specify the type
C		of computation required
C     		If IFLG=-1, both elimination and solution are calculated
C     			without pivoting and IFLG is set to 2
C		If IFLG=0, both elimination and solution are computed
C     			with partial pivoting and IFLG is set to 2
C     		If IFLG=1, only elimination is done with pivoting and IFLG is set to 2
C     		If IFLG.GE.2 only solution is calculated, the triangular
C     			decomposition should have been calculated earlier
C     WK : Real array of length 3*KB+1 used as scratch space
C
C	Required routines : None
 
      SUBROUTINE GAUBND(N,KB,NUM,A,X,DET,IDET,INC,LJ,IER,IFLG,WK)
      IMPLICIT REAL*8(A-H,O-Z)
C     For complex matrices use the following statements instead
C      IMPLICIT REAL*8(R)
C      IMPLICIT COMPLEX*16(A-H,S-Z)
 
      DIMENSION A(LJ,3*KB+1),INC(N),X(LJ,NUM),WK(3*KB+1)
 
      IF(N.LE.0.OR.N.GT.LJ.OR.KB.GT.N) THEN
        IER=104
        RETURN
      ENDIF
 
      KB1=KB+1
      IER=124
      IF(IFLG.LE.1) THEN
C     Perform elimination
        DO 2000 I=1,N
          DO 2000 J=2*KB+2,3*KB+1
            A(I,J)=0.0
2000    CONTINUE
 
        DET=1.0
        IDET=0
        DO 2600 K=1,N-1
C     Find the maximum element in the Kth column
          R1=0.0
          KM=K
          IF(IFLG.GE.0) THEN
            DO 2200 L=K,MIN(N,K+KB)
              IF(ABS(A(L,K-L+KB1)).GT.R1) THEN
                R1=ABS(A(L,K-L+KB1))
                KM=L
              ENDIF
2200        CONTINUE
          ENDIF
 
          INC(K)=KM
          IF(KM.NE.K) THEN
C     Interchange the rows if needed
            DO 2300 L=K,MIN(N,2*KB+K)
              WK(L-K+1)=A(K,L-K+KB1)
2300        CONTINUE
            DO 2400 L=K,MIN(N,2*KB+K)
              A(K,L-K+KB1)=A(KM,L-KM+KB1)
2400        A(KM,L-KM+KB1)=WK(L-K+1)
            DET=-DET
          ENDIF
 
          DET=DET*A(K,KB1)
          IF(A(K,KB1).EQ.0.0) RETURN
C     To check for singular or nearly singular matrices replace this
C     statement by, where REPS is approximately \hcross*Max(A(I,J))
C         IF(ABS(A(K,KB1)).LT.REPS) RETURN
          IF(DET.NE.0.0) THEN
 
C     Scale the value of the determinant DET
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
 
          DO 2500 L=K+1,MIN(N,K+KB)
            A(L,K-L+KB1)=A(L,K-L+KB1)/A(K,KB1)
            DO 2500 L1=K+1,MIN(N,2*KB+K)
2500      A(L,L1-L+KB1)=A(L,L1-L+KB1)-A(L,K-L+KB1)*A(K,L1-K+KB1)
2600    CONTINUE
        DET=DET*A(N,KB1)
        INC(N)=N
C     If pivot is zero then return, IER has been set to 124
        IF(A(N,KB1).EQ.0.0) RETURN
C     To check for singular or nearly singular matrices replace this
C     statement by, where REPS is approximately \hcross*Max(A(I,J))
C         IF(ABS(A(N,KB1)).LT.REPS) RETURN
 
        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF
 
      IER=0
C     Solution for the NUM different right-hand sides
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
C     Forward substitution
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF
          DO 3000 L=K+1,MIN(N,K+KB)
3000    X(L,J)=X(L,J)-A(L,K-L+KB1)*X(K,J)
 
C     back-substitution
        X(N,J)=X(N,J)/A(N,KB1)
        DO 3300 K=N-1,1,-1
          DO 3200 L=MIN(N,K+2*KB),K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L-K+KB1)
3300    X(K,J)=X(K,J)/A(K,KB1)
5000  CONTINUE
      END
 
C     ---------------------------------------------------------------
 
C	To generate uniformly distributed random numbers in interval (0,1)
C
C	ISEED : (input/output) is an integer value used as the seed
C		It should be initialised to negative value before first call
C		and should not be modified between successive calls.
C
C	Required routines : None

      FUNCTION RANF(ISEED)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(M1=714025,IA1=1366,IC1=150889)
      PARAMETER(M2=214326,IA2=3613,IC2=45289)
      PARAMETER(M3=139968,IA3=3877,IC3=29573,ISH=43)
      DIMENSION RAN(ISH)
      SAVE
      DATA IFLG/0/

C	Initialise on first call or when ISEED<0
      IF(ISEED.LT.0.OR.IFLG.EQ.0) THEN
        IFLG=1
        RM1=1.D0/M1
        RM2=1.D0/M2

C	Seeds for the three random number generators
        IS1=MOD(-ISEED,M1)
        IS2=MOD(IA1*IS1+IC1,M1)
        IS3=MOD(IA2*IS2+IC2,M2)
        ISEED=1

C	Store ISH random numbers in the array RAN
        DO 1000 J=1,ISH
          IS1=MOD(IA1*IS1+IC1,M1)
          IS2=MOD(IA2*IS2+IC2,M2)
          RAN(J)=(FLOAT(IS1)+FLOAT(IS2)*RM2)*RM1
1000    CONTINUE
      ENDIF

      IS1=MOD(IA1*IS1+IC1,M1)
      IS2=MOD(IA2*IS2+IC2,M2)
      IS3=MOD(IA3*IS3+IC3,M3)
C	Select a random entry from RAN and store a new number in its place
      I=1+(ISH*IS3)/M3
      RANF=RAN(I)
      RAN(I)=(FLOAT(IS1)+FLOAT(IS2)*RM2)*RM1
      END
