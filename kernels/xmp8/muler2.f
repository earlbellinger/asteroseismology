C     PROGRAM TO FIND COMPLEX ROOTS OF A NONLINEAR EQUATION USING
C     MULLER'S METHOD WITH DEFLATION FOR KNOWN ROOTS OR SECANT METHOD
C     THE FUNCTION IS CALCULATED IN A SCALED FORM

C	If we try to start from larger zeros, then Muller fails to
C	find them as iteration always tends to go towards small X.
C	But if we start with small starting values then all zeros
C	can be determined without much difficulty as the smaller
C	zeros are then deflated.

      PROGRAM ZERO
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      EXTERNAL CDET
      DIMENSION CZERO(20),CX(20)

C     EXERCISE 7.44:  EIGENVALUES OF HILBERT MATRIX

51    FORMAT('   IER =',I4,5X,I3,' ZEROS CALCULATED')
52    FORMAT(1P2D14.6,2X,2D14.6)
53    FORMAT('  STARTING VALUES'/(1P2D14.6,2X,2D14.6))
54    FORMAT('   IER =',I4,5X,'ZERO =',1P2D14.6)

      REPS=1.D-14
      AEPS=1.D-27
      RMX=1.D2
      NZ=0

C     THE ZEROS WILL KEEP ACCUMULATING IN THE ARRAY CZERO

100   PRINT *,'TYPE NUM=NO. OF ZEROS TO BE TRIED  (QUITS WHEN NUM.LE.0)' 
      READ *,NUM
      IF(NUM.LE.0) STOP
      PRINT *,'TYPE NUM COMPLEX STARTING VALUES'
      READ *,(CX(I),I=1,NUM)
      WRITE(6,53) (CX(I),I=1,NUM)
      CALL ZROOT2(NUM,CX,CZERO,NZ,REPS,AEPS,IER,RMX,CDET)
      WRITE(6,51) IER,NZ
      WRITE(6,52) (CZERO(I),I=1,NZ)

C     Use secant iteration with same starting values, SECANC_2 does
C     not use deflation

      DO 1000 I=1,NUM
        CX0=CX(I)
        CALL SECANC_2(CX0,RMX,CZ,REPS,AEPS,IER,CDET)
        WRITE(6,54) IER,CZ
1000  CONTINUE

      GO TO 100
      END

C     -------------------------------------------------------

C	Solution of a system of linear equations with complex coefficients
C	using Crout's algorithm with partial pivoting
C
C	N : (input) Number of equations to be solved
C	NUM : (input) Number of different sets (each with N equations) of
C	         equations to be solved
C	A : (input/output) The matrix of coefficient of size LJ*N
C	         A(I,J) is the coefficient of x_J in Ith equation
C	     at output it will contain the triangular decomposition
C	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
C	         X(I,J) is the Ith element of Jth right hand side
C	     at output it will contain the solutions
C	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
C	INC : (output) Integer array of length N containing information about
C		interchanges performed during elimination
C	LJ : (input) First dimension of arrays A and X in calling program
C	IER : (output) Error flag, IER=0 signifies successful execution
C		IER=102 implies (N.LE.0 or N.GT.LJ) 
C		IER=122 implies some pivot turned out to be zero and hence
C			matrix must be nearly singular
C	IFLG : (input) Integer parameter which determines the type of computation
C		required.
C		If IFLG.LE.0, both elimination and solution are calculated
C			and IFLG is set to 2
C		If IFLG=1, only elimination is done and IFLG is set to 2
C		If IFLG.GE.2 only solution is calculated, the triangular
C		    	decomposition should have been calculated earlier
C	WK : scratch array of length N
C
C	Required routines : None

      SUBROUTINE CROUT_C(N,NUM,A,X,DET,IDET,INT,LJ,IER,IFLG,WK)
      IMPLICIT COMPLEX*16(A-H,S-Z)
      IMPLICIT REAL*8(O-R)

      DIMENSION A(LJ,N),INT(N),X(LJ,NUM),WK(N)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=102
        RETURN
      ENDIF
      IER=122

      IF(IFLG.LE.1) THEN
C	Perform LU decomposition
        DO 1500 I=1,N
          R1=0.0
          DO 1200 J=1,N
            IF(ABS(A(I,J)).GT.R1) R1=ABS(A(I,J))
1200      CONTINUE
          WK(I)=R1
C	If any row is zero, then quit
          IF(R1.EQ.0.0) RETURN
1500    CONTINUE

        DET=1.0
        IDET=0
        DO 2600 K=1,N
          R1=0.0
          KM=K
C	Generate the Kth column of L
          DO 2200 L=K,N
            D1=A(L,K)
            DO 2100 L1=1,K-1
2100        D1=D1-A(L,L1)*A(L1,K)
            A(L,K)=D1

C	Finding the pivot
            R2=ABS(D1/WK(L))
            IF(R2.GT.R1) THEN
              R1=R2
              KM=L
            ENDIF
2200      CONTINUE

          INT(K)=KM
C	Interchange the rows if needed
          IF(KM.NE.K) THEN
            DET=-DET
            DO 2300 L=1,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
2300        A(KM,L)=T1
            T1=WK(K)
            WK(K)=WK(KM)
            WK(KM)=T1
          ENDIF

          DET=DET*A(K,K)
C	If the pivot is zero, then quit
          IF(A(K,K).EQ.0.0) RETURN
C	To check for singular or nearly singular matrices replace this statement by
C         IF(ABS(A(K,K)).LT.REPS) RETURN

          IF(DET.NE.0.0) THEN
C	Scale the value of the determinant DET
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

C	Generate the Kth row of U
          DO 2500 L=K+1,N
            D1=A(K,L)
            DO 2400 L1=1,K-1
2400        D1=D1-A(K,L1)*A(L1,L)
            A(K,L)=D1/A(K,K)
2500      CONTINUE
2600    CONTINUE
        IER=0

        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
C	Solution for NUM different right-hand sides
      DO 5000 J=1,NUM
C	Forward substitution
        DO 3000 K=1,N
          IF(K.NE.INT(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INT(K),J)
            X(INT(K),J)=T1
          ENDIF

          D1=X(K,J)
          DO 2800 L=1,K-1
2800      D1=D1-A(K,L)*X(L,J)
          X(K,J)=D1/A(K,K)
3000    CONTINUE

C	Back-substitution
        DO 3300 K=N-1,1,-1
          D1=X(K,J)
          DO 3200 L=N,K+1,-1
3200      D1=D1-X(L,J)*A(K,L)
3300    X(K,J)=D1
5000  CONTINUE
      END

C     ------------------------------------------------------------

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

C      ---------------------------------------------------------

C	Complex zero of a given function using secant iteration 
C	Function is calculated as FX*2**JF
C
C	X0 : (input) Initial guess for the zero
C	R : (input) limiting magnitude for the zero
C	X : (output) Computed value of the zero
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=40 implies that function value is equal at two points
C			and it is not possible to continue the iteration
C		IER=402 implies ABS(X0)>R in which case no calculations are done
C		IER=422 implies that iteration goes outside the specified limits
C		IER=423 implies that iteration failed to converge to specified accuracy
C	FUN : (input) Name of the function routine to calculate the function
C		FUNCTION FUN(X,JF) must be supplied by the user.
C		The function value should be FUN*2**JF
C
C	Required routines : FUN

      SUBROUTINE SECANC_2(X0,R,X,REPS,AEPS,IER,FUN)
      IMPLICIT REAL*8(A-C,O-T)
      IMPLICIT COMPLEX*16(D-H,U-Z)
      PARAMETER(NIS=75)

      IER=0
      IF(ABS(X0).GT.R) THEN
C	If X0 is outside the specified interval (XL,XU) then quit
        IER=402
        RETURN
      ENDIF

      X=X0
C	Select the increment for the next point X+DX
      DX=X*0.001D0
      F1=0.0
      JF1=0

      DO 1000 L=1,NIS
        F=FUN(X,JF)
        DX1=DX
        F1=F1*2.D0**(JF1-JF)

        IF(F1-F.EQ.0.0) THEN
          IF(F.EQ.0.0) RETURN
C	If F1=F and F.NE.0, then quit
          IER=40
          RETURN
        ENDIF

C	The secant iteration
        IF(L.GT.1) DX=DX1*F/(F1-F)
        X=X+DX
        F1=F
        JF1=JF

        IF(ABS(DX).LT.MAX(REPS*ABS(X),AEPS).AND.L.GT.2) RETURN
        IF(ABS(X).GT.R) THEN
C	If iteration goes outside the specified limits (XL,XU), then quit
          IER=422
          RETURN
        ENDIF
1000  CONTINUE

C	The iteration fails to converge
      IER=423
      END

C      ---------------------------------------------------------

C	Complex zeros of a given function using Muller iteration with deflation
C	For use with MULER2 which requires function in form CF*2**IDET
C
C	N : (input) Number of zeros to be determined
C	CX : (input) Complex array of length N containing the initial guesses
C		for the zeros
C	CZERO : (input/output) Complex array of length N+NZ containing the
C		computed values of the zeros. The first NZ zeros which are
C		already known must be supplied as input while other zeros will be added
C	NZ : (input/output) Number of known zeros. At input it should
C		contain the number already known which are included in
C		array CZERO. This number will be incremented as more zeros
C		are determined successfully.
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=41 implies that iteration did not converge to satisfactory accuracy for
C			at least one zero, but value is acceptable at lower accuracy.
C		IER=429 implies that iteration failed to find at least one zero.
C			The number of zeros found can be checked from the value of NZ.
C	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
C		when ABS(CZ) > RMAX
C	CF : (input) Name of the function routine to calculate the function
C		FUNCTION CF(CZ,JF) must be supplied by the user, where CF and CZ
C		are complex. To avoid overflow and underflow the function is
C		calculated in the form CF*2**JF.
C
C	Required routines : MULER2, CF

      SUBROUTINE ZROOT2(N,CX,CZERO,NZ,REPS,AEPS,IER,RMAX,CF)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      EXTERNAL CF
      DIMENSION CZERO(NZ+N),CX(N)

      IF(N.LE.0) RETURN
      IF(NZ.LT.0) NZ=0
      IER=0

      DO 1000 I=1,N
C	The starting values for Muller's iteration
        CX3=CX(I)
        CDX=0.01D0*CX3
        IF(I.LT.N.AND.CDX.EQ.0.0) CDX=0.1D0*(CX(I+1)-CX(I))
        IF(I.GT.1.AND.CDX.EQ.0.0) CDX=0.1D0*(CX(I-1)-CX(I))
        IF(CDX.EQ.0.0) CDX=1.D-3
C	These values may need to be changed in some cases
        CX2=CX3+CDX
        CX1=CX3-CDX

C	Find the next zero
C        CALL MULLER(CX1,CX2,CX3,REPS,AEPS,IER1,CF,NZ,CZERO,RMAX)
C	For MULER use the preceding statements instead
        IER1=0
500     CALL MULER2(CX1,CX2,CX3,REPS,AEPS,IER1,CF0,CX0,IX,NZ,CZERO,RMAX)
        IF(IER1.LT.0) THEN
          CF0=CF(CX0,IX)
          GO TO 500
        ENDIF

        IF(IER1.LT.100) THEN
C	The zero is accepted
          NZ=NZ+1
          CZERO(NZ)=CX3
          IF(IER.EQ.0.AND.IER1.NE.0) IER=41
        ELSE
          IER=429
        ENDIF
1000  CONTINUE
      END

C      ---------------------------------------------------------

      FUNCTION CDET(CX,IDET)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      DIMENSION CR(20),CW(20),CA(20,20),INT(20)

C     FUNCTION ROUTINE TO CALCULATE THE DETERMINANT OF HILBERT MATRIX
C     THE FUNCTION IS REAL, BUT FOR MULLER'S METHOD IT IS NECESSARY TO
C     USE COMPLEX ARITHMETIC

C     THE ORDER OF MATRIX
      N=10
      NUM=1
      LJ=20
      IFLG=1

C     SETTING UP THE HILBERT MATRIX  (H-X I)
      DO 1000 I=1,N
        DO 800 J=1,N
          CA(I,J)=(1.D0,0.D0)/(I+J-1.)
800     CONTINUE
        CA(I,I)=CA(I,I)-CX
1000  CONTINUE

      CALL CROUT_C(N,NUM,CA,CR,CDET,IDET,INT,LJ,IER,IFLG,CW)
      END

