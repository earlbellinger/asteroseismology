C     PROGRAM TO SOLVE A SYSTEM OF NONLINEAR EQUATIONS USING NEWTON'S METHOD
C     IT USES DAVIDENKO'S METHOD TO IMPROVE CHANCES OF CONVERGENCE

      PROGRAM NONLIN
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FCN
      DIMENSION X(20),F(20),WK(500),IWK(20)
      DATA (X(I),I=1,8)/-0.25D0,0.1D0,-0.25D0,0.35D0,-0.25D0,0.6D0,
     1   -0.25D0,0.85D0/

C     EXAMPLE 7.14

51    FORMAT(/'  F =',1P5D14.6/(5X,5D14.6))
52    FORMAT('   IER =',I4,'    X =',1P4D14.6/(6X,5D14.6))

      NP=8
      REPS=1.D-7
      AEPS=1.D-6
      CALL DAVIDN(FCN,NP,X,F,REPS,AEPS,IER,WK,IWK)
      WRITE(6,52) IER,(X(I),I=1,NP)
      WRITE(6,51) (F(I),I=1,NP)
      END

C     --------------------------------------------

C	Real root of a system of nonlinear equations using Davidenko's method
C
C	FCN : (input) Name of the subroutine to calculate the vector function
C		A sample routine can be found at the end of this file.
C	NP : (input) Number of variables and equations to be solved
C	X : (input/output) Real array of length NP containing the starting
C		values. After execution it should contain the computed roots
C	F : (output) Real array of length NP containing the function values
C		at the point specified by array X
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(X(I)))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=407 implies NP > NMAX or NP < 1, in which case no
C			calculations are done
C		IER=440 implies that THETA is nonzero even after 100 steps
C		Other values may be set by subroutines NEWTON or BROYDN
C	WK : Scratch real array of length NP*(NP+1) (for use with NEWTON) or
C		length MAX(2*NP*NP,NP*(NP+4)) (for use with BROYDN) 
C	IWK : integer array of length NP used as scratch space
C	THETA : which is read by the subroutine is the parameter introduced
C		into the equations. Normally THETA=1 to start with and
C		initial guess should satisfy the equations with this parameter
C		THETA=0 corresponds to the equations that need to be solved.
C		The value of THETA should be gradually reduced from 1 to 0
C
C	SUBROUTINE FCN(NP,X,F,DF) or SUBROUTINE FCN(NP,X,F)
C	must be supplied by the user. Here X and F are real arrays of length
C	NP with F containing the value of vector function at point X. DF is
C	real array of length NP*NP containing the Jacobian which is required
C	by subroutine NEWTON. If subroutine BROYDN is used then DF is not
C	required. The subroutine FCN should parameterise the equations in
C	terms of THETA as explained in section 7.16. A sample routine is
C	appended at the end of this file.
C	DAVIDN_B is the version to be used with BROYDN
C
C	Required routines : NEWTON (or BROYDN), GAUELM, FCN
C
      SUBROUTINE DAVIDN(FCN,NP,X,F,REPS,AEPS,IER,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=200)
      EXTERNAL FCN
      DIMENSION X(NP),F(NP),WK(*),IWK(NP)
C	To pass on parameters to FCN this common block must be included
C	in subroutine FCN
      COMMON/DNFUN/THETA,X0(NMAX),F0(NMAX)

51    FORMAT(/5X,'THETA =',F6.2/)

      IF(NP.GT.NMAX.OR.NP.LT.1) THEN
        IER=407
        RETURN
      ENDIF

      PRINT *,' TYPE IN THE INITIAL GUESS'
      READ *,(X(I),I=1,NP)

C	To initialise the parameterised function
      THETA=0.0
      CALL FCN(NP,X,F0,WK)
C	For BROYDN use the following form
C     CALL FCN(NP,X,F0)
      DO 500 I=1,NP
500   X0(I)=X(I)
      REPS0=MAX(1.D-3,REPS)
      AEPS0=MAX(1.D-4,AEPS)

      DO 1000 I=1,100
        PRINT *,' TYPE THETA'
        READ *,THETA
        PRINT 51, THETA
        IF(THETA.EQ.0.0) THEN
          REPS0=REPS
          AEPS0=AEPS
        ENDIF

C	Use either NEWTON or BROYDN
        CALL NEWTON(FCN,NP,X,F,REPS0,AEPS0,IER,WK,IWK)
C       CALL BROYDN(FCN,NP,X,F,REPS0,AEPS0,IER,WK,IWK)
        IF(IER.GT.0) RETURN
        PRINT *,' X= ',(X(J),J=1,NP)
        IF(THETA.EQ.0.0) RETURN
1000  CONTINUE

      IER=440
      END

C	-------------------
C
C	For use with BROYDN omit DF.
C
C	SUBROUTINE FCN(NP,X,F,DF)
C	IMPLICIT REAL*8(A-H,O-Z)
C	PARAMETER (NMAX=200)
C	DIMENSION X(NP),F(NP),DF(NP,NP)
C	COMMON/DNFUN/THETA,X0(NMAX),F0(NMAX)
C
C	call the required function to calculate the vector function
C	There is no provision to pass the name FUN, it has to be put explicitly
C	CALL FUN(NP,X,F,DF)
C
C	APPLY PARAMETRISATION IN TERMS OF THETA
C
C	DO 2222 I=1,NP
C	F(I)=F(I)-THETA*F0(I)
C2222	CONTINUE
C	END

C     --------------------------------------------------

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

C	Real root of a system of nonlinear equations using Newton's method
C
C	FCN : (input) Name of the subroutine to calculate the vector function
C	NP : (input) Number of variables and equations to be solved
C	X : (input/output) Real array of length NP containing the starting
C		values. After execution it should contain the computed roots
C	F : (output) Real array of length NP containing the function values
C		at the point specified by array X
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(X(I)))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=441 implies that GAUELM failed to solve the required
C			system of linear equations
C		IER=442 implies that iteration did not converge to required
C			accuracy
C	WK : Real array of length NP*(NP+1) used as scratch space
C	IWK : integer array of length NP used as scratch space
C
C	SUBROUTINE FCN(NP,X,F,DF) must be supplied by the user.
C	Here X and F are real arrays of length NP containing the value of
C	vector function at point X. DF is a real array of length NP*NP 
C	containing the Jacobian. The first dimension of DF in FCN should be NP.
C	The system of equations to be solved are F(I)=0
C
C	Required routines : GAUELM, FCN
C
      SUBROUTINE NEWTON(FCN,NP,X,F,REPS,AEPS,IER,WK,IWK)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION X(NP),F(NP),WK(NP*(NP+1)),IWK(NP)

      LJ=NP
      NUM=1
      IER=0
      NP2=NP*NP

C	Loop for Newton's iteration
      DO 5000 IPAS=1,200
        CALL FCN(NP,X,F,WK)
C	The right hand side of system of equations for Newton's iteration
        DO 2000 I=1,NP
2000    WK(NP2+I)=F(I)

        IFLG=0
        CALL GAUELM(NP,NUM,WK,WK(NP2+1),DET,IWK,LJ,IER1,IFLG)
        IF(IER1.GT.0) THEN
C	If GAUELM fails, then quit
          IER=441
          RETURN
        ENDIF

C	Convergence check
        QCHK=.TRUE.
        DO 3500 J=1,NP
          IF(ABS(WK(NP2+J)).GT.MAX(REPS*ABS(X(J)),AEPS)) QCHK=.FALSE.
          X(J)=X(J)-WK(NP2+J)
3500    CONTINUE
        IF(QCHK) RETURN

5000  CONTINUE

C	Iteration fails to converge
      IER=442
      END

C  ------------------------------------------------------------

      SUBROUTINE FCN(N,X,F,DF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(*),F(*),DF(N,N)
      COMMON/DNFUN/THETA,X0(200),F0(200)

C	The required system of nonlinear equations, F(I)=0

      DO 2000 I=1,N
        FI=1.D0/I**2
        DO 1500 J=1,N,2
          IF(I.EQ.1) THEN
            FI=FI+X(J)
            DF(I,J)=1.
            DF(I,J+1)=0.0
          ELSE
            FI=FI+X(J)*X(J+1)**(I-1)
            DF(I,J)=X(J+1)**(I-1)
            IF(I.GT.2) THEN
              DF(I,J+1)=(I-1)*X(J)*X(J+1)**(I-2)
            ELSE
              DF(I,J+1)=X(J)
            ENDIF
          ENDIF
1500    CONTINUE

C     PARAMETERISATION FOR DAVIDENKO'S METHOD

        F(I)=FI-THETA*F0(I)
2000  CONTINUE
      END
