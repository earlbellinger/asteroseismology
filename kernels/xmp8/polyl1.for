C     PROGRAM TO OBTAIN POLYNOMIAL L1-APPROXIMATION FOR DISCRETE DATA
C     BOTH POLYL1 AND LINL1 ARE USED TO CALCULATE THE SAME APPROXIMATION

      PROGRAM L1POLY
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION B(20),WK(9000),IWK(230),F(200),X(200),g(100,200)

C     EXAMPLE 9.13
C     TO GENERATE DATA SET WITH RANDOM ERROR
      FM(Y)=((231*Y*Y-315)*Y*Y+105)*Y*Y-5+RANGAU(IS)*1.D-3

51    FORMAT('   IER =',I4,5X,'DEGREE =',I3,5X,'NO. OF PTS =',I4,
     1       5X,'ERROR SUM =',1PD14.6/'   COEF. :',4D14.6/(2X,5D14.6))
52    FORMAT('   USING LINL1 :   IER =',I4,5X,
     1       5X,'ERROR SUM =',1PD14.6/'   COEF. :',4D14.6/(2X,5D14.6))

      EPS=1.D-9
      IG=100
      XL=0.0
      IS=2

100   PRINT *,'TYPE M=DEGREE,  N=NO. OF DATA PTS'
      PRINT *,'         (QUITS WHEN N.EQ.0)'
      READ *,M,N
      IF(N.LE.0) STOP

C     GENERATING INPUT DATA SET
      H=1.D0/(N-1.)
      DO 2000 I=1,N
        XI=XL+(I-1)*H
        X(I)=XI
        F(I)=FM(XI)

C	THE BASIS FUNCTION FOR USE WITH LINL1
        G(1,I)=1.0
        DO 2000 J=1,M
          G(J+1,I)=XI**J
2000  CONTINUE

      CALL POLYL1(M,N,B,X,F,EPS,ESUM,IER,WK,IWK)
      WRITE(6,51) IER,M,N,ESUM,(B(I),I=1,M+1)

C	NO. OF BASIS FUNCTIONS
      M1=M+1
      CALL LINL1(M1,N,B,F,G,IG,EPS,ESUM,IER,WK,IWK)
      WRITE(6,52) IER,ESUM,(B(I),I=1,M+1)
      GO TO 100
      END

C     -------------------------------------------

C	To calculate coefficients of linear L_1 approximation in terms of
C	specified basis functions for a tabulated function
C
C	M : (input) Number of basis functions in the required approximation
C	N : (input) Number of points in the table of values. N> M
C	A : (output) Real array of length M+1 containing the coefficients
C		of approximation. A(I) is the coefficient of Phi_I(x) in the
C		approximation
C	F : (input) Real array of length N containing the function values
C	G : (input) Real array of length IG*N, containing the values of
C		basis functions at each point in table. G(I,J) is the
C		value of Phi_I(x) at Jth tabular point.
C	IG : (input) First dimension of array G as declared in the calling
C		program, IG .GE. M
C	EPS : (input) Estimate of roundoff error, passed on to SIMPL1
C	ESUM : (output) L_1 norm of the residual at the calculated approximation
C	IER : (output) error parameter, IER=0 implies successful execution
C		IER=616 implies that M<1  or N<M+1
C			in this case no calculations are done
C		Other values of IER may be set by SIMPL1
C	WK : Real array of length (N+2)*(M+3) used as scratch space
C	IWK : Integer array of length N+M+3 used as scratch space
C
C	Required routines : SIMPL1
C
 
      SUBROUTINE LINL1(M,N,A,F,G,IG,EPS,ESUM,IER,WK,IWK)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION A(M+2),F(N),IWK(N+M+3),WK(N+2,M+3),X(N),G(IG,N)
 
      IF(M.LE.0.OR.N.LT.M+1) THEN
        IER=616
        RETURN
      ENDIF
 
      IER=0
      LJ=N+2
      NV=M+1
      M3=N
 
C	Setting up the tableau for simplex algorithm
      DO 2100 J=1,M+2
2100  WK(1,J)=0.0
      DO 2300 I=1,N
        FI=F(I)
        SI=SIGN(1.D0,FI)
        IWK(I+1)=I*SI
        WK(I+1,1)=FI*SI
        WK(1,1)=WK(1,1)-FI*SI
        S=0.0
        DO 2200 J=1,M
          TI=G(J,I)*SI
          S=S+TI
          WK(I+1,J+1)=TI
          WK(1,J+1)=WK(1,J+1)-TI
2200    CONTINUE
        WK(I+1,M+2)=-S
        WK(1,M+2)=WK(1,M+2)+S
2300  CONTINUE
 
      DO 2400 J=1,M+1
        IWK(N+1+J)=N+J
2400  CONTINUE
 
      NC=NV+M3
      CALL SIMPL1(WK,LJ,NC,M3,IWK,IWK(M3+1),IER,EPS)

C	L_1 norm of the residual
      ESUM=-WK(1,1)
C	Finding the coefficients from the tableau
      DO 3200 I=M3+2,NC+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=0.0
3200  CONTINUE
      DO 3400 I=2,M3+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=WK(I,1)
3400  CONTINUE
      DO 3600 I=1,M
3600  A(I)=A(I)-A(M+1)
      END

C     -------------------------------------------------------------

C	To calculate coefficients of polynomial L_1 approximation
C	for a tabulated function
C
C	M : (input) Required degree of polynomial
C	N : (input) Number of points in the table of values. N> M+1
C	A : (output) Real array of length M+2 containing the coefficients
C		of approximation. A(I+1) is the coefficient of x**I
C	X : (input) Real array of length N containing the points at which
C		function value is available
C	F : (input) Real array of length N containing the function values at X(I)
C	EPS : (input) Estimate of roundoff error, passed on to SIMPL1
C	ESUM : (output) L_1 norm of the residual at the calculated approximation
C	IER : (output) error parameter, IER=0 implies successful execution
C		IER=616 implies that M<0 or N<M+2
C			in this case no calculations are done
C		Other values of IER may be set by SIMPL1
C	WK : Real array of length (N+2)*(M+3) used as scratch space
C	IWK : Integer array of length N+M+3 used as scratch space
C
C	Required routines : SIMPL1
C
      SUBROUTINE POLYL1(M,N,A,X,F,EPS,ESUM,IER,WK,IWK)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION A(M+2),F(N),IWK(N+M+3),WK(N+2,M+3),X(N)

      IF(M.LT.0.OR.N.LT.M+2) THEN
        IER=616
        RETURN
      ENDIF

      IER=0
      LJ=N+2
      NV=M+2
C	The number of constraints
      M3=N

C	Setting up the tableau for simplex algorithm
      DO 2100 J=1,M+3
2100  WK(1,J)=0.0
      DO 2300 I=1,N
        FI=F(I)
        SI=SIGN(1.D0,FI)
        IWK(I+1)=I*SI
        WK(I+1,1)=FI*SI
        WK(1,1)=WK(1,1)-FI*SI
        WK(I+1,2)=SI
        WK(1,2)=WK(1,2)-SI
        S=SI
        DO 2200 J=1,M
          TI=X(I)**J*SI
          S=S+TI
          WK(I+1,J+2)=TI
          WK(1,J+2)=WK(1,J+2)-TI
2200    CONTINUE
        WK(I+1,M+3)=-S
        WK(1,M+3)=WK(1,M+3)+S
2300  CONTINUE

      DO 2400 J=1,M+2
        IWK(N+1+J)=N+J
2400  CONTINUE

      NC=NV+M3
      CALL SIMPL1(WK,LJ,NC,M3,IWK,IWK(M3+1),IER,EPS)

C	L_1 norm of the residual
      ESUM=-WK(1,1)
C	Finding the coefficients from the tableau
      DO 3200 I=M3+2,NC+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=0.0
3200  CONTINUE
      DO 3400 I=2,M3+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=WK(I,1)
3400  CONTINUE
      DO 3600 I=1,M+1
3600  A(I)=A(I)-A(M+2)
      END

C     --------------------------------------------------------------

C	To solve a linear programming problem in the standard form arising
C		 in L_1 minimisation problems using the simplex method
C
C	A : (input/output) Real array of length IA*(N-M+1) containing
C		the tableau of simplex algorithm
C		A(1,I+1)=c_i, the cost coefficients
C		Rows 2 to M+1 contain constraints with A(j,1)=b_j and A(j,i+1)=a_i
C	IA : (input) First dimension of array A as declared in the calling
C		program. IA .GE. M+2
C	N : (input) Number of variables, each is constraint to be .GE.0
C	M : (input) Number of constraints of form a^T X = b_i .GE. 0
C	ID : (input/output) integer array of length M+1 which contains
C		information about interchange of variables on LHS
C	IV : (input/output) integer array of length N-M+1 which contains
C		information about interchange of variables on RHS
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=63 implies that the objective function is unbounded from below
C		IER=635 implies that the simplex algorithm failed to find
C			the optimal feasible vector
C	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
C		assumed to be zero
C
C	Required routines : None
C
      SUBROUTINE SIMPL1(A,IA,N,M,ID,IV,IER,AEPS)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION A(IA,N-M+1),ID(M+1),IV(N-M+1)
      PARAMETER(NIT=20)

      JF=1
      M1=M+1
      N1=N-M+1
      IER=0

C	The simplex iteration
      DO 4000 IT=1,NIT*(N+M)
C	Finding the minimum reduced cost coefficient
        RMIN=0.0
        K=0
        DO 2000 J=2,N1
          IF(A(JF,J).LT.RMIN) THEN
            RMIN=A(JF,J)
            K=J
          ELSE IF(IV(J).LE.M.AND.2-A(JF,J).LT.RMIN) THEN
            RMIN=2-A(JF,J)
            K=-J
          ENDIF
2000    CONTINUE
        IF(RMIN.GE.-AEPS) RETURN

C	Finding the pivot element
        K1=ABS(K)
        RMIN=0.0
        L=0
        DO 2400 J=2,M+1
          AJ=A(J,K1)
          IF(K.LT.0) AJ=-AJ
          IF(AJ.GT.AEPS) THEN
            R1=A(J,1)/AJ
            IF(R1.LT.RMIN.OR.L.EQ.0) THEN
              RMIN=R1
              L=J
            ENDIF
          ENDIF
2400    CONTINUE
        IF(L.EQ.0) THEN
C	The objective function is unbounded from below
          IER=63
          RETURN
        ENDIF

        IF(K.LT.0) THEN
          A(JF,K1)=-2+A(JF,K1)
          IV(K1)=-IV(K1)
        ENDIF

C	Exchange the variables
        L1=ID(L)
        ID(L)=IV(K1)
        IV(K1)=L1
        DO 3000 J=1,N1
          IF(J.NE.K1) THEN
            R1=A(L,J)/A(L,K1)
            DO 2800 I=1,M1
              IF(I.NE.L) A(I,J)=A(I,J)-A(I,K1)*R1
2800        CONTINUE
          ENDIF
3000    CONTINUE

        R1=ABS(A(L,K1))
        DO 3200 J=1,N1
          IF(J.NE.K1) A(L,J)=A(L,J)/R1
3200    CONTINUE
        DO 3400 I=1,M1
          IF(I.NE.L) A(I,K1)=-A(I,K1)/A(L,K1)
3400    CONTINUE
        A(L,K1)=1./R1
4000  CONTINUE

C	Iteration fails to converge
      IER=635
      END

C     ---------------------------------------------------

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
