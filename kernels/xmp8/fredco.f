C     PROGRAM TO SOLVE FREDHOLM EQUATION USING COLLOCATION METHOD

      PROGRAM INTEQ
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION X(65),WK(65,65),INT(65),F(65)
      EXTERNAL FG,FKER,PHI,PSI

C     EXAMPLE 12.4 : FREDHOLM EQUATION OF THE FIRST KIND

51    FORMAT('   IER =',I4,5X,'N =',I3,5X,'COEF. =',1P3D14.6/
     1       (2X,5D14.6))
52    FORMAT('    X =',1PD14.6,5X,'F(X) =',2D14.6)

      REPS=1.D-7
      AEPS=1.D-8
      A=0.0
      B=1.0
      IQ=0
      IT=1

C     For collocation method the number of points is equal to the
C     number of basis functions.
100   PRINT *,'TYPE N=NO. OF PTS     (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0) STOP

C     CHOOSE UNIFORMLY SPACED POINTS AS THE COLLOCATION POINTS

      DO 1000 I=1,N
        X(I)=DFLOAT(I)/N
1000  CONTINUE

      CALL FREDCO(N,A,B,F,X,REPS,AEPS,WK,INT,IQ,IT,IER)
      WRITE(6,51) IER,N,(F(J),J=1,N)

C     CALCULATING THE SOLUTION AT SOME SELECTED POINTS

      DO 2000 I=1,5
        XI=(I-1)*0.25D0
        FI=0.0
        DO 1500 J=1,N
1500    FI=FI+PHI(J,XI)*F(J)
        WRITE(6,52) XI,FI
2000  CONTINUE

      GO TO 100
      END

C     -------------------------------------------------

C	To integrate a function over finite interval using adaptive control
C	of step size
C
C	RINT : (output) Calculated value of the integral
C	XL : (input) The lower limit
C	XU : (input) The upper limit
C	REPS : (input) The required relative accuracy
C	AEPS : (input) The required absolute accuracy
C		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
C	DIF : (output) estimated (absolute) error achieved by the subroutine
C	F : (input) Name of the function routine to calculate the integrand
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=31 implies specified accuracy was not achieved on
C			at least one subinterval
C		IER=32 implies that this failure occurred more than IFMAX (=5) times
C		IER=325 implies that subroutine failed to attain required
C			accuracy using NMAX function evaluations
C		In all cases DIF will contain the estimated accuracy
C	NPT : (output) Number of function evaluations used by subroutine
C	NMAX : (input/output) Maximum number of function evaluations to be tried
C		If NMAX.LE.0 it is set to MAXPT (=100000)
C
C		FUNCTION F(X) must be supplied by the user.
C
C	Required routines : KRONRD (or GAUS16), F

      SUBROUTINE ADPINT(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(IPMAX=100,IFMAX=5,MAXPT=100000)
      EXTERNAL F
      DIMENSION XU1(IPMAX)

      IER=0
      IFAIL=0
      RINT=0.0
      DIF=0.0
      IF(XL.EQ.XU) RETURN
      IF(NMAX.LE.0) NMAX=MAXPT
      AEPSL=AEPS
      NPT=0
      RL=XL
      RU=XU
      IU=0

C	To evaluate the integral over [RL,RU]
1000  CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
C1000  CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RM=0.5*(RL+RU)
C	Q=.TRUE. if the interval cannot be divided further
      Q=IU.GE.IPMAX.OR.RM.EQ.RL.OR.RM.EQ.RU

      IF(DIF0.LT.MAX(ABS(FINT)*REPS,AEPSL).OR.Q) THEN
C	Accept the value of FINT if adequate convergence or if the interval
C	cannot be subdivided further
        RINT=RINT+FINT
        DIF=DIF+DIF0
        IF(Q.AND.DIF0.GT.MAX(ABS(RINT)*REPS,AEPSL)) THEN
C	Integration fails to converge on this subinterval. Go to the next subinterval
          IER=31
          IFAIL=IFAIL+1
          IF(IFAIL.GT.IFMAX) THEN
C	If failure is frequent then adjust the convergence criterion.
            IER=32
            AEPSL=DIF*0.5
          ENDIF
        ENDIF

C	If all subintervals are exhausted then return
        IF(IU.LE.0) RETURN

C	otherwise try next subinterval
        RL=RU
        RU=XU1(IU)
        IU=IU-1
      ELSE

C	Subdivide the current interval and try again
        IU=IU+1
        XU1(IU)=RU
        RU=RM
      ENDIF

      IF(NPT.LT.NMAX) GO TO 1000
C	If the number of function evaluations has exceeded the limit then
C	try a last call to estimate the integral over the remaining interval
      IER=325
      RU=XU
      CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
C      CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END

C     --------------------------------------------------

C	To solve linear Fredholm equation of first or second kind using
C	collocation method
C	It can be used to solve Volterra equations by defining the kernel
C	to be zero for t>x
C
C	N : (input) Number of points to be used in the collocation method
C	A : (input) Lower limit of the integral
C	B : (input) Upper limit of the integral
C	F : (output) Real array of length N containing the calculated coefficients
C		of expansion
C		solution = SUM_I F(I)*PHI(I,X)
C	X : (input) Real array of length N containing the points to
C		be used for collocation.
C	REPS : (input) Required relative accuracy to which integrals are to
C		be calculated
C	AEPS : (input) Required absolute accuracy to which integrals are to
C		be calculated
C		REPS and AEPS are passed on to subroutine ADPINT for calculating
C		the integrals when IQ=0. Otherwise these variables are not used.
C	WK : Real array of length N*N used as scratch space.
C	IWK : Integer array of length N used as a scratch space
C	IQ : (input) Integer variable to specify the treatment for integrals
C		PSI(I,X)=Integral[FKER(X,T)*PHI(I,T) dT] over [A,B]
C		If IQ=0, then the integrals are evaluated using subroutine
C			ADPINT, using function routine FUNK to calculate the
C			integrand, which in turn requires, PHI(I,T) and FKER(X,T). 
C		Otherwise the integrals are calculated using a user supplied
C			routine PSI(I,X).
C	IT : (input) Integer variable to specify the type of integral equation
C		If IT=1 Fredholm equation of the first kind is solved
C		If IT=2 Fredholm equation of the second kind is solved
C	IER : (output) The error parameter, IER=0 implies successful execution
C		IER=708 implies N<1, IT>2 or IT.LE.0, No calculations are done.
C		Other values of IER may be set by GAUELM and ADPINT
C
C	FUNCTION FG(X), FUNCTION PHI(I,X) and FUNCTION FKER(X,T) (for IQ=0)
C	or FUNCTION PSI(I,X) (for IQ.NE.0) must be supplied by the user.
C	Names of these function routines are fixed. FG(X) is the right
C	hand side function g(x), FKER(X,T) is the kernel, PHI(I,X) calculates
C	the basis functions phi_i(x), while PSI(I,X) calculates the integrals
C	defined above. The common block ZZFRED is used to pass on the variables
C	to FUNK. XI and II are the values of X and I.
C
C	Required routines : GAUELM, ADPINT, KRONRD, FUNK, FG, FKER, PHI, PSI
C	
      SUBROUTINE FREDCO(N,A,B,F,X,REPS,AEPS,WK,IWK,IQ,IT,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNK
C	To pass arguments to FUNK
      COMMON/ZZFRED/XI,II
      DIMENSION F(N),X(N),WK(N,N),IWK(N)

      IER=0
      IF(N.LT.1.OR.IT.GT.2.OR.IT.LE.0) THEN
        IER=708
        RETURN
      ENDIF
      NMAX=10000

C	Setting up the system of linear equations
      DO 2000 J=1,N
        F(J)=FG(X(J))
        DO 2000 I=1,N
          IF(IQ.EQ.0) THEN
C	Evaluate the integrals numerically, split the range into two 
C	to tackle possible discontinuity at t=X
            XI=X(I)
            II=J
            CALL ADPINT(RI,A,XI,REPS,AEPS,DIF,FUNK,IER,NPT,NMAX)
            IF(IER.GT.100) RETURN
            CALL ADPINT(RI1,XI,B,REPS,AEPS,DIF1,FUNK,IER,NPT1,NMAX)
            WK(I,J)=RI+RI1
            IF(IER.GT.100) RETURN
          ELSE
C	Calculate the integrals PSI(I,X) using user supplied routine
            WK(I,J)=PSI(J,X(I))
          ENDIF
          IF(IT.EQ.2) WK(I,J)=WK(I,J)-PHI(J,X(I))
2000  CONTINUE

C	Solve the resulting system of linear equations
      NUM=1
      LJ=N
      IFLG=0
      CALL GAUELM(N,NUM,WK,F,DET,IWK,LJ,IER,IFLG)
      END

C     ---------------------------------------------------------

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

C   -----------------------------------------------------

C	To integrate a function over a finite interval using Gauss-Kronrod formula
C	For use with ADPINT
C
C	RI : (output) Calculated value of the integral
C	A : (input) The lower limit
C	B : (input) The upper limit
C	DIF : (output) estimated (absolute) error achieved by the subroutine
C	N : (output) Number of function evaluations used by subroutine
C	F : (input) Name of the function routine to calculate the integrand
C
C	FUNCTION F(X) must be supplied by the user
C
C	Required routines : F

      SUBROUTINE KRONRD(RI,A,B,DIF,N,F)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  W7(4),A7(4),WK7(4),WK15(4),AK15(4)

C	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
C	WK7 are the weights for these points in Kronrod formula
C	WK15 and AK15 are the weights and abscissas for the remaining points
C	in Kronrod formula.
C	Because of symmetry only half the points are given.

      DATA W7  /0.12948496616886969327D0, 0.27970539148927666790D0,
     *          0.38183005050511894495D0, 0.41795918367346938775D0/
      DATA A7  /0.94910791234275852452D0, 0.74153118559939443986D0,
     *          0.40584515137739716690D0, 0.0/
      DATA WK7 /0.06309209262997855329D0, 0.14065325971552591874D0,
     *          0.19035057806478540991D0, 0.20948214108472782801D0/
      DATA WK15/0.02293532201052922496D0, 0.10479001032225018383D0,
     *          0.16900472663926790282D0, 0.20443294007529889241D0/
      DATA AK15/0.99145537112081263920D0, 0.86486442335976907278D0,
     *          0.58608723546769113029D0, 0.20778495500789846760D0/

      AT=(B-A)/2.
      BT=(B+A)/2.
      FBT=F(BT)
      R1=W7(4)*FBT
      RI=WK7(4)*FBT
      DO 2000 K=1,3
        F1=F(AT*A7(K)+BT)
        F2=F(BT-AT*A7(K))
C	7-point Gauss-Legendre formula
        R1=R1+W7(K)*(F1+F2)
C	15-point Kronrod formula
        RI=RI+WK7(K)*(F1+F2)
2000  CONTINUE

      DO 2500 K=1,4
2500  RI=RI+WK15(K)*(F(AT*AK15(K)+BT)+F(BT-AT*AK15(K)))

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=15
      END

C     -----------------------------

C	Function routine to calculate the integrand for calculating
C	PSI(I,X) as required by subroutine FREDCO.
C
C	Function FKER(X,T) is the kernel K(x,t) and PHI(I,T) is the Ith
C	basis function, phi_i(t). The argument X and I are passed through
C	the common block.
C
C	FUNCTION FKER(X,T) and FUNCTION PHI(I,T) must be supplied by the user
C
C	Required routines : FKER, PHI

      FUNCTION FUNK(T)
      IMPLICIT REAL*8(A-H,O-Z)
C	To pass parameters from subroutine FREDCO
      COMMON/ZZFRED/X,I

      FUNK=FKER(X,T)*PHI(I,T)
      END

C     ------------------------------------------

      FUNCTION FG(X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.1415926535D0)

C     RHS FUNCTION

      FG=(EXP(1.+X)-1.)/(X+1)
      END

C     ---------------------------------------

      FUNCTION FKER(X,T)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE KERNEL

      FKER=EXP(X*T)
      END

C     ------------------------------

      FUNCTION PSI(I,X)
      IMPLICIT REAL*8(A-H,O-Z)

C     DUMMY FUNCTION, SINCE INTEGRALS ARE EVALUATED NUMERICALLY

      PSI=0.0
      END

C     -------------------------------------

      FUNCTION PHI(I,X)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE BASIS FUNCTIONS

      PHI=0.0
      IF(X.EQ.0.0) RETURN
      PHI=X**(I-1)
      END
