C     PROGRAM TO SOLVE BOUNDARY VALUE PROBLEMS IN
C     ORDINARY DIFFERENTIAL EQUATIONS USING EXPANSION METHOD WITH B-SPLINES
 
      PROGRAM BVP
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL EQN,BCS
      DIMENSION PAR(10),X(2,801),WK(90000),TX(800),T(801),A(1000)
 
C     EXAMPLE 11.13
 
51    FORMAT('   IER =',I4,5X,'NO. OF KNOTS =',I5,5X,
     1       'ORDER OF B-SPLINES =',I4,5X,'NO. OF PTS =',I5/5X,
     2       'LAMDA =',1PD14.6/11X,'T',20X,'SOLUTION')
52    FORMAT(I4,2X,1PD14.6,2X,2D14.6)
 
      M=2
      ML=1
      IFLAG=0
      REPS=1.D-8
 
C     PASS THE PARAMETER LAMDA IN THE EQUATION VIA THE ARRAY PAR
 
100   PRINT *,'TYPE NK=NO. OF KNOTS, K=ORDER OF B-SPLINE'
      PRINT *,'N = NO. OF POINTS,  LAMDA        (QUITS WHEN N.LE.0)'
      READ *,NK,K,N,PAR(1)
      IF(N.LE.0) STOP
 
C     SET UP THE KNOTS WITH UNIFORM SPACING AND THE INITIAL VALUES
 
      H=1.D0/(NK-1)
      DO 500 I=1,NK
500   T(I)=(I-1)*H
      DO 800 I=1,NK+K-2
        A(I)=1.0
        A(I+NK+K-2)=1.0
800   CONTINUE
      CALL BSPODE(NK,K,M,ML,PAR,X,A,T,N,TX,EQN,BCS,WK,IFLAG,REPS,IER)
      WRITE(6,51)IER,NK,K,N,PAR(1)
      DO 1000 I=1,N,10
        WRITE(6,52) I,TX(I),X(1,I),X(2,I)
1000  CONTINUE
      GO TO 100
      END
 
C     ------------------------------------------------------
 
C	To calculate function value using B-spline expansion
C
C	N : (input) Number of knots to define B-splines
C	X : (input) Real array of length N+2K+1 containing the knots.
C		The knots must be distinct and in ascending order.
C	K : (input) Order of B-splines, K=4 for cubic B-splines
C	NDERIV : (input) Number of derivatives required
C		For NDERIV.LE.0 only function value is calculated
C		For NDERIV=1 first derivative is also calculated
C		For NDERIV>1 both first and second derivatives are calculated
C	WT : (input) Coefficients of B-spline expansion
C	X0 : (input) The point at which expansion has to be evaluated
C	DF : (output) First derivative of function at X0
C	DDF : (output) Second derivative of function at X0
C	WK : Scratch array of length 4N+5K+2
C	IER : (output) Error parameter, IER=0 implies successful execution
C		Nonzero values of IER may be set by BSPLIN which is called
C
C	BSPEVL = SUM_{i=1}^{N+K-2} WT(I) \phi_i(X0)
C	where \phi_i(x) are B-spline basis functions on knots X
C
C	Required routines : BSPLIN

      FUNCTION BSPEVL(N,X,K,NDERIV,WT,X0,DF,DDF,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N+K),WT(N+K-2),WK(4*N+5*K+2)
 
      BSPEVL=0.0
      NK=(N+K)
      CALL BSPLIN(X,N,K,X0,NDERIV,WK,WK(NK),WK(2*NK),LEFT,IER,WK(3*NK))
      IF(IER.GT.100) RETURN
 
      F=0.0
      DF=0.0
      DDF=0.0
      N1=N+K-1
      N2=2*(N+K)-1
      DO 2000 I=LEFT,LEFT+K-1
        F=F+WT(I)*WK(I)
        DF=DF+WT(I)*WK(N1+I)
        DDF=DDF+WT(I)*WK(N2+I)
2000  CONTINUE
      BSPEVL=F
      END
 
C     ----------------------------------------------------------
 
C	To calculate the B-spline basis functions at a specified point
C
C	X : (input) Real array of length NX containing the knots.
C		The knots must be distinct and in ascending order.
C	NX : (input) Number of knots
C	K : (input) Order of B-spline, 0< K <KMAX+1
C		K=4 gives cubic B-splines
C	XB : (input) The point at which B-spline basis functions are to be evaluated
C	NDERIV : (input) Number of derivatives required
C		NDERIV.LE.0 only B-splines are calculated
C		NDERIV=1 first derivative is also calculated
C		NDERIV>1 first and second derivatives are also calculated
C	B : (output) Array of length NX+K-2 containing the value of
C		B-spline basis functions
C	DB : (output) Array of length NX+K-2 containing the value of
C		the first derivative of B-spline basis functions (if NDERIV>0)
C	DDB : (output) Array of length NX+K-2 containing the value of
C		the second derivative of B-spline basis functions (if NDERIV>1)
C	LEFT : (output) XB is located between X(LEFT) and X(LEFT+1)
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=26 implies XB > X(NX)
C		IER=27 implies XB < X(1)
C		IER=203 implies NX<2, K<1 or K>KMAX
C	WK : Real array of length NX+2K+1 used as scratch space
C
C	Required routines : None

      SUBROUTINE BSPLIN(X,NX,K,XB,NDERIV,B,DB,DDB,LEFT,IER,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (KMAX=20)
      DIMENSION X(NX),B(NX+K-2),DR(KMAX),DL(KMAX),DB(NX+K-2),
     1       DDB(NX+K-2),WK(-K:NX+K)
 
      SAVE
      DATA LOW/0/
 
      IF(NX.LE.1.OR.K.LT.1.OR.K.GT.KMAX) THEN
        IER=203
        RETURN
      ENDIF
 
      IER=0
 
      IF(LOW.LT.1.OR.LOW.GE.NX) THEN
C	If the previous value of LOW is inadmissible, set the range to (1,N)
        LOW=1
        IGH=NX
      ELSE
        IGH=LOW+1
      ENDIF
 
1000  IF((XB.LT.X(LOW).AND.XB.LT.X(IGH)).OR.
     1  (XB.GT.X(LOW).AND.XB.GT.X(IGH))) THEN
C	Extend the range
        IF(XB.GT.X(LOW)) THEN
C	Extend the range on higher side
          IF(IGH.GE.NX) THEN
            IER=26
            LOW=NX-1
          ELSE
            NIGH=MIN(NX,IGH+2*(IGH-LOW))
            LOW=IGH
            IGH=NIGH
            GO TO 1000
          ENDIF
        ELSE
C	Extend the range on lower side
          IF(LOW.LE.1) THEN
            IER=27
          ELSE
            NIGH=LOW
            LOW=MAX(1,LOW-2*(IGH-LOW))
            IGH=NIGH
            GO TO 1000
          ENDIF
        ENDIF
      ELSE
 
C	Once the point is bracketed between two tabular points locate it by bisection
1500    IF(IGH-LOW.GT.1.AND.XB.NE.X(LOW)) THEN
          MID=(LOW+IGH)/2
          IF(XB.LE.X(MID).EQV.XB.LE.X(LOW)) THEN
            LOW=MID
          ELSE
            IGH=MID
          ENDIF
          GO TO 1500
        ENDIF
      ENDIF
 
C	Evaluate the B-spline basis functions

C	Define the extra knots on either side of table
C	Note that the program assumes knots from -K+2 to NX+K-1
C	and the B-splines B_{i,k}, i ranges from 1 to NX+K-2
C	The knots are stored in scratch array WK.
      DO 1700 I=1,NX
        WK(I)=X(I)
1700  CONTINUE
      DO 1800 I=1,K
        WK(1-I)=X(1)
        WK(NX+I)=X(NX)
1800  CONTINUE

      DO 1900 I=1,NX+K-2
        B(I)=0.0
        DB(I)=0.0
        DDB(I)=0.0
1900  CONTINUE
      LEFT=LOW
      LX=LOW-1
      J=1
      B(LX+1)=1.
 
C	The recurrence relation for B-splines
      DO 3000 J=1,K-1
        DR(J)=WK(LOW+J)-XB
        DL(J)=XB-WK(LOW+1-J)
        T1=0.0
        DO 2000 I=1,J
          T2=B(LX+I)/(DR(I)+DL(J+1-I))
          B(LX+I)=T1+T2*DR(I)
          T1=T2*DL(J+1-I)
2000    CONTINUE
        B(LX+J+1)=T1

C	Calculate the first derivative using recurrence relations
        IF(J.EQ.K-2.AND.NDERIV.GT.0) THEN
          T1=0.0
          DO 2200 I=1,J+1
            T2=B(LX+I)/(WK(LOW+I)-WK(LOW+I+1-K))
            DB(LX+I)=(K-1)*(T1-T2)
            T1=T2
2200      CONTINUE
          DB(LX+J+2)=(K-1)*T1
        ENDIF

C	Calculate the second derivative using recurrence relations
        IF(J.EQ.K-3.AND.NDERIV.GT.1) THEN
          T2=0.0
          P1=0.0
          DO 2400 I=1,J+1
            T3=B(LX+I)/(WK(LOW+I)-WK(LOW+I+2-K))
            P2=(T2-T3)/(WK(LOW+I)-WK(LOW+I-K+1))
            DDB(LX+I)=(K-2)*(K-1)*(P1-P2)
            T2=T3
            P1=P2
2400      CONTINUE
          P2=T2/(WK(LOW+J+2)-WK(LOW+J+3-K))
          DDB(LX+J+2)=(K-2)*(K-1)*(P1-P2)
          DDB(LX+J+3)=(K-2)*(K-1)*P2
        ENDIF
3000  CONTINUE
 
C	For K=2 the first derivative has to be calculated outside the loop
      IF(K.EQ.2.AND.NDERIV.GT.0) THEN
        T2=1/(WK(LOW+1)-WK(LOW+2-K))
        DB(LX+1)=-T2
        DB(LX+2)=T2
      ENDIF
 
C	For K=3 the second derivative has to be calculated outside the loop
      IF(K.EQ.3.AND.NDERIV.GT.1) THEN
        T3=1./(WK(LOW+1)-WK(LOW+3-K))
        P2= -T3/(WK(LOW+1)-WK(LOW-K+2))
        DDB(LX+1)=-2.*P2
        P1=P2
        P2=T3/(WK(LOW+2)-WK(LOW+3-K))
        DDB(LX+2)=2.*(P1-P2)
        DDB(LX+3)=2.*P2
      ENDIF
 
      END
 
C     -------------------------------------------------------
 
C	To solve two-point boundary value problem in ordinary differential
C	equations using expansion method with B-spline basis functions
C
C	NK : (input) Number of knots to be used for calculating B-splines.
C	K : (input) Order of B-splines to be used, K=4 for cubic B-splines
C	M : (input) Number of first order differential equations in the system
C	ML : (input) Number of boundary conditions at the first boundary t=T(1)
C	PAR : (input) Real array to be passed on to EQN and BCS for calculating
C		the equations and boundary conditions. This array is not used
C		by BSPODE, but is merely used to pass on any extra parameters
C		that may be required to specify the equations
C	X : (output) Real array of length M*N containing the solution.
C		X(i,j) is the ith component of solution at jth mesh point
C		TX(j). First dimension of X in calling program must be M
C	A : (input/output) Real array of length (NK+K-2)*M containing the
C		coefficients of expansion. At the time of calling it should
C		contain the initial guess. After execution it will contain
C		the final values. The first dimension of A must be NK+K-2.
C	T : (input) Real array of length NK containing the knots.
C		These points must be in ascending order with T(1) and T(NK)
C		as the two boundaries. 
C	N : (input) The number of mesh points to be used for calculating
C		the coefficients. The solution will be calculated at
C		all these points. N.GE.NK+K-2
C	TX : (input/output) Real array of length N containing the mesh points
C		For IFLAG>1 the mesh points must be supplied, while for
C		other values of IFLAG the routine calculates the mesh point
C		assuming uniform spacing. 
C	EQN : (input) Name of subroutine to specify the differential equation
C		to be solved.
C	BCS : (input) Name of subroutine to calculate the boundary conditions
C		at t=T(1) and T(NK)
C	WK : Real array of length M*(N+1)*(M*(NK+K-2)+7)+(M*(NK+K-2))**2
C		used as scratch space
C	IFLAG : (input) Integer variable used as a flag to decide the type of
C		computation required.
C		IFLAG=0 implies that equations are nonlinear and mesh points
C			are not supplied in array TX. These would be calculated
C			assuming uniform spacing.
C		IFLAG=1 implies that equations are linear and mesh points
C			are not supplied in array TX. These would be calculated
C			assuming uniform spacing.
C		IFLAG=2 implies that equations are nonlinear and mesh points
C			are supplied in array TX.
C		IFLAG=3 implies that equations are linear and mesh points
C			are supplied in array TX.
C	REPS : (input) Required accuracy. This is only used to check convergence
C		of Newton's method for solving the resulting system of equations.
C		The truncation error depends on NK and K and is not
C		controlled in this routine.
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=705 implies that NK<3, M.LE.ML, ML.LE.0, or N.LT.NK+K-2
C			in which case no calculations are done
C		IER=741 implies that Newton's iteration for solving the
C			system of equations failed to converge.
C		Other values may be set by BSPLIN or SVD
C
C	SUBROUTINE EQN and BCS must be supplied by the user.
C	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) calculates the right hand
C		sides for differential equations y'_i=f_i(t,y,par)
C		J is the serial number of mesh point at which calculation
C		is required. M is the number of first order differential
C		equations in the system, ML the number of boundary conditions
C		at the first point, PAR is a real array which can be used
C		to pass on any required parameters to define the equations.
C		A and B are real arrays of length M*M defining the
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
C		M*M which should contain the coefficients of boundary
C		conditions. First ML rows will specify the boundary conditions
C		at t=T1, while remaining rows will specify those at t=TN.
C		G is a real array of length M specifying the boundary conditions
C		G(I)=g_i(T1,PAR,Y1) (I.LE.ML) are the boundary conditions
C		at T1 (g_i=0), while G(I)=g_i(TN,PAR,YN) (I.GT.ML) are the
C		boundary conditions at TN. BC is the Jacobian matrix dg_i/dY(K)
C		Y1 and YN are real arrays of length M specifying the
C		solution at t=T1 and TN.
C
C	Required routines : BSPLIN, BSPEVL, SVD, SVDEVL, EQN, BCS
C
 
      SUBROUTINE BSPODE(NK,K,M,ML,PAR,X,A,T,N,TX,EQN,BCS,WK,IFLAG,
     1            REPS,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=50,EPS=1.D-30)
      DIMENSION PAR(*),X(M,N),T(NK),TX(N),WK(M,N+1,*),A(NK+K-2,M)
 
      NB=NK+K-2
      IER=0
      IF(NK.LT.3.OR.M.LE.ML.OR.ML.LE.0.OR.N.LT.NB) THEN
        IER=705
        RETURN
      ENDIF
 
      IF(IFLAG.LE.1) THEN
C	Setup the mesh assuming uniform spacing
        H=(T(NK)-T(1))/(N-1)
        DO 2000 I=1,N
          TX(I)=T(1)+(I-1)*H
2000    CONTINUE
      ENDIF
      NE=N*M+M
      ME=NB*M
 
C	The iteration loop
      NDERIV=1
      DO 5000 IT=1,NIT
 
C	Setup the equation matrix
        DO 3000 I=1,N
          XB=TX(I)
          CALL BSPLIN(T,NK,K,XB,NDERIV,WK(1,1,ME+1),WK(1,1,ME+2),
     1                WK(1,1,ME+3),LEFT,IER,WK(1,1,ME+4))
          IF(IER.GT.100) RETURN
          DO 2500 J=1,M
            X(J,I)=BSPEVL(NK,T,K,NDERIV,A(1,J),XB,DF,DDF,WK(1,1,ME+3)
     1            ,IER)
            IF(IER.GT.100) RETURN
            X(J,I+1)=DF
2500      CONTINUE
          JI=I
          CALL EQN(JI,M,ML,PAR,WK(1,1,ME+3),WK(1,1,ME+4),X(1,I),
     1                WK(1,1,ME+5),XB)
 
          DO 2800 J=1,M
            S=0
            DO 2600 J1=1,M
              S=S+WK(J,J1,ME+4)*X(J1,I+1)
2600        CONTINUE
            WK(J,I,ME+6)=-S+WK(J,1,ME+5)
            DO 2700 J1=1,M*NB
              JK=(J1-1)/NB+1
              JI=J1-(JK-1)*NB
              IF(JI.EQ.0) JI=NB
              WK(J,I,J1)=WK(J,JK,ME+4)*WK(JI,1,ME+2)-WK(J,JK,ME+3)*
     1                   WK(JI,1,ME+1)
2700        CONTINUE
 
2800      CONTINUE
3000    CONTINUE
 
        T1=T(1)
        CALL BSPLIN(T,NK,K,T1,NDERIV,WK(1,1,ME+1),WK(1,1,ME+2),
     1              WK(1,1,ME+3),LEFT,IER,WK(1,1,ME+4))
        IF(IER.GT.100) RETURN
        T2=T(NK)
        CALL BSPLIN(T,NK,K,T2,NDERIV,WK(1,1,ME+2),WK(1,1,ME+3),
     1              WK(1,1,ME+4),LEFT,IER,WK(1,1,ME+5))
        IF(IER.GT.100) RETURN
 
        CALL BCS(M,ML,PAR,WK(1,1,ME+3),WK(1,1,ME+4),T1,T2,X(1,1),X(1,N))
 
        DO 3500 I=1,M
          WK(I,N+1,ME+6)=-WK(I,1,ME+4)
          DO 3200 J1=1,M*NB
            JK=(J1-1)/NB+1
            JI=J1-(JK-1)*NB
            IF(JI.EQ.0) JI=NB
            IF(I.LE.ML) THEN
              WK(I,N+1,J1)=WK(I,JK,ME+3)*WK(JI,1,ME+1)
            ELSE
              WK(I,N+1,J1)=WK(I,JK,ME+3)*WK(JI,1,ME+2)
            ENDIF
3200      CONTINUE
3500    CONTINUE
 
C	Solve the system of equations using SVD
        CALL SVD(ME,NE,WK,WK(1,1,ME+7),WK(1,1,ME+1),NE,ME,WK(1,1,ME+2),
     1          IER)
        IF(IER.GT.100) RETURN
 
        CALL SVDEVL(ME,NE,WK,WK(1,1,ME+7),WK(1,1,ME+1),NE,ME,
     1              WK(1,1,ME+6),WK(1,1,ME+2),REPS)
 
C	Convergence check
        RERR=0.0
        DO 4000 J=1,NB
          J1=J+1
          IF(J.EQ.NB) J1=J-1
          DO 4000 I=1,M
            JK=J+(I-1)*NB
            XJ=A(J,I)+WK(JK,1,ME+6)
            R2=ABS(XJ)+ABS(X(J,I)-X(J1,I))+EPS
            RE=ABS(WK(JK,1,ME+6)/R2)
            IF(RE.GT.RERR) RERR=RE
            A(J,I)=XJ
4000    CONTINUE
        IF(RERR.LT.REPS.OR.IFLAG.EQ.1.OR.IFLAG.EQ.3) RETURN
 
5000  CONTINUE
      IER=741
      RETURN
 
      END

C	--------------------------------
C	
C	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T)
C       IMPLICIT REAL*8(A,B,D-H,O-Z)
C	DIMENSION A(M,M),B(M,M),Y(M),PAR(*),F(M)
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
C	DIMENSION PAR(*),BC(M,M),G(M),Y1(M),YN(M)
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
 
C     -------------------------------------------
 
C	To calculate the Singular Value Decomposition of a matrix A=U D V-transpose
C
C	N : (input) Number of variables
C	M : (input) Number of equations
C	A : (input/output) Matrix of coefficients of size LA*N
C		After execution it will contain the matrix U
C	V : (output) The matrix V of size LV*N
C	SIGMA : (output) Array of length N, containing the singular values
C	LA : (input) Actual value of first dimension of A in the calling program
C	LV : (input) Actual value of first dimension of V in the calling program
C	E : Scratch array of length N
C	IER : (output) Error parameter, IER=0 if execution is successful
C		IER=12 QR iteration failed to converge to required accuracy
C		IER=105 implies N.LE.0, N.GT.LV, M.LE.0, M.GT.LA, N.GT.M
C
C	Required routines : None

      SUBROUTINE SVD(N,M,A,V,SIGMA,LA,LV,E,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ITMAX=30,REPS=1.D-16)
C	For REAL*4 use REPS=6.E-8
C      PARAMETER(ITMAX=30,REPS=6.E-8)
      DIMENSION A(LA,N),V(LV,N),SIGMA(N),E(N)

      IF(N.GT.M.OR.N.LE.0.OR.M.LE.0.OR.M.GT.LA.OR.N.GT.LV) THEN
        IER=105
        RETURN
      ENDIF

      IER=0
C	Reduction to Bidiagonal form using Householder transformations
      G=0
      RMAX=0

      DO 3000 I=1,N
C	Off-diagonal elements of bidiagonal form
        E(I)=G
        S=0
        DO 1200 J=I,M
1200    S=S+A(J,I)**2

        IF(S.LE.0.0) THEN
C	transformation not required
          G=0
        ELSE
          F=A(I,I)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I)=F-G

          DO 1800 J=I+1,N
            S=0
            DO 1400 K=I,M
1400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 1600 K=I,M
1600        A(K,J)=A(K,J)+F*A(K,I)
1800      CONTINUE
        ENDIF

C	Diagonal elements of bidiagonal form
        SIGMA(I)=G
        S=0
        DO 2000 J=I+1,N
2000    S=S+A(I,J)**2

        IF(S.LE.0.0) THEN
C	Transformation not required
          G=0
        ELSE
          F=A(I,I+1)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I+1)=F-G
          DO 2200 J=I+1,N
C	Temporary storage of intermediate results
2200      E(J)=A(I,J)/H

          DO 2800 J=I+1,M
            S=0
            DO 2400 K=I+1,N
2400        S=S+A(J,K)*A(I,K)
            DO 2600 K=I+1,N
2600        A(J,K)=A(J,K)+S*E(K)
2800      CONTINUE
        ENDIF
        R1=ABS(SIGMA(I))+ABS(E(I))
        IF(R1.GT.RMAX) RMAX=R1
3000  CONTINUE

C	Accumulation of right hand transformation in array V
      DO 4000 I=N,1,-1
        IF(G.NE.0.0) THEN
          H=A(I,I+1)*G
          DO 3200 J=I+1,N
3200      V(J,I)=A(I,J)/H

          DO 3800 J=I+1,N
            S=0
            DO 3400 K=I+1,N
3400        S=S+A(I,K)*V(K,J)
            DO 3600 K=I+1,N
3600        V(K,J)=V(K,J)+S*V(K,I)
3800      CONTINUE
        ENDIF

        DO 3900 J=I+1,N
          V(I,J)=0.0
          V(J,I)=0.0
3900    CONTINUE
        V(I,I)=1
        G=E(I)
4000  CONTINUE

C	Accumulation of left hand transformation overwritten on matrix A
      DO 5000 I=N,1,-1
        G=SIGMA(I)
        DO 4200 J=I+1,N
4200    A(I,J)=0
        IF(G.NE.0.0) THEN
          H=A(I,I)*G

          DO 4700 J=I+1,N
            S=0
            DO 4400 K=I+1,M
4400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 4600 K=I,M
4600        A(K,J)=A(K,J)+F*A(K,I)
4700      CONTINUE

          DO 4800 J=I,M
4800      A(J,I)=A(J,I)/G
        ELSE
          DO 4900 J=I,M
4900      A(J,I)=0.0
        ENDIF
        A(I,I)=A(I,I)+1
5000  CONTINUE

C	Diagonalisation of the bidiagonal form
      AEPS=REPS*RMAX
C	Loop over the singular values
      DO 8000 K=N,1,-1
C	The QR transformation
        DO 7500 ITR=1,ITMAX

C	Test for splitting
          DO 5200 L=K,1,-1
            IF(ABS(E(L)).LT.AEPS) GO TO 6000
            IF(ABS(SIGMA(L-1)).LT.AEPS) GO TO 5400
5200      CONTINUE

C	cancellation of E(L) if L>1
5400      C=0.0
          S=1.0
          DO 5800 I=L,K
            F=S*E(I)
            E(I)=C*E(I)
            IF(ABS(F).LT.AEPS) GO TO 6000
            G=SIGMA(I)
            SIGMA(I)=SQRT(F*F+G*G)
            C=G/SIGMA(I)
            S=-F/SIGMA(I)

            DO 5600 J=1,M
              R1=A(J,L-1)
              R2=A(J,I)
              A(J,L-1)=R1*C+R2*S
              A(J,I)=C*R2-S*R1
5600        CONTINUE
5800      CONTINUE

6000      Z=SIGMA(K)
          IF(L.EQ.K) THEN
C	QR iteration has converged
            IF(Z.LT.0.0) THEN
              SIGMA(K)=-Z
              DO 6200 J=1,N
6200          V(J,K)=-V(J,K)
            ENDIF
            GO TO 8000
          ENDIF

          IF(ITR.EQ.ITMAX) THEN
            IER=12
            GO TO 7500
          ENDIF

C	calculating shift from bottom 2x2 minor
          X=SIGMA(L)
          Y=SIGMA(K-1)
          G=E(K-1)
          H=E(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.*H*Y)
          G=SQRT(1.+F*F)
          IF(F.LT.0.0) G=-G
          F=((X-Z)*(X+Z)+H*(Y/(F+G)-H))/X

C	next QR transformation
          C=1.0
          S=1.0
C	Given's rotation
          DO 7000 I=L+1,K
            G=E(I)
            Y=SIGMA(I)
            H=S*G
            G=C*G
            E(I-1)=SQRT(F*F+H*H)
            C=F/E(I-1)
            S=H/E(I-1)
            F=C*X+S*G
            G=C*G-S*X
            H=S*Y
            Y=C*Y

            DO 6400 J=1,N
              X=V(J,I-1)
              Z=V(J,I)
              V(J,I-1)=C*X+S*Z
              V(J,I)=C*Z-S*X
6400        CONTINUE

            SIGMA(I-1)=SQRT(F*F+H*H)
            IF(SIGMA(I-1).NE.0.0) THEN
              C=F/SIGMA(I-1)
              S=H/SIGMA(I-1)
            ENDIF
            F=C*G+S*Y
            X=C*Y-S*G
            DO 6600 J=1,M
              Y=A(J,I-1)
              Z=A(J,I)
              A(J,I-1)=C*Y+S*Z
              A(J,I)=C*Z-S*Y
6600        CONTINUE
7000      CONTINUE

          E(L)=0
          E(K)=F
          SIGMA(K)=X
7500    CONTINUE
8000  CONTINUE
      END
 
C     -------------------------------------------------------------
 
C	To evaluate the solution of a system of linear equations using SVD
C
C	N : (input) Number of variables
C	M : (input) Number of equations
C	U : (input) array of size LU*N containing the left-hand transformation
C	V : (input) array of size LV*N containing the right-hand transformation
C	SIGMA : (input) array of size N containing the singular values
C	LU : (input) First dimension of array U in the calling program
C	LV : (input) First dimension of array V in the calling program
C	B : (input/output) Array of length M containing the RHS
C		after execution it will contain the solution
C	WK : Scratch array of length N
C	REPS : Relative accuracy. All singular values < REPS*(Max of singular values)
C		will be reduced to zero
C
C	Required routines : None

      SUBROUTINE SVDEVL(N,M,U,V,SIGMA,LU,LV,B,WK,REPS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(LU,N),V(LV,N),SIGMA(N),B(*),WK(N)

C	Finding the largest singular value
      SMAX=0.0
      DO 2000 I=1,N
        IF(SIGMA(I).GT.SMAX) SMAX=SIGMA(I)
2000  CONTINUE

      AEPS=SMAX*REPS
      DO 3000 I=1,N
        S=0.0
C	Only SIGMA(I) > AEPS contribute to the solution
        IF(SIGMA(I).GT.AEPS) THEN
          DO 2400 J=1,M
2400      S=S+U(J,I)*B(J)
          S=S/SIGMA(I)
        ENDIF
        WK(I)=S
3000  CONTINUE

      DO 4000 I=1,N
        S=0.0
        DO 3400 J=1,N
3400    S=S+V(I,J)*WK(J)
        B(I)=S
4000  CONTINUE
      END
 
C     -------------------------------------------------------------
 
      SUBROUTINE EQN(J,M,ML,PAR,A,B,X,F,T)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M,M),B(M,M),X(M),PAR(*),F(*)
 
C     THE DIFFERENTIAL EQ.  X1'=X2,  X2'=LAMDA*SINH(LAMDA*X1)
C     WITH LAMDA=PAR(1)
 
      DO 1000 I=1,M
        F(I)=0.0
        DO 1000 K=1,M
          A(K,I)=0.0
1000  B(K,I)=0.0
 
      B(1,1)=1.
      A(1,2)=1.
      B(2,2)=1.
      A(2,1)=PAR(1)**2*COSH(PAR(1)*X(1))
 
      F(1)=X(2)
      F(2)=PAR(1)*SINH(PAR(1)*X(1))
      END
 
C     ----------------------------------------------
 
      SUBROUTINE BCS(M,ML,PAR,BC,G,T1,T2,X1,X2)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PAR(*),BC(M,M),G(M),X1(*),X2(*)
 
C     THE BOUNDARY CONDITIONS  X1(T=T1)=0,  X1(T=T2)=1
 
      BC(1,1)=1.0
      BC(1,2)=0.0
      BC(2,1)=1.0
      BC(2,2)=0.0
 
      G(1)=X1(1)
      G(2)=X2(1)-1.
      END
