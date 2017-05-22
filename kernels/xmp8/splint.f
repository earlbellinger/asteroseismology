C     PROGRAM TO INTEGRATE A TABULATED FUNCTION USING CUBIC SPLINE/ B-SPLINE

      PROGRAM SPL
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(20),F(20),C(3,20)
      DIMENSION AB(20,20),BC(20),XN(20),INT(20),WK(200)

C     EXERCISE 6.4

C     SPECIFIED TABLE OF VALUES
      DATA (X(I),I=1,11)/0.0,0.1D0,0.2D0,0.3D0,0.4D0,0.5D0,0.6D0,
     1  0.7D0,0.8D0,0.9D0,1.0D0/
      DATA (F(I),I=1,11)/0.0,0.3162D0,0.4472D0,0.5477D0,0.6325D0,
     1  0.7071D0,0.7746D0,0.8367D0,0.8944D0,0.9487D0,1.0D0/

51    FORMAT('  IER =',I4,5X,'A =',1PD14.6,5X,'B =',D14.6/
     1 5X,'CUBIC SPLINE:',1PD14.6,5X,'TRAPEZOIDAL RULE:',D14.6) 
52    FORMAT('  IER =',I4,5X,'B-SPLINE OF ORDER K =',I4,5X,
     1  'INTEGRAL =',1PD14.6,5X,'EXACT VALUE =',D14.6)
53    FORMAT('   SPLINE :    IER =',I4)
54    FORMAT('   BSPINT :    IER =',I4,'   K =',I4)
 
C     CALCULATE THE CUBIC SPLINE COEFFICIENTS

      N=11
      CALL SPLINE(X,F,N,C,IER)
      WRITE(6,53) IER

      PRINT *,'TYPE K= ORDER OF B-SPLINE'
      READ *,K

C     CALCULATE THE B-SPLINE COEFFICIENTS

      LA=20
      IFLG=0
      CALL BSPINT(N,X,F,K,AB,LA,BC,XN,NO,IFLG,INT,WK,IER)
      WRITE(6,54) IER,K
 
100   PRINT *,'TYPE  A=LOWER LIMIT,  B=UPPER LIMIT  (QUITS WHEN A.EQ.B)'
      READ *,A,B
      IF(A.EQ.B) STOP
      REX=(B**1.5-A**1.5)/1.5
      CALL SPLINT(A,B,SINT,TINT,N,X,F,C,IER)
      WRITE(6,51) IER,A,B,SINT,TINT

      BINT=BSPQD(NO,XN,K,BC,A,B,WK,IER)
      WRITE(6,52) IER,K,BINT,REX

      GO TO 100
      END
 
C     ------------------------------------------------
 
C	To calculate coefficients for B-spline interpolation
C
C	N : (input) Number of entries in the table
C	X : (input) Array of length N containing the abscissas
C	F : (input) Array of length N containing the function values
C		F(I) is the tabulated function value at X(I).
C	K : (input) Order of B-spline required. K=4 gives cubic B-splines
C	A : (input/output) Real array of length LA*3K containing the
C		triangular decomposition of equation matrix in band form
C		For IFLG=2, this array must be supplied, for other values
C		of IFLG it is calculated by the subroutine
C	LA : (input) The first dimension of A as specified in calling program
C		LA.GE.N
C	C : (output) Coefficients of expansion, which will be calculated
C		provided IFLG.NE.1
C	XF : (input/output) Real array of size NO, containing
C		the knots used for B-spline calculations.
C		The knots must be distinct and in ascending order.
C		For IFLG=2, this array must be supplied, for other values
C		of IFLG it is calculated by the subroutine
C	NO : (input/output) Number of knots for B-splines
C		For IFLG=2, this number must be supplied, for other values
C		of IFLG it is calculated by the subroutine
C	IFLG : (input/output) Integer specifying the type of calculation required
C		IFLG=0 The matrix will be calculated and solved for coefficients
C		IFLG=1 The matrix will be calculated and triangular decomposition
C			is obtained, but coefficients are not calculated
C		IFLG=2 The triangular decomposition of matrix is assumed
C			to be available in A and coefficients C are calculated
C		IFLG=-1 same as 0, except that no pivoting will be used
C	INC : (input/output) Integer array containing information about
C		pivoting during solution of system of linear equations
C		For IFLG=2, this array must be supplied, for other values
C		of IFLG it is calculated by the subroutine
C	WK : Scratch array of length 3*N+K+7
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=204 implies N<K or K<2
C		other values may be set by BSPLIN or GAUBND
C
C	Required routines : BSPLIN, GAUBND

      SUBROUTINE BSPINT(N,X,F,K,A,LA,C,XF,NO,IFLG,INC,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),A(LA,3*K),WK(3*N+K+7),C(N),F(N),XF(NO),INC(N)
 
      IF(N.LT.K.OR.K.LT.2) THEN
        IER=204
        RETURN
      ENDIF
 
      IF(IFLG.LE.1) THEN
C	set up the knots for B-splines by dropping points near the ends
        XF(1)=X(1)
        KL=(K-2)/2
        KU=(K-1)/2
        DO 2000 I=2+KL,N-1-KU
          XF(I-KL)=X(I)
2000    CONTINUE
        XF(N-KL-KU)=X(N)
        NO=N-KL-KU
        NDB=N+2
        NDERIV=0
 
C	Set up the equation matrix for calculating coefficients of expansion
C	The matrix is in band form A_{i,j} is stored in A(I,J-I+K)
        DO 2500 I=1,N
          XB=X(I)
          CALL BSPLIN(XF,NO,K,XB,NDERIV,WK,WK(NDB),WK(NDB+2),LEFT,IER,
     1                WK(2*NDB+2))
          IF(IER.GT.100) RETURN
          DO 2200 J=MAX(1,I-K+1),MIN(N,I+K-1)
            A(I,J-I+K)=WK(J)
2200      CONTINUE
2500    CONTINUE
      ENDIF

C	Solve the system of equations for a band matrix
      NUM=1
      KB=K-1
      DO 3000 I=1,N
        C(I)=F(I)
3000  CONTINUE
      CALL GAUBND(N,KB,NUM,A,C,DET,IDET,INC,LA,IER,IFLG,WK)
      END
 
C   ----------------------------------------------------------
 
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
 
C     ---------------------------------------
 
C	To calculate the integral of a function defined by B-spline expansion
C
C	N : (input) Number of knots for B-spline expansion
C	X : (input) Real array of length N containing the knots.
C		The knots must be distinct and in ascending order.
C       K : (input) Order of B-splines, K=4 for cubic B-splines
C	WT : (input) Real array of length N+K-2 containing the
C		coefficients of B-spline expansion
C	XL : (input) Lower limit of integration
C	XU : (input) Upper limit of integration
C	WK : Scratch array of length 5N+6K+9
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=31 implies that XL is outside the range of table
C		IER=32 implies that XU is outside the range of table
C
C	BSPQD = Integral of \sum WT(I)\phi_I(x) over [XL,XU]
C
C	Required routines : BSPLIN

      FUNCTION BSPQD(N,X,K,WT,XL,XU,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),WT(N+K-2),WK(5*N+6*K+9)
 
      IER=0
      BSPQD=0.0
      K1=K+1
      IF(XL.EQ.XU) RETURN
 
C	Calculate the B-spline basis of order K+1 at XL and XU
      NDERIV=0
      NK=N+K1
      CALL BSPLIN(X,N,K1,XL,NDERIV,WK,WK(NK),WK(2*NK),LF1,IER,WK(3*NK))
      IF(IER.GT.100) RETURN
      CALL BSPLIN(X,N,K1,XU,NDERIV,WK(NK),WK(2*NK),WK(3*NK),LF2,IER,
     1            WK(4*NK))
      IF(IER.GT.100) RETURN
      IF(XL.LT.X(1).OR.XL.GT.X(N)) IER=31
      IF(XU.LT.X(1).OR.XU.GT.X(N)) IER=32

C	Define the extra knots
      K2=2*NK+K
      DO 1200 I=1,N
        WK(K2+I)=X(I)
1200  CONTINUE
      DO 1400 I=1,K
        WK(K2+1-I)=X(1)
        WK(K2+N+I)=X(N)
1400  CONTINUE
 
C	The sum for x=XL
      F1=0.0
      N1=N+K
      DO 2000 I=LF1,LF1+K
        F=0
        DO 1500 J=2,I
          F=F+WT(J-1)*(WK(J+K2)-WK(J-K+K2))
1500    CONTINUE
        F1=F1+F*WK(I)/K
2000  CONTINUE
 
C	The sum for x=XU
      F2=0.0
      DO 2500 I=LF2,LF2+K
        F=0
        DO 2200 J=2,I
          F=F+WT(J-1)*(WK(J+K2)-WK(J-K+K2))
2200    CONTINUE
        F2=F2+F*WK(N1+I)/K
2500  CONTINUE

      BSPQD=F2-F1
      END
 
C     ---------------------------------------
 
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
 
C     ---------------------------------------
 
C	To evaluate the cubic spline interpolant at a specified point
C
C	XB : (input) point at which interpolation is required
C	N : (input) Number of points in the table
C	X : (input) Array of length N, containing the abscissas
C	F : (input) Array of length N, containing the function values at X(I)
C	C : (input) Array of length 3*N containing the spline coefficients
C		which should have been calculated using SPLINE
C	DFB : (output) First derivative of spline at x=XB
C	DDFB : (output) Second derivative of spline at x=XB
C	IER : (output) error parameter, IER=0 if execution is successful
C		IER=24 implies XB is outside the range of table on higher side
C		IER=25 implies XB is outside the range of table on lower side
C		IER=201 implies N<2
C	SPLEVL will be the interpolated value at x=XB
C
C	Required routines : None

      FUNCTION SPLEVL(XB,N,X,F,C,DFB,DDFB,IER)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION X(N),F(N),C(3,N)
      SAVE
      DATA LOW/0/

      SPLEVL=0.0
      IF(N.LE.1) THEN
        IER=201
        RETURN
      ENDIF

C	QASCND is true if table is in ascending order
      QASCND=X(N).GT.X(1)
      IER=0

      IF(LOW.LT.1.OR.LOW.GE.N) THEN
C	If the previous value of LOW is inadmissible, set the range to (1,N)
        LOW=1
        IGH=N
      ELSE
        IGH=LOW+1
      ENDIF

1000  IF((XB.LT.X(LOW).AND.XB.LT.X(IGH)).OR.
     1    (XB.GT.X(LOW).AND.XB.GT.X(IGH))) THEN
C	Extend the range
        IF(XB.GT.X(LOW).EQV.QASCND) THEN
C	Extend the range on higher side
          IF(IGH.GE.N) THEN
            IER=24
            LOW=N-1
          ELSE
            NIGH=MIN(N,IGH+2*(IGH-LOW))
            LOW=IGH
            IGH=NIGH
            GO TO 1000
          ENDIF
        ELSE
C	Extend the range on lower side
          IF(LOW.LE.1) THEN
            IER=25
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

C	Calculate the interpolant and its derivatives
      DX=XB-X(LOW)
      SPLEVL=((C(3,LOW)*DX+C(2,LOW))*DX+C(1,LOW))*DX+F(LOW)
      DFB=(3.*C(3,LOW)*DX+2.*C(2,LOW))*DX+C(1,LOW)
      DDFB=6.*C(3,LOW)*DX+2.*C(2,LOW)
      END
 
C     ------------------------------------------------------
 
C	To calculate coefficients of cubic spline interpolation with
C		not-a-knot boundary conditions
C
C	X : (input) Real array of length N containing x values
C	F : (input) Real array of length N containing values of function at X(I)
C		F(I) is the tabulated function value at X(I).
C	N : (input) Length of table X,F
C	C : (output) Real array of length 3*N containing the spline coefficients
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=201 implies that N<2
C
C	Required routines : None

      SUBROUTINE SPLINE(X,F,N,C,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N),C(3,N)

      IF(N.LE.1) THEN
        IER=201
        RETURN

      ELSE IF(N.EQ.2) THEN
C	Use linear interpolation
        C(1,1)=(F(2)-F(1))/(X(2)-X(1))
        C(2,1)=0.0
        C(3,1)=0.0
        RETURN

      ELSE IF(N.EQ.3) THEN
C	Use quadratic interpolation
        DIV12=(F(2)-F(1))/(X(2)-X(1))
        DIV23=(F(3)-F(2))/(X(3)-X(2))
        C(3,1)=0.0
        C(3,2)=0.0
        C(2,1)=(DIV23-DIV12)/(X(3)-X(1))
        C(2,2)=C(2,1)
        C(1,1)=DIV12+C(2,1)*(X(1)-X(2))
        C(1,2)=DIV23+C(2,1)*(X(2)-X(3))
        RETURN

      ELSE
C	Use cubic splines 

C	Setting up the coefficients of tridiagonal matrix
        C(3,N)=(F(N)-F(N-1))/(X(N)-X(N-1))
        DO 1000 I=N-1,2,-1
          C(3,I)=(F(I)-F(I-1))/(X(I)-X(I-1))
          C(2,I)=2.*(X(I+1)-X(I-1))
C	The right hand sides
1000    C(1,I)=3.*(C(3,I)*(X(I+1)-X(I))+C(3,I+1)*(X(I)-X(I-1)))

C	The not-a-knot boundary conditions
        C1=X(3)-X(1)
        C(2,1)=X(3)-X(2)
        C(1,1)=C(3,2)*C(2,1)*(2.*C1+X(2)-X(1))+C(3,3)*(X(2)-X(1))**2
        C(1,1)=C(1,1)/C1
        CN=X(N)-X(N-2)
        C(2,N)=X(N-1)-X(N-2)
        C(1,N)=C(3,N)*C(2,N)*(2.*CN+X(N)-X(N-1))
        C(1,N)=(C(1,N)+C(3,N-1)*(X(N)-X(N-1))**2)/CN

C	Solving the equation by Gaussian elimination
        G=(X(3)-X(2))/C(2,1)
        C(2,2)=C(2,2)-G*C1
        C(1,2)=C(1,2)-G*C(1,1)
        DO 2000 J=2,N-2
          G=(X(J+2)-X(J+1))/C(2,J)
          C(2,J+1)=C(2,J+1)-G*(X(J)-X(J-1))
2000    C(1,J+1)=C(1,J+1)-G*C(1,J)
        G=CN/C(2,N-1)
        C(2,N)=C(2,N)-G*(X(N-1)-X(N-2))
        C(1,N)=C(1,N)-G*C(1,N-1)

C	The back-substitution
        C(1,N)=C(1,N)/C(2,N)
        DO 3000 I=N-1,2,-1
3000    C(1,I)=(C(1,I)-C(1,I+1)*(X(I)-X(I-1)))/C(2,I)
        C(1,1)=(C(1,1)-C(1,2)*C1)/C(2,1)

C	Calculating the coefficients of cubic spline
        DO 4000 I=1,N-1
          C(2,I)=(3.*C(3,I+1)-2.*C(1,I)-C(1,I+1))/(X(I+1)-X(I))
          C(3,I)=(C(1,I)+C(1,I+1)-2.*C(3,I+1))/(X(I+1)-X(I))**2
4000    CONTINUE
C	Set the coefficients for interval beyond X(N) using continuity
C	of second derivative, although they may not be used.
        C(2,N)=C(2,N-1)+3*(X(N)-X(N-1))*C(3,N-1)
        C(3,N)=0.0
      ENDIF
      END
 
C     ------------------------------------------------------
 
C	To compute integral of a tabulated function using cubic splines
C
C	A : (input) Lower limit of integration
C	B : (input) Upper limit of integration (B > A)
C	SINT : (output) Computed value of integral using cubic splines
C	TINT : (output) Computed value of integral using trapezoidal rule
C	N : (input) Number of tabular points
C	X : (input) Array of length N containing abscissas (in ascending order)
C	F : (input) Array of length N containing the function values
C	C : (input) Array of length 3*N containing the coefficients of splines
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=31 implies lower limit A is outside the table
C		IER=32 implies upper limit B is outside the table
C		IER=301 implies A>B or X(1)>X(N) and no calculations are done
C		Other values may be set by SPLEVL
C
C	Required routines : SPLEVL
C
      SUBROUTINE SPLINT(A,B,SINT,TINT,N,X,F,C,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N),C(3,N)

      IER=0
      SINT=0.0
      TINT=0.0
      IF(A.EQ.B) RETURN
      IER=301
      IF(A.GT.B.OR.X(1).GE.X(N)) RETURN

      A1=A
C	Evaluate the function at x=A for trapezoidal rule
      F1=SPLEVL(A1,N,X,F,C,DFB,DDFB,IER)
      IF(IER.GT.100) RETURN
      IER=0
      IF(A.LT.X(1).OR.A.GT.X(N)) IER=31
      IF(B.LT.X(1).OR.B.GT.X(N)) IER=32

C	Integrating over the N-1 subintervals
      DO 2000 I=1,N-1
        IF(A1.LT.X(I+1).OR.I.EQ.N-1) THEN
          B1=MIN(B,X(I+1))
          IF(I.EQ.N-1) B1=B
          D1=A1-X(I)
          D2=B1-X(I)
C	Add integral of cubic spline over [A1,B1]
          SINT=SINT+F(I)*(D2-D1)+C(1,I)*(D2*D2-D1*D1)/2.+
     1         C(2,I)*(D2**3-D1**3)/3.+C(3,I)*(D2**4-D1**4)/4.
          F2=F(I+1)
          IF(B1.NE.X(I+1)) F2=SPLEVL(B1,N,X,F,C,DFB,DDFB,IER1)
C	Trapezoidal rule approximation to integral
          TINT=TINT+0.5*(B1-A1)*(F1+F2)
          IF(B.LE.X(I+1)) RETURN
          A1=B1
          F1=F2
        ENDIF
2000  CONTINUE
      END
