C     PROGRAM TO PERFORM INTERPOLATION IN 4 DIMENSIONS USING B-SPLINES
 
      PROGRAM BSPN
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL RAN
      DIMENSION B(10,10,10,10),X0(90),DF(90)
      DIMENSION AX(20*20,10),INTX(20,10),WK(20000),XF(20,10)
      DIMENSION NK(20),IWK(90),F(10,10,10,10),X(20,10),MK(20),DDF(4,4)
 
C     Sample function for generating table of values and testing
C     interpolated values and derivatives
 
      FUN(XX,Y,Z,Z1)=SIN(XX)*SIN(Y)*SIN(Z)*SIN(Z1)
      FX1(XX,Y,Z,Z1)=COS(XX)*SIN(Y)*SIN(Z)*SIN(Z1)
      FX2(XX,Y,Z,Z1)=SIN(XX)*COS(Y)*SIN(Z)*SIN(Z1)
      FX3(XX,Y,Z,Z1)=SIN(XX)*SIN(Y)*COS(Z)*SIN(Z1)
      FX4(XX,Y,Z,Z1)=SIN(XX)*SIN(Y)*SIN(Z)*COS(Z1)
      FX11(XX,Y,Z,Z1)=-SIN(XX)*SIN(Y)*SIN(Z)*SIN(Z1)
      FX22(XX,Y,Z,Z1)=-SIN(XX)*SIN(Y)*SIN(Z)*SIN(Z1)
      FX33(XX,Y,Z,Z1)=-SIN(XX)*SIN(Y)*SIN(Z)*SIN(Z1)
      FX44(XX,Y,Z,Z1)=-SIN(XX)*SIN(Y)*SIN(Z)*SIN(Z1)
      FX12(XX,Y,Z,Z1)=COS(XX)*COS(Y)*SIN(Z)*SIN(Z1)
      FX13(XX,Y,Z,Z1)=COS(XX)*SIN(Y)*COS(Z)*SIN(Z1)
      FX14(XX,Y,Z,Z1)=COS(XX)*SIN(Y)*SIN(Z)*COS(Z1)
      FX23(XX,Y,Z,Z1)=SIN(XX)*COS(Y)*COS(Z)*SIN(Z1)
      FX24(XX,Y,Z,Z1)=SIN(XX)*COS(Y)*SIN(Z)*COS(Z1)
      FX34(XX,Y,Z,Z1)=SIN(XX)*SIN(Y)*COS(Z)*COS(Z1)
 
51    FORMAT('   IER =',I4,'    ORDER OF B-SPLINES =',I4,
     1       '    NO. OF KNOTS : ', 4I3)
52    FORMAT('  MAX. ERROR IN FUNCTION VALUE =', 1PD12.4)
53    FORMAT('  MAX. ERROR IN FIRST DERIVATIVES =', 1P5D12.4)
54    FORMAT('  MAX. ERROR IN SECOND DERIVATIVES =', 1P10D12.4)
 
C     The number of points along each axis
C     To change the number of points along any axis dimensions too will
C     need to be changed
      NX=10
      N=4
100   PRINT *,'TYPE K= ORDER OF B-SPLINE  (QUITS WHEN K<1)'
      READ *,K
      IF(K.LE.0) STOP
 
C     Initialise the mesh and function values for interpolation
      H=1.D0/(NX-1.)
      DO 1000 I=1,NX
        X(I,1)=(I-1)*H
        X(I,2)=(I-1)*H
        X(I,3)=(I-1)*H
        X(I,4)=(I-1)*H
1000  CONTINUE
      DO 1500 I2=1,NX
        DO 1500 I1=1,NX
          DO 1500 I=1,NX
            DO 1500 J=1,NX
              F(J,I,I1,I2)=FUN(X(J,1),X(I,2),X(I1,3),X(I2,4))
1500  CONTINUE
      DO 1600 I=1,N
        NK(I)=NX
1600  CONTINUE
 
C     la should be first dimension of X, AX, XF, INTX
      LA=20
 
      CALL BSPINTN(N,NK,X,LA,F,K,AX,B,XF,MK,INTX,WK,IER)
      WRITE(6,51)IER,K,(MK(I),I=1,N)
 
C     Test the interpolation by evaluating the function at
C     a set of random  points and comparing with exact values.
C     Testing the function at 100000 points will take some time
C     If needed the number of trials can be reduced.
 
      S1=0.0
      S2=0.0
      S3=0.0
      S4=0.0
      S11=0.0
      S12=0.0
      S22=0.0
      S23=0.0
      S33=0.0
      S13=0.0
      S14=0.0
      S24=0.0
      S34=0.0
      S44=0.0
      SEED=123.0
      DO 3000 J=1,100000
        XB=RAN(SEED)
        X0(1)=XB
        YB=RAN(SEED)
        X0(2)=YB
        ZB=RAN(SEED)
        X0(3)=ZB
        TB=RAN(SEED)
        X0(4)=TB
C	For only function value use BSPEVN
C              FX=BSPEVN(N,MK,XF,LA,K,B,X0,WK,IWK,IER)
C	For only function value and first derivative use BSPEVN1
C              FX=BSPEVN1(N,MK,XF,LA,K,B,X0,DF,WK,IWK,IER)
C	For calculating second derivatives use BSPEVN2
        FX=BSPEVN2(N,MK,XF,LA,K,B,X0,DF,DDF,WK,IWK,IER)
        FX0=FUN(XB,YB,ZB,TB)
        S1=MAX(S1,ABS(FX-FX0))
        S2=MAX(S2,ABS(DF(1)-FX1(XB,YB,ZB,TB)))
        S3=MAX(S3,ABS(DF(2)-FX2(XB,YB,ZB,TB)))
        S4=MAX(S4,ABS(DF(3)-FX3(XB,YB,ZB,TB)))
        S5=MAX(S5,ABS(DF(4)-FX4(XB,YB,ZB,TB)))
        S11=MAX(S11,ABS(DDF(1,1)-FX11(XB,YB,ZB,TB)))
        S12=MAX(S12,ABS(DDF(1,2)-FX12(XB,YB,ZB,TB)))
        S22=MAX(S22,ABS(DDF(2,2)-FX22(XB,YB,ZB,TB)))
        S13=MAX(S13,ABS(DDF(1,3)-FX13(XB,YB,ZB,TB)))
        S33=MAX(S33,ABS(DDF(3,3)-FX33(XB,YB,ZB,TB)))
        S23=MAX(S23,ABS(DDF(2,3)-FX23(XB,YB,ZB,TB)))
        S14=MAX(S14,ABS(DDF(1,4)-FX14(XB,YB,ZB,TB)))
        S24=MAX(S24,ABS(DDF(2,4)-FX24(XB,YB,ZB,TB)))
        S34=MAX(S34,ABS(DDF(3,4)-FX34(XB,YB,ZB,TB)))
        S44=MAX(S44,ABS(DDF(4,4)-FX44(XB,YB,ZB,TB)))
3000  CONTINUE
      WRITE(6,52) S1
      WRITE(6,53) S2,S3,S4,S5
      WRITE(6,54) S11,S12,S13,S14,S22,S23,S24,S33,S34,S44
      GO TO 100
 
 
      END
 
 
C     -------------------------------------------------
 
C	To calculate function value using B-spline expansion in n-dimensions
C
C	N : (input) Number of dimensions
C	NK : (input) Array of length N containing the number of knots
C		along each dimension
C	X : (input) Real array of length NXD*N containing the knots
C		X(I,J) is the Ith knot along Jth dimension
C	NXD : (input) First dimension of array X in calling program.
C		NXD must be greater than or equal to maximum of NK(I).
C	K : (input) Order of B-splines, K=4 for cubic B-splines
C	WT : (input) Coefficients of B-spline expansion, the coefficients
C		are assumed to be stored in natural FORTRAN order with
C		no gaps in data. The calling program should have
C	 	dimension  WT(NK(1)+K-2,NK(2)+K-2,...,NK(N)+K-2)
C	X0 : (input) Real array of length N, containing the coordinates
C		of the point at which function needs to be evaluated
C	WK : Real scratch array of length (NXD+K)*(3N+1)+K+2
C	IWK : Integer scratch array of length 2N
C	IER : (output) Error parameter, IER=0 implies successful execution
C		Nonzero values of IER may be set by BSPLIN which is called
C
C	BSPEVN will give the value of the B-spline expansion at X0.
C	This routine does not calculate the derivatives, for that
C	use BSPEVN1 or BSPEVN2.
C
C	Required routines : BSPLIN
 
      FUNCTION BSPEVN(N,NK,X,NXD,K,WT,X0,WK,IWK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),X(NXD,N),WT(*),X0(N),WK(*),IWK(2*N)
 
      K1=(NXD+K)*3*N+1
      NDERIV=0
      BSPEVN=0.0
C	Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*(NXD+K)*(I-1)+1
        N2=(NXD+K)*(3*I-2)+1
        N3=(NXD+K)*(3*I-1)+1
        CALL BSPLIN(X(1,I),NK(I),K,X0(I),NDERIV,WK(N1),WK(N2),WK(N3),
     1              IWK(I),IER,WK(K1))
        IF(IER.GT.100) RETURN
        IWK(N+I)=0
1000  CONTINUE
 
C	Calculate the summation over n dimensions
      F=0.0
2000  INDEX=IWK(1)+IWK(N+1)
      TERM=WK(INDEX)
      NDP=NK(1)+K-2
      DO 2200 I=2,N
        N1=(I-1)*3*(NXD+K)
        TERM=TERM*WK(N1+IWK(I)+IWK(N+I))
        INDEX=INDEX+(IWK(I)+IWK(I+N)-1)*NDP
        NDP=NDP*(NK(I)+K-2)
2200  CONTINUE
      F=F+TERM*WT(INDEX)
 
      J=1
C	Choose the next point
2400  IF(IWK(J+N).GE.K-1) GO TO 2600
      IWK(J+N)=IWK(J+N)+1
      GO TO 2000
 
C	If Jth dimension is exhausted go to next one
2600  IWK(J+N)=0
      J=J+1
      IF(J.LE.N) GO TO 2400
 
      BSPEVN=F
 
      END
 
C     -------------------------------------------------
 
C	To calculate function value and first derivative using B-spline
C		expansion in n-dimensions
C
C	N : (input) Number of dimensions
C	NK : (input) Array of length N containing the number of knots
C		along each dimension
C	X : (input) Real array of length NXD*N containing the knots
C		X(I,J) is the Ith knot along Jth dimension
C	NXD : (input) First dimension of array X in calling program.
C		NXD must be greater than or equal to maximum of NK(I).
C	K : (input) Order of B-splines, K=4 for cubic B-splines
C	WT : (input) Coefficients of B-spline expansion, the coefficients
C		are assumed to be stored in natural FORTRAN order with
C		no gaps in data. The calling program should have
C	 	dimension  WT(NK(1)+K-2,NK(2)+K-2,...,NK(N)+K-2)
C	X0 : (input) Real array of length N, containing the coordinates
C		of the point at which function needs to be evaluated
C	DF : (output) Real array of length N which will contain the
C		calculated values of first derivatives.
C		DF(I) will be derivative w.r.t. X_I
C	WK : Real scratch array of length (NXD+K)*(3N+1)+K+N+2
C	IWK : Integer scratch array of length 2N
C	IER : (output) Error parameter, IER=0 implies successful execution
C		Nonzero values of IER may be set by BSPLIN which is called
C
C	BSPEVN1 will give the value of the B-spline expansion at X0.
C	This routine does not calculate the second derivative, for that
C	use BSPEVN2. To improve efficiency if derivatives are not
C	required then use BSPEVN.
C
C	Required routines : BSPLIN
 
 
      FUNCTION BSPEVN1(N,NK,X,NXD,K,WT,X0,DF,WK,IWK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),X(NXD,N),WT(*),X0(N),DF(N),WK(*),IWK(2*N)
 
      K1=(NXD+K)*3*N+1
      NDERIV=1
      BSPEVN1=0.0
C	Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*(NXD+K)*(I-1)+1
        N2=(NXD+K)*(3*I-2)+1
        N3=(NXD+K)*(3*I-1)+1
        CALL BSPLIN(X(1,I),NK(I),K,X0(I),NDERIV,WK(N1),WK(N2),WK(N3),
     1              IWK(I),IER,WK(K1))
        IF(IER.GT.100) RETURN
        IWK(N+I)=0
1000  CONTINUE
 
C	calculate the sum over all points
      F=0.0
      DO 1800 I=1,N
        DF(I)=0.0
1800  CONTINUE
      N0=(NXD+K)*3*N

2000  INDEX=IWK(1)+IWK(N+1)
      TERM=WK(INDEX)
      WK(N0+1)=WK(NXD+K+INDEX)
      DO 2100 I=2,N
        WK(N0+I)=WK(INDEX)
2100  CONTINUE
      NDP=NK(1)+K-2
      DO 2300 I=2,N
        N1=(I-1)*3*(NXD+K)
        N2=IWK(I)+IWK(N+I)
        TERM=TERM*WK(N1+N2)

C	terms for first derivatives
        DO 2200 J=1,N
          IF(J.EQ.I) THEN
            WK(N0+J)=WK(N0+J)*WK(N1+N2+NXD+K)
          ELSE
            WK(N0+J)=WK(N0+J)*WK(N1+N2)
          ENDIF
2200    CONTINUE
        INDEX=INDEX+(N2-1)*NDP
        NDP=NDP*(NK(I)+K-2)
2300  CONTINUE
      F=F+TERM*WT(INDEX)
      DO 2400 I=1,N
        DF(I)=DF(I)+WK(N0+I)*WT(INDEX)
2400  CONTINUE
 
      J=1
C	Go to the next point
2500  IF(IWK(J+N).GE.K-1) GO TO 2600
      IWK(J+N)=IWK(J+N)+1
      GO TO 2000
 
C	If Jth dimension is exhausted go to the next one
2600  IWK(J+N)=0
      J=J+1
      IF(J.LE.N) GO TO 2500
 
      BSPEVN1=F
 
      END
 
C     -------------------------------------------------------
 
C	To calculate function value and first and second derivatives using
C		B-spline expansion in n-dimensions
C
C	N : (input) Number of dimensions
C	NK : (input) Array of length N containing the number of knots
C		along each dimension
C	X : (input) Real array of length NXD*N containing the knots
C		X(I,J) is the Ith knot along Jth dimension
C	NXD : (input) First dimension of array X in calling program.
C		NXD must be greater than or equal to maximum of NK(I).
C	K : (input) Order of B-splines, K=4 for cubic B-splines
C	WT : (input) Coefficients of B-spline expansion, the coefficients
C		are assumed to be stored in natural FORTRAN order with
C		no gaps in data. The calling program should have
C	 	dimension  WT(NK(1)+K-2,NK(2)+K-2,...,NK(N)+K-2)
C	X0 : (input) Real array of length N, containing the coordinates
C		of the point at which function needs to be evaluated
C	DF : (output) Real array of length N which will contain the
C		calculated values of first derivatives.
C		DF(I) will be derivative w.r.t. X_I
C	DDF : (output) Real array of length N*N which will contain the
C		calculated values of second derivatives. It must have
C		a dimension of (N,N) in the calling program.
C		DDF(I,J) will be derivative w.r.t. X_I and X_J
C	WK : Real scratch array of length (NXD+K)*(3N+1)+K+N+N*N+2
C	IWK : Integer scratch array of length 2N
C	IER : (output) Error parameter, IER=0 implies successful execution
C		Nonzero values of IER may be set by BSPLIN which is called
C
C	BSPEVN2 will give the value of the B-spline expansion at X0.
C	To improve efficiency if derivatives are not required then use BSPEVN,
C	while for only first derivatives use BSPEVN1
C
C	Required routines : BSPLIN
 
      FUNCTION BSPEVN2(N,NK,X,NXD,K,WT,X0,DF,DDF,WK,IWK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),X(NXD,N),WT(*),X0(N),DF(N),DDF(N,N),WK(*),IWK(2*N)
 
      K1=(NXD+K)*3*N+1
      NDERIV=2
      BSPEVN2=0.0
C	Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*(NXD+K)*(I-1)+1
        N2=(NXD+K)*(3*I-2)+1
        N3=(NXD+K)*(3*I-1)+1
        CALL BSPLIN(X(1,I),NK(I),K,X0(I),NDERIV,WK(N1),WK(N2),WK(N3),
     1              IWK(I),IER,WK(K1))
        IF(IER.GT.100) RETURN
        IWK(N+I)=0
1000  CONTINUE
 
C	Calculate the sum over all points
      F=0.0
      DO 1800 I=1,N
        DF(I)=0.0
        DO 1800 J=1,N
          DDF(J,I)=0.0
1800  CONTINUE
      N0=(NXD+K)*3*N
2000  INDEX=IWK(1)+IWK(N+1)
      TERM=WK(INDEX)

C	Terms for the first and second derivatives
      WK(N0+1)=WK(NXD+K+INDEX)
      DO 2100 I=2,N+N*N
        WK(N0+I)=WK(INDEX)
2100  CONTINUE
      WK(N0+N+1)=WK(2*NXD+2*K+INDEX)
      DO 2200 I=2,N
        WK(N0+N+I)=WK(INDEX+NXD+K)
2200  CONTINUE
      NDP=NK(1)+K-2
      DO 2500 I=2,N
        N1=(I-1)*3*(NXD+K)
        N2=IWK(I)+IWK(N+I)
        TERM=TERM*WK(N1+N2)
        DO 2300 J=1,N
          IF(J.EQ.I) THEN
            WK(N0+J)=WK(N0+J)*WK(N1+N2+NXD+K)
          ELSE
            WK(N0+J)=WK(N0+J)*WK(N1+N2)
          ENDIF
          DO 2300 J1=1,J
            IF(J.EQ.I.AND.J1.EQ.I) THEN
             WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2+2*NXD+2*K)
            ELSE IF(J.EQ.I.OR.J1.EQ.I) THEN
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2+NXD+K)
            ELSE
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2)
            ENDIF
2300    CONTINUE
        INDEX=INDEX+(N2-1)*NDP
        NDP=NDP*(NK(I)+K-2)
2500  CONTINUE
      F=F+TERM*WT(INDEX)
      DO 2600 I=1,N
        DF(I)=DF(I)+WK(N0+I)*WT(INDEX)
        DO 2600 J=1,I
          DDF(I,J)=DDF(I,J)+WK(N0+N+(J-1)*N+I)*WT(INDEX)
2600  CONTINUE
 
      J=1
C	Go to the next point
3000  IF(IWK(J+N).GE.K-1) GO TO 3200
      IWK(J+N)=IWK(J+N)+1
      GO TO 2000
 
C	If Jth dimension is exhausted go to the next one
3200  IWK(J+N)=0
      J=J+1
      IF(J.LE.N) GO TO 3000
 
      BSPEVN2=F
      DO 3500 I=1,N
        DO 3500 J=I+1,N
          DDF(I,J)=DDF(J,I)
3500  CONTINUE
 
      END
 
C     -------------------------------------------------------
 
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
 
C     -------------------------------------------------------
 
C	To calculate coefficients for B-spline interpolation in n-dimensions
C
C	N : (input) Number of dimensions
C	NK : (input) Integer array of length N giving the number of
C		tabular points in each dimension
C	X : (input) Array of length NXD*N containing the abscissas
C		X(I,J) is the Ith abscissa in Jth dimension
C	NXD : (input) First dimension of arrays X, XF, INTX as specified
C		in the calling program. NXD must be greater than or equal
C		to the maximum of NK(I)
C	F : (input) Array of length NK(1)*NK(2)* ... *NK(N) containing
C		the table of values. The dimensions of F in the calling
C		program must exactly match the size of table so that
C		there are no gaps in memory allocation.
C	K : (input) Order of B-spline required. K=4 gives cubic B-splines
C	AX : (output) Real array of length NXD*3*K*N containing the
C		triangular decomposition of equation matrix in band form
C		for each dimension.
C	C : (output) Real array of length NK(1)*NK(2)*... *NK(N)
C		containing the coefficients of expansion.
C	XF : (output) Real array of size NXD*N, containing
C		the knots used for B-spline calculations.
C		XF(I,J) is the Ith knot along Jth dimension.
C	MK : (output) Integer array of length N containing number of
C		knots for B-splines in each dimension.
C	INTX : (output) Integer array of length NXD*N containing information
C		about pivoting during solution of system of linear equations
C		for each dimension.
C	WK : Scratch array of length NK(1)*NK(2)*...*NK(N)+3*K
C	IER : (output) Error parameter, IER=0 implies successful execution
C		Nonzero values may be set by BSPINT, BSPLIN or GAUBND
C
C	Required routines : BSPINT, BSPLIN, GAUBND

 
      SUBROUTINE BSPINTN(N,NK,X,NXD,F,K,AX,C,XF,MK,INTX,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NXD,N),NK(N),F(*),AX(NXD*3*K,N),XF(NXD,N),C(*),
     1        MK(N),WK(*),INTX(NXD,N)
 
C	Calculate the triangular decomposition of matrices for interpolation
C	along each dimension
      DO 2000 I=1,N
        IFLG1=1
        LJ=NK(I)
        CALL BSPINT(NK(I),X(1,I),F,K,AX(1,I),LJ,C,XF(1,I),MK(I),
     1         IFLG1,INTX(1,I),WK,IER)
        IF(IER.GT.100) RETURN
2000  CONTINUE
 
      KB=K-1
      NX=1
      DO 2200 I=1,N-1
        NX=NX*NK(I)
2200  CONTINUE
      N0=NX*NK(N)+1
      NY=1
      JU=N
 
C	If N is odd interpolate along last dimension outside the loop
      IF(MOD(N,2).EQ.1) THEN
        LJ=NK(N)
C	Set up the RHS
        DO 2500 I=1,NX
          DO 2500 J=1,NK(N)
            WK(J+(I-1)*LJ)=F(I+(J-1)*NX)
2500    CONTINUE
        NUM=NX
        IFLG1=2
        CALL GAUBND(NK(N),KB,NUM,AX(1,N),WK,DET,IDET,INTX(1,N),LJ,
     1         IER,IFLG1,WK(N0))
        IF(IER.GT.100) RETURN
 
C	Set up the RHS for next interpolation along N-1 th dimension
        NX1=NX/NK(N-1)
        NY=NK(N)
        LJ1=NK(N-1)
        DO 3000 I1=1,NK(N-1)
          DO 3000 I=1,NX1
            DO 3000 J=1,NY
              C(I1+(I-1)*LJ1+(J-1)*NX1*LJ1)=WK(J+(I-1)*LJ+(I1-1)*NX1*LJ)
3000    CONTINUE
        NX=NX1
        JU=N-1
      ELSE
C	Set up the RHS for interpolation along N th dimension
        LJ=NK(N)
        DO 3200 I=1,NX
          DO 3200 J=1,NK(N)
            C(J+(I-1)*LJ)=F(I+(J-1)*NX)
3200    CONTINUE
 
      ENDIF
 
C	Loop for interpolation in each dimension, each pass
C	interpolates along 2 dimensions
      DO 5000 J1=JU,1,-2
        NUM=NX*NY
        IFLG1=2
        LJ=NK(J1)
        CALL GAUBND(NK(J1),KB,NUM,AX(1,J1),C,DET,IDET,INTX(1,J1),LJ,
     1         IER,IFLG1,WK(N0))
        IF(IER.GT.100) RETURN
 
C	Set up the RHS for interpolation along the next dimension
        NX1=NX/NK(J1-1)
        NY1=NY*NK(J1)
        LJ1=NK(J1-1)
        DO 3500 I1=1,NY
          DO 3500 I2=1,NK(J1)
            DO 3500 I=1,NK(J1-1)
              DO 3500 J=1,NX1
                WK(I+(J-1)*LJ1+(I2-1)*NX+(I1-1)*NX*NK(J1))=
     1          C(I2+(J-1)*LJ+(I-1)*NX1*LJ+(I1-1)*NX*LJ)
3500    CONTINUE
        NX=NX1
        NY=NY1
        NUM=NY*NX
        IFLG1=2
        LJ=NK(J1-1)
        CALL GAUBND(NK(J1-1),KB,NUM,AX(1,J1-1),WK,DET,IDET,
     1         INTX(1,J1-1),LJ,IER,IFLG1,WK(N0))
        IF(IER.GT.100) RETURN
        IF(J1.EQ.2) THEN
C	Store the coefficients in array C
          DO 3800 I=1,NK(1)
            DO 3800 J=1,NY
              C(I+(J-1)*NK(1))=WK(I+(J-1)*LJ)
3800      CONTINUE
        ELSE
 
C	Set up the RHS for interpolation along the next dimension
          LJ1=NK(J1-2)
          NX1=NX/NK(J1-2)
          NY1=NY*NK(J1-1)
          DO 4000 I1=1,NK(J1-1)
            DO 4000 I2=1,NK(J1-2)
              DO 4000 I=1,NX1
                DO 4000 J=1,NY
                  C(I2+(I-1)*LJ1+(I1-1)*NX+(J-1)*NX*NK(J1-1))=
     1            WK(I1+(I-1)*LJ+(I2-1)*LJ*NX1+(J-1)*LJ*NX)
4000      CONTINUE
          NX=NX1
          NY=NY1
        ENDIF
5000  CONTINUE
      END
 
C     -------------------------------------------------------------
 
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
 
C     ---------------------------------------------------------
 
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

C	---------------------------------------------------------

C	To generate uniformly distributed random numbers in interval (0,1)
C
C	SEED : (input/output) is a real value used as the seed
C		It should be positive during initial call and
C		should not be modified between different calls.
C
C	Required routines : None

      FUNCTION RAN(SEED)
      IMPLICIT REAL*8(A-H,O-Z)
C	Retain the following declaration even for REAL*4 version
C	otherwise AM and AC will be rounded
      REAL*8 AM,A,AC
      PARAMETER(AM=2147483647D0,A=16807D0,AC=2147483648D0)

      SEED=MOD(SEED*A,AM)
      RAN=SEED/AC
      END
