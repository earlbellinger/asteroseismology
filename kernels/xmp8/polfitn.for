C     PROGRAM TO OBTAIN LEAST SQUARES POLYNOMIAL FIT TO A GIVEN DATA
C     USING ORTHOGONAL POLYNOMIALS IN 3 DIMENSIONS
 
      PROGRAM POLYLS
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WK(99000),X(100,9),FB(20,20,20),NK(9),MK(9),XB(9),
     1          IWK(90),DFF(9),D2FF(3,3)
      DIMENSION Y(100),AX(100,16),FXY(20,20,20),C(20020)
 
C     FUNCTION WITH RANDOM ERROR TO GENERATE DATA
 
      F(X1,Y1,Y2)=((231*X1*Y1-315)*Y2*X1+105)*Y2*Y1-5+RANGAU(IS)*1.D-3
 
C     THE TRUE FUNCTION AND ITS DERIVATIVES
 
      F1(X1,Y1,Y2)=((231*X1*Y1-315)*Y2*X1+105)*Y2*Y1-5
      DF1(X1,Y1,Y2)=(462*X1*Y1-315)*Y2**2*Y1
      DF2(X1,Y1,Y2)=((462*X1*Y1-315)*X1*Y2+105)*Y2
      DF3(X1,Y1,Y2)=((462*X1*Y1-630)*X1*Y2+105)*Y1
      DF11(X1,Y1,Y2)=462*Y2**2*Y1**2
      DF12(X1,Y1,Y2)=(924*X1*Y1-315)*Y2*Y2
      DF22(X1,Y1,Y2)=462*X1**2*Y2**2
      DF13(X1,Y1,Y2)=(924*X1*Y1-630)*Y1*Y2
      DF23(X1,Y1,Y2)=(924*X1*Y1-630)*X1*Y2+105
      DF33(X1,Y1,Y2)=(462*X1*Y1-630)*Y1*X1
 
51    FORMAT('   NO. OF PTS =',3I4,5X,'DEGREE OF POLYNOMIAL =',3I3/
     1       5X,'IER =',I4,5X,'CHI SQUARE =',1PD14.7)
52    FORMAT('  X =',1P3D14.6,3X,'F =',D14.6)
53    FORMAT(3X,'F =',1PD14.6,5X,4HF' =,3D14.6)
54    FORMAT('   COEFFICIENTS :'/(1P5D14.6))
56    FORMAT(3X,'F =',1PD14.6,5X,4HF' =,3D14.6/3X,
     1       4HF''=,6D14.6)
57    FORMAT('   TRUE VALUES :',6X,'F =',1PD14.6/3X,4HF' =,3D14.6/3X,
     1       4HF''=,6D14.6)
 
      IS=2
      ND=3
C	No. of data points is fixed by dimension
      N=20
 
C     GENERATING THE DATA WITH RANDOM ERROR
 
      HH=1.D0/(N-1.)
      DO 1000 I=1,N
        X(I,1)=(I-1)*HH
        X(I,2)=(I-1)*HH
        X(I,3)=(I-1)*HH
1000  CONTINUE
      DO 1200 I1=1,N
        DO 1200 I=1,N
          DO 1200 J=1,N
            FXY(I,J,I1)=F(X(I,1),X(J,2),X(I1,3))
1200  CONTINUE
      NK(1)=N
      NK(2)=N
      NK(3)=N
      LA=100
      IC=100
 
      PRINT *,'TYPE M=DEGREE OF POLYNOMIAL TO BE FITTED'
      READ *,(MK(I),I=1,ND)
 
      CALL POLFITN(ND,NK,X,FXY,AX,LA,C,MK,FB,WK,IWK,CHISQ,IER)
 
      WRITE(6,51) (NK(I),I=1,ND),(MK(I),I=1,ND),IER,CHISQ
      WRITE(6,54) (C(I),I=1,(MK(1)+1)*(MK(2)+1)*(MK(3)+1))
 
2000  PRINT *,'TYPE THE X VALUE WHERE FITTED VALUE IS REQUIRED'
      PRINT *,'                 (QUITS WHEN XB(1)<-100)'
      READ *,(XB(I),I=1,ND)
      IF(XB(1).LT.-100) STOP

C	Evaluate the fitted polynomial and derivatives using all three versions of POLEVN
      CALL POLEVN(ND,MK,AX,LA,C,XB,FF,WK,IWK)
      WRITE(6,52)(XB(I),I=1,ND),FF

      CALL POLEVN1(ND,MK,AX,LA,C,XB,FF,DFF,WK,IWK)
      WRITE(6,53) FF,(DFF(I),I=1,ND)

      CALL POLEVN2(ND,MK,AX,LA,C,XB,FF,DFF,D2FF,WK,IWK)
      WRITE(6,56) FF,(DFF(I),I=1,ND),((D2FF(I,J),I=1,J),J=1,ND)

C	The exact value of function and its derivatives
      WRITE(6,57) F1(XB(1),XB(2),XB(3)),DF1(XB(1),XB(2),XB(3)),
     1           DF2(XB(1),XB(2),XB(3)),DF3(XB(1),XB(2),XB(3)),
     1           DF11(XB(1),XB(2),XB(3)),DF12(XB(1),XB(2),XB(3)),
     1           DF22(XB(1),XB(2),XB(3)),DF13(XB(1),XB(2),XB(3)),
     1           DF23(XB(1),XB(2),XB(3)),DF33(XB(1),XB(2),XB(3))
      GO TO 2000
      END
 
C     --------------------------------------------------
 
C	Evaluating the fitted polynomial at any point using known
C	coefficients of orthogonal polynomials in N dimensions.
C	Should be used to evaluate the polynomial using coefficients calculated
C	by POLFITN. It does not calculate the derivatives, for derivatives
C	use POLEVN1 or POLEVN2
C
C	N : (input) Number of dimensions
C	NK : (input) Integer array of length N containing the degree of
C		polynomial in each direction
C	AX : (input) Real array of length LA*(3*N+3) containing the coefficients
C		alpha and beta for orthogonal polynomials in each direction
C		AX(I,3*J-2) contains alpha and AX(I,3*J-1) contains beta
C		for polynomials in Jth dimension
C	LA : (input) First dimension of array AX in the calling program.
C		It must be same as what was used in call to POLFITN while
C		calculating the coefficients.
C	WT : (input) Real array of length (MK(1)+1)(MK(2)+1)...(MK(N)+1)
C		containing the coefficients of the fit. The dimensions of WT in
C		the calling program must match the size along each dimension,
C		WT(MK(1)+1,MK(2)+1,...,MK(N)+1)
C	X0 : (input) Real array of length N containing the coordinates of
C		the point at which polynomial needs to be evaluated
C	F : (output) Calculated value of the fitted polynomial at X0
C	WK : Real array of length 3*N*LA used as scratch space
C	IWK : Integer array of length N used as scratch space
C	
C	Required routines : POLORT
C
      SUBROUTINE POLEVN(N,NK,AX,LA,WT,X0,F,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),AX(LA,3*N),WT(*),X0(N),WK(LA*3*N),IWK(N)
 
C     Calculate the orthogonal polynomials along each dimension
      DO 1000 I=1,N
        N1=3*LA*(I-1)+1
        N2=3*LA*(I-1)+LA+1
        N3=3*LA*(I-1)+2*LA+1
        NJ=3*(I-1)+1
        XB=X0(I)
        CALL POLORT(NK(I),AX(1,NJ),AX(1,NJ+1),XB,WK(N1),WK(N2),WK(N3))
        IWK(I)=1
1000  CONTINUE
 
C     Calculate the summation over n dimensions
      F=0.0
2000  INDEX=IWK(1)
      TERM=WK(INDEX)
      NDP=NK(1)+1
      DO 2200 I=2,N
        N1=(I-1)*3*LA
        TERM=TERM*WK(N1+IWK(I))
        INDEX=INDEX+(IWK(I)-1)*NDP
        NDP=NDP*(NK(I)+1)
2200  CONTINUE
      F=F+TERM*WT(INDEX)
 
      J=1
C     Choose the next point
2400  IF(IWK(J).GE.NK(J)+1) GO TO 2600
      IWK(J)=IWK(J)+1
      GO TO 2000
 
C     If Jth dimension is exhausted go to next one
2600  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 2400
 
      END
 
C     -------------------------------------------------------
 
C	Evaluating the fitted polynomial and its first derivative at any point
C	using known coefficients of orthogonal polynomials in N dimensions.
C	Should be used to evaluate the polynomial using coefficients calculated
C	by POLFITN. It does not calculate the second derivatives, for that
C	use POLEVN2, for no derivatives use POLEVN
C
C	N : (input) Number of dimensions
C	NK : (input) Integer array of length N containing the degree of
C		polynomial in each direction
C	AX : (input) Real array of length LA*(3*N+3) containing the coefficients
C		alpha and beta for orthogonal polynomials in each direction
C		AX(I,3*J-2) contains alpha and AX(I,3*J-1) contains beta
C		for polynomials in Jth dimension
C	LA : (input) First dimension of array AX in the calling program.
C		It must be same as what was used in call to POLFITN while
C		calculating the coefficients.
C	WT : (input) Real array of length (MK(1)+1)(MK(2)+1)...(MK(N)+1)
C		containing the coefficients of the fit. The dimensions of WT in
C		the calling program must match the size along each dimension,
C		WT(MK(1)+1,MK(2)+1,...,MK(N)+1)
C	X0 : (input) Real array of length N containing the coordinates of
C		the point at which polynomial needs to be evaluated
C	F : (output) Calculated value of the fitted polynomial at X0
C	DF : (output) Real array of length N containing the first derivatives
C		of F at X0
C	WK : Real array of length 3*N*LA+N used as scratch space
C	IWK : Integer array of length N used as scratch space
C	
C	Required routines : POLORT
C
      SUBROUTINE POLEVN1(N,NK,AX,LA,WT,X0,F,DF,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),AX(LA,3*N),WT(*),X0(N),WK(LA*3*N+N),IWK(N),DF(N)
 
C     Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*LA*(I-1)+1
        N2=3*LA*(I-1)+LA+1
        N3=3*LA*(I-1)+2*LA+1
        NJ=3*(I-1)+1
        XB=X0(I)
        CALL POLORT(NK(I),AX(1,NJ),AX(1,NJ+1),XB,WK(N1),WK(N2),WK(N3))
        IWK(I)=1
1000  CONTINUE
 
C     calculate the sum over all points
      F=0.0
      DO 1800 I=1,N
        DF(I)=0.0
1800  CONTINUE
      N0=LA*3*N
 
2000  INDEX=IWK(1)
      TERM=WK(INDEX)
      WK(N0+1)=WK(LA+INDEX)
      DO 2100 I=2,N
        WK(N0+I)=WK(INDEX)
2100  CONTINUE
      NDP=NK(1)+1
      DO 2300 I=2,N
        N1=(I-1)*3*LA
        N2=IWK(I)
        TERM=TERM*WK(N1+N2)
 
C     terms for first derivatives
        DO 2200 J=1,N
          IF(J.EQ.I) THEN
            WK(N0+J)=WK(N0+J)*WK(N1+N2+LA)
          ELSE
            WK(N0+J)=WK(N0+J)*WK(N1+N2)
          ENDIF
2200    CONTINUE
        INDEX=INDEX+(N2-1)*NDP
        NDP=NDP*(NK(I)+1)
2300  CONTINUE
      F=F+TERM*WT(INDEX)
      DO 2400 I=1,N
        DF(I)=DF(I)+WK(N0+I)*WT(INDEX)
2400  CONTINUE
 
      J=1
C     Go to the next point
2500  IF(IWK(J).GE.NK(J)+1) GO TO 2600
      IWK(J)=IWK(J)+1
      GO TO 2000
 
C     If Jth dimension is exhausted go to the next one
2600  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 2500
 
 
      END
 
C     ---------------------------------------------------
 
C	Evaluating the fitted polynomial and its derivatives at any point
C	using known coefficients of orthogonal polynomials in N dimensions.
C	Should be used to evaluate the polynomial using coefficients calculated
C	by POLFITN. If second derivative is not required use POLEVN1, if
C	no derivatives are required then use POLEVN
C
C	N : (input) Number of dimensions
C	NK : (input) Integer array of length N containing the degree of
C		polynomial in each direction
C	AX : (input) Real array of length LA*(3*N+3) containing the coefficients
C		alpha and beta for orthogonal polynomials in each direction
C		AX(I,3*J-2) contains alpha and AX(I,3*J-1) contains beta
C		for polynomials in Jth dimension
C	LA : (input) First dimension of array AX in the calling program.
C		It must be same as what was used in call to POLFITN while
C		calculating the coefficients.
C	WT : (input) Real array of length (MK(1)+1)(MK(2)+1)...(MK(N)+1)
C		containing the coefficients of the fit. The dimensions of WT in
C		the calling program must match the size along each dimension,
C		WT(MK(1)+1,MK(2)+1,...,MK(N)+1)
C	X0 : (input) Real array of length N containing the coordinates of
C		the point at which polynomial needs to be evaluated
C	F : (output) Calculated value of the fitted polynomial at X0
C	DF : (output) Real array of length N containing the first derivatives
C		of F at X0
C	DDF : (output) Real array of length N*N containing the second derivatives
C		of F at X0, DDF(I,J)=d^2F/(dX(I) dX(J))
C	WK : Real array of length 3*N*LA+N+N*N used as scratch space
C	IWK : Integer array of length N used as scratch space
C	
C	Required routines : POLORT
C
      SUBROUTINE POLEVN2(N,NK,AX,LA,WT,X0,F,DF,DDF,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),AX(LA,3*N),WT(*),X0(N),WK(*),IWK(N),DF(N),DDF(N,N)
 
C     Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*LA*(I-1)+1
        N2=3*LA*(I-1)+LA+1
        N3=3*LA*(I-1)+2*LA+1
        NJ=3*(I-1)+1
        XB=X0(I)
        CALL POLORT(NK(I),AX(1,NJ),AX(1,NJ+1),XB,WK(N1),WK(N2),WK(N3))
        IWK(I)=1
1000  CONTINUE
 
C     calculate the sum over all points
 
      F=0.0
      DO 1800 I=1,N
        DF(I)=0.0
        DO 1800 J=1,N
          DDF(J,I)=0.0
1800  CONTINUE
      N0=LA*3*N
2000  INDEX=IWK(1)
      TERM=WK(INDEX)

C	Terms for the derivatives
      WK(N0+1)=WK(LA+INDEX)
      DO 2100 I=2,N+N*N
        WK(N0+I)=WK(INDEX)
2100  CONTINUE
      WK(N0+N+1)=WK(2*LA+INDEX)
      DO 2200 I=2,N
        WK(N0+N+I)=WK(INDEX+LA)
2200  CONTINUE
      NDP=NK(1)+1
      DO 2500 I=2,N
        N1=(I-1)*3*LA
        N2=IWK(I)
        TERM=TERM*WK(N1+N2)
        DO 2300 J=1,N
          IF(J.EQ.I) THEN
            WK(N0+J)=WK(N0+J)*WK(N1+N2+LA)
          ELSE
            WK(N0+J)=WK(N0+J)*WK(N1+N2)
          ENDIF
          DO 2300 J1=1,J
            IF(J.EQ.I.AND.J1.EQ.I) THEN
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2+2*LA)
            ELSE IF(J.EQ.I.OR.J1.EQ.I) THEN
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2+LA)
            ELSE
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2)
            ENDIF
2300    CONTINUE
        INDEX=INDEX+(N2-1)*NDP
        NDP=NDP*(NK(I)+1)
2500  CONTINUE
      F=F+TERM*WT(INDEX)
      DO 2600 I=1,N
        DF(I)=DF(I)+WK(N0+I)*WT(INDEX)
        DO 2600 J=1,I
          DDF(I,J)=DDF(I,J)+WK(N0+N+(J-1)*N+I)*WT(INDEX)
2600  CONTINUE
 
      J=1
C	Go to the next point
3000  IF(IWK(J).GE.NK(J)+1) GO TO 3200
      IWK(J)=IWK(J)+1
      GO TO 2000
 
C	If Jth dimension is exhausted, go to the next one
3200  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 3000
 
      DO 3500 I=1,N
        DO 3500 J=I+1,N
          DDF(I,J)=DDF(J,I)
3500  CONTINUE
 
      END
 
C     ---------------------------------------------------
 
C	Least squares polynomial fit using orthogonal polynomials in 1 dimension
C	Modified version of POLFIT to fit multiple sets of function values
C	This routine is meant to be used for fit in multiple dimensions
C
C	N : (input) Number of data points
C	M : (input) Required degree of polynomial
C	NUM : (input) Number of different RHS (function values) to be fitted
C		Each set must be defined over the same abscissas and
C		with same weights.
C	X : (input) Real array of length N containing the abscissas
C	F : (input) Real array of length N*NUM containing the function values
C		for each RHS, F(I,J) is the value at X(I) for Jth set
C		The first dimension of F is assumed to be exactly equal
C		to N to minimise storage requirement.
C	W : (input) Real array of length N containing the weights associated
C		with each point W(I) is the weight for F(I,J). Normally,
C		W(I) should be 1/err(I)**2, where err(I) is error in F(I,J)
C		The weights are assumed to be the same for all NUM data sets.
C	A : (output) Real array of length (M+1)*NUM containing the coefficients
C		for the fit for each RHS. The first dimension of A
C		is assumed to be M+1 to minimise storage requirements
C		A(I,J) is the Ith coefficient for Jth RHS.
C	ALP, BETA : (output) Real arrays of length M+1, containing the coefficients
C		required for defining the orthogonal polynomials
C	GAM : (output) Real array of length M+1, containing the quantities
C		\gamma_i for the orthogonal polynomials
C	WK : Real array of length 2N used as scratch space
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=601 implies that N<M+1 or M<0 or N < 1
C		IER=621 implies that GAM(I) vanishes at some I and calculations
C			are abandoned
C	
C	The fitted polynomial can be calculated at any value of x using POLEVL
C
C	Required routines : None
 
      SUBROUTINE POLFIT1(N,M,NUM,X,F,W,A,ALP,BETA,GAM,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N,NUM),W(N),A(M+1,NUM),ALP(M+1),BETA(M+1),
     1      WK(N,2),GAM(M+1)
 
      IF(M.GE.N.OR.M.LT.0.OR.N.LT.1) THEN
        IER=601
        RETURN
      ENDIF
 
C	initialisation
      I1=1
      I2=2
      GAM0=0.0
      DO 2000 I=1,N
        WK(I,I1)=0.0
        WK(I,I2)=1.0
        GAM0=GAM0+W(I)
2000  CONTINUE
      GAM(1)=GAM0
      BETA(1)=0
      IER=0
 
C	Loop over the degree
      DO 5000 J=0,M
        IF(GAM(J+1).LE.0) THEN
          IER=621
          RETURN
        ENDIF
 
        DO 3200 J1=1,NUM
          S=0.0
          DO 3000 I=1,N
3000      S=S+W(I)*F(I,J1)*WK(I,I2)
          A(J+1,J1)=S/GAM(J+1)
3200    CONTINUE
        IF(J.EQ.M) RETURN

        S=0
        DO 3400 I=1,N
3400    S=S+W(I)*X(I)*WK(I,I2)**2
        ALP(J+1)=S/GAM(J+1)
        GAM0=0.0
        DO 3600 I=1,N
          WK(I,I1)=(X(I)-ALP(J+1))*WK(I,I2)-BETA(J+1)*WK(I,I1)
          GAM0=GAM0+W(I)*WK(I,I1)**2
3600    CONTINUE
        BETA(J+2)=GAM0/GAM(J+1)
        GAM(J+2)=GAM0
        IT=I1
        I1=I2
        I2=IT
5000  CONTINUE
 
      END
 
C     ---------------------------------------------------
 
C	Least squares polynomial fit using orthogonal polynomials in n dimensions
C	Weights are assumed to be equal for all points and points are
C	assumed to be on a hyper-rectangular mesh
C
C	N : (input) Number of dimensions
C	NK : (input) Integer array of length N containing the number of
C		data points along each direction
C	X : (input) Real array of length LA*N containing the coordinates
C		of tabular points, X(I,J) contains the Ith point along
C		Jth dimension
C	F : (input) Real array of length NK(1)*NK(2)*...*NK(N) containing
C		the function values. The dimension of F in the calling program
C		must match the size along each dimension, F(NK(1),...,NK(N))
C	AX : (output) Real array of length LA*(3*N+3) containing information about
C		fit along each direction. AX(I,3*J-2), AX(I,3*J-1), AX(I,3*J) will
C		respectively contain the coefficients, alpha, beta, gamma
C		for fit along Jth dimension.
C		The rest of the array is used as scratch space
C	LA : (input) First dimension of arrays X and AX as declared
C		in the calling program. LA .GE. MAX(NK(I))
C	C : (output) Real array of length (MK(1)+1)(MK(2)+1)...(MK(N)+1) containing
C		the fitted coefficients of product of orthogonal polynomials 
C		The dimension of F in the calling program must match the size
C		along each dimension, C(MK(1)+1,MK(2)+1,...,MK(N)+1)
C	MK : (input) Integer array of length N containing the required
C		degree of polynomial in each dimension
C	FY : (output) Real array of same size and shape as F containing the
C		fitted function values at each tabular point	
C	WK : Real array of length 2*NK(1)*NK(2)*...*NK(N)  used as scratch space
C	IWK : Integer array of length 2*N used as scratch space
C	CHISQ : (output) the Chi square value for the fit
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=605 implies that LA < MAX(NK(I))
C		In  this case no calculations are done
C		other values may be set by POLFIT1
C	
C	The fitted polynomial can be calculated at any value of x using POLEVN
C
C	Required routines : POLFIT1, POLEVN, POLORT
 
      SUBROUTINE POLFITN(N,NK,X,F,AX,LA,C,MK,FY,WK,IWK,CHISQ,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(LA,N),AX(LA,3*N+3),C(*),F(*),FY(*),WK(*),NK(N),
     1    MK(N),IWK(2*N)
 
      N1=3*N+1
C     Set the weights to one for fits along each dimension
      DO 1000 I=1,LA
        AX(I,N1)=1.0
1000  CONTINUE
 
      NX=1
      NMAX=NK(N)
      DO 2200 I=1,N-1
        NX=NX*NK(I)
        NMAX=MAX(NMAX,NK(I))
2200  CONTINUE
      IF(NMAX.GT.LA) THEN
        IER=605
        RETURN
      ENDIF
 
      N0=NX*NK(N)+1
      NY=1
 
C     Set up the RHS for fit along Nth dimension
      LJ=NK(N)
      DO 2500 I=1,NX
        DO 2500 J=1,NK(N)
          WK(J+(I-1)*LJ)=F(I+(J-1)*NX)
2500  CONTINUE
 
C     Loop for fit along each dimension
      DO 5000 J1=N,1,-1
        NUM=NX*NY
        LJ=NK(J1)
        M=MK(J1)+1
        NI=1+(J1-1)*3
        CALL POLFIT1(LJ,MK(J1),NUM,X(1,J1),WK,AX(1,N1),WK(N0),
     1          AX(1,NI),AX(1,NI+1),AX(1,NI+2),AX(1,N1+1),IER)
        IF(IER.GT.100) RETURN
 
        IF(J1.GT.1) THEN
C     Set up the RHS for fit along next dimension
          NX1=NX/NK(J1-1)
          NY1=NY*M
          LJ1=NK(J1-1)
          DO 3500 I1=1,NY
            DO 3500 I2=1,M
              DO 3500 I=1,NK(J1-1)
                DO 3500 J=1,NX1
                  WK(I+(J-1)*LJ1+(I2-1)*LJ1*NX1+(I1-1)*NX1*LJ1*M)=
     1            WK(N0+I2-1+(J-1)*M+(I-1)*NX1*M+(I1-1)*NX*M)
3500      CONTINUE
          NX=NX1
          NY=NY1
 
        ELSE
C     Store the fitted coefficients in array C
          DO 3800 I=1,M
            DO 3800 J=1,NY
              C(I+(J-1)*M)=WK(N0+I-1+(J-1)*M)
3800      CONTINUE
        ENDIF
5000  CONTINUE
 
C     Calculate the Chi Square
      CHISQ=0.0
      DO 6000 I=1,N
        IWK(I)=1
6000  CONTINUE
      N0=N+1
 
C     Loop over all points
6200  IJ=IWK(1)
      WK(1)=X(IJ,1)
      ID=NK(1)
      DO 6500 I=2,N
        IJ=IJ+ID*(IWK(I)-1)
        ID=ID*NK(I)
        WK(I)=X(IWK(I),I)
6500  CONTINUE
 
      CALL POLEVN(N,MK,AX,LA,C,WK,FY(IJ),WK(N0),IWK(N0))
      CHISQ=CHISQ+(F(IJ)-FY(IJ))**2
 
C     Choose the next point
      J=1
6800  IF(IWK(J).GE.NK(J)) GO TO 7000
      IWK(J)=IWK(J)+1
      GO TO 6200
 
C     If Jth dimension is exhausted go to next one
7000  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 6800
 
      END
 
C     ---------------------------------------------------
 
C	Evaluating the orthogonal polynomial basis functions at any value
C	of x using known coefficients
C	Should be used to evaluate the basis using coefficients calculated
C	by POLFIT or POLFIT1.
C
C	M : (input) Degree of polynomial
C	ALP, BETA : (input) Real arrays of length M+1, containing the coefficients
C		required for defining the orthogonal polynomials
C		ALP, BETA could be calculated using POLFIT
C	X : (input) Value of x at which polynomials needs to be evaluated
C	F : (output) Real array of length M+1 containing the value of
C		orthogonal polynomials at X
C	DF : (output) Real array of length M+1 containing first
C		derivative of F at X
C	DDF : (output) Real array of length M+1 containing second
C		derivative of F at X
C	
C	Required routines : None

 
      SUBROUTINE POLORT(M,ALP,BETA,X,F,DF,DDF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ALP(M+1),BETA(M+1),F(M+1),DF(M+1),DDF(M+1)

      F(1)=1.0
      F(2)=(X-ALP(1))
      DF(1)=0.0
      DF(2)=1.0
      DDF(1)=0.0
      DDF(2)=0.0
 
C	The recurrence relations
      DO 1000 J=3,M+1
        DDF(J)=2.*DF(J-1)+(X-ALP(J-1))*DDF(J-1)-BETA(J-1)*DDF(J-2)
        DF(J)=F(J-1)+(X-ALP(J-1))*DF(J-1)-BETA(J-1)*DF(J-2)
        F(J)=(X-ALP(J-1))*F(J-1)-BETA(J-1)*F(J-2)
1000  CONTINUE
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
