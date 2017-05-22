C     PROGRAM TO OBTAIN LEAST SQUARES POLYNOMIAL FIT TO A GIVEN DATA
C     USING ORTHOGONAL POLYNOMIALS IN TWO DIMENSIONS
 
      PROGRAM POLYLS
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WK(90000),X(100),A(20),FB(100,100),IWK(90),AW(500)
      DIMENSION Y(100),AX(50,6),AY(50,3),FXY(100,100),C(50,100)
 
 
C     FUNCTION WITH RANDOM ERROR TO GENERATE DATA
 
      F(X1,Y1)=((231*X1*Y1-315)*X1*X1+105)*X1*Y1-5+RANGAU(IS)*5.D-3
 
C     THE TRUE FUNCTION AND ITS DERIVATIVES
 
      F1(X1,Y1)=((231*X1*Y1-315)*X1*X1+105)*X1*Y1-5
      DF1(X1,Y1)=((924*X1*Y1-945)*X1*X1+105)*Y1
      DF2(X1,Y1)=((462*X1*Y1-315)*X1*X1+105)*X1
      DF11(X1,Y1)=(2772*X1*Y1-1890)*X1*Y1
      DF12(X1,Y1)=(1848*X1*Y1-945)*X1*X1+105
      DF22(X1,Y1)=462*X1**4
 
51    FORMAT('   NO. OF PTS =',2I4,5X,'DEGREE OF POLYNOMIAL =',2I3,5X,
     1       'IER =',I4/5X,'CHI SQUARE=',1PD14.7)
52    FORMAT('   COEFFICIENTS :'/(1P5D14.6))
56    FORMAT('  X =',1P2D14.6,5X,'F =',D14.6/5X,4HF' =,2D14.6/
     1       5X,4HF''=,3D14.6)
57    FORMAT('   TRUE VALUES :',6X,'F =',1PD14.6/5X,4HF' =,2D14.6/5X,
     1        4HF''=,3D14.6)
 
      IS=2
100   PRINT *, 'TYPE N=NO OF DATA POINTS ALONG EACH AXIS'
      PRINT *,'                                     (QUITS WHEN N.LE.1)'
      READ *,N
      IF(N.LE.1) STOP
 
C     GENERATING THE DATA WITH KNOWN FUNCTION PLUS RANDOM ERROR
 
      HH=1.D0/(N-1.)
      DO 1000 I=1,N
        X(I)=(I-1)*HH
        Y(I)=(I-1)*HH
1000  CONTINUE
      DO 1200 I=1,N
        DO 1200 J=1,N
          FXY(I,J)=F(X(I),Y(J))
1200  CONTINUE
      NX=N
      NY=N
      LA=100
      IC=50
      ND=2
 
      PRINT *,'TYPE MX, MY =DEGREE OF POLYNOMIAL TO BE FITTED ALONG X,Y'
      READ *,MX,MY
 
      CALL POLFIT2(NX,NY,X,Y,FXY,AX,AY,LA,C,IC,MX,MY,FB,WK,AW,CHISQ,IER)
 
      WRITE(6,51) NX,NY,MX,MY,IER,CHISQ
      WRITE(6,52) ((C(I,J),I=1,MX+1),J=1,MY+1)
 
2000  PRINT *,'TYPE THE X0,Y0 VALUES WHERE FITTED VALUE IS REQUIRED'
      PRINT *,'                 (QUITS WHEN X0<-100)'
      READ *,X0,Y0
      IF(X0.LT.-100) STOP
      CALL POLEV2(MX,MY,AX,AY,IC,C,X0,Y0,FF,DFX,DFY,
     1            DFXX,DFXY,DFYY,WK,IER)
      WRITE(6,56) X0,Y0,FF,DFX,DFY,DFXX,DFXY,DFYY
      WRITE(6,57) F1(X0,Y0),DF1(X0,Y0),DF2(X0,Y0),DF11(X0,Y0),
     1           DF12(X0,Y0),DF22(X0,Y0)
      GO TO 2000
      END
 
C     --------------------------------------------------
 
C	Evaluating the fitted polynomial and its derivatives at any value
C	of x using known coefficients of orthogonal polynomials in 2 dimensions
C	Should be used to evaluate the polynomial using coefficients calculated
C	by POLFIT2.
C
C	NX : (input) Degree of polynomial in X
C	NY : (input) Degree of polynomial in Y
C	AX : (input) Real array of length LA*2 containing the coefficients
C		alpha and beta for orthogonal polynomials in X
C		AX(I,1) contains alpha and AX(I,2) contains beta
C	AY : (input) Real array of length LA*2 containing the coefficients
C		alpha and beta for orthogonal polynomials in Y
C		AY(I,1) contains alpha and AY(I,2) contains beta
C		The arrays AX and AY can be calculated using POLFIT2
C	LA : (input) First dimension of arrays AX, AY and WT in the calling
C		program. LA > MAX(NX,NY)
C	WT : (input) Real array of length LA*(MY+1) containing the coefficients
C		of the fit. WT(I,J) is the coefficient of
C		PHI_I(X)PSI_J(Y), where PHI_I and PSI_J are orthogonal
C		polynomials in X and Y
C	X0,Y0 : (input) Coordinates of the point at which polynomial
C		needs to be evaluated
C	F : (output) Calculated value of the fitted polynomial at (X0,Y0)
C	DFX : (output) First derivative  dF/dX at X0,Y0
C	DFY : (output) First derivative dF/dY at X0,Y0
C	DFXX : (output) Second derivative d^2F/dXdX at X0,Y0
C	DFXY : (output) Second derivative d^2F/dXdY at X0,Y0
C	DFYY : (output) Second derivative d^2F/dYdY at X0,Y0
C	WK : Real array of length 6(MAX(NX,NY)+2) used as scratch space
C	IER : (output) error parameter, IER=0 implies successful execution
C		IER=604 implies that LA.LE.MAX(NX,NY), in which case no
C			calculations are done.
C	
C	Required routines : POLORT
C
      SUBROUTINE POLEV2(NX,NY,AX,AY,LA,WT,X0,Y0,F,DFX,DFY,
     1      DFXX,DFXY,DFYY,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AX(LA,2),AY(LA,2),WT(LA,NY+1),WK(*)
 
      NK=MAX(NX,NY)+2
      IF(NK-2.GE.LA) THEN
        IER=604
        RETURN
      ENDIF
C	Calculate the orthogonal polynomials along each dimension
      CALL POLORT(NX,AX(1,1),AX(1,2),X0,WK,WK(NK),WK(2*NK))
      CALL POLORT(NY,AY(1,1),AY(1,2),Y0,WK(3*NK),WK(4*NK),WK(5*NK))
 
C	Calculate the fitted polynomial and its derivatives
      F=0.0
      DFX=0.0
      DFY=0.0
      DFXX=0.0
      DFYY=0.0
      DFXY=0.0
      N1=NK-1
      N2=2*NK-1
      N3=3*NK-1
      N4=4*NK-1
      N5=5*NK-1
      DO 2000 I=1,NX+1
        DO 2000 J=1,NY+1
          F=F+WT(I,J)*WK(I)*WK(N3+J)
          DFX=DFX+WT(I,J)*WK(N1+I)*WK(N3+J)
          DFY=DFY+WT(I,J)*WK(I)*WK(N4+J)
          DFXX=DFXX+WT(I,J)*WK(N2+I)*WK(N3+J)
          DFXY=DFXY+WT(I,J)*WK(N1+I)*WK(N4+J)
          DFYY=DFYY+WT(I,J)*WK(I)*WK(N5+J)
2000  CONTINUE
      END
 
C     -------------------------------------------------------
 
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

C	Least squares polynomial fit using orthogonal polynomials in 2 dimensions
C	Weights are assumed to be equal for all points
C
C	NX : (input) Number of data points along X direction
C	NY : (input) Number of data points along Y direction
C	X,Y : (input) Real arrays of length NX,NY containing the coordinates
C		of tabular points
C	F : (input) Real array of length LA*NY containing the function values
C		F(I,J) is the value at X(I),Y(J)
C	AX : (output) Real array of length IC*3 containing information about
C		fit along X direction. AX(I,1), AX(I,2), AX(I,3) will
C		respectively contain the coefficients, alpha, beta, gamma
C	AY : (output) Real array of length IC*3 containing information about
C		fit along Y direction. AY(I,1), AY(I,2), AY(I,3) will
C		respectively contain the coefficients, alpha, beta, gamma
C	LA : (input) First dimension of arrays  F, FY as declared
C		in the calling program. LA .GE. MAX(NX,NY)
C	C : (output) Real array of length IC*(MY+1) containing the fitted
C		coefficients of product of orthogonal polynomials in X & Y
C	IC : (input) First dimension of arrays C, AX, AY as declared in the calling
C		program. IC > MAX(MX,MY)+1
C	MX : (input) Required degree of polynomial in X
C	MY : (input) Required degree of polynomial in Y
C	FY : (output) Real array of length LA*NY containing the values of
C		fitted function at each tabular point	
C	WK : Real array of length MAX[ NX*(NY+MY+1), 6(MAX(MX,MY)+2)]
C		 used as scratch space
C	AW : Real array of length LA*3 used as scratch space
C	CHISQ : (output) the Chi square value for the fit
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=602 implies that IC<MX+1 or IC<MY+1 or LA<NX or LA<NY
C		IER=603 implies that NX<MX+1 or MX<0 or NY < MY+1 or MY<0
C		In both these cases no calculations are done
C		Other values may be set by POLFIT1.
C	
C	The fitted polynomial can be calculated at any value of x using POLEV2
C
C	Required routines : POLFIT1, POLEV2, POLORT
 
      SUBROUTINE POLFIT2(NX,NY,X,Y,F,AX,AY,LA,C,IC,MX,MY,FY,WK,AW,
     1         CHISQ,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NX),Y(NY),AX(IC,3),AY(IC,3),C(IC,MY+1),F(LA,NY),
     1         FY(LA,NY),WK(*),AW(LA,3)
 
      IF(MAX(NX,NY).GT.LA.OR.IC.LT.MX+1.OR.IC.LT.MY+1) THEN
        IER=602
        RETURN
      ENDIF
 
      IF(MX.LT.0.OR.MY.LT.0.OR.NX.LE.MX.OR.NY.LE.MY) THEN
        IER=603
        RETURN
      ENDIF
 
C     Set the weights to 1
      DO 1000 I=1,MAX(NX,NY)
        AW(I,1)=1.0
1000  CONTINUE
 
C     Set up the RHS for calculating the fits along y-axis
      LJ=NY
      DO 2500 I=1,NX
        DO 2500 J=1,NY
          WK(J+(I-1)*LJ)=F(I,J)
2500  CONTINUE
 
      LN1=LJ*NX+1
      M=MY+1
      NUM=NX
      CALL POLFIT1(NY,MY,NUM,Y,WK,AW(1,1),WK(LN1),AY,AY(1,2),
     1      AY(1,3),AW(1,2),IER)
      IF(IER.GT.100) RETURN
 
C     Set up the RHS for calculating the fits along x-axis
      DO 3000 J=1,M
        DO 3000 I=1,NX
          WK(I+(J-1)*NX)=WK(LN1+J-1+(I-1)*M)
3000  CONTINUE
      NUM=M
      LN1=M*NX+1
      CALL POLFIT1(NX,MX,NUM,X,WK,AW(1,1),WK(LN1),AX,AX(1,2),
     1      AX(1,3),AW(1,2),IER)
      IF(IER.GT.100) RETURN
 
C	Store the calculated coefficients in array C
      M=MX+1
      DO 3500 I=1,MY+1
        DO 3500 J=1,MX+1
          C(J,I)=WK(LN1+J-1+(I-1)*M)
3500  CONTINUE
 
C     Calculate the CHI square
      CHISQ=0.0
      DO 4000 I=1,NY
        DO 4000 J=1,NX
          CALL POLEV2(MX,MY,AX,AY,IC,C,X(J),Y(I),FY(J,I),DFX,DFY,
     1                DFXX,DFXY,DFYY,WK,IER1)
          CHISQ=CHISQ+(F(J,I)-FY(J,I))**2
4000  CONTINUE
 
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
