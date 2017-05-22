C     PROGRAM TO OBTAIN LEAST SQUARES POLYNOMIAL FIT TO A GIVEN DATA
C     USING ORTHOGONAL POLYNOMIALS

      PROGRAM POLYLS
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WK(200),X(100),A(20),FX(100),W(100),BETA(20),ALP(20)
     1          ,FA(100),H(20),GAM(20)

C     EXAMPLE 9.3

C     FUNCTION WITH RANDOM ERROR TO GENERATE DATA

      F(Y)=((231*Y*Y-315)*Y*Y+105)*Y*Y-5+RANGAU(IS)*5.D-3

C     THE TRUE FUNCTION AND ITS DERIVATIVES

      F1(Y)=((231*Y*Y-315)*Y*Y+105)*Y*Y-5
      DF(Y)=((1386*Y*Y-1260)*Y*Y+210)*Y
      DDF(Y)=(6930*Y*Y-3780)*Y*Y+210

51    FORMAT('   NO. OF PTS =',I4,5X,'DEGREE =',I3,5X,'IER =',I4/
     1       '  DEGREE',4X,'COEFFICIENT',3X,'STND. DEV.',6X,'ALPHA',9X,
     2       'BETA',10X,'CHI-SQ.')
52    FORMAT(I6,2X,1P5D14.6)
56    FORMAT('  X =',1PD14.6,3X,'F =',D14.6,3X,4HF' =,D14.6,3X,
     1       4HF''=,D14.6)
57    FORMAT('   TRUE VALUES :',6X,'F =',1PD14.6,3X,4HF' =,D14.6,3X,
     1       4HF''=,D14.6)

      IS=2
100   PRINT *, 'TYPE N=NO OF DATA POINTS    (QUITS WHEN N.LE.1)'
      READ *,N
      IF(N.LE.1) STOP

C     GENERATING THE DATA WITH RANDOM ERROR

      HH=1.D0/(N-1.)
      DO 1000 I=1,N
        X(I)=(I-1)*HH
        FX(I)=F(X(I))
        W(I)=4.0E4
1000  CONTINUE

      PRINT *,'TYPE M=DEGREE OF POLYNOMIAL TO BE FITTED'
      READ *,M

      CALL POLFIT(N,M,X,FX,W,A,ALP,BETA,FA,H,GAM,WK,IER)

      WRITE(6,51) N,M,IER
      DO 1500 I=1,M+1
1500  WRITE(6,52) I-1,A(I),1./SQRT(GAM(I)),ALP(I),BETA(I),H(I)

2000  PRINT *,'TYPE THE X VALUE WHERE FITTED VALUE IS REQUIRED'
      PRINT *,'                 (QUITS WHEN XX<-100)'
      READ *,XX
      IF(XX.LT.-100) STOP
      CALL POLEVL(M,A,ALP,BETA,XX,FF,DA,DDA)
      WRITE(6,56) XX,FF,DA,DDA
      WRITE(6,57) F1(XX),DF(XX),DDF(XX)
      GO TO 2000
      END

C     --------------------------------------------------

C	Evaluating the fitted polynomial and its derivatives at any value
C	of x using known coefficients of orthogonal polynomials
C	Should be used to evaluate the polynomial using coefficients calculated
C	by POLFIT.
C
C	M : (input) Degree of polynomial
C	A : (input) Real array of length M+1 containing the coefficients
C		of the fit
C	ALP, BETA : (input) Real arrays of length M+1, containing the coefficients
C		required for defining the orthogonal polynomials
C		A, ALP, BETA could be calculated using POLFIT
C	X : (input) Value of x at which polynomial needs to be evaluated
C	F : (output) Calculated value of the polynomial at X
C	DF : (output) First derivative of F(X) at X
C	DDF : (output) Second derivative of F(X) at X
C	
C	Required routines : None

      SUBROUTINE POLEVL(M,A,ALP,BETA,X,F,DF,DDF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+1),ALP(M+1),BETA(M+1)

      F=A(M)+(X-ALP(M))*A(M+1)
      F1=A(M+1)
      DF=A(M+1)
      DF1=0.0
      DDF=0.0
      DDF1=0.0
      IF(M.LE.1) RETURN

C	Clenshaw's recurrence for F, DF and DDF
      DO 1000 J=M-2,0,-1
        DD=2.*DF+(X-ALP(J+1))*DDF-BETA(J+2)*DDF1
        D=F+(X-ALP(J+1))*DF-BETA(J+2)*DF1
        FF=A(J+1)+(X-ALP(J+1))*F-BETA(J+2)*F1
        F1=F
        F=FF
        DF1=DF
        DF=D
        DDF1=DDF
        DDF=DD
1000  CONTINUE
      END

C     -----------------------------------------------------------

C	Least squares polynomial fit using orthogonal polynomials in 1 dimension
C
C	N : (input) Number of data points
C	M : (input) Required degree of polynomial
C	X : (input) Real array of length N containing the abscissas
C	F : (input) Real array of length N containing the function values
C		F(I) is the value at X(I)
C	W : (input) Real array of length N containing the weights associated
C		with each point. W(I) is the weight for F(I), normally,
C		W(I) should be 1/err(I)**2, where err(I) is error in F(I)
C	A : (output) Real array of length M+1 containing the coefficients
C		for the fit
C	ALP, BETA : (output) Real arrays of length M+1, containing the coefficients
C		required for defining the orthogonal polynomials
C	Y : (output) Real array of length N containing the fitted values at X(I)
C	H : (output) Real array of length M+1 containing the \chi^2 values
C		for residuals using polynomial fit of degrees 0,..,M
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

      SUBROUTINE POLFIT(N,M,X,F,W,A,ALP,BETA,Y,H,GAM,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N),W(N),A(M+1),ALP(M+1),BETA(M+1),Y(N),WK(N,2),
     1          H(M+1),GAM(M+1)

      IF(M.GE.N.OR.M.LT.0.OR.N.LT.1) THEN
        IER=601
        RETURN
      ENDIF

C	Initialisation
      I1=1
      I2=2
      GAM0=0.0
      H0=0.0
      DO 2000 I=1,N
        WK(I,I1)=0.0
        WK(I,I2)=1.0
        Y(I)=0.0
        GAM0=GAM0+W(I)
        H0=H0+W(I)*F(I)**2
2000  CONTINUE
      GAM(1)=GAM0
      BETA(1)=0
      IER=0

C	Loop over the degree of polynomial
      DO 5000 J=0,M
        S=0.0
        DO 3000 I=1,N
3000    S=S+W(I)*F(I)*WK(I,I2)
        IF(GAM(J+1).LE.0) THEN
          IER=621
          RETURN
        ENDIF

C	The coefficient a_j
        A(J+1)=S/GAM(J+1)
        H0=H0-A(J+1)**2*GAM(J+1)
        H(J+1)=H0
        DO 3200 I=1,N
3200    Y(I)=Y(I)+A(J+1)*WK(I,I2)
        IF(J.EQ.M) RETURN
        S=0
        DO 3400 I=1,N
3400    S=S+W(I)*X(I)*WK(I,I2)**2
C	The coefficient \alpha_{j+1}
        ALP(J+1)=S/GAM(J+1)
        GAM0=0.0
        DO 3600 I=1,N
          WK(I,I1)=(X(I)-ALP(J+1))*WK(I,I2)-BETA(J+1)*WK(I,I1)
          GAM0=GAM0+W(I)*WK(I,I1)**2
3600    CONTINUE
C	The coefficient \beta_{j+1}
        BETA(J+2)=GAM0/GAM(J+1)
C	The coefficient \gamma_{j+1}
        GAM(J+2)=GAM0
C	Interchange indices I1, I2 so that only last two columns of WK are stored
        IT=I1
        I1=I2
        I2=IT
5000  CONTINUE

      END

C     --------------------------------------------------------

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
