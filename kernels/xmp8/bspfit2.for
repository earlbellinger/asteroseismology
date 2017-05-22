C     LEAST SQUARES FIT USING B-SPLINES IN 2 DIMENSIONS

      PROGRAM BSP2
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION YO(90),XO(90),T(900),B(100,90),B1(90),DB(100,52)
      DIMENSION AX(100,90),WK(90000),F(100,90),XF(90),YF(90)
      DIMENSION AY(100,90),INTY(90),ZO(90),ZF(90),VX(90,90),VY(90,90)
      DIMENSION EF(100,90),Y0(100,90),A2(9000,1000),V2(1000,1000)
 
      FUN(X,Y)=SIN(X)*SIN(Y)+1.D-5*RANGAU(IS)
      F0(X,Y)=SIN(X)*SIN(Y)
      DF0X(X,Y)=COS(X)*SIN(Y)
      DF0Y(X,Y)=SIN(X)*COS(Y)
      DF0XX(X,Y)=-SIN(X)*SIN(Y)
      DF0XY(X,Y)=COS(X)*COS(Y)
      DF0YY(X,Y)=-SIN(X)*SIN(Y)
 
51    FORMAT('   IER =',I4,3X,'NO. OF POINTS =',I5,3X,'NO. OF KNOTS =',
     1       I3,4X,'ORDER OF B-SPLINE =',I3/5X,
     2       'REGULARISATION PARAMETER =',1PD12.4,5X,
     3       'ORDER OF DERIVATIVE =',I2/5X,'CHI SQUARE =',E12.4,5X,
     4       'IT =',I3)
52    FORMAT('   IER =',I4,3X,'X,Y =',1P2D14.6,5X,'F(X) =',D14.6/
     1       5X,4HF' =,2D14.6/5X,5HF'' =,3D14.6)
53    FORMAT('   EXACT VALUES :   F(X) =',1PD14.6/
     1       5X,4HF' =,2D14.6/5X,5HF'' =,3D14.6)
54    FORMAT('   COEFFICIENTS : '/(1P5D14.6))
 
      NX=6
C	For simplicity assume that no. of points and knots are the same
C	along both axes

100   PRINT *,'TYPE NX=NO. OF POINTS, K=ORDER OF B-SPLINE,'
      PRINT *,'     NO=NO. OF KNOTS, RLM=REGULARISATION PARAMETER,'
	  PRINT *,'     IDE=ORDER OF DERIVATIVE FOR REGULARISATION'
	  READ *,NX,K,NO,RLM,IDE

	  PRINT *,'TYPE IT =1/2 FOR BSPFIT2/BSPFITW2'
      PRINT *,'              (QUITS WHEN IT<1)'
	  READ *,IT
      IF(IT.LT.1.OR.IT.GT.2) STOP
 
C   Set the seed so that same random numbers are used every time
      IS=2
C	Set up the table of values with random error added
      H=1.0D0/(NX-1.)
      NY=NX
      DO 1000 I=1,NX
        XO(I)=(I-1)*H
        YO(I)=(I-1)*H
1000  CONTINUE
      H=1.D0/(NO-1.)
      DO 1200 I=1,NO
        XF(I)=(I-1)*H
        YF(I)=(I-1)*H
1200  CONTINUE
        DO 1500 I=1,NY
          DO 1500 J=1,NX
            F(J,I)=FUN(XO(J),YO(I))
            EF(J,I)=1.D-5
1500  CONTINUE
      IC=4
      LA=100
      L2=9000
      IV2=1000
      IV=90
      IFLG=0
      REPS=1.D-12
	  IF(IT.EQ.1) THEN
        CALL BSPFIT2(NX,NY,XO,YO,F,K,AX,AY,LA,VX,VY,IV,T,
     1     B1,B,XF,YF,NO,NO,Y0,WK,REPS,RLM,IDE,CHISQ,IER)
      ELSE
        CALL BSPFITW2(NX,NY,XO,YO,F,EF,K,A2,L2,V2,IV2,T,
     1     B,LA,XF,YF,NO,NO,Y0,WK,REPS,RLM,IDE,CHISQ,IER)
	  ENDIF
      WRITE(6,51) IER,NX,NO,K,RLM,IDE,CHISQ,IT
      WRITE(6,54) ((B(I,J),I=1,NO+K-2),J=1,NO+K-2)
      IDE=2
 
2000  PRINT *,'TYPE X, Y VALUE AT WHICH FUNCTION IS REQUIRED'
      PRINT *,'                 (QUITS WHEN XX<-100)'
      READ *,XX,YY
      IF(XX.LT.-100) GO TO 100
      F1= BSPEV2(NO,NO,XF,YF,K,IDE,B,LA,XX,YY,DFX,DFY,DFXX,DFXY,
     1           DFYY,WK,IER)
      WRITE(6,52) IER,XX,YY,F1,DFX,DFY,DFXX,DFXY,DFYY
      WRITE(6,53) F0(XX,YY),DF0X(XX,YY),DF0Y(XX,YY),DF0XX(XX,YY),
     1            DF0XY(XX,YY),DF0YY(XX,YY)
      GO TO 2000
 
 
      END
 
 
C     -------------------------------------------------
 
C	To evaluate a B-spline expansion in 2 dimensions
C
C	NX : (input) Number of knots to define B-splines along 1st dimension
C	NY : (input) Number of knots to define B-splines along 2nd dimension
C	X : (input) Real array of length NX containing the knots.
C		The knots must be distinct and in ascending order.
C	Y : (input) Real array of length NY containing the knots.
C		The knots must be distinct and in ascending order.
C	K : (input) Order of B-splines, K=4 for cubic B-splines
C	NDERIV : (input) Number of derivatives required
C		For NDERIV.LE.0 only function value is calculated
C		For NDERIV=1 first derivative is also calculated
C		For NDERIV>1 both first and second derivatives are calculated
C	WT : (input) Array of length IW*(NY+K-2) containing coefficients
C		of B-spline expansion,
C	IW : (input) First dimension of array WT as defined in the calling program
C	X0,Y0 : (input) The point at which expansion has to be evaluated
C	DFX : (output) First derivative of function w.r.t. X at X0, Y0
C	DFY : (output) First derivative of function w.r.t. Y at X0, Y0
C	DFXX : (output) Second derivative of function w.r.t X,X at X0, Y0
C	DFXY : (output) Second derivative of function w.r.t X,Y at X0, Y0
C	DFYY : (output) Second derivative of function w.r.t Y,Y at X0, Y0
C	WK : Scratch array of length 7*MAX(NX,NY)+8K+2
C	IER : (output) Error parameter, IER=0 implies successful execution
C		Nonzero values of IER may be set by BSPLIN which is called
C
C	BSPEV2 =SUM_{i=1}^{NX+K-2} SUM_{j=1}^{NY+K-2} WT(i,j)\phi_i(X0)\psi_j(Y0)
C	where \phi_i(x) are B-spline basis functions on knots X
C	and \psi_j(y) are B-spline basis functions on knots Y
C
C	Required routines : BSPLIN

      FUNCTION BSPEV2(NX,NY,X,Y,K,NDERIV,WT,IW,X0,Y0,DFX,DFY,
     1          DFXX,DFXY,DFYY,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NX),Y(NY),WT(IW,NY+K-2),WK(*)
 
      BSPEV2=0.0
      NK=MAX(NX,NY)+K
      CALL BSPLIN(X,NX,K,X0,NDERIV,WK,WK(NK),WK(2*NK),LX,IER,WK(3*NK))
      IF(IER.GT.100) RETURN
      CALL BSPLIN(Y,NY,K,Y0,NDERIV,WK(3*NK),WK(4*NK),WK(5*NK),LY,IER,
     1            WK(6*NK))
      IF(IER.GT.100) RETURN
 
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
      DO 2000 I=LX,LX+K-1
        DO 2000 J=LY,LY+K-1
          F=F+WT(I,J)*WK(I)*WK(N3+J)
          DFX=DFX+WT(I,J)*WK(N1+I)*WK(N3+J)
          DFY=DFY+WT(I,J)*WK(I)*WK(N4+J)
          DFXX=DFXX+WT(I,J)*WK(N2+I)*WK(N3+J)
          DFXY=DFXY+WT(I,J)*WK(N1+I)*WK(N4+J)
          DFYY=DFYY+WT(I,J)*WK(I)*WK(N5+J)
2000  CONTINUE
      BSPEV2=F
      END
 
C     -------------------------------------------------
 
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
 
C     --------------------------------------------------------
 
C	To calculate linear least squares fit to B-spline basis functions in
C	1 dimension
C
C	N : (input) Number of data points to be fitted
C	X : (input) Real array of length N containing the coordinates
C		of points at which function values is available
C	F : (input) Real array of length N containing the function values
C		F(I) should be function value at X(I)
C	EF : (input) Real array of length N containing the estimated error
C		in F(I). 
C	K : (input) Order of B-splines required, K=4 gives cubic B-splines
C	A : (output) Real array of length LA*(NO+K-2) containing the matrix
C		U of SVD of the design matrix
C	LA : (input) First dimension of A in the calling program (LA.GE.2N)
C	V : (output) Real array of length IV*(NO+K-2) containing the matrix
C		V of SVD of the design matrix
C	IV : (input) First dimension of V in the calling program (IV.GE.NO+K-2)
C	SIGMA : (output) Real array of length NO+K-2 containing the singular
C		values of the design matrix
C	C : (output) Real array of length 2N containing the fitted coefficients
C		Note that although the number of coefficients is NO+K-2, the
C		rest of array is used as scratch space
C	XF : (input) Real array of size NO, containing
C		the knots used for defining B-spline basis functions.
C		The knots must be distinct and in ascending order.
C	NO : (input) Number of knots for B-splines, the number of basis
C		functions would be NO+K-2
C	Y : (output) Real array of length N containing the values of fitted
C		function at each of the tabular points
C	IFLG : (input/output) Integer specifying the type of calculation required
C		IFLG=0 The matrix will be calculated and solved for coefficients
C			the fitted values Y and CHISQ are also calculated
C		IFLG=1 The matrix will be calculated and SVD
C			is obtained, but coefficients are not calculated
C		IFLG=2 The SVD of matrix is assumed to be available in
C			arrays A, V, SIGMA and coefficients C are calculated
C		IFLG=3 The SVD of matrix is assumed to be available in arrays
C			A, V, SIGMA and coefficients C are calculated and in
C			addition fitted values Y and CHISQ are also calculated
C	WK : Real array of length 4*NO+5*K+2 used as scratch space
C	REPS : (input) Required accuracy for solution of equations using SVD
C		singular values less than REPS times maximum will be set to zero
C	RLM : (input) Parameter lambda for smoothing. If RLM.LE.0 no smoothing
C		is applied
C	IDE : (input) Order of derivative to be used for smoothing
C		This is used only when RLM>0. IDE=1 for first derivative
C		and IDE=2 for second derivative smoothing
C	CHISQ : (output) The value of Chi square at minimum
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=608 implies that NO+K-2>N or K<2
C		IER=609 implies that RLM>0 and IDE is not acceptable
C		IER=610 implies that EF(I).LE.0 for some I
C		No calculations are done in all these cases
C		Other values of IER may be set by SVD or BSPLIN
C
C	Required routines : BSPLIN, BSPEVL, SVD, SVDEVL
C
      SUBROUTINE BSPFIT(N,X,F,EF,K,A,LA,V,IV,SIGMA,C,XF,NO,Y,IFLG,WK,
     1       REPS,RLM,IDE,CHISQ,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),A(LA,NO+K-2),WK(4*NO+5*K+2),C(N),F(N),XF(NO),EF(N),
     1   Y(N),SIGMA(NO+K-2),V(IV,NO+K-2)
 
      IF(N.LT.NO+K-2.OR.K.LT.2) THEN
        IER=608
        RETURN
      ENDIF
 
      IF(RLM.GT.0.0.AND.(IDE.LT.1.OR.IDE.GT.2)) THEN
        IER=609
        RETURN
      ENDIF
C	N1 is the number of equations to be solved
      N1=N
      IF(RLM.GT.0.0) N1=2*N
 
      IF(IFLG.LE.1) THEN
C	Set up the matrix equation and obtain SVD
C	M is the number of coefficients to be determined
        M=NO+K-2
        NDB=M+2
        NDERIV=0
        IF(RLM.GT.0.0) THEN
          NDERIV=IDE
          NB=NDB-1
          IF(IDE.EQ.2) NB=2*NDB-1
        ENDIF
 
C	Set up the matrix for equations
        DO 2500 I=1,N
          XB=X(I)
          IF(EF(I).LE.0) THEN
            IER=610
            RETURN
          ENDIF
          CALL BSPLIN(XF,NO,K,XB,NDERIV,WK,WK(NDB),WK(NDB*2),LEFT,IER,
     1                WK(3*NDB))
          IF(IER.GT.100) RETURN
          DO 2200 J=1,M
            A(I,J)=WK(J)/EF(I)
            IF(RLM.GT.0.0) A(I+N,J)=RLM*WK(NB+J)
2200      CONTINUE
2500    CONTINUE
        CALL SVD(M,N1,A,V,SIGMA,LA,IV,WK,IER)
        IF(IER.GT.100) RETURN
 
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
 
      ENDIF
 
C	Setup the RHS and solve the equations
      DO 3000 I=1,N
        C(I)=F(I)/EF(I)
        IF(RLM.GT.0.0) C(I+N)=0.0
3000  CONTINUE
 
      CALL SVDEVL(M,N1,A,V,SIGMA,LA,IV,C,WK,REPS)
 
      IF(IFLG.EQ.2) RETURN
      IFLG=2
 
C	Calculate the \chi^2
      CHISQ=0.0
      NDERIV=0
      DO 4000 I=1,N
        Y(I)=BSPEVL(NO,XF,K,NDERIV,C,X(I),DF,DDF,WK,IER)
        CHISQ=CHISQ+((F(I)-Y(I))/EF(I))**2
4000  CONTINUE
 
      END
 
C     --------------------------------------------------------
 
C	To calculate linear least squares fit to B-spline basis functions in 2 dimension
C	Weights are assumed to be unity
C
C	NX : (input) Number of data points along x-axis
C	NY : (input) Number of data points along y-axis
C	X,Y : (input) Real arrays of length NX,NY containing the coordinates
C		of points at which function values are available
C	F : (input) Real array of length LA*NY containing the function values
C		F(I,J) should be function value at X(I),Y(J)
C	K : (input) Order of B-splines required, K=4 gives cubic B-splines
C	AX : (output) Real array of length LA*(MX+K-2) containing the matrix
C		U of SVD of the design matrix for fit along x-axis
C	AY : (output) Real array of length LA*(MY+K-2) containing the matrix
C		U of SVD of the design matrix for fit along y-axis
C	LA : (input) First dimension of arrays F, AX, AY, C, FY as declared
C		in the calling program (LA .GE. 2*MAX(NX,NY))
C	VX : (output) Real array of length IV*(MX+K-2) containing the matrix
C		V of SVD of the design matrix for fit along x-axis
C	VY : (output) Real array of length IV*(MY+K-2) containing the matrix
C		V of SVD of the design matrix for fit along y-axis
C	IV : (input) First dimension of VX, VY in the calling program
C		IV .GE. MAX(MX,MY)+K-2
C	SIGMAX : (output) Real array of length MX+K-2 containing the singular
C		values of the design matrix for fit along x-axis
C	SIGMAY : (output) Real array of length MY+K-2 containing the singular
C		values of the design matrix for fit along y-axis
C	C : (output) Real array of length LA*NY containing the fitted coefficients
C		Note that although the number of coefficients is (MX+K-2)*(MY+K-2)
C		the rest of array is used as scratch space
C	XF : (input) Real array of size MX, containing the knots
C		along x-axis used for defining B-spline basis functions.
C		The knots must be distinct and in ascending order.
C	YF : (input) Real array of size MY, containing the knots
C		along y-axis used for defining B-spline basis functions.
C		The knots must be distinct and in ascending order.
C	MX : (input) Number of knots for B-splines along x-axis,
C		the number of basis functions would be MX+K-2
C	MY : (input) Number of knots for B-splines along y-axis,
C		the number of basis functions would be MY+K-2
C	FY : (output) Real array of length LA*NY containing the values of fitted
C		function at each of the tabular points
C	WK : Real array of length LA*NX+NX+NY used as scratch space
C	REPS : (input) Required accuracy for solution of equations using SVD
C		singular values less than REPS times maximum will be set to zero
C	RLM : (input) Parameter lambda for smoothing. If RLM.LE.0 no smoothing
C		is applied
C	IDE : (input) Order of derivative to be used for smoothing
C		This is used only when RLM>0. IDE=1 for first derivative
C		and IDE=2 for second derivative smoothing
C	CHISQ : (output) The value of Chi square at minimum
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=608 implies that MX+K-2>NX, MY+K-2>NY or K<2
C		IER=609 implies that RLM>0 and IDE is not acceptable
C		No calculations are done in all these cases
C		Other values of IER may be set by SVD or BSPLIN
C
C	Required routines : BSPFIT, BSPLIN, BSPEVL, BSPEV2, SVD, SVDEVL
C
      SUBROUTINE BSPFIT2(NX,NY,X,Y,F,K,AX,AY,LA,VX,VY,IV,SIGMAX,
     1        SIGMAY,C,XF,YF,MX,MY,FY,WK,REPS,RLM,IDE,CHISQ,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NX),Y(NY),AX(LA,MX+K-2),AY(LA,MY+K-2),XF(MX),YF(MY),
     1        C(LA,NY),F(LA,NY),VX(IV,MX+K-2),VY(IV,MY+K-2),
     2        SIGMAX(MX+K-2),SIGMAY(MY+K-2),FY(LA,NY),WK(NX*LA+NX+NY)
 
      N1=4*MAX(NX,NY)+K+12
C	Set the weights to 1
      DO 1000 I=0,MAX(NX,NY)
        WK(N1+I)=1.0
1000  CONTINUE

C	Calculate the SVD of matrix for fit along x-axis
      IFLG1=1
      CALL BSPFIT(NX,X,F,WK(N1),K,AX,LA,VX,IV,SIGMAX,C,XF,MX,FY,
     1        IFLG1,WK,REPS,RLM,IDE,CHISQ,IER)
      IF(IER.GT.100) RETURN

C	Calculate the SVD of matrix for fit along y-axis
      IFLG1=1
      CALL BSPFIT(NY,Y,F,WK(N1),K,AY,LA,VY,IV,SIGMAY,C,YF,MY,FY,
     1        IFLG1,WK,REPS,RLM,IDE,CHISQ,IER)
      IF(IER.GT.100) RETURN
 
C	Set up the RHS for calculating the fits along y-axis
      DO 2500 I=1,NX
        DO 2500 J=1,NY
          WK(J+(I-1)*LA)=F(I,J)
          IF(RLM.GT.0.0) WK(J+NY+(I-1)*LA)=0.0
2500  CONTINUE

C	N1 is the number of equations for fit along y-axis
      N1=NY
      IF(RLM.GT.0.0) N1=2*NY
      LN1=LA*NX+1
      M=MY+K-2
      DO 2800 I=1,NX
        NI=1+(I-1)*LA
        CALL SVDEVL(M,N1,AY,VY,SIGMAY,LA,IV,WK(NI),WK(LN1),REPS)
2800  CONTINUE
 
C	Set up the RHS for calculating the fits along x-axis
      DO 3000 J=1,M
        DO 3000 I=1,NX
          C(I,J)=WK(J+(I-1)*LA)
          IF(RLM.GT.0.0) C(I+NX,J)=0.0
3000  CONTINUE
      M1=MX+K-2
      N1=NX
      IF(RLM.GT.0.0) N1=2*NX
      DO 3200 I=1,M
        CALL SVDEVL(M1,N1,AX,VX,SIGMAX,LA,IV,C(1,I),WK,REPS)
3200  CONTINUE
 
C	Calculate the CHI square
      CHISQ=0.0
      NDERIV=0
      DO 4000 I=1,NY
        DO 4000 J=1,NX
          FY(J,I)=BSPEV2(MX,MY,XF,YF,K,NDERIV,C,LA,X(J),Y(I),DFX,DFY,
     1            DFXX,DFXY,DFYY,WK,IER1)
          CHISQ=CHISQ+(F(J,I)-FY(J,I))**2
4000  CONTINUE
 
      END
 
C     ------------------------------------------------------------
 
C	To calculate linear least squares fit to B-spline basis functions in 2 dimension
C	version for general weights, but much slower than BSPFIT2
C
C	NX : (input) Number of data points along x-axis
C	NY : (input) Number of data points along y-axis
C	X,Y : (input) Real array of length NX,NY containing the coordinates
C		of points at which function values is available
C	F : (input) Real array of length IC*NY containing the function values
C		F(I,J) should be function value at X(I),Y(J)
C	EF : (input) Real array of length IC*NY containing the error estimates in F(I,J)
C	K : (input) Order of B-splines required, K=4 gives cubic B-splines
C	A : (output) Real array of length LA*(MX+K-2)*(MY+K-2) containing
C		the matrix U of SVD of the design matrix
C	LA : (input) First dimension of array A as declared
C		in the calling program (LA .GE. 3*NX*NY)
C	V : (output) Real array of length IV*(MX+K-2)*(MY+K-2) containing
C		the matrix V of SVD of the design matrix
C	IV : (input) First dimension of V in the calling program
C		IV .GE. (MX+K-2)*(MY+K-2)
C	SIGMA : (output) Real array of length (MX+K-2)*(MY+K-2)
C		containing the singular values of the design matrix
C	C : (output) Real array of length IC*(MY+K-2) containing the
C		fitted coefficients
C	IC : (input) First dimension of arrays C, F, EF, FY as declared
C		in the calling program (IC .GE. NX)
C	XF : (input) Real array of size MX, containing the knots
C		along x-axis used for defining B-spline basis functions.
C	YF : (input) Real array of size MY, containing the knots
C		along y-axis used for defining B-spline basis functions.
C	MX : (input) Number of knots for B-splines along x-axis,
C		the number of basis functions would be MX+K-2
C	MY : (input) Number of knots for B-splines along y-axis,
C		the number of basis functions would be MY+K-2
C	FY : (output) Real array of length IC*NY containing the values of fitted
C		function at each of the tabular points
C	WK : Real array of length 3*NX*NY+(MX+K)*(MY+K) used as scratch space
C	REPS : (input) Required accuracy for solution of equations using SVD
C		singular values less than REPS times maximum will be set to zero
C	RLM : (input) Parameter lambda for smoothing. If RLM.LE.0 no smoothing
C		is applied
C	IDE : (input) Order of derivative to be used for smoothing
C		This is used only when RLM>0. IDE=1 for first derivative
C		and IDE=2 for second derivative smoothing
C	CHISQ : (output) The value of Chi square at minimum
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=608 implies that MX+K-2>NX, MY+K-2>NY or K<2
C		IER=609 implies that RLM>0 and IDE is not acceptable
C		IER=610 implies that EF(I,J).LE.0 for some I,J
C		No calculations are done in all these cases
C		Other values of IER may be set by SVD or BSPLIN
C
C	Required routines : BSPLIN, BSPEV2, SVD, SVDEVL
C
      SUBROUTINE BSPFITW2(NX,NY,X,Y,F,EF,K,A,LA,V,IV,SIGMA,
     1      C,IC,XF,YF,MX,MY,FY,WK,REPS,RLM,IDE,CHISQ,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NX),Y(NY),A(LA,*),XF(MX),YF(MY),C(IC,*),F(IC,NY),
     1      EF(IC,NY),V(IV,*),SIGMA(*),FY(IC,NY),WK(*)
 
      IF(NX.LT.K.OR.NY.LT.K.OR.K.LT.2) THEN
        IER=608
        RETURN
      ENDIF
 
      IF(RLM.GT.0.0.AND.(IDE.LT.1.OR.IDE.GT.2)) THEN
        IER=609
        RETURN
      ENDIF
C     N1 is the number of equations to be solved
      N1=NX*NY
      IF(RLM.GT.0.0) N1=3*N1
 
C     Set up the matrix equation and obtain SVD
C     M is the number of coefficients to be determined
      M1=MX+K-2
      M2=MY+K-2
      M=M1*M2
      NDB=MAX(MX,MY)+5
      NDERIV=0
      IF(RLM.GT.0.0) THEN
        NDERIV=IDE
        NB=NDB-1
        NB1=NDB*4-1
        IF(IDE.EQ.2) NB=2*NDB-1
        IF(IDE.EQ.2) NB1=5*NDB-1
      ENDIF
 
C     Set up the matrix for equations
      DO 2500 I=1,NX
        XB=X(I)
        CALL BSPLIN(XF,MX,K,XB,NDERIV,WK,WK(NDB),WK(NDB*2),LEFT,IER,
     1              WK(NDB*3))
        IF(IER.GT.100) RETURN
        DO 2500 J=1,NY
          YB=Y(J)
          CALL BSPLIN(YF,MY,K,YB,NDERIV,WK(NDB*3),WK(NDB*4),WK(NDB*5),
     1                LEFT,IER,WK(NDB*6))
          IF(IER.GT.100) RETURN
          IF(EF(I,J).LE.0) THEN
            IER=610
            RETURN
          ENDIF
          IJ=I+(J-1)*NX
          DO 2200 K1=1,M1
            DO 2200 K2=1,M2
              IJ1=K1+(K2-1)*M1
              A(IJ,IJ1)=WK(K1)*WK(NDB*3+K2-1)/EF(I,J)
              IF(RLM.GT.0) A(IJ+NX*NY,IJ1)=RLM*WK(NB+K1)*WK(NDB*3+K2-1)
              IF(RLM.GT.0.0) A(IJ+2*NX*NY,IJ1)=RLM*WK(K1)*WK(NB1+K2)
2200      CONTINUE
2500  CONTINUE
      CALL SVD(M,N1,A,V,SIGMA,LA,IV,WK,IER)
      IF(IER.GT.100) RETURN
 
 
 
C     Setup the RHS and solve the equations
      DO 3000 I=1,NX
        DO 3000 J=1,NY
          IJ=I+(J-1)*NX
          WK(IJ)=F(I,J)/EF(I,J)
          IF(RLM.GT.0.0) WK(IJ+NX*NY)=0.0
          IF(RLM.GT.0.0) WK(IJ+2*NX*NY)=0.0
3000  CONTINUE
 
      CALL SVDEVL(M,N1,A,V,SIGMA,LA,IV,WK,WK(N1+2),REPS)
 
      DO 3400 I=1,M1
        DO 3400 J=1,M2
          C(I,J)=WK(I+M1*(J-1))
3400  CONTINUE
 
 
C     Calculate the \chi^2
      CHISQ=0.0
      NDERIV=0
      DO 4000 I=1,NY
        DO 4000 J=1,NX
          FY(J,I)=BSPEV2(MX,MY,XF,YF,K,NDERIV,C,IC,X(J),Y(I),DFX,DFY,
     1                   DFXX,DFXY,DFYY,WK,IER1)
          CHISQ=CHISQ+((F(J,I)-FY(J,I))/EF(J,I))**2
4000  CONTINUE
 
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
 
C     ---------------------------------------------------------------
 
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
 
C     ------------------------------------------------------
 
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
 
