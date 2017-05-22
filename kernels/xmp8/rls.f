C	To solve an inverse problem with specified kernels.
C	It first generates artificial data using a specified solution using
C	FORW and then these data are used for inversion

      PROGRAM INVERS
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUN
      PARAMETER(NPT=1000,NMD=2000,NV=100,PI=3.1415926535D0)
      DIMENSION XO(90),A(NMD,NV),AV(NV,NV),B(NMD),SIGMA(NV),F(NPT)
      DIMENSION R(NPT),AC(NMD,NV),RKER(NMD,NPT),DI(NMD),DF(NMD),DE(NMD)
      DIMENSION F0(NPT),FE(NPT),WK(NPT*70)

C	EXERCISE 12.23
 
51    FORMAT(I6,1P9D14.6)
52    FORMAT('  IER =',I4,5X,'CHI SQUARE =',1PD12.4,3X,
     1       'REGULARISATION TERM =',D12.4/2X,'COEFFICIENTS :',5D14.6
     2       /(5D14.6))
53    FORMAT('   NO. OF POINTS IN KERNEL =',I5,5X,
     1       'NO. OF DATA POINTS =',I5)
54    FORMAT('   ORDER OF B-SPLINES =',I3,5X,'NO. OF KNOTS =',I4,5X,
     1       'IFLG =',I2/5X,'REGULARISATION PARAMETER =',1PD14.6,5X,
     2       'ORDER OF DERIVATIVE FOR SMOOTHING =',I2)
55    FORMAT('# S. NO.',2X,'INPUT DATA',2X,'ERROR',8X,'RESIDUALS')
56    FORMAT('#  IDE',2X,'ALP',9X,'CHI SQUARE',4X,'SUMD')
57    FORMAT(1P9D14.6)
58    FORMAT('#  R',12X,'INVERTED SOLN    ERROR',7X,'EXACT VALUE')
 
      OPEN(UNIT=12,FILE='kernel',STATUS='OLD')
      OPEN(UNIT=21,FILE='invers.out',STATUS='UNKNOWN')
      OPEN(UNIT=22,FILE='lcurv.out',STATUS='UNKNOWN')
 
 
C     SEED FOR RANDOM NUMBER
      IDUM=1
      ES=0.2

C     Read in the kernels
      READ(12,*)NUMPT
      READ(12,*)(R(I),I=1,NUMPT)
      DO 3000 I=1,NMD
        READ(12,*,END=3200) N
        READ(12,*)(RKER(I,J),J=1,NUMPT)
        NMODE=I
        DE(NMODE)=ES
3000  CONTINUE
 
3200  CLOSE(12)
      WRITE(6,53)NUMPT,NMODE
      WRITE(22,56)
 
      IFLG=0
      IDE=2
      IK=NMD
      IV=NV
C	No. of points to be used for smoothing
      NS=400
C	No. of sets to be used for simulations
      NSIM=20

C	Calculate the data points using FORW
      CALL FORW(NUMPT,NMODE,R,RKER,IK,DI,F0,FUN,IER,IFLG)
      REPS=1.D-6

C     Add random errors to calculated data points
      DO 3500 I=1,NMODE
        DI(I)=DI(I)+DE(I)*RANGAU(IDUM)
3500  CONTINUE

4000  PRINT *,'TYPE K=ORDER OF B-SPLINE, NW=NO. OF KNOTS,'
      PRINT *,'IDE = ORDER OF DERIVATIVE FOR SMOOTHING,'
      PRINT *,'IFLG = FLAG,  ALP = REGULARISATION PARAMETER'
      PRINT *,'                             (QUITS WHEN K.LE.1)'
      READ *,K,NW,IDE,IFLG,ALP
      IF(K.LE.1) STOP
      WRITE(6,54) K,NW,IFLG,ALP,IDE
 
C     Set up the knots, with uniform spacing
      H=(R(NUMPT)-R(1))/(NW-1)
      DO 4200 I=1,NW
        XO(I)=R(1)+H*(I-1)
4200  CONTINUE
      REWIND 21
 
      CALL RLS(NW,XO,K,NUMPT,R,RKER,IK,AC,NMODE,NS,ALP,IDE,DI,DE,
     1    DF,F,B,IFLG,IER,REPS,CHISQ,SUMD,A,AV,IV,SIGMA,NSIM,FE,WK)
 
      WRITE(6,52) IER,CHISQ,SUMD,(B(I),I=1,NW+K-2)
 
C     Write out chisq and sumd in lcurv.out for plotting L-curve
      WRITE(22,51) IDE,ALP,CHISQ,SUMD
 
C     Write out the fitted data points in invers.out for comparison
C     only the last fit is preserved
C     First NMODE+1 lines contain the residuals in fitted points and
C     subsequent lines contain the inverted solution compared with
C     exact values.

      WRITE(21,55)
      DO 4800 I=1,NMODE
        WRITE(21,51) I,DI(I),DE(I),DF(I)
4800  CONTINUE

      WRITE(21,58)
      DO 5000 I=1,NUMPT
        WRITE(21,57) R(I),F(I),FE(I),F0(I)
5000  CONTINUE
      GO TO 4000
      END
 
C     ------------------------------------------
 
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
 
C     ------------------------------------------
 
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
 
C     ------------------------------------------------------------
 
C	To solve the forward problem corresponding to a linear inverse
C	problem. This routine may be used to generate artificial data
C	for testing inversion techniques.
C
C	NP : (input) Number of points used in defining the kernels.
C	NM : (input) Number of data points in the inversion problem
C		which should be same as the number of kernels that are
C		supplied in array RKER.
C	R : (input) Real array of length NP containing the coordinates
C		of points at which kernels are available.
C	RKER : (input) Real array of length IK*NP containing the kernels
C		for the inverse problem. RKER(I,J) should contain the
C		value at R(J) for the Ith kernel.
C	IK : (input) The first dimension of RKER as declared in the calling
C		program
C	DI : (output) Real array of length NM containing the calculated
C		data points using the kernel.
C	F : (input/output) Real array of length NP containing the function
C		value at points in R. F(I) should contain the function
C		value at R(I). If IFLG=0, the function values are
C		calculated using user supplied function routine FUN,
C		otherwise, these values must be supplied while calling
C		the routine.
C	FUN : (input) Name of function routine to calculate the given
C		function. This is used only if IFLG=0, otherwise the
C		function values are to be supplied in array F.
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=711 implies that IK<NM and no calculations are done
C	IFLG : (input/output) Integer parameter used as a flag to decide
C		the type of computation required.
C		If IFLG=0, then the function values are calculated using
C			a user supplied routine FUN. These values are stored
C			in array F and IFLG is set to 1 so that next time
C			the values need not be calculated.
C		For other values of IFLG the function values must be
C			supplied in array F.
C
C	FUNCTION FUN(X) must be supplied by the user
C
C	Required routines : FUN

      SUBROUTINE FORW(NP,NM,R,RKER,IK,DI,F,FUN,IER,IFLG)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(NP),RKER(IK,NP),DI(NM),F(NP)
 
      IF(IK.LT.NM) THEN
        IER=711
        RETURN
      ENDIF
      IER=0

      IF(IFLG.EQ.0) THEN
C     Calculate the function value using supplied routine
        DO 2000 I=1,NP
          F(I)=FUN(R(I))
2000    CONTINUE
        IFLG=1
      ENDIF
 
C     Calculate the integrals
      DO 3000 I=1,NM
        S1=0.0
        H=(R(2)-R(1))/2.
        DO 2500 IR=1,NP
          S1=S1+H*F(IR)*RKER(I,IR)
          IF(IR.LT.NP-1) THEN
            H=(R(IR+2)-R(IR))/2.0
          ELSE IF(IR.EQ.NP-1) THEN
            H=(R(IR+1)-R(IR))/2.0
          ENDIF
2500    CONTINUE
        DI(I)=S1
3000  CONTINUE
 
 
      END
 
C     --------------------------------------------------------------
 
C	To solve a linear inversion problem in one dimension using RLS
C	technique with B-spline basis functions
C
C	NK : (input) Number of knots for defining B-splines, the number
C		of basis functions would be NK+K-2
C	XO : (input) Real array of length NK containing the knots
C		used for defining B-spline basis functions.
C		The knots must be distinct and in ascending order.
C	K : (input) Order of B-splines required, K=4 gives cubic B-splines
C	NR : (input) Number of points used in defining the kernels
C	R : (input) Real array of length NR containing the coordinates
C		of points at which kernels are available.
C	RKER : (input) Real array of length IK*NR containing the kernels
C		for the inverse problem. RKER(I,J) should contain the
C		value at R(J) for the Ith kernel. This array
C		must be supplied if IFLG<2, otherwise it is not required
C	IK : (input) First dimension of arrays RKER, AC and A, as specified
C		in the calling program. IK.GE.NM+NS
C	AC : (input/output) Real array of length IK*(NK+K-2) containing
C		the coefficients of matrix defining the inversion problem
C		If IFLG<2, these coefficients are calculating by integrating
C		the kernels with appropriate weights. For IFLG=2,3 these
C		coefficients must be supplied.
C	NM : (input) Number of data points in the inversion problem
C	NS : (input) Number of points to be used for applying regularisation
C		The routine chooses a uniform mesh covering the full interval
C		for applying smoothing.
C	ALP : (input) Regularisation parameter, ALP>0.
C	IDE : (input) Order of derivative to be used for regularisation,
C		IDE should be 1 or 2 for first or second derivative smoothing
C	DI : (input) Real array of length NM, containing the data points
C		for inversion
C	DE : (input) Real array of length NM, containing the estimated
C		error in DI. 
C	DF : (output) Real array of length NM, containing the normalised
C		residuals (DI-DI(fit))/DE for each data point.
C	F : (output) Real array of length NR which will contain the
C		calculated solution at each point in array R.
C	B : (output) Real array of length NM+NS containing the coefficients
C		of basis functions in fitted solution. Although the
C		the number of coefficients is only NK+K-2, the rest
C		of array is used as scratch space
C	IFLG : (input/output) Integer specifying the type of calculation required.
C		IFLG=0 : The matrix coefficients will be calculated using
C			the kernels and then the equations are solved to
C			find the solution for given data points DI.
C			IFLG is set to 4 after calculations.
C		IFLG=1 : The matrix coefficients will be calculated using
C			the kernels and then the SVD of the full matrix is
C			computed, but the solution is not computed.
C			IFLG is set to 4 after calculations.
C		IFLG=2 : The matrix coefficients are assumed to be available
C			in array AC and the matrix is setup and solved
C			to find the solution for given data points DI.
C			IFLG is set to 4 after calculations.
C		IFLG=3 : The matrix coefficients are assumed to be available
C			in array AC and the matrix is setup and the SVD
C			is computed, but the solution is not computed.
C			IFLG is set to 4 after calculations.
C		IFLG=4 : The SVD of matrix is assumed to be available
C			from previous calculations and only the solution
C			for given DI is computed.
C		Since IFLG is set to 4 every-time, it should be reset to
C		0 or 2 before next call when the data or error estimates
C		or smoothing are changed.
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=709 implies that NM.LE.NK+K-2 or IK<NM+NS or IV<NK+K-2
C		IER=710 implies that ALP<0 or IDE<1 or IDE>2
C		other values of IER may be set by BSPLIN or SVD or BSPEVL
C	REPS : (input) Required accuracy for solution of equations using
C		SVD. Singular values less than REPS times maximum will be
C		set to zero.
C	CHISQ : (output) The computed value of Chi square for the solution
C	SUMD : (output) The computed value of the smoothing term
C	A : (input/output) Real array of length IK*(NK+K-2) containing
C		the SVD of the matrix of equations. If IFLG<4 this matrix
C		will be calculated, otherwise it must be supplied.
C	AV : (input/output) Real array of length IV*(NK+K-2) containing
C		the matrix V or SVD of the matrix of equations. if IFLG<4
C		this matrix will be calculated, otherwise it must be supplied.
C	IV : (input) The first dimension of AV as declared in the calling
C		program. IV.GE.NK+K-2
C	SIGMA : (input/output) Real array of length NK+K-2 containing the
C		singular values of the matrix A. If IFLG<4 this array will
C		be calculated, otherwise it must be supplied.
C	NSIM : (input) Number of sets to be tried for simulations to
C		calculate the error estimates. If NSIM.LE.1 error estimates
C		are not calculated.
C	FE : (output) Real array of length NR containing the estimated
C		error in F(I). This is calculated only if NSIM>1.
C	WK : Real array of length NR*NSIM+NM+NS+4NK+5K+2, used as scratch space.
C	
C
C	Required routines : BSPLIN, BSPEVL, SVD, SVDEVL, RANGAU
C	
C	
      SUBROUTINE RLS(NK,XO,K,NR,R,RKER,IK,AC,NM,NS,ALP,IDE,DI,DE,DF,F,
     1           B,IFLG,IER,REPS,CHISQ,SUMD,A,AV,IV,SIGMA,NSIM,FE,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XO(NK),R(NR),RKER(IK,NR),AC(IK,NK+K),DI(NM),DE(NM),
     1   DF(NM),F(NR),B(NK+K),A(IK,NK+K),AV(IV,NK+K),SIGMA(NK+K),D0(2),
     2   FE(NR),WK(NR,*)
 
      NV=NK+K-2
      NE=NM+NS
      IF(NM.LE.NV.OR.IK.LT.NE.OR.IV.LT.NV) THEN
        IER=709
        RETURN
      ENDIF
      IF(ALP.LT.0.0.OR.IDE.LT.1.OR.IDE.GT.2) THEN
        IER=710
        RETURN
      ENDIF
      IER=0
 
      IF(IFLG.LT.2) THEN
C     Setting up the system of linear equations
        NDERIV=0
        DO 3000 I=1,NM
 
          DO 2000 J=1,NV
            SIGMA(J)=0.0
2000      CONTINUE
          H=(R(2)-R(1))/2.
          DO 2500 IR=1,NR
            CALL BSPLIN(XO,NK,K,R(IR),NDERIV,DF,AV,AV(1,2),LEFT,IER,WK)
            IF(IER.GT.100) RETURN
 
            DO 2200 J=1,NV
              SIGMA(J)=SIGMA(J)+H*DF(J)*RKER(I,IR)
2200        CONTINUE
 
            IF(IR.LT.NR-1) THEN
              H=(R(IR+2)-R(IR))/2.0
            ELSE IF(IR.EQ.NR-1) THEN
              H=(R(IR+1)-R(IR))/2.0
            ENDIF
2500      CONTINUE
 
          DO 2800 J=1,NV
            A(I,J)=SIGMA(J)/DE(I)
            AC(I,J)=SIGMA(J)
2800      CONTINUE
3000    CONTINUE
 
      ELSE IF(IFLG.LT.4) THEN
C	The coefficients of matrix are available
        DO 3200 I=1,NM
          DO 3200 J=1,NV
            A(I,J)=AC(I,J)/DE(I)
3200    CONTINUE
 
      ENDIF
 
      IF(IFLG.LT.4) THEN
C     The equations arising from regularisation term
        H=(R(NR)-R(1))/(NS-1)
        NDERIV=IDE
        FA=ALP*SQRT(H)
 
        DO 4000 I=1,NS
          XI=R(1)+H*(I-1)
          CALL BSPLIN(XO,NK,K,XI,NDERIV,DF,AV,AV(1,2),LEFT,IER,WK)
          IF(IER.GT.100) RETURN
          DO 3500 J=1,NV
            A(I+NM,J)=AV(J,IDE)*FA
3500      CONTINUE
4000    CONTINUE
 
        CALL SVD(NV,NE,A,AV,SIGMA,IK,IV,DF,IER)
        IF(IER.GT.0) RETURN
      ENDIF
      IF(IFLG.EQ.1.OR.IFLG.EQ.3) THEN
        IFLG=4
        RETURN
      ENDIF
      IFLG=4
 
C     Set up the RHS of equations
      DO 4200 I=1,NM
        B(I)=DI(I)/DE(I)
4200  CONTINUE
      DO 4300 I=NM+1,NM+NS
        B(I)=0.0
4300  CONTINUE

C	Solve the system of equations using SVD
      CALL SVDEVL(NV,NE,A,AV,SIGMA,IK,IV,B,DF,REPS)
 
      NDERIV=0
      DO 4400 I=1,NR
        F(I)=BSPEVL(NK,XO,K,NDERIV,B,R(I),DF0,DDF,WK,IER)
4400  CONTINUE
 
C	Calculate the smoothing term
      SUMD=0.0
      NDERIV=IDE
      H=(R(NR)-R(1))/(NS-1)
      DO 4600 I=1,NS
        XI=R(1)+(I-1)*H
        FI=BSPEVL(NK,XO,K,NDERIV,B,XI,D0(1),D0(2),WK,IER)
        SUMD=SUMD+D0(IDE)**2
4600  CONTINUE
      SUMD=SUMD*H
 
C	Calculate the Chi square
      CHISQ=0.0
      DO 5000 I=1,NM
        S1=0.0
        DO 4800 J=1,NV
          S1=S1+AC(I,J)*B(J)
4800    CONTINUE
        DF(I)=(DI(I)-S1)/DE(I)
        CHISQ=CHISQ+DF(I)**2
5000  CONTINUE
 
C	Calculate the error estimates in the solution */
      IF(NSIM.LE.1) RETURN
C     Seed for random number generator
      ISD=123
      NDERIV=0
      DO 5100 I=1,NR
        FE(I)=0.0
5100  CONTINUE
 
      DO 6000 IS=1,NSIM
        DO 5200 I=1,NM
          WK(I,NSIM+1)=DI(I)/DE(I)+RANGAU(ISD)
5200    CONTINUE
        DO 5300 I=NM+1,NM+NS
          WK(I,NSIM+1)=0.0
5300    CONTINUE
 
        CALL SVDEVL(NV,NE,A,AV,SIGMA,IK,IV,WK(1,NSIM+1),DF,REPS)
 
        DO 5400 I=1,NR
          WK(I,IS)=BSPEVL(NK,XO,K,NDERIV,WK(1,NSIM+1),R(I),DF0,DDF,
     1            WK(NM+NS+2,NSIM+1),IER)
          FE(I)=FE(I)+WK(I,IS)
5400    CONTINUE
6000  CONTINUE
 
      DO 6500 I=1,NR
        A1=FE(I)/NSIM
        S1=0.0
        DO 6200 IS=1,NSIM
          S1=S1+(WK(I,IS)-A1)**2
6200    CONTINUE
        FE(I)=SQRT(S1/NSIM)
6500  CONTINUE
 
      END
 
C     --------------------------------------------------------------
 
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
 
C     ----------------------------------------------------------

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

C     --------------------------------------------------------------
 
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
 
C     --------------------------
 
C	Function for generating data points using FORW

      FUNCTION FUN(X)
      IMPLICIT REAL*8(A-H,O-Z)
 
      FUN=430+20*SIN(30*X)
      END
 
