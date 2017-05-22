C     PROGRAM TO CALCULATE WEIGHTS AND ABSCISSAS OF GAUSSIAN QUADRATURE
C     FORMULAS WITH A SPECIFIED WEIGHT FUNCTION.

      PROGRAM GAUSSL
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      EXTERNAL FMOM
      DIMENSION X(65),W(65)

C     EXERCISE  6.17: QUADRATURE FORMULAS FOR LOGARITHMIC SINGULARITY
 
51    FORMAT(5X,'IER =',I4/9X,'ABSCISSAS',18X,'WEIGHTS'/
     1       (5X,1PD22.14,5X,D22.14))
52    FORMAT(/'   MAXIMUM RELATIVE ERROR USING ',I4,' POINTS =',1PD12.4)
 

100   PRINT *,'TYPE  NT=NO. OF POINTS IN QUADRATURE FORMULA, QGAUS=T/F'
      PRINT *,'QGAUS=T  FOR GAUSSIAN FORMULA'
      PRINT *,'QGAUS=F  FOR NEWTON-COTES FORMULA  (QUITS WHEN NT.LE.0)'
      READ *, NT,QGAUS
      IF(NT.LE.0) STOP

      IF(.NOT.QGAUS) THEN

C     CHOOSE UNIFORMLY SPACED ABSCISSAS IN [0,1]

        DO 150 I=1,NT
150     X(I)=(I-1.D0)/(NT-1.)
      ENDIF 
      CALL GAUSWT(NT,W,X,FMOM,IER,QGAUS)
      WRITE(6,51) IER,(X(I),W(I),I=1,NT)

C     USE THE WEIGHTS AND ABSCISSAS TO EVALUATE THE INTEGRAL X**N FOR TESTING
C     The weights and abscissas may not agree perfectly with those given in
C     Appendix A, unless quadruple precision is used.
 
      ERR=0.0
      NT2=NT-1
      IF(QGAUS) NT2=2*NT-1
      RI=0.0
      DO 800 I=1,NT
        RI=RI+W(I)
800   CONTINUE
      REX=FMOM(0)
      ERR=MAX(ERR,ABS(RI-REX)/REX)

      DO 2000 J=1,NT2
        RI=0.0
        DO 1000 I=1,NT
          RI=RI+W(I)*X(I)**J
1000    CONTINUE
    	REX=FMOM(J)
        ERR=MAX(ERR,ABS(RI-REX)/REX)
2000  CONTINUE
 
      WRITE(6,52) NT,ERR
      GO TO 100
      END
 
C     --------------------------------------------------
 
C	To calculate weights and abscissas of a quadrature formula with
C	specified weight function.
C
C	N : (input) Number of points in the required quadrature formula
C	W : (output) Array of length N, which will contain the weights
C	AB : (input/output) Array of length N containing the abscissas
C		For Gaussian formulas (QGAUS=.TRUE.) AB will be calculated, 
C		while for QGAUS=.FALSE. abscissas must be supplied
C	FMOM : (input) Name of the function routine to calculate the moments
C		FUNCTION FMOM(I) should calculate integral of w(x)x**I
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=303 implies N.LE.0 or N.GE.NPMAX
C		IER=322 implies GAUELM failed to find coefficients of polynomial
C		IER=323 implies POLYR failed to find roots of polynomial
C		IER=324 implies GAUELM failed to find weights 
C	QGAUS : (input) Logical parameter to decide type of formula to be obtained
C		If QGAUS=.TRUE. a Gaussian formula is calculated. In this
C		case both abscissas and weights are calculated.
C		If QGAUS=.FALSE. an interpolatory formula is calculated.
C		In this case only weights are calculated, while abscissas
C		must be supplied.
C
C	Required routines : GAUELM, POLYR, LAGITR, FMOM

      SUBROUTINE GAUSWT(N,W,AB,FMOM,IER,QGAUS)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER (NPMAX=65)
      EXTERNAL FMOM
      DIMENSION A(NPMAX,NPMAX),INC(NPMAX),COF(NPMAX),W(N),AB(N)
      COMPLEX*16  ZERO(NPMAX)

      IER=303
      IF(N.LE.0.OR.N.GE.NPMAX) RETURN
      LJ=NPMAX
      IER=0
      IF(QGAUS) THEN

C	Calculating the coefficients of the orthogonal polynomial
        DO 2000 I=1,N
          DO 1500 J=1,N
1500      A(I,J)=FMOM(I+J-2)
          COF(I)=-FMOM(N+I-1)
2000    CONTINUE
        IFLG=0
        CALL GAUELM(N,1,A,COF,DET,INC,LJ,IER1,IFLG)
        IF(IER1.GT.100) THEN
          IER=322
          RETURN
        ENDIF

C	Find the roots of polynomial, which will be the abscissas
        COF(N+1)=1.
        QREFIN=.TRUE.
C	Array A is used as scratch space
        CALL POLYR(N,COF,ZERO,IER1,QREFIN,A)
        IF(IER1.GT.100) THEN
          IER=323
          RETURN
        ENDIF
        DO 2800 I=1,N
2800    AB(I)=ZERO(I)

      ENDIF

C	Calculate the weights
      DO 4000 I=1,N
        DO 3000 J=1,N
          IF(I.EQ.1) A(I,J)=1.0
          IF(I.GT.1) A(I,J)=AB(J)**(I-1)
3000    CONTINUE
4000  W(I)=FMOM(I-1)
      IFLG=0
      CALL GAUELM(N,1,A,W,DET,INC,LJ,IER1,IFLG)
      IF(IER1.GT.100) IER=324

      END
 
C     --------------------------------------------------
 
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
 
C  ------------------------------------------------------------
 
C	Root of a polynomial with real coefficients using Laguerre iteration
C
C	N : (input) The degree of the polynomial
C	A : (input) Real array of length N+1 containing the coefficients of
C		the polynomial. A(1) is the constant term and A(N+1) is the
C		coefficient of X**N
C	CXI : (input/output) Complex variable containing the initial guess,
C		 after execution it will contain the computed root
C	IER : (output) Error parameter, IER=0 for successful execution
C		IER=438 implies that denominator is zero and iteration cannot
C			be continued further
C		IER=439 implies that iteration has failed to converge
C
C	Required routines : None

      SUBROUTINE LAGITR(N,A,CXI,IER)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      PARAMETER(ITMAX=50,REPS=1.D-4,AEPS=1.D-6)
      DIMENSION A(N+1)

      IER=0
      CDX1=ABS(CXI)+1.0
      QC=.FALSE.
      IC=ITMAX

C	The Laguerre iteration
      DO 2000 I=1,ITMAX
C	Calculate the polynomial and its derivatives
        CF=A(N+1)
        CFP=0.0
        CFPP=0.0
        DO 1000 J=N,1,-1
          CFPP=CXI*CFPP+2.*CFP
          CFP=CXI*CFP+CF
          CF=CXI*CF+A(J)
1000    CONTINUE

        CH=(N-1)*((N-1)*CFP*CFP-N*CF*CFPP)
        CH=SQRT(CH)
        CDEN=CFP+CH
        IF(ABS(CFP-CH).GT.ABS(CDEN)) CDEN=CFP-CH

        IF(CDEN.NE.0.0) THEN
C	Laguerre's iteration
          CDX=-N*CF/CDEN
          IF(ABS(CDX).LT.MAX(REPS*ABS(CXI),AEPS).AND.I.GT.1.AND.
     1        (.NOT.QC)) THEN
            QC=.TRUE.
            IC=I
          ENDIF
          IF(QC.AND.ABS(CDX/CDX1).GT.1.0) RETURN
          IF(I-IC.GT.5.AND.ABS(CDX/CDX1).GT.0.99D0) RETURN
          CDX1=CDX
          IF(CDX.EQ.0.0) RETURN
          CXI=CXI+CDX
        ELSE
          IF(CF.EQ.0.0) RETURN
C	If denominator vanishes, then quit
          IER=438
          RETURN
        ENDIF
2000  CONTINUE

C	Iteration fails to converge
      IF(.NOT.QC) IER=439
      END
 
C     ----------------------------------------
 
C	Roots of a polynomial with real coefficients using Laguerre iteration
C
C	N : (input) The degree of the polynomial
C	A : (input) Real array of length N+1 containing the coefficients of
C		the polynomial. A(1) is the constant term and A(N+1) is the
C		coefficient of X**N
C	CX : (output) Complex array of length N, containing the computed roots
C	IER : (output) Error parameter, IER=0 for successful execution
C		IER=k*11 implies that iteration for refining the roots failed
C			to converge for k of the roots
C		IER=406 implies that N.LE.0 and no calculations are done
C		IER=408 implies that A(N+1)=0 and no calculations are done
C		IER=430 implies that iteration failed to converge for some root
C	QREFIN : (input) Logical parameter to decide if roots need to be refined
C		If QREFIN=.TRUE. the roots are refined using original polynomial
C		otherwise no refinement is done.
C	WK : Real array of length N+1 used as a scratch space
C
C	Required routines : LAGITR
C
      SUBROUTINE POLYR(N,A,CX,IER,QREFIN,WK)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      PARAMETER(EPS=1.D-16)
C	For REAL*4 use the following value
C      PARAMETER(EPS=1.E-7)
      DIMENSION A(N+1),CX(N),WK(N+1)

      IF(N.LE.0) THEN
        IER=406
        RETURN
      ENDIF
      IF(A(N+1).EQ.0.0) THEN
C	If the leading coefficient is zero, then quit
        IER=408
        RETURN
      ENDIF

C	Make a copy of coefficients to be used for deflation
      DO 1000 I=1,N+1
1000  WK(I)=A(I)
      NP=N
C	Starting value for iteration to allow the roots to be found
C	in ascending order of magnitude
      CXR=0.0
      IER=0

2000  CALL LAGITR(NP,WK,CXR,IER1)
      IF(IER1.NE.0) THEN
C	If the iteration fails to converge, then try once more
        CXR=1.123456
        CALL LAGITR(NP,WK,CXR,IER1)
        IF(IER1.NE.0) THEN
C	if it fails again, then quit
          IER=430
          RETURN
        ENDIF
      ENDIF

C	Intrinsic functions REAL and IMAG are not treated as generic
C	names on many compilers, and their interpretation may be
C	ambiguous as they may return a value of type REAL*4, even
C	though the argument is COMPLEX*16. But it does not matter
C	as far as the use in this subroutine is concerned.
C	In all other cases where it matters we have not used these
C	functions. If the compiler objects then replace REAL by DREAL.

      IF(ABS(IMAG(CXR)).LE.10.*EPS*ABS(REAL(CXR))) THEN
C	Perform deflation for a real root
        XRT=CXR
        IF(NP.GT.1) THEN
          BN=WK(NP)
          WK(NP)=WK(NP+1)
          DO 2400 I=NP-1,1,-1
            AN=XRT*WK(I+1)+BN
            BN=WK(I)
            WK(I)=AN
2400      CONTINUE
        ENDIF
        NP=NP-1
        CX(N-NP)=CXR

      ELSE
C	Perform deflation for a pair of complex conjugate roots
        XR=CXR
        XI=ABS(IMAG(CXR))
        IF(NP.GT.2) THEN
          P=2.*XR
          S=XR**2+XI**2
          BN=WK(NP-1)
          BN1=WK(NP-2)
          WK(NP-1)=WK(NP+1)
          WK(NP-2)=WK(NP)+WK(NP-1)*P
          DO 2600 I=NP-3,1,-1
            AN=BN+P*WK(I+1)-WK(I+2)*S
            BN=BN1
            BN1=WK(I)
            WK(I)=AN
2600      CONTINUE
        ENDIF
        NP=NP-2
        CX(N-NP-1)=DCMPLX(XR,XI)
        CX(N-NP)=DCMPLX(XR,-XI)
      ENDIF
C	If some roots are remaining, find next
      IF(NP.GT.0) GO TO 2000

      IF(QREFIN) THEN
C	Refine the roots by iterating on original polynomial
        DO 3000 I=2,N
          CXR=CX(I)
          CALL LAGITR(N,A,CXR,IER1)
          IF(IER1.EQ.0) THEN
C	If the iteration has converged, accept the root
            CX(I)=CXR
          ELSE
C	else retain the old approximation to the root
            IER=IER+11
          ENDIF
3000    CONTINUE
      ENDIF

C	Sort the roots in ascending order of real part
      DO 4000 I=2,N
        CXR=CX(I)
        DO 3500 J=I-1,1,-1
          IF(REAL(CX(J)).LE.REAL(CXR)) GO TO 3600
          CX(J+1)=CX(J)
3500    CONTINUE
        J=0
3600    CX(J+1)=CXR
4000  CONTINUE
      END

C     ---------------------------------------------------------
 
      FUNCTION FMOM(N)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE INTEGRALS OF W(X)X**N OVER THE REQUIRED INTERVAL

      FMOM=1.D0/(1.+N)**2
      END
