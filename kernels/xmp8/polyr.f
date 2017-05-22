C     PROGRAM TO FIND ALL ROOTS OF A POLYNOMIAL WITH REAL COEFFICIENTS
C     USING LAGUERRE'S METHOD

      PROGRAM POLY
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      DIMENSION CZERO(51),XCOF(51),WK(51)

51    FORMAT('    THE COEFFICIENTS OF POLYNOMIAL ARE :'/(1P5D15.7))
55    FORMAT('  IER =',I4,4X,'ZEROS =',1P2D14.6,2X,2D14.6/
     1       (2X,1P2D14.6,2X,2D14.6))

100   PRINT *,'TYPE NDEG=DEGREE OF POLYNOMIAL  (QUITS WHEN NDEG.LT.0)'
      READ *,NDEG
      IF(NDEG.LE.0) STOP
      PRINT *,'TYPE THE COEFFICIENTS STARTING FROM HIGHEST DEGREE TERM'
      READ *,(XCOF(I),I=NDEG+1,1,-1)

      WRITE(6,51) (XCOF(I),I=NDEG+1,1,-1)
      QREFIN=.TRUE.
      CALL POLYR(NDEG,XCOF,CZERO,IER,QREFIN,WK)
      WRITE(6,55) IER,(CZERO(I),I=1,NDEG)

      GO TO 100
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
