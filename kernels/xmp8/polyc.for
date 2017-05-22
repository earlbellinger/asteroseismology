C     PROGRAM TO FIND ALL ROOTS OF A POLYNOMIAL WITH COMPLEX COEFFICIENTS
C     USING LAGUERRE'S METHOD

      PROGRAM POLY
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      DIMENSION CZERO(51),COF(51),CWK(51)

51    FORMAT('    THE COEFFICIENTS OF POLYNOMIAL ARE :'/
     1       (2X,1P2D15.7,2X,2D15.7))
55    FORMAT('  IER =',I4,4X,'ZEROS =',1P2D14.6,2X,2D14.6/
     1       (2X,1P2D14.6,2X,2D14.6))

100   PRINT *,'TYPE NDEG=DEGREE OF POLYNOMIAL  (QUITS WHEN NDEG.LT.0)'
      READ *,NDEG
      IF(NDEG.LE.0) STOP
      PRINT *,'TYPE THE COEFFICIENTS STARTING FROM HIGHEST DEGREE TERM'
      READ *,(COF(I),I=NDEG+1,1,-1)

      PRINT 51,(COF(I),I=NDEG+1,1,-1)
      QREFIN=.TRUE.
      CALL POLYC(NDEG,COF,CZERO,IER,QREFIN,CWK)
      PRINT 55,IER,(CZERO(I),I=1,NDEG)
      GO TO 100
      END

C  ------------------------------------------------------------

C	Root of a polynomial with complex coefficients using Laguerre iteration
C
C	N : (input) The degree of the polynomial
C	COF : (input) Complex array of length N+1 containing the coefficients of
C		the polynomial. COF(I) is the coefficient of X**(I-1)
C	CXI : (input/output) Complex variable containing the initial guess,
C		 after execution it will contain the computed root
C	IER : (output) Error parameter, IER=0 for successful execution
C		IER=438 implies that denominator is zero and iteration cannot
C			be continued further
C		IER=439 implies that iteration has failed to converge
C
C	Required routines : None

      SUBROUTINE LAGITC(N,COF,CXI,IER)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      PARAMETER(ITMAX=50,REPS=1.D-4,AEPS=1.D-6)
      DIMENSION COF(N+1)

      IER=0
      CDX1=ABS(CXI)+1
      QC=.FALSE.
      IC=ITMAX

C	The Laguerre's iteration
      DO 2000 I=1,ITMAX
        CF=COF(N+1)
        CFP=0.0
        CFPP=0.0
        DO 1000 J=N,1,-1
          CFPP=CXI*CFPP+2.*CFP
          CFP=CXI*CFP+CF
          CF=CXI*CF+COF(J)
1000    CONTINUE

        CH=(N-1)*((N-1)*CFP*CFP-N*CF*CFPP)
        CH=SQRT(CH)
        CDEN=CFP+CH
        IF(ABS(CFP-CH).GT.ABS(CDEN)) CDEN=CFP-CH

        IF(CDEN.NE.0.0) THEN
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
C	If the denominator vanishes, then quit
          IER=438
          RETURN
        ENDIF
2000  CONTINUE

C	Iteration fails to converge
      IF(.NOT.QC) IER=439
      END

C     ----------------------------------------

C	Roots of a polynomial with complex coefficients using Laguerre iteration
C
C	N : (input) The degree of the polynomial
C	COF : (input) Complex array of length N+1 containing the coefficients of
C		the polynomial. COF(I) is the coefficient of X**(I-1)
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
C	CWK : Complex array of length N+1 used as a scratch space
C
C	Required routines : LAGITC
C
      SUBROUTINE POLYC(N,COF,CX,IER,QREFIN,CWK)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      DIMENSION COF(N+1),CX(N),CWK(N+1)

      IF(N.LE.0) THEN
        IER=406
        RETURN
      ENDIF
      IF(COF(N+1).EQ.0.0) THEN
        IER=408
        RETURN
      ENDIF

      DO 1000 I=1,N+1
1000  CWK(I)=COF(I)
      NP=N
      CXR=0.0
      IER=0

C	Find the next root
2000  CALL LAGITC(NP,CWK,CXR,IER1)
      IF(IER1.NE.0) THEN
C	If iteration fails to converge, try once more
        CXR=1.123456
        CALL LAGITC(NP,CWK,CXR,IER1)
        IF(IER1.NE.0) THEN
C	If iteration fails again, then quit
          IER=430
          RETURN
        ENDIF
      ENDIF

      CXRT=CXR
      IF(NP.LT.N.AND.QREFIN) THEN
C	Refine the roots with original polynomial
        CALL LAGITC(N,COF,CXRT,IER1)
        IF(IER1.NE.0) THEN
C	If iteration fails to converge then retain the old value
          IER=IER+11
          CXRT=CXR
        ENDIF
      ENDIF

      IF(NP.GT.1) THEN
C	Perform deflation using unrefined root
        CN0=CWK(NP)
        CWK(NP)=CWK(NP+1)
        DO 2400 I=NP-1,1,-1
          CN=CXR*CWK(I+1)+CN0
          CN0=CWK(I)
          CWK(I)=CN
2400    CONTINUE
      ENDIF

      NP=NP-1
      CX(N-NP)=CXRT
C	If any more roots are left, find them
      IF(NP.GT.0) GO TO 2000
      END
