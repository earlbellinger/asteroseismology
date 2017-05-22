C     PROGRAM TO FIND COMPLEX ROOTS OF A NONLINEAR EQUATION
C     USING MULLER'S METHOD WITH DEFLATION TO REMOVE KNOWN ZEROS

      PROGRAM ZERO
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      EXTERNAL CF
      DIMENSION CZERO(20),CX(20)

C     EXAMPLE 7.6:  COMPLEX ZEROS OF Z+SIN(Z)

51    FORMAT('   IER =',I4,5X,I3,' ZEROS CALCULATED')
52    FORMAT(2X,1P2D14.6,2X,2D14.6)
53    FORMAT('  STARTING VALUES'/(2X,1P2D14.6,2X,2D14.6))

      REPS=1.D-7
      AEPS=1.D-7
      RMX=1.D2
      NZ=0

C     THE ZEROS WILL KEEP ACCUMULATING IN THE ARRAY CZERO

100   PRINT *,'TYPE NUM=NO. OF ZEROS TO BE TRIED  (QUITS WHEN NUM.EQ.0)' 
      READ *,NUM
      IF(NUM.EQ.0) STOP
      PRINT *,'TYPE NUM COMPLEX STARTING VALUES'
      READ *,(CX(I),I=1,NUM)
      WRITE(6,53) (CX(I),I=1,NUM)
      CALL ZROOT(NUM,CX,CZERO,NZ,REPS,AEPS,IER,RMX,CF)
      WRITE(6,51) IER,NZ
      WRITE(6,52) (CZERO(I),I=1,NZ)
      GO TO 100
      END

C     -------------------------------------------------------

C	Complex zero of a given function using Muller iteration with deflation
C
C	CX1,CX2,CX3 : (input/output) Complex starting values for iteration
C		These will be updated during execution and CX3 should be
C		the best estimate for the zero.
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(CX3))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=42 implies that roundoff errors appear to be dominating
C			and calculations are terminated.
C		IER=43 implies that iteration failed to converge to specified accuracy
C			but a weaker convergence criterion is satisfied
C		IER=404 implies 2 of the 3 starting values are identical and
C			no calculations are done
C		IER=431 implies that iteration goes outside the specified limits
C		IER=432 implies that iteration failed to converge to specified accuracy
C		IER=433 implies that denominator for calculating the iteration
C			vanishes and it is not possible to continue further
C	CF : (input) Name of the function routine to calculate the function value
C	NZ : (input) Number of zeros already known (for deflation)
C	CZERO : (input) Complex array of length NZ containing the known
C		zeros for deflation.
C	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
C		when ABS(CX) > RMAX
C	
C	FUNCTION CF(CX) must be supplied by the user. Here CF and CX are
C		both complex variables.
C
C	Required routines : CF

      SUBROUTINE MULLER(CX1,CX2,CX3,REPS,AEPS,IER,CF,NZ,CZERO,RMAX)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      PARAMETER(REPS0=1.D-4,NIT=50)
      DIMENSION CZERO(NZ)

      IF(CX2.EQ.CX1.OR.CX2.EQ.CX3.OR.CX3.EQ.CX1) THEN
C	If two of the three starting values are equal then quit
        IER=404
        RETURN
      ENDIF

      CF1=CF(CX1)
      CF2=CF(CX2)
C	Perform deflation
      DO 2000 J=1,NZ
        CF1=CF1/(CX1-CZERO(J))
2000  CF2=CF2/(CX2-CZERO(J))

C	The divided difference f[CX1,CX2]
      CH1=(CF2-CF1)/(CX2-CX1)
      IER=0
      DX1=ABS(CX3-CX2)

C	Loop for Muller iteration
      DO 5000 I=1,NIT
        CF3=CF(CX3)
C	Perform deflation
        DO 3000 J=1,NZ
3000    CF3=CF3/(CX3-CZERO(J))

        IF(CX3.EQ.CX1.OR.CX3.EQ.CX2) RETURN
        IF(CF3.EQ.0.0) RETURN
        CH2=(CF3-CF2)/(CX3-CX2)
        C2=(CH2-CH1)/(CX3-CX1)
        CG=CH2+(CX3-CX2)*C2
        CI=SQRT(CG*CG-4.*CF3*C2)
        CD=CG+CI
        CD1=CG-CI
        IF(ABS(CD1).GT.ABS(CD)) CD=CD1

        IF(CD.EQ.0.0) THEN
C	If denominator is zero, then quit
          IER=433
          RETURN
        ENDIF

        CDX=-2.*CF3/CD
        CX1=CX2
        CX2=CX3
C	The new approximation to zero
        CX3=CX3+CDX
        CF1=CF2
        CF2=CF3
        CH1=CH2

        DX=ABS(CDX)
        IF(I.GT.2.AND.DX.LT.MAX(REPS*ABS(CX3),AEPS)) RETURN
        IF(I.GT.5.AND.DX.LT.MAX(REPS0*ABS(CX3),AEPS)
     1      .AND.DX.GT.DX1) THEN
C	Roundoff errors appear to be dominating, hence quit
          IER=42
          RETURN
        ENDIF

        DX1=DX
        IF(ABS(CX3).GT.RMAX) THEN
C	Iteration goes outside the specified range
          IER=431
          RETURN
        ENDIF
5000  CONTINUE

C	Iteration fails to converge.
      IER=43
      IF(DX.LT.MAX(REPS0*ABS(CX3),AEPS)) RETURN
      IER=432
      END

C     ---------------------------------------------

C	Complex zeros of a given function using Muller iteration with deflation
C
C	N : (input) Number of zeros to be determined
C	CX : (input) Complex array of length N containing the initial guesses
C		for the zeros
C	CZERO : (input/output) Complex array of length N+NZ containing the
C		computed values of the zeros. The first NZ zeros which are
C		already known must be supplied as input while other zeros will be added
C	NZ : (input/output) Number of known zeros. At input it should
C		contain the number already known which are included in
C		array CZERO. This number will be incremented as more zeros
C		are determined successfully.
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=41 implies that iteration did not converge to satisfactory accuracy for
C			at least one zero, but value is acceptable at lower accuracy.
C		IER=429 implies that iteration failed to find at least one zero.
C			The number of zeros found can be checked from the value of NZ.
C	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
C		when ABS(CZ) > RMAX
C	CF : (input) Name of the function routine to calculate the function
C		FUNCTION CF(CZ) must be supplied by the user, where CF and CZ
C		are complex. For use with MULER2 the function is calculated
C		in the form CF*2**JF and FUNCTION CF(CZ,JF) should be provided.
C
C	Required routines : MULLER (or MULER2), CF

      SUBROUTINE ZROOT(N,CX,CZERO,NZ,REPS,AEPS,IER,RMAX,CF)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      EXTERNAL CF
      DIMENSION CZERO(NZ+N),CX(N)

      IER=0
      IF(N.LE.0) RETURN
      IF(NZ.LT.0) NZ=0

      DO 1000 I=1,N
C	The starting values for Muller's iteration
        CX3=CX(I)
        CDX=0.01D0*CX3
        IF(I.LT.N.AND.CDX.EQ.0.0) CDX=0.1D0*(CX(I+1)-CX(I))
        IF(I.GT.1.AND.CDX.EQ.0.0) CDX=0.1D0*(CX(I-1)-CX(I))
        IF(CDX.EQ.0.0) CDX=1.D-3
C	These values may need to be changed in some cases
        CX2=CX3+CDX
        CX1=CX3-CDX

C	Find the next zero
        CALL MULLER(CX1,CX2,CX3,REPS,AEPS,IER1,CF,NZ,CZERO,RMAX)
C	For MULER2 use the following statements instead
C         IER1=0
C 500     CALL MULER2(CX1,CX2,CX3,REPS,AEPS,IER1,CF0,CX,IX,NZ,CZERO,RMAX)
C         IF(IER1.LT.0) THEN
C           CF0=CF(CX,IX)
C           GO TO 500
C         ENDIF

        IF(IER1.LT.100) THEN
C	The zero is accepted
          NZ=NZ+1
          CZERO(NZ)=CX3
          IF(IER.EQ.0.AND.IER1.NE.0) IER=41
        ELSE
          IER=429
        ENDIF
1000  CONTINUE
      END

C     ----------------------------------------------

      FUNCTION CF(Z)
      IMPLICIT COMPLEX*16(C,Z)

C	The required function for finding zeros
      CF=Z+SIN(Z)
      END
