C     TO EVALUATE A POLYNOMIAL AND ITS DERIVATIVES

      PROGRAM POLY
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(100),PD(20)
 
51    FORMAT('   N =',I4,4X,'NO. OF DERIVATIVES =',I3,5X,'X =',1PD14.6,
     1       3X,'P(X) =',D14.6)
52    FORMAT('   DERIVATIVES : ',1P5D14.6/(5D14.6))
53    FORMAT('   COEFFICIENTS :',1P5D14.6/(5D14.6))

      PRINT *,'TYPE N = DEGREE OF POLYNOMIAL (N<100)'
      READ *,N
      PRINT *,'TYPE COEFFICIENTS OF POLYNOMIAL, STARTING FROM HIGHEST'
     1      ,' DEGREE TERM'
      READ *,(A(I),I=N+1,1,-1)
      WRITE(6,53) (A(I),I=N+1,1,-1)
 
100   PRINT *,'TYPE X, ND= NO. OF DERIVATIVES REQUIRED'
      PRINT *,'                 (QUITS WHEN X<-100)'
      READ *,X,ND
      IF(X.LT.-100) STOP
      PF=POLD(N,A,X,ND,PD)
      WRITE(6,51) N,ND,X,PF
      WRITE(6,52) (PD(I),I=1,ND)
      GO TO 100
      END
 
C     ---------------------------------------------------------
 
C	To evaluate a polynomial and its derivatives at any point
C
C	N : (input) Degree of the polynomial
C	A : (input) Real array of length N+1 containing the coefficients
C		of the polynomial. A(1) is the constant term and A(N+1) is
C		the coefficient of X**N
C	X : (input) The value of x at which the polynomial is to be evaluated
C	ND : (input) Number of derivatives required. The subroutine always
C		calculates the first derivative irrespective of value of ND
C	PD : (output) Real array of length PD containing the calculated
C		values of the derivatives. PD(I) will contain the Ith derivative
C	The polynomial value is returned through POLD
C
C	Required routines : None

      FUNCTION POLD(N,A,X,ND,PD)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N+1),PD(ND)
 
      POLD=A(N+1)
      DO 2000 I=1,ND
        PD(I)=0.0
2000  CONTINUE
 
      DO 2500 J=N,1,-1
 
        DO 2200 I=ND,2,-1
          PD(I)=PD(I)*X+I*PD(I-1)
2200    CONTINUE
        PD(1)=PD(1)*X+POLD
        POLD=POLD*X+A(J)
2500  CONTINUE
      END
