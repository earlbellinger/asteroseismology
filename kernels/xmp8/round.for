C     PROGRAM FOR ROUNDING A GIVEN FLOATING POINT NUMBER TO SPECIFIED
C     NO. OF DIGITS USING THE SPECIFIED BASE. ALL NUMBERS ARE PRINTED
C     OUT IN DECIMAL SYSTEM AND HENCE IT MAY BE DIFFICULT TO VERIFY THE
C     RESULT UNLESS THE BASE B=10.

      PROGRAM FROUND
      IMPLICIT REAL*8(A-H,O-Z)

51    FORMAT(1PD15.6,'  ROUNDED TO',D15.6)

100   PRINT *,'TYPE X=NO. TO BE ROUNDED,   N=NO. OF DIGITS,   B=BASE'
      PRINT *,'        (QUITS WHEN N.LE.0)'
      READ *,X,N,B
      IF(N.LE.0) STOP
      XR=ROUND(X,N,B)
      WRITE(6,51) X,XR
      GO TO 100
      END
 
C     -----------------------------------------------------
 
C	To round a number to N digits using base B.
C
C	X : (input) The number to be rounded
C	N : (input) The number of digits required in rounded number
C	B : (input) Base of number system to be used for rounding
C
C	Required routines : None

      FUNCTION ROUND(X,N,B)
      IMPLICIT REAL*8(A-H,O-Z)

      ROUND=0.0
      IF(X.EQ.0.0) RETURN
      XA=ABS(X)
      LGX=LOG(XA)/LOG(B)
      IF(XA.LT.1.0) LGX=LGX-1
      FX=B**(N-LGX-1)
      NX=XA*FX+0.5
      XA=NX/FX
      ROUND=SIGN(XA,X)
      END
