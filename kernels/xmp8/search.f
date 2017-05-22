C     PROGRAM TO LOCATE COMPLEX ROOTS OF A NONLINEAR EQUATION BY LOOKING
C     FOR SIGN CHANGES IN THE COMPLEX PLANE

      PROGRAM CROOT
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL CF

C     EXAMPLE 7.6

100   PRINT *,'TYPE RX1,RX2=RANGE ON REAL AXIS'
      PRINT *,'      (QUITS WHEN RX1.EQ.RX2)'
      READ *,RX1,RX2
      IF(RX1.EQ.RX2) STOP
      PRINT *,'TYPE RY1,RY2=RANGE ON IMAGINARY AXIS'
      READ *,RY1,RY2
      PRINT *,'TYPE NX,NY=NO. OF POINTS ALONG X & Y AXES'
      READ *,NX,NY
      CALL SEARCH(RX1,RX2,RY1,RY2,NX,NY,CF)
      GO TO 100
      END

C     -----------------------------------------------------

C	To search for complex zeros in a rectangular region in complex plane
C
C	RX1 : (input) Lower limit of real part for the region to be searched
C	RX2 : (input) Upper limit of real part for the region to be searched
C	RY1 : (input) Lower limit of imaginary part for the region to be searched
C	RY2 : (input) Upper limit of imaginary part for the region to be searched
C		The region with real part between RX1 & RX2 and 
C		imaginary part between RY1 & RY2 will be searched
C	NX : (input/output) number of points to be used along real axis
C	NY : (input/output) number of points to be used along imaginary axis
C	CFUN : (input) Name of the function routine to calculate the complex
C		function. Function CFUN(CZ) must be supplied by the user.
C		CZ and CFUN are both complex variables.
C	The zero will be located near the point where all 4 quadrants
C	in function values meet in the figure that is produced.
C
C	Required routines : CFUN

      SUBROUTINE SEARCH(RX1,RX2,RY1,RY2,NX,NY,CFUN)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16(C)
      EXTERNAL CFUN
      PARAMETER(IMAX=41)
      DIMENSION IQ(IMAX),X(IMAX)

51    FORMAT(/3X,1P18D15.3)
52    FORMAT(1PD12.3,2X,81I3)

      CI=(0.0D0,1.D0)
      IF(RX1.EQ.RX2.OR.RY1.EQ.RY2) RETURN
      IF(NX.GT.IMAX) NX=IMAX
      IF(NX.LE.1) NX=21
      IF(NY.LE.1) NY=21

      HX=(RX2-RX1)/(NX-1)
      HY=(RY2-RY1)/(NY-1)
      DO 1000 I=1,NX
1000  X(I)=RX1+HX*(I-1)

      DO 2000 I=1,NY
        Y=RY2-HY*(I-1)
        DO 1500 J=1,NX
          CZ=X(J)+CI*Y
          CF=CFUN(CZ)

C	Determine the quadrant in which the function value lies
          RX=CF
          RY=IMAG(CF)
          IF(RX.GE.0.0.AND.RY.GE.0.0) IQ(J)=1
          IF(RX.LT.0.0.AND.RY.GE.0.0) IQ(J)=2
          IF(RX.LT.0.0.AND.RY.LT.0.0) IQ(J)=3
          IF(RX.GE.0.0.AND.RY.LT.0.0) IQ(J)=4

1500    CONTINUE
        WRITE(6,52) Y,(IQ(J),J=1,NX)
2000  CONTINUE

      WRITE(6,51)(X(I),I=1,NX,5)
      END

C     ------------------------------------------

      FUNCTION CF(Z)
      COMPLEX*16 Z,CF

C	The required function

      CF=Z+SIN(Z)
      END
