C     PROGRAM FOR DRAWING A SMOOTH CURVE THROUGH A TABLE OF POINTS

      PROGRAM SMDRAW
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(20),F(20),FP(500),C(3,20),XP(500)

C     EXERCISE 4.22
C     THE GIVEN TABLE OF FUNCTION VALUES
      DATA (X(I),I=1,9)/0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0/
      DATA (F(I),I=1,9)/4.579D0,6.543D0,9.209D0,12.788D0,17.535D0,
     1                  23.756D0,31.824D0,41.175D0,55.324D0/

51    FORMAT('    IER =',I4,5X,'NP =',I4)
52    FORMAT(1P2D14.6)

      NTAB=9
      PRINT *, 'TYPE NP = NO. OF POINTS REQUIRED'
      READ *,NP
      IF(NP.GT.500) STOP
      CALL SMOOTH(NTAB,X,F,C,NP,XP,FP,IER)
      WRITE(6,51) IER,NP

C    	The interpolated points are written out in file smooth.out
C    	This file can be used to plot the required smooth curve
      OPEN(UNIT=16,FILE='smooth.out',STATUS='UNKNOWN')
      WRITE(16,52)(XP(I),FP(I),I=1,NP)

      END

C	---------------------------------------------------

C	To draw a smooth curve passing through a set of data points
C	    using cubic spline interpolation
C
C	NTAB : (input) Number of points in the table
C	X : (input) Array of length NTAB containing X values
C	F : (input) Array of length NTAB containing function values at X(I)
C	C : (output) Array of length 3*NTAB which will contain the spline coefficients
C	NP : (input) Number of points at which interpolation is to be calculated
C	XP : (output) Real array of length NP containing the x values at
C	           NP uniformly spaced points for use in plotting
C	FP : (output) Real array of length NP containing interpolated
C		function values at XP(I)
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=202 implies NP.LE.1
C		other values may be set by SPLINE
C
C	Arrays XP and FP can be used to draw a smooth curve through the
C	tabulated points.
C
C	Required routines : SPLINE, SPLEVL

      SUBROUTINE SMOOTH(NTAB,X,F,C,NP,XP,FP,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NTAB),F(NTAB),C(3,NTAB),XP(NP),FP(NP)

      CALL SPLINE(X,F,NTAB,C,IER)
      IF(IER.GT.100) RETURN

      IF(NP.LE.1) THEN
        IER=202
        RETURN
      ENDIF

      DX=(X(NTAB)-X(1))/(NP-1)
      DO 1000 I=1,NP
        XP(I)=X(1)+DX*(I-1)
        FP(I)=SPLEVL(XP(I),NTAB,X,F,C,DFB,DDFB,IER)
1000  CONTINUE
      END

C	------------------------------------------------------------

C	To evaluate the cubic spline interpolant at a specified point
C
C	XB : (input) point at which interpolation is required
C	N : (input) Number of points in the table
C	X : (input) Array of length N, containing the abscissas
C	F : (input) Array of length N, containing the function values at X(I)
C	C : (input) Array of length 3*N containing the spline coefficients
C		which should have been calculated using SPLINE
C	DFB : (output) First derivative of spline at x=XB
C	DDFB : (output) Second derivative of spline at x=XB
C	IER : (output) error parameter, IER=0 if execution is successful
C		IER=24 implies XB is outside the range of table on higher side
C		IER=25 implies XB is outside the range of table on lower side
C		IER=201 implies N<2
C	SPLEVL will be the interpolated value at x=XB
C
C	Required routines : None

      FUNCTION SPLEVL(XB,N,X,F,C,DFB,DDFB,IER)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION X(N),F(N),C(3,N)
      SAVE
      DATA LOW/0/

      SPLEVL=0.0
      IF(N.LE.1) THEN
        IER=201
        RETURN
      ENDIF

C	QASCND is true if table is in ascending order
      QASCND=X(N).GT.X(1)
      IER=0

      IF(LOW.LT.1.OR.LOW.GE.N) THEN
C	If the previous value of LOW is inadmissible, set the range to (1,N)
        LOW=1
        IGH=N
      ELSE
        IGH=LOW+1
      ENDIF

1000  IF((XB.LT.X(LOW).AND.XB.LT.X(IGH)).OR.
     1    (XB.GT.X(LOW).AND.XB.GT.X(IGH))) THEN
C	Extend the range
        IF(XB.GT.X(LOW).EQV.QASCND) THEN
C	Extend the range on higher side
          IF(IGH.GE.N) THEN
            IER=24
            LOW=N-1
          ELSE
            NIGH=MIN(N,IGH+2*(IGH-LOW))
            LOW=IGH
            IGH=NIGH
            GO TO 1000
          ENDIF
        ELSE
C	Extend the range on lower side
          IF(LOW.LE.1) THEN
            IER=25
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

C	Calculate the interpolant and its derivatives
      DX=XB-X(LOW)
      SPLEVL=((C(3,LOW)*DX+C(2,LOW))*DX+C(1,LOW))*DX+F(LOW)
      DFB=(3.*C(3,LOW)*DX+2.*C(2,LOW))*DX+C(1,LOW)
      DDFB=6.*C(3,LOW)*DX+2.*C(2,LOW)
      END

C	------------------------------------------------------------

C	To calculate coefficients of cubic spline interpolation with
C		not-a-knot boundary conditions
C
C	X : (input) Real array of length N containing x values
C	F : (input) Real array of length N containing values of function at X(I)
C		F(I) is the tabulated function value at X(I).
C	N : (input) Length of table X,F
C	C : (output) Real array of length 3*N containing the spline coefficients
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=201 implies that N<2
C
C	Required routines : None

      SUBROUTINE SPLINE(X,F,N,C,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N),C(3,N)

      IF(N.LE.1) THEN
        IER=201
        RETURN

      ELSE IF(N.EQ.2) THEN
C	Use linear interpolation
        C(1,1)=(F(2)-F(1))/(X(2)-X(1))
        C(2,1)=0.0
        C(3,1)=0.0
        RETURN

      ELSE IF(N.EQ.3) THEN
C	Use quadratic interpolation
        DIV12=(F(2)-F(1))/(X(2)-X(1))
        DIV23=(F(3)-F(2))/(X(3)-X(2))
        C(3,1)=0.0
        C(3,2)=0.0
        C(2,1)=(DIV23-DIV12)/(X(3)-X(1))
        C(2,2)=C(2,1)
        C(1,1)=DIV12+C(2,1)*(X(1)-X(2))
        C(1,2)=DIV23+C(2,1)*(X(2)-X(3))
        RETURN

      ELSE
C	Use cubic splines 

C	Setting up the coefficients of tridiagonal matrix
        C(3,N)=(F(N)-F(N-1))/(X(N)-X(N-1))
        DO 1000 I=N-1,2,-1
          C(3,I)=(F(I)-F(I-1))/(X(I)-X(I-1))
          C(2,I)=2.*(X(I+1)-X(I-1))
C	The right hand sides
1000    C(1,I)=3.*(C(3,I)*(X(I+1)-X(I))+C(3,I+1)*(X(I)-X(I-1)))

C	The not-a-knot boundary conditions
        C1=X(3)-X(1)
        C(2,1)=X(3)-X(2)
        C(1,1)=C(3,2)*C(2,1)*(2.*C1+X(2)-X(1))+C(3,3)*(X(2)-X(1))**2
        C(1,1)=C(1,1)/C1
        CN=X(N)-X(N-2)
        C(2,N)=X(N-1)-X(N-2)
        C(1,N)=C(3,N)*C(2,N)*(2.*CN+X(N)-X(N-1))
        C(1,N)=(C(1,N)+C(3,N-1)*(X(N)-X(N-1))**2)/CN

C	Solving the equation by Gaussian elimination
        G=(X(3)-X(2))/C(2,1)
        C(2,2)=C(2,2)-G*C1
        C(1,2)=C(1,2)-G*C(1,1)
        DO 2000 J=2,N-2
          G=(X(J+2)-X(J+1))/C(2,J)
          C(2,J+1)=C(2,J+1)-G*(X(J)-X(J-1))
2000    C(1,J+1)=C(1,J+1)-G*C(1,J)
        G=CN/C(2,N-1)
        C(2,N)=C(2,N)-G*(X(N-1)-X(N-2))
        C(1,N)=C(1,N)-G*C(1,N-1)

C	The back-substitution
        C(1,N)=C(1,N)/C(2,N)
        DO 3000 I=N-1,2,-1
3000    C(1,I)=(C(1,I)-C(1,I+1)*(X(I)-X(I-1)))/C(2,I)
        C(1,1)=(C(1,1)-C(1,2)*C1)/C(2,1)

C	Calculating the coefficients of cubic spline
        DO 4000 I=1,N-1
          C(2,I)=(3.*C(3,I+1)-2.*C(1,I)-C(1,I+1))/(X(I+1)-X(I))
          C(3,I)=(C(1,I)+C(1,I+1)-2.*C(3,I+1))/(X(I+1)-X(I))**2
4000    CONTINUE
C	Set the coefficients for interval beyond X(N) using continuity
C	of second derivative, although they may not be used.
        C(2,N)=C(2,N-1)+3*(X(N)-X(N-1))*C(3,N-1)
        C(3,N)=0.0
      ENDIF
      END
