C     PROGRAM TO SOLVE SECOND ORDER LINEAR ELLIPTIC EQUATIONS IN TWO
C     DIMENSIONS USING THE SUCCESSIVE OVER-RELAXATION METHOD

      PROGRAM ELIPTC
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(60),Y(60),U(51,51),WK(35000)
      EXTERNAL COF,BC

C     EXAMPLE 13.4

51    FORMAT('   IER =',I4,5X,'NX =',I4,5X,'NY =',I4/5X,
     1    'NO. OF ITERATIONS =',I6,5X,'RELAXATION PARAMETER =',1PD14.6)
52    FORMAT('  Y =',1PD14.6,2X,4D14.6/(21X,4D14.6))

      IU=51
      X0=0
      XN=1.0
      Y0=0
      YN=1.0
      AEPS=1.D-6

100   PRINT *,'TYPE NX,NY=NO. OF PTS ALONG X,Y AXES,'
     1       ,'    OMEGA=RELAXATION PARA'
      PRINT *,'      (QUITS WHEN NX.LE.0.OR.NY.LE.0)'
      READ *,NX,NY,OMEGA
      IF(NX.LE.0.OR.NY.LE.0) STOP

C     INITIAL GUESS FOR SOLUTION

      DO 111 I=1,NX
      DO 111 J=1,NY
        U(I,J)=0.0
111   CONTINUE

      CALL SOR(X0,XN,Y0,YN,NX,NY,X,Y,U,IU,COF,BC,OMEGA,IER,
     1         AEPS,NIT,WK)
      WRITE(6,51) IER,NX,NY,NIT,OMEGA
      DO 1000 I=1,NY
1000  WRITE(6,52) Y(I),(U(J,I),J=1,NX)
      GO TO 100

      END

C     ----------------------------------------------

C	To solve linear second order elliptic differential equation using the
C	successive over-relaxation (SOR) method
C	The differential equation is assumed to be of the form
C
C	Axx(x,y)d^2u/dx^2 + Axy(x,y)d^2u/dydx + Ayy(x,y)d^2u/dy^2 +
C		Ax(x,y)du/dx + Ay(x,y)du/dy + A0(x,y)u + F(x,y)=0
C
C	with Dirichlet boundary conditions on a rectangular region
C
C	X0 : (input) Lower limit on X where the solution is required
C	XN : (input) Upper limit on X where the solution is required.
C		Solution is computed in the interval (X0,XN)
C	Y0 : (input) Lower limit on Y where the solution is required
C	YN : (input) Upper limit on Y where the solution is required.
C		Solution is computed in the interval (Y0,YN)
C	NX : (input) Number of mesh points in the X direction.
C	NY : (input) Number of mesh points in the Y direction.
C	X : (output) Real array of length NX containing the mesh points used
C		in X direction. These are calculated by the routine by
C		assuming uniform spacing.
C	Y : (output) Real array of length NY containing the mesh points used
C		in Y direction. These are calculated by the routine by
C		assuming uniform spacing.
C	U : (input/output) Real array of length IU*NY containing the solution
C		It should contain the initial values at the time of calling.
C		After execution it will contain	the computed solution.
C		U(I,J) is the solution at (x_i,y_j)
C	IU : (input) The first dimension of U as declared in the calling
C		program, IU.GE.NX
C	COF : (input) Name of the subroutine to calculate the coefficients
C		in the equation
C	BC : (input) Name of the subroutine to calculate the boundary conditions
C	OMEGA : (input/output) Value of the relaxation parameter, 1<OMEGA<2
C		If OMEGA .LE. 0 then the routine sets it to the
C		optimal value for Poisson's equation.
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=716 implies that YN=Y0, XN=X0, NX<3, NY<3 or IU<NX
C			in which case no calculations are done
C		IER=764 implies that the diagonal term in the difference
C			equation vanishes and calculations have to be abandoned.
C		IER=765 implies that SOR iteration failed to converge to
C			specified accuracy
C	AEPS : (input) Required absolute accuracy. The SOR iteration is
C			continued until the change in all elements is less than AEPS
C	NIT : (output) Number of SOR iterations required by the subroutine.
C	WK : Real array of length 9*NX*NY used as scratch space
C
C	SUBROUTINE COF(X,Y,AXX,AXY,AYY,AX,AY,A0,F) and FUNCTION BC(IB,X,Y)
C	must be supplied by the user 
C	Subroutine COF should calculate the coefficients AXX, AXY, AYY, 
C	AX, AY, A0, F as defined above for given values of X,Y.
C	Subroutine BC should calculate the Boundary values at each boundary.
C	Here IB is an integer denoting which boundary is being considered.
C	The boundary conditions are assumed to be
C	u(X0,Y)=BC(1,X0,Y);	u(XN,Y)=BC(2,XN,Y);
C	u(x,Y0)=BC(3,x,Y0);	u(x,YN)=BC(4,x,YN)
C
C	Required routines : COF, BC

      SUBROUTINE SOR(X0,XN,Y0,YN,NX,NY,X,Y,U,IU,COF,BC,OMEGA,IER,
     1               AEPS,NIT,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.14159265358979324D0,MAXIT=1000)
      DIMENSION X(NX),Y(NY),U(IU,NY),WK(9,NX,NY)

      IF(XN.EQ.X0.OR.YN.EQ.Y0.OR.NX.LE.2.OR.NY.LE.2.OR.IU.LT.NX) THEN
        IER=716
        RETURN
      ENDIF
      IER=764

C	Setting up the mesh points
      DX=(XN-X0)/(NX-1)
      DO 1000 I=1,NX
1000  X(I)=X0+(I-1)*DX
      DY=(YN-Y0)/(NY-1)
      DO 1200 I=1,NY
1200  Y(I)=Y0+(I-1)*DY

      DX2=DX*DX
      DY2=DY*DY
      DXY=4.*DX*DY
      IF(OMEGA.LE.0.0) THEN
C	Estimate the optimal value of OMEGA
        RJ=(DY2*COS(PI/(NX-1.))+DX2*COS(PI/(NY-1.)))/(DX2+DY2)
        OMEGA=2./(1.+SQRT(1.-RJ**2))
      ENDIF

C	Calculate the coefficients of the difference equations
      DO 2000 I=1,NX
        DO 2000 J=1,NY
          CALL COF(X(I),Y(J),AXX,AXY,AYY,AX,AY,A0,F)
          AD=2.*AXX/DX2+2.*AYY/DY2-A0
          IF(AD.EQ.0.0) RETURN
          WK(1,I,J)=-AXY/(DXY*AD)
          WK(2,I,J)=(-AYY/DY2+0.5*AY/DY)/AD
          WK(3,I,J)=-WK(1,I,J)
          WK(4,I,J)=(-AXX/DX2+0.5*AX/DX)/AD
          WK(5,I,J)=(-AXX/DX2-0.5*AX/DX)/AD
          WK(6,I,J)=WK(3,I,J)
          WK(7,I,J)=(-AYY/DY2-0.5*AY/DY)/AD
          WK(8,I,J)=WK(1,I,J)
          WK(9,I,J)=F/AD
2000  CONTINUE

C	Calculate the boundary values
      DO 2200 K=1,NY
        U(1,K)=BC(1,X0,Y(K))
        U(NX,K)=BC(2,XN,Y(K))
2200  CONTINUE
      DO 2600 J=1,NX
        U(J,1)=BC(3,X(J),Y0)
        U(J,NY)=BC(4,X(J),YN)
2600  CONTINUE

C	Loop for the SOR iteration
      DO 6000 IT=1,MAXIT
        ERR=0.0
        DO 3500 K=2,NY-1
          DO 3500 J=2,NX-1
            RES=WK(9,J,K)-U(J,K)-WK(1,J,K)*U(J-1,K-1)-WK(2,J,K)*U(J,K-1)
     1      -WK(3,J,K)*U(J+1,K-1)-WK(4,J,K)*U(J-1,K)-WK(5,J,K)*U(J+1,K)-
     2      WK(6,J,K)*U(J-1,K+1)-WK(7,J,K)*U(J,K+1)-WK(8,J,K)*U(J+1,K+1)
            ERR=MAX(ERR,ABS(RES))
            U(J,K)=U(J,K)+OMEGA*RES
3500    CONTINUE

        IF(ERR.LT.AEPS) THEN
          IER=0
          NIT=IT
          RETURN
        ENDIF

6000  CONTINUE

C	Iteration fails to converge
      NIT=MAXIT
      IER=765
      END

C     --------------------------------------------------

      SUBROUTINE COF(X,Y,AXX,AXY,AYY,AX,AY,AU,A0)
      IMPLICIT REAL*8(A-H,O-Z)

C     THE COEFFICIENTS OF ELLIPTIC EQUATIONS

      AXX=1.0
      AXY=0.0
      AYY=1.0
      AX=0.0
      AY=0.0
      AU=0.0
      A0=0.0
      END

C     ---------------------------------------------

      FUNCTION BC(IB,X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.1415926535D0)

C     DIRICHLET BOUNDARY CONDITION  U(0,Y)=U(1,Y)=U(X,0)=0, U(X,1)=SIN(PI*X)

      IF(IB.LE.3) THEN
        BC=0
      ELSE
        BC=SIN(PI*X)
      ENDIF

      END
