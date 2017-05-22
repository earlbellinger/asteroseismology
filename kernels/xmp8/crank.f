C     PROGRAM TO SOLVE LINEAR SECOND ORDER PARABOLIC EQUATIONS USING
C     CRANK-NICOLSON METHOD

      PROGRAM DIFFUS
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(100),U(100),WK(800)
      EXTERNAL COF,BC,FIC

C     EXAMPLE 13.1 : DIFFUSION EQUATION

51    FORMAT('   IER =',I4,5X,'T =',1PD14.6,5X,'TIME STEP =',D14.6,5X,
     1       'NX =',I4/5X,'SOLUTION =',4D14.6/(2X,5D14.6))
52    FORMAT('   EXACT SOLUTION =',1P4D14.6/(2X,5D14.6))

      X0=0.0
      XN=0.5D0
      IFLG=0
      PRINT *,'TYPE T=INITIAL TIME,  NX=NO. OF PTS'
      READ *,T,NX

100   PRINT *,'TYPE DT=TIME STEP,  NT=NO. OF TIME STEPS'
      PRINT *,'           (QUITS WHEN NT.LE.0)'
      READ *,DT,NT
      IF(NT.LE.0) STOP
      CALL CRANK(T,DT,X0,XN,NT,NX,X,U,COF,BC,FIC,IER,IFLG,WK)
      WRITE(6,51) IER,T,DT,NX,(U(I),I=1,NX)
      WRITE(6,52) (FIC(X(I),T),I=1,NX)
      GO TO 100
      END

C     ----------------------------------------------

C	To solve linear parabolic differential equation using the
C	Crank-Nicolson difference scheme
C	The differential equation is assumed to be of the form
C
C	du/dt = A(x,t)d^2u/dx^2 + B(x,t)du/dx + C(x,t)u + D(x,t)
C
C	with boundary conditions
C
C	A0(t) u(x0,t)+B0(t) du(x0,t)/dx = F0(t)
C	AN(t) u(xN,t)+BN(t) du(xN,t)/dx = FN(t)
C
C	T : (input/output) Initial value of "time" where the initial conditions
C		are specified. After execution, it will be replaced by the
C		value of T at the last point where execution is successful
C	DT : (input) The time step to be used for computations. It is kept fixed.
C	X0 : (input) Lower limit on X where the solution is required
C	XN : (input) Upper limit on X where the solution is required.
C		Solution is computed in the interval (X0,XN)
C	NT : (input) Number of time steps each of length DT  to be executed.
C	NX : (input) Number of mesh points in the X direction.
C	X : (output) Real array of length NX containing the mesh points used
C		in X direction. These are calculated by the routine by
C		assuming uniform spacing.
C	U : (input/output) Real array of length NX containing the solution
C		at the current time step. It should contain the initial
C		values at the time of calling, if IFLG.NE.0. Otherwise it
C		is computed using Function FIC. After execution it will contain
C		the computed solution at t=T.
C		It may be noted that solution at intermediate time steps
C		is not preserved. Hence, if solution is required at number
C		of time steps, multiple calls to CRANK will be needed.
C	COF : (input) Name of the subroutine to calculate the coefficients
C		in the equation
C	BC : (input) Name of the subroutine to calculate the coefficients
C		in the boundary conditions
C	FIC : (input) Name of the function routine to calculate the initial
C		values when IFLG=0. For other values of IFLG this routine
C		is not used, but a dummy routine may be required by the compiler.
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=713 implies that DT=0, XN=X0 or NX<3
C			in which case no calculations are done
C		IER=761 implies that the difference equations are singular
C			and solution cannot be continued further.
C	IFLG : (input/output) Integer variable used as a flag to denote how the
C			initial values are calculated.
C		If IFLG=0 then the initial values are calculated using the
C			function FIC to be supplied by the user. IFLG is set to 1
C		Otherwise the initial values must be supplied in array U.
C	WK : Real array of length 8*NX used as scratch space
C
C	SUBROUTINE COF(X,T,A,B,C,D), SUBROUTINE BC(T,X0,XN,A0,B0,F0,AN,BN,FN)
C	and FUNCTION FIC(X,T) must be supplied by the user 
C	Subroutine COF should calculate the coefficients A, B, C, D as
C	defined above for given values of X,T.
C	Subroutine BC should calculate the coefficients A0, B0, F0, AN, BN, FN
C	for given values of T, X0 and XN.
C	Function FIC is required only if IFLG=0, otherwise a dummy routine
C	with this name will suffice. If IFLG=0, function FIC must calculate
C	the initial values at required X,T.
C
C	Required routines : COF, BC, FIC

      SUBROUTINE CRANK(T,DT,X0,XN,NT,NX,X,U,COF,BC,FIC,IER,IFLG,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NX),U(NX),WK(NX,8)

      IF(DT.EQ.0.OR.XN.EQ.X0.OR.NX.LE.2) THEN
        IER=713
        RETURN
      ENDIF
      IER=761

C	Set up the step size in X as well as the initial values
      DX=(XN-X0)/(NX-1)
      DO 1000 I=1,NX
        X(I)=X0+(I-1)*DX
        CALL COF(X(I),T,WK(I,5),WK(I,6),WK(I,7),WK(I,8))
        IF(IFLG.EQ.0) U(I)=FIC(X(I),T)
1000  CONTINUE
      IFLG=1
      R=0.5*DT/DX**2
      R1=0.25*DT/DX

C	Boundary condition at X0
      CALL BC(T,X0,XN,A0,B0,F0,AN,BN,FN)
      JL=2
      IF(B0.NE.0.0) THEN
        JL=1
        U0=U(2)+2.*DX*(A0*U(1)-F0)/B0
      ENDIF

C	boundary condition at XN
      JU=NX-1
      IF(BN.NE.0.0) THEN
        JU=NX
        UN=U(NX-1)-2.*DX*(AN*U(NX)-FN)/BN
      ENDIF

C	Loop for time steps
      DO 5000 IT=1,NT
        T1=T+DT
        CALL BC(T1,X0,XN,A0,B0,F0,AN,BN,FN)

        DO 2000 J=JL,JU
          IF(J.EQ.1) THEN
            UM=U0
          ELSE
            UM=U(J-1)
          ENDIF
          IF(J.EQ.NX) THEN
            UP=UN
          ELSE
            UP=U(J+1)
          ENDIF

C	The Crank-Nicolson difference scheme
          CALL COF(X(J),T1,A,B,C,D)
          WK(J,1)=-R*A+R1*B
          WK(J,2)=1+2.*R*A-0.5*DT*C
          WK(J,3)=-R*A-R1*B
          WK(J,4)=U(J)+R*WK(J,5)*(UP-2*U(J)+UM)+
     1            R1*WK(J,6)*(UP-UM)+0.5*DT*(WK(J,7)*U(J)+WK(J,8)+D)
          WK(J,5)=A
          WK(J,6)=B
          WK(J,7)=C
          WK(J,8)=D
2000    CONTINUE

C	Boundary condition at X0
        IF(JL.EQ.2) THEN
          IF(A0.EQ.0) RETURN
          U(1)=F0/A0
          WK(2,4)=WK(2,4)-U(1)*WK(2,1)
        ELSE
          WK(1,2)=WK(1,2)+2.*WK(1,1)*A0*DX/B0
          WK(1,3)=WK(1,3)+WK(1,1)
          WK(1,4)=WK(1,4)+2.*WK(1,1)*F0*DX/B0
        ENDIF

C	Boundary condition at XN
        IF(JU.EQ.NX) THEN
          WK(JU,1)=WK(JU,1)+WK(JU,3)
          WK(JU,2)=WK(JU,2)-2.*WK(JU,3)*AN*DX/BN
          WK(JU,4)=WK(JU,4)-2.*WK(JU,3)*FN*DX/BN
        ELSE
          IF(AN.EQ.0.0) RETURN
          U(NX)=FN/AN
          WK(NX-1,4)=WK(NX-1,4)-U(NX)*WK(NX-1,3)
        ENDIF

C	Gaussian elimination for tridiagonal matrix
        DO 2500 J=JL+1,JU
          IF(WK(J-1,2).EQ.0) RETURN
          RP=-WK(J,1)/WK(J-1,2)
          WK(J,2)=WK(J,2)+RP*WK(J-1,3)
          WK(J,4)=WK(J,4)+RP*WK(J-1,4)
2500    CONTINUE
        IF(WK(JU,2).EQ.0) RETURN

C	Back-substitution for tridiagonal system
        U(JU)=WK(JU,4)/WK(JU,2)
        DO 3000 J=JU-1,JL,-1
          U(J)=(WK(J,4)-WK(J,3)*U(J+1))/WK(J,2)
3000    CONTINUE

        IF(JL.EQ.1) U0=U(2)+2.*DX*(A0*U(1)-F0)/B0
        IF(JU.EQ.NX) UN=U(NX-1)-2.*DX*(AN*U(NX)-FN)/BN
        T=T1
5000  CONTINUE
      IER=0
      END

C     --------------------------------------------------

      SUBROUTINE COF(X,T,A,B,C,D)
      IMPLICIT REAL*8(A-H,O-Z)

C     COEF. OF DIFFUSION EQ.

      A=1
      B=0.0
      C=0.0
      D=0.0
      END

C     ---------------------------------------------

      SUBROUTINE BC(T,X0,XN,A0,B0,F0,AN,BN,FN)
      IMPLICIT REAL*8(A-H,O-Z)

C     BOUNDARY CONDITIONS  U(X0,T)=0,  DU/DX(XN,T)=0

      A0=1.0
      B0=0.0
      F0=0.0

      AN=0.0
      BN=1.0
      FN=0.0
      END

C     ---------------------------------------

      FUNCTION FIC(X,T)
      IMPLICIT REAL*8(A-H,O-Z)

C     EXACT SOLUTION USED TO GENERATE INITIAL CONDITIONS

      PARAMETER(PI=3.1415926535D0)

      FIC=SIN(PI*X)*EXP(-PI*PI*T)
      END
