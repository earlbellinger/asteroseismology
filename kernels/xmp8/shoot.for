C     PROGRAM TO SOLVE A BOUNDARY VALUE PROBLEM IN ORDINARY DIFFERENTIAL
C     EQUATION USING SHOOTING METHOD

      PROGRAM SHOOT
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FS,DIF
      DIMENSION X(20),DX(20),WK(100)

C     EXAMPLE 11.11

51    FORMAT('   IER =',I4,5X,'INITIAL SLOPE =',1PD14.6)
52    FORMAT('   IER =',I4,5X,'T =',1PD14.6,5X,'Y(T) =',D14.6)
53	  FORMAT('   LOWER LIMIT ON X =',1PD14.6,5X,'UPPER LIMIT =',D14.6)

      REPS=1.D-7
      AEPS=1.D-8

C     THE INITIAL VALUE X FOR F'(0) IS DETERMINED BY SOLVING A NONLINEAR
C     EQUATION DEFINED BY THE SECOND BOUNDARY CONDITION
C     SUPPLY THE EXPECTED LIMITS ON THIS VALUE

100   PRINT *,'TYPE XL=LOWER LIMIT ON X,  XU=UPPER LIMIT ON X'
      PRINT *,'           (QUITS WHEN XL.GE.XU)'
      READ *,XL,XU
      IF(XL.GE.XU) STOP
      WRITE(6,53) XL,XU

      CALL BRENT(XL,XU,X1,REPS,AEPS,IER,FS)
      WRITE(6,51) IER,X1
      IF(IER.GT.0) GO TO 100

C     ONCE THE INITIAL VALUES ARE KNOWN GENERATE THE SOLUTION AT
C     REQUIRED POINTS

      T0=0
      X(1)=0
      X(2)=X1
      N=2
      H=0.01D0
      NMAX=10000
      DO 1000 I=1,5
        TN=I*0.2D0
        CALL RKM(N,X,DX,DIF,H,T0,TN,REPS,NSTEP,NMAX,IER,WK)
        WRITE(6,52) IER,T0,X(1)
1000  CONTINUE
      GO TO 100
      END

C     -------------------------------------------

C	Real zero of a given function using Brent's method
C
C	A : (input/output) Lower limit of interval where zero is located
C	B : (input/output) Upper limit of interval where zero is located
C		The limits A, B are updated by the subroutine
C	X : (output) Computed value of the zero
C	REPS : (input) Required relative accuracy
C	AEPS : (input) Required absolute accuracy
C		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=427 implies that function has same sign at x=A and B
C			in which case no calculations are done
C		IER=428 implies that iteration failed to converge to specified
C			accuracy
C	F : (input) Name of the function routine to calculate the function
C		FUNCTION F(X) must be supplied by the user.
C
C	Required routines : F

      SUBROUTINE BRENT(A,B,X,REPS,AEPS,IER,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=75)

      IER=0
      FA=F(A)
      FB=F(B)
      IF(FA.GT.0.0.EQV.FB.GT.0.0) THEN
C	If the function has the same sign at both the end points, then quit
        IER=427
        RETURN
      ENDIF

      C=A
      FC=FA
      D=B-A
      E=D

C	The iteration loop
      DO 2000 I=1,NIT
        IF(ABS(FC).LT.ABS(FB)) THEN
C	If F(C) < F(B), then interchange B and C
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF

        EPS=MAX(REPS*ABS(B),AEPS)
        DX=0.5*(C-B)
        IF(ABS(DX).LT.EPS.OR.FB.EQ.0) THEN
C	Iteration has converged
          IER=0
          X=B
          RETURN
        ENDIF

        IF(ABS(E).LT.EPS.OR.ABS(FA).LE.ABS(FB)) THEN
C	If the previous step is too small or if F(A).LE.F(B), perform bisection
          D=DX
          E=DX
        ELSE
          R3=FB/FA
          IF(A.EQ.C) THEN
C	Try linear interpolation
            P=2.*DX*R3
            P1=1.-R3
          ELSE
C	Try inverse quadratic interpolation
            R1=FA/FC
            R2=FB/FC
            P=R3*(2.*DX*R1*(R1-R2)-(B-A)*(R2-1))
            P1=(R1-1.)*(R2-1.)*(R3-1.)
          ENDIF
          IF(P.GT.0) THEN
            P1=-P1
          ELSE
            P=-P
          ENDIF

          E=D
          IF(2.*P.LT.3.*DX*P1-ABS(EPS*P1).AND.P.LT.ABS(0.5*E*P1)) THEN
C	Accept the interpolated value
            D=P/P1
          ELSE
C	otherwise, perform bisection
            D=DX
            E=DX
          ENDIF
        ENDIF

        A=B
        FA=FB
        IF(ABS(D).GT.EPS) THEN
          B=B+D
        ELSE
C	If the change is too small, shift by EPS
          B=B+SIGN(EPS,DX)
        ENDIF
        FB=F(B)
        IF(FB.GT.0.0.EQV.FC.GT.0.0) THEN
C	If F(B) and F(C) have the same sign, then replace C by A
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
2000  CONTINUE

C	Iteration fails to converge
      IER=428
      X=B
      END

C     ----------------------------------------------------

C	To perform one step of integration of ordinary differential equations
C	using a fourth-order Runge-Kutta method 
C
C	N : (input) Number of first order differential equations to be solved
C	T : (input) Initial value of independent variable t, at
C		which initial values are specified. This value is not updated.
C	Y0 : (input) Real array of length N containing the initial
C		values of variables.
C	DY0 : (input) Real array of length N containing the derivatives
C		of Y at the initial point Y0
C	H : (input) The step size to be used for integration
C	Y1 : (output) Real array of length N containing the solution at t=T+H
C	DIF : (input) Name of subroutine to calculate the right hand side
C		of differential equation y'=f(t,y)
C	WK : Real array of length 2N used as scratch space
C
C	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
C	the differential equation. T is the value of independent variable,
C	N is the number of variables, Y is a real array of length N containing
C	the values of variables. DY is a real array of length N which should
C	contain the calculated values of derivatives at (T,Y).
C
C	Required routines : DIF
C	
C	
      SUBROUTINE RK4(N,T,Y0,DY0,H,Y1,DIF,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y0(N),DY0(N),Y1(N),WK(N,2)

      H2=H/2.
      DO 1000 I=1,N
        Y1(I)=H2*DY0(I)
C	y_0+0.5k_1
1000  WK(I,1)=Y0(I)+Y1(I)

      T1=T+H2
      CALL DIF(T1,N,WK,WK(1,2))
      DO 1200 I=1,N
        Y1(I)=Y1(I)+H*WK(I,2)
C	y_0+0.5k_2
1200  WK(I,1)=Y0(I)+H2*WK(I,2)

      CALL DIF(T1,N,WK,WK(1,2))
      DO 1400 I=1,N
        Y1(I)=Y1(I)+H*WK(I,2)
C	y_0+k_3
1400  WK(I,1)=Y0(I)+H*WK(I,2)

      T1=T+H
      CALL DIF(T1,N,WK,WK(1,2))
      DO 1600 I=1,N
1600  Y1(I)=Y0(I)+(Y1(I)+H2*WK(I,2))/3.
      END

C     ----------------------------------------------------

C	To solve initial value problems in ordinary differential equations
C	using a second or fourth-order Runge-Kutta method with adaptive
C	step size control
C
C	N : (input) Number of first order differential equations to be solved
C	Y : (input/output) Real array of length N containing the initial
C		values of variables. After execution it will contain the
C		values at the last point where the integration
C		has been successful.
C	DY : (output) Real array of length N containing the derivatives
C		of Y at the last point
C	DIF : (input) Name of subroutine to calculate the right hand side
C		of differential equation y'=f(t,y)
C	H : (input/output) Initial guess for the step size. After execution
C		it will contain the step size used by the program
C	T0 : (input/output) Initial value of independent variable t, at
C		which initial values are specified. After execution it will
C		be set to the point up to which integration has been successful
C	TN : (input) The final value of t at which the solution is required.
C		If integration is successful T0 will be set equal to TN.
C		Intermediate values will not be preserved so if solution
C		is required at intermediate points, TN must be set to first
C		such value and multiple calls will be needed to calculate
C		all required values. For each subsequent call only TN needs
C		to be updated.
C	REPS : (input) Required accuracy in each component of the solution.
C		The subroutine only controls local truncation error and hence
C		actual error could be larger
C	NSTEP : (output) Number of steps required to complete the integration
C		Each step requires 10 or 11 calls to DIF (using RK4)
C	NMAX : (input/output) Maximum number of steps to be used. If NMAX.LE.0
C		it will be set to a default value of NMX=10000.
C	IER : (output) Error parameter; IER=0 implies successful execution
C		IER=701 implies N.LE.0, no calculations are done
C		IER=721 implies that step-size has become smaller than
C			REPS*|TN-T0|
C		IER=722 implies that step size is too small for arithmetic used
C		IER=723 implies that integration could not be completed in
C			the specified number of steps.
C	WK : Real array of length 5N used as scratch space
C
C	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
C	the differential equation. T is the value of independent variable,
C	N is the number of variables, Y is a real array of length N containing
C	the values of variables. DY is a real array of length N which should
C	contain the calculated values of derivatives at (T,Y).
C
C	Required routines : RK4 (or RK2), DIF
C	
C	
      SUBROUTINE RKM(N,Y,DY,DIF,H,T0,TN,REPS,NSTEP,NMAX,IER,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=10000,SFAC=0.9D0,EPS=1.D-30)
      PARAMETER(E1=0.2D0,E2=0.25D0,CE=15.D0)
C	For use with RK2 use the following statement instead of the preceding one
C      PARAMETER(E1=0.33D0,E2=0.5D0,CE=3.D0)
      EXTERNAL DIF
      DIMENSION Y(N),DY(N),WK(N,5)

      IF(N.LE.0) THEN
        IER=701
        RETURN
      ENDIF
      IF(NMAX.LE.0) NMAX=NMX
      NSTEP=0
      TSTEP=TN-T0
      IER=0
      IF(TSTEP.EQ.0.0) RETURN
C	Adjust the initial value of H if needed
      IF(H.EQ.0.0.OR.ABS(H).GT.ABS(TSTEP)) H=TSTEP
      IF(H.LT.0.0.EQV.TN.GT.T0) H=-H
C	Calculate the derivatives at the initial point
      CALL DIF(T0,N,Y,DY)

C	Loop for integration
1000  NSTEP=NSTEP+1
C	Use two steps of h/2
      H2=H/2.
      CALL RK4(N,T0,Y,DY,H2,WK(1,4),DIF,WK)
      T1=T0+H2
      CALL DIF(T1,N,WK(1,4),WK(1,5))
      CALL RK4(N,T1,WK(1,4),WK(1,5),H2,WK(1,3),DIF,WK)

C	Use single step of h
      CALL RK4(N,T0,Y,DY,H,WK(1,4),DIF,WK)

C	Estimate the truncation error
      ERR=0
      DO 2000 I=1,N
        R2=ABS(Y(I))+ABS(WK(I,3)-Y(I))+EPS
        R1=ABS(WK(I,4)-WK(I,3))/R2
        IF(R1.GT.ERR) ERR=R1
2000  CONTINUE
      ERR=ABS(ERR*TSTEP/(CE*H))

      IF(T1.EQ.T0) THEN
C	Step size is too small for arithmetic used
        IER=722
        RETURN
      ENDIF

      IF(ERR.LT.REPS) THEN
C	Integration at this step is successful, update T0 and Y
        T0=T0+H
        DO 3000 I=1,N
3000    Y(I)=WK(I,3)-(WK(I,4)-WK(I,3))/CE
        CALL DIF(T0,N,Y,DY)
        IF(ABS((TN-T0)/TSTEP).LT.REPS) RETURN

C	Adjust the step size
        IF(ERR.EQ.0.0) THEN
          H=2.*H
        ELSE
          H=SFAC*H*(REPS/ERR)**E1
        ENDIF
        IF(T0+H.GT.TN.EQV.H.GT.0) H=TN-T0
      ELSE

C	If the integration is not successful, try again with smaller step size
        H=SFAC*H*(REPS/ERR)**E2
      ENDIF

      IF(ABS(H/TSTEP).LT.REPS) THEN
C	Step size is too small
        IER=721
        RETURN
      ENDIF
      IF(NSTEP.LT.NMAX) GO TO 1000
      IER=723
      END

C     ----------------------------------------------------

      FUNCTION FS(XS)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL DIF
      DIMENSION X(100),DX(100),WK(500)

C     FUNCTION ROUTINE SPECIFYING THE SECOND BOUNDARY CONDITION
C     WHOSE ZERO WILL GIVE THE REQUIRED INITIAL CONDITIONS

C     SETUP THE INITIAL VALUES
      T0=0.0
      X(1)=0.0
      X(2)=XS
      N=2
      H=0.01D0
      NMAX=10000
      TN=1
      REPS=1.D-7

C     INTEGRATE THE EQUATION UP TO THE SECOND BOUNDARY
      CALL RKM(N,X,DX,DIF,H,T0,TN,REPS,NSTEP,NMAX,IER,WK)
      IF(IER.GT.0) THEN
        FS=1000.
      ELSE
        FS=X(1)-1.
      ENDIF
      END

C   --------------------------------------------------


      SUBROUTINE DIF(T,N,X,DX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(*),DX(*)
      DATA AL/5.0/

C     THE DIFFERENTIAL EQUATION  Y1'=Y2,  Y2'=AL*SINH(AL*Y1)

      DX(1)=X(2)
      DX(2)=AL*SINH(AL*X(1))
      END
