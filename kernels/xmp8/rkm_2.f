C     PROGRAM FOR SOLVING INITIAL VALUE PROBLEM IN ORDINARY DIFFERENTIAL
C     EQUATIONS USING SECOND ORDER RUNGE-KUTTA METHOD WITH ADAPTIVE
C     STEP SIZE CONTROL

      PROGRAM RUNGE
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL DIF
      DIMENSION X(100),DX(100),WK(500)

C     EXAMPLE 11.6
C     AL IS THE PARAMETER LAMBDA IN THE EQUATION, WHICH IS PASSED
C     ON TO SUBROUTINE DIF VIA THE COMMON BLOCK

      COMMON/FN/AL

51    FORMAT('   IER =',I4,5X,'T =',1PD14.6,5X,'NO. OF STEPS =',I7/
     1      5X,'STEP LENGTH =',D14.6,5X,'SOLUTION =',2D14.6/(2X,5D14.6))
52    FORMAT('    INITIAL VALUES =',1P4D14.6/(2X,5D14.6))
53    FORMAT('   LAMDA =',1PD13.5,4X,'INITIAL STEP SIZE =',D13.5,4X,
     1       'T0 =',D13.5)

      N=1
      REPS=1.D-4
      NMAX=100000

C     ACCUMULATE THE NUMBER OF TIME STEPS USED IN NM

      NM=0
      PRINT *,'TYPE T0=INITIAL TIME,  H=INITIAL STEP LENGTH, LAMDA' 
      READ *,T0,H,AL
      PRINT *,'TYPE INITIAL VALUES'
      READ *,(X(I),I=1,N)
      WRITE(6,53) AL,H,T0
      WRITE(6,52) (X(I),I=1,N)

C     GO ON TYPING DIFFERENT INTERMEDIATE VALUES OF T
C     AT WHICH THE SOLUTION IS REQUIRED

100   PRINT *,'TYPE TN=FINAL TIME    (QUITS WHEN TN.LT.-10)'
      READ *,TN
      IF(TN.LT.-10.0) STOP
      CALL RKM_2(N,X,DX,DIF,H,T0,TN,REPS,NSTEP,NMAX,IER,WK)
      NM=NM+NSTEP
      WRITE(6,51) IER,T0,NM,H,(X(I),I=1,N)
      IF(IER.EQ.0) GO TO 100
      END

C     ---------------------------------------

C	To perform one step of integration of ordinary differential equations
C	using a second-order Runge-Kutta method 
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
C	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
C	the differential equation. T is the value of independent variable,
C	N is the number of variables, Y is a real array of length N containing
C	the values of variables. DY is a real array of length N which should
C	contain the calculated values of derivatives at (T,Y).
C
C	Required routines : DIF
C	
C	
      SUBROUTINE RK2(N,T,Y0,DY0,H,Y1,DIF,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y0(N),DY0(N),Y1(N),WK(N,2)

      DO 1000 I=1,N
        Y1(I)=H*DY0(I)
1000  WK(I,1)=Y0(I)+Y1(I)

      T1=T+H
      CALL DIF(T1,N,WK,WK(1,2))
      DO 1200 I=1,N
1200  Y1(I)=Y0(I)+(Y1(I)+H*WK(I,2))/2.
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
C		Each step requires 4 or 5 calls to DIF (using RK2)
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
C	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
C	the differential equation. T is the value of independent variable,
C	N is the number of variables, Y is a real array of length N containing
C	the values of variables. DY is a real array of length N which should
C	contain the calculated values of derivatives at (T,Y).
C
C	Required routines : RK2, DIF
C	
C	
      SUBROUTINE RKM_2(N,Y,DY,DIF,H,T0,TN,REPS,NSTEP,NMAX,IER,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=10000,SFAC=0.9D0,EPS=1.D-30)
C      PARAMETER(E1=0.2D0,E2=0.25D0,CE=15.D0)
C	For use with RK4 use the preceding statement instead of the following one
      PARAMETER(E1=0.33D0,E2=0.5D0,CE=3.D0)
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
      CALL RK2(N,T0,Y,DY,H2,WK(1,4),DIF,WK)
      T1=T0+H2
      CALL DIF(T1,N,WK(1,4),WK(1,5))
      CALL RK2(N,T1,WK(1,4),WK(1,5),H2,WK(1,3),DIF,WK)

C	Use single step of h
      CALL RK2(N,T0,Y,DY,H,WK(1,4),DIF,WK)

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


C   --------------------------------------------------

      SUBROUTINE DIF(T,N,X,DX)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/FN/AL
      DIMENSION X(*),DX(*)

C	The differential equation

      DX(1)=AL*(X(1)-T**3)+3.*T**2
      END
