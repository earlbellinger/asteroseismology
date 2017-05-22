C     PROGRAM TO SOLVE INITIAL VALUE PROBLEMS IN 
C     ORDINARY DIFFERENTIAL EQUATIONS USING EXTRAPOLATION METHOD

      PROGRAM BULIR
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL DIF
      COMMON/FN/AL
      DIMENSION X(20),DX(20),WK(800)
      
C     EXAMPLE 11.4

51    FORMAT('    INITIAL VALUES =',1P4D14.6/(2X,5D14.6))
52    FORMAT('   IER =',I4,5X,'T =',1PD14.6,5X,
     1     'NO. OF FUNCTION EVALUATIONS =',I8/5X,'STEP SIZE =',D14.6,
     2      5X,'SOLUTION =',2D14.6/(2X,5D14.6))
53    FORMAT('   T0 =',1PD12.4,5X,'INITIAL STEP SIZE =',D12.4,5X,
     1      'LAMDA =',D12.4/5X,'IFLG =',I3)

      N=1
      REPS=1.D-8
      IFLG=0
      NMAX=100000
      NM=0
100   PRINT *,'TYPE T0=INITIAL TIME,  H=INITIAL STEP SIZE, LAMDA,',
     1        'IFLG =0/1'
      PRINT *,'IFLG =0/1 FOR POLYNOMIAL/RATIONAL FUNCTION EXTRAPOLATION'
      READ *,T0,H,AL,IFLG
      WRITE(6,53) T0,H,AL,IFLG
      PRINT *,'TYPE INITIAL VALUES'
      READ *,(X(I),I=1,N)
      WRITE(6,51) (X(I),I=1,N)

200   PRINT *,'TYPE TN=FINAL TIME,   (QUITS WHEN TN.LT.-10)'
      READ *,TN
      IF(TN.LT.-10.0) STOP
      CALL EXTP(N,X,DX,DIF,H,T0,TN,REPS,NSTEP,NMAX,IER,WK,IFLG)
      NM=NM+NSTEP
      WRITE(6,52) IER,T0,NM,H,(X(I),I=1,N)
      GO TO 200
      END

C     ----------------------------------------------------------

C	To solve initial value problems in ordinary differential equations
C	using extrapolation method
C
C	N : (input) Number of first order differential equations to be solved
C	Y : (input/output) Real array of length N containing the initial
C		values of variables. After execution it will contain the
C		values of variable at the last point where the integration
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
C		to be changed.
C	REPS : (input) Required accuracy in each component of the solution.
C		The subroutine only controls local truncation error and hence
C		actual error could be larger
C	NSTEP : (output) Number of calls to DIF made by the subroutine to
C		complete the integration during each call to EXTP.
C	NMAX : (input/output) Maximum number of function evaluations to be used.
C		If NMAX.LE.0 it will be set to a default value of NMX=100000.
C	IER : (output) Error parameter; IER=0 implies successful execution
C		IER=703 implies N.LE.0, no calculations are done
C		IER=730 implies that step-size has become smaller than
C			REPS*|TN-T0|
C		IER=731 implies that step size is too small for arithmetic used
C		IER=732 implies that integration could not be completed in
C			the specified number of steps.
C		IER=733 implies that denominator for evaluating rational
C			function extrapolation vanished
C	WK : Real array of length 39N used as scratch space
C	IFLG : (input) Integer variable used as a flag to decide the type
C		of extrapolation to be used.
C		If IFLG=0 polynomial extrapolation is used
C		otherwise rational function extrapolation is used
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
      SUBROUTINE EXTP(N,Y,DY,DIF,H,T0,TN,REPS,NSTEP,NMAX,IER,WK,IFLG)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=100000,SFAC=0.8D0,EPS=1.D-30)
      DIMENSION Y(N),DY(N),WK(N,39),NSEQ(12),HI(12)
      DATA NSEQ/2,4,6,8,12,16,24,32,48,64,96,128/

      IF(N.LE.0) THEN
        IER=703
        RETURN
      ENDIF
      IF(NMAX.LE.0) NMAX=NMX
      NSTEP=0
      TSTEP=TN-T0
      IER=0
      IF(TSTEP.EQ.0.0) RETURN
      IF(H.EQ.0.0.OR.ABS(H).GT.ABS(TSTEP)) H=TSTEP
      IF(H.LT.0.0.EQV.TN.GT.T0) H=-H

C	Loop for integration
1000  NSTEP=NSTEP+1
      CALL DIF(T0,N,Y,DY)

C	Loop for extrapolation
      DO 2500 IE=1,12
C	Step size for the midpoint method
        H2=H/NSEQ(IE)
        HI(IE)=H2**2
        DO 1200 I=1,N
          WK(I,1)=Y(I)
1200    WK(I,2)=Y(I)+H2*DY(I)

        T1=T0
        DO 1500 J=1,NSEQ(IE)-1
          T1=T1+H2
          CALL DIF(T1,N,WK(1,2),WK(1,3))
          DO 1400 I=1,N
C	The midpoint rule
            C1=WK(I,1)+2.*H2*WK(I,3)
            WK(I,1)=WK(I,2)
1400      WK(I,2)=C1
1500    CONTINUE
        T1=T1+H2
        CALL DIF(T1,N,WK(1,2),WK(1,3))
        NSTEP=NSTEP+NSEQ(IE)
        DO 1800 I=1,N
C	Modified midpoint method
          WK(I,3+IE)=0.5*(WK(I,1)+WK(I,2)+H2*WK(I,3))
          WK(I,15+IE)=WK(I,3+IE)
          WK(I,27+IE)=WK(I,3+IE)
1800    CONTINUE

        IF(IE.GT.1) THEN
C	Perform extrapolation
          DO 2000 J=IE-1,1,-1
            DO 2000 I=1,N
              IF(IFLG.EQ.0) THEN
C	Use polynomial extrapolation
                CN=(WK(I,4+J)-WK(I,15+J))/(HI(J)-HI(IE))
C	D(J,IE-J)
                WK(I,15+J)=HI(IE)*CN
C	C(J,IE-J)
                WK(I,3+J)=HI(J)*CN
              ELSE

C	Use rational function extrapolation
                DEN=HI(J)*WK(I,15+J)/HI(IE)-WK(I,4+J)
                IF(DEN.EQ.0.0) THEN
                  IER=733
                  GO TO 2800
                ENDIF
                CN=(WK(I,4+J)-WK(I,15+J))/DEN
C	C(J,IE-J)
                WK(I,3+J)=HI(J)*WK(I,15+J)*CN/HI(IE)
C	D(J,IE-J)
                WK(I,15+J)=WK(I,4+J)*CN
              ENDIF
              WK(I,27+J)=WK(I,3+J)+WK(I,27+J)
2000      CONTINUE

C	Estimating the truncation error
          ERR=0
          DO 2200 I=1,N
            R2=ABS(Y(I))+ABS(WK(I,28)-Y(I))+EPS
            R1=ABS(WK(I,28)-WK(I,29))/R2
            IF(R1.GT.ERR) ERR=R1
2200      CONTINUE
          ERR=ABS(ERR*TSTEP/H)
          IF(ERR.LT.REPS) GO TO 2800
        ENDIF
2500  CONTINUE

2800  IF(T1.EQ.T0) THEN
C	The step size is too small for the arithmetic used
        IF(IER.EQ.0) IER=731
        RETURN
      ENDIF

      IF(ERR.LT.REPS) THEN
C	Integration is successful
        T0=T0+H
        DO 3000 I=1,N
3000    Y(I)=WK(I,28)
        IF(ABS((TN-T0)/TSTEP).LT.REPS) RETURN

        IF(IE.LT.7) THEN
C	Increase the step size
          H=SFAC*NSEQ(7)*H/NSEQ(IE)
        ELSE IF(IE.GT.7) THEN
C	Decrease the step size
          H=SFAC*H
        ENDIF
        IF(T0+H.GT.TN.EQV.H.GT.0) H=TN-T0
      ELSE
C	If the integration has failed, then decrease the step size
        H=H/32.
      ENDIF

      IF(ABS(H/TSTEP).LT.REPS) THEN
C	The step size is too small
        IF(IER.EQ.0) IER=730
        RETURN
      ENDIF
      IF(NSTEP.LT.NMAX) GO TO 1000
      IER=732
      END

C   --------------------------------------------------

      SUBROUTINE DIF(T,N,X,DX)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/FN/AL
      DIMENSION X(*),DX(*)

C     THE DIFFERENTIAL EQUATION

      DX(1)=AL*(X(1)-T**3)+3.*T**2
      END
