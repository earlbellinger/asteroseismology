C     PROGRAM TO FIND INVERSE LAPLACE TRANSFORM

      PROGRAM LAPLAC
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      EXTERNAL CFS
      DIMENSION F(200),T(200)

C     EXAMPLE 9.7

C	The exact solution
      FX(Y)=2.-EXP(-Y)

51    FORMAT('  IER =',I4,5X,'EXPONENTIAL ORDER =',1PD14.6/
     1       (3X,'F(',D14.6,') =',D14.6,5X,'EXACT VALUE =',D14.6))

      REPS=1.D-6

C     THE PROGRAM CHOOSES N UNIFORMLY SPACED POINTS IN THE INTERVAL
C     [0,TMAX] TO EVALUATE THE RESULT

100   PRINT *,'TYPE TMAX, N=NO. OF PTS, ALPHA=EXPONENTIAL ORDER'
      PRINT *,'                 (QUITS WHEN N.LE.0)'
      READ *,TMAX,N,ALPHA
      IF(N.LE.0) STOP
      H=TMAX/(N-1)

      DO 1000 I=1,N
        T(I)=H*(I-1)
1000  CONTINUE

      CALL LAPINV(N,T,F,CFS,ALPHA,REPS,IER)
      WRITE(6,51) IER,ALPHA,(T(I),F(I),FX(T(I)),I=1,N)
      GO TO 100
      END

C     ---------------------------------------

C	To calculate inverse Laplace transform
C
C	N : (input) Number of points at which inverse transform needs to
C		be evaluated
C	T : (input) Real array of length N, containing the points at
C		which inverse transform needs to be calculated. These
C		elements can be in any order but the last element must
C		be the maximum.
C	F : (output) Real array of length N, containing the calculated
C		function values. F(I) is the inverse Laplace transform at T(I)
C	CFS : (input) Name of the function routine to calculate the
C		Laplace transform, which is to be inverted
C	ALPHA : (input) The estimated exponential order of the function
C	REPS : (input/output) The required relative accuracy
C	IER : (output) The error parameter, IER=0 implies successful execution
C		IER=61 implies that epsilon-algorithm failed for at least
C			one value of T(I)
C		IER=62 implies that denominator is zero for some T(I),
C			in this case the last value is accepted
C
C	FUNCTION CFS(CS) must be supplied by the user. Here both CFS and
C		and CS are complex variables
C
C	Required routines : CFS
C
      SUBROUTINE LAPINV(N,T,F,CFS,ALPHA,REPS,IER)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      PARAMETER(PI=3.14159265358979324D0,NMIN=3,NMAX=100)
      DIMENSION T(N),F(N),CWK(NMAX+1),ET(NMAX,2)

      TN=T(N)
      IF(REPS.LE.0.0) REPS=1.D-6
      A=ALPHA-0.5*LOG(0.5*REPS)/TN
      PIT=PI/TN
      CS=A
      CWK(1)=CFS(CS)
      NF=1
      IER=0
      CI=(0.D0,1.D0)

C	Loop over the N points
      DO 5000 IT=1,N
        PI0=PIT*T(IT)
        F1=EXP(A*T(IT))/TN
        S1=0.5*CWK(1)
        ND=1
        NT=2
        ET(1,1)=0.0
        I1=1
        I2=2

        DO 4000 I=1,NMAX-1
          DO 2200 J=ND,NT
            CS=A+CI*J*PIT
            IF(J.GT.NF-1) THEN
C	Evaluate the function
              CF=CFS(CS)
              NF=NF+1
              CWK(NF)=CF
            ELSE
C	else use the old value of function
              CF=CWK(J+1)
            ENDIF

            T2=CF
            T1=T2*COS(J*PI0)-IMAG(CF)*SIN(J*PI0)
2200      S1=S1+T1
C	The second column of epsilon-table
          ET(2,I2)=S1
          ET(1,I2)=0.0
          IF(I.EQ.2) THEN
            DIF=ABS(ET(2,I2)-ET(2,I1))
            RI=ET(2,I2)
          ENDIF

C	The epsilon-algorithm (ET(I,J)=epsilon(j-i+2,i-2))
          DO 2400 J=3,I+1
            DEN=ET(J-1,I2)-ET(J-1,I1)
            IF(DEN.NE.0.0) THEN
              ET(J,I2)=ET(J-2,I1)+1./DEN
            ELSE
C	denominator is zero in epsilon-algorithm
              IER=62
              ET(J,I2)=ET(J-2,I1)
            ENDIF
2400      CONTINUE

          IF(I.GT.2) THEN
C	Perform convergence check
            DO 2600 J=2,I-1,2
              DIF1=ABS(ET(J,I2)-ET(J,I1))
              IF(DIF1.LT.DIF) THEN
                DIF=DIF1
                RI=ET(J,I2)
              ENDIF
2600        CONTINUE
          ENDIF

          ND=NT+1
          IF(100.*DIF.LT.REPS*ABS(RI).AND.I.GT.NMIN) GO TO 4500

          II=I1
          I1=I2
          I2=II
4000    NT=ND
C	epsilon-algorithm fails to converge
        IER=61

4500    F(IT)=F1*RI
5000  CONTINUE
      END

C     --------------------------------------

      FUNCTION CFS(CS)
      IMPLICIT COMPLEX*16(C)

C     THE REQUIRED FUNCTION

      CFS=2./CS-1./(CS+1)
      END

