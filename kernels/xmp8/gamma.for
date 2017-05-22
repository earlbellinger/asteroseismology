C     TO EVALUATE GAMMA FUNCTION AND ITS NATURAL LOGARITHM

      PROGRAM GAM
      IMPLICIT REAL*8(A-H,O-Z)
 
51    FORMAT('   X =',1PD14.6,5X,'GAMMA(X) =',D14.6,5X,
     1       'LOG(|GAMMA(X)|) =',D14.6)

 
100   PRINT *,'TYPE X    (QUITS WHEN X<-1000)'
      READ *,X
      IF(X.LT.-1000) STOP
      F1=GAMMA(X)
      F2=GAMLN(X)
      WRITE(6,51) X,F1,F2
      GO TO 100
      END
 
C     ----------------------------------------------------
 
C	To calculate the natural logarithm of Gamma function
C	For negative GAMMA(XG) it gives LOG(ABS(GAMMA(XG)))
C
C	Required routines : None
 
      FUNCTION GAMLN(XG)
      IMPLICIT REAL*8(A-H,O-Z)
C	PIS is LOG(SQRT(2*PI))
      PARAMETER(PI=3.1415926535897932385D0,PIS=0.91893853320467274178D0)
      DIMENSION A(2),B(3),A1(6),B1(7)

C	The coefficients of rational function approximations
      DATA A/1.767971449569122937D+00,  2.909421117928672645D-01/
      DATA B/8.333333333333231537D-02,  1.445531763554246280D-01,
     1       2.012779361583001035D-02/
      DATA A1/3.905731686764559737D+03,  2.204952264401381785D+03,
     1       -1.932467485468849660D+03,  4.643360871045442213D+02,
     1       -4.818088806916028754D+01,  1.896853765546068169D+00/
      DATA B1/3.918055655523400310D+03, -1.088116266563809683D+02,
     1        8.203258626193993149D+02, -9.289402000761705906D+01,
     1        6.521113026294866877D+01, -6.090618615608719044D+00,
     1        1.475909104740280784D+00/
 
      X=ABS(XG)
      IF(X.GT.1000.0) THEN
C	Use asymptotic formula (Stirling formula)
        GX=(1+1.D0/(12*X)+1./(288*X*X)
     1      -139/(51840*X**3)-571./(2488320D0*X**4))
        GAMLN=(X-0.5)*LOG(X)-X+PIS+LOG(GX)

      ELSE IF(X.GT.8.0) THEN
C	Use rational function approximation to Log(Gamma)
        Y=1./X**2
        RMK=((B(3)*Y+B(2))*Y+B(1))/((A(2)*Y+A(1))*Y+1)
        GAMLN=(X-0.5)*LOG(X)-X+PIS+RMK/X

      ELSE IF(X.GE.2.0) THEN
C	Use rational function approximation to Gamma over [2,3]
C	after translating the range if necessary
        F1=1.0
        X1=X
2500    IF(X1.LE.3) GO TO 3000
        F1=F1*(X1-1)
        X1=X1-1
        GO TO 2500
3000    IF(X1.EQ.3) THEN
          GAMLN=LOG(2.*F1)
        ELSE IF(X1.EQ.2) THEN
          GAMLN=LOG(F1)
        ENDIF

        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1      B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1      A1(1))*X1+1
        GAMLN=LOG(F1*FN/FD)

      ELSE IF(X.GT.0.0) THEN
C	Use rational function approximation to Gamma over [2,3]
C	after translating the range
        F1=1./X
        X1=X+1
        IF(X.LT.1) THEN
          F1=F1/X1
          X1=X1+1
        ENDIF
        IF(X1.EQ.2) THEN
          GAMLN=LOG(F1)
        ENDIF
 
        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1      B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1      A1(1))*X1+1
        GAMLN=LOG(F1*FN/FD)
 
      ENDIF

      IF(XG.GT.0.0) RETURN
      GAMLN=LOG(PI/(X*ABS(SIN(PI*X))))-GAMLN
 
      END
 
C     ----------------------------------------------------
 
C	To calculate Gamma function for any real value of XG
C	Use GAMLN for calculating the logarithm of Gamma function
C	which may be useful for large arguments or when argument is
C	close to a negative integer.
C
C	Required routines : None
 
      FUNCTION GAMMA(XG)
      IMPLICIT REAL*8(A-H,O-Z)
C	PIS is SQRT(2*PI)
      PARAMETER(PI=3.14159265358979323846D0,PIS=2.5066282746310005024D0)
      DIMENSION A(2),B(3),A1(6),B1(7)

C	The coefficients for rational function approximations
      DATA A/1.767971449569122937D+00,  2.909421117928672645D-01/
      DATA B/8.333333333333231537D-02,  1.445531763554246280D-01,
     1       2.012779361583001035D-02/
      DATA A1/3.905731686764559737D+03,  2.204952264401381785D+03,
     1       -1.932467485468849660D+03,  4.643360871045442213D+02,
     1       -4.818088806916028754D+01,  1.896853765546068169D+00/
      DATA B1/3.918055655523400310D+03, -1.088116266563809683D+02,
     1        8.203258626193993149D+02, -9.289402000761705906D+01,
     1        6.521113026294866877D+01, -6.090618615608719044D+00,
     1        1.475909104740280784D+00/
 
      X=ABS(XG)
      IF(X.GT.1000.0) THEN
C	Use asymptotic approximation (Stirling formula)
        GX=(1+1.D0/(12*X)+1./(288*X*X)
     1         -139/(51840*X**3)-571./(2488320D0*X**4))
        GAMMA=X**(X-0.5)*EXP(-X)*PIS*GX

      ELSE IF(X.GT.8.0) THEN
C	Use rational function approximation for Log(Gamma) 
        Y=1./X**2
        RMK=((B(3)*Y+B(2))*Y+B(1))/((A(2)*Y+A(1))*Y+1)
        GAMMA=X**(X-0.5)*EXP(-X)*PIS*EXP(RMK/X)

      ELSE IF(X.GE.2.0) THEN
C	Use rational function approximation for (Gamma) over [2,3]
C	after translating the range if necessary
        F1=1.0
        X1=X
2500    IF(X1.LE.3) GO TO 3000
        F1=F1*(X1-1)
        X1=X1-1
        GO TO 2500
3000    IF(X1.EQ.3) THEN
          GAMMA=F1*2
        ELSE IF(X1.EQ.2) THEN
          GAMMA=F1
        ENDIF

        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1          B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1          A1(1))*X1+1
        GAMMA=F1*FN/FD
      ELSE IF(X.GT.0.0) THEN
C	Use rational function approximation for (Gamma) over [2,3]
C	after translating the range if necessary
        F1=1./X
        X1=X+1
        IF(X.LT.1) THEN
          F1=F1/X1
          X1=X1+1
        ENDIF
        IF(X1.EQ.2) GAMMA=F1

        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1          B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1          A1(1))*X1+1
        GAMMA=F1*FN/FD
 
      ENDIF

      IF(XG.GT.0.0) RETURN
      GAMMA=PI/(XG*SIN(PI*X)*GAMMA)
 
      END
