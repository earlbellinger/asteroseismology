C     PROGRAM TO CALCULATE WEIGHTS AND ABSCISSAS OF GAUSS-LAGUERRE QUADRATURE FORMULAS

      PROGRAM GAUSSL
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(65),W(65),COF(3,65),WK(9000)
 
51    FORMAT(5X,'IER =',I4,'    N =',I4/9X,'ABSCISSAS',10X,'WEIGHTS'/
     1       (1PD22.14,5X,D22.14))
52    FORMAT('   MAXIMUM RELATIVE ERROR USING ',I4,' POINTS =',1PD12.4)
 
C	For Laguerre polynomials
      RI0=1.D0

100   PRINT *,'TYPE  NT=NO. OF POINTS IN QUADRATURE FORMULA'
      PRINT *,'           (QUITS WHEN NT.LE.0)'
      READ *, NT
      IF(NT.LE.0) STOP

      DO 150 I=1,NT+1
C	For Laguerre polynomials
        COF(1,I)=1.0
        COF(2,I)=-(2*I-1)
        COF(3,I)=(I-1.D0)**2
150   CONTINUE
      CALL GAUSRC(NT,W,X,COF,RI0,IER,WK)
      WRITE(6,51) IER,NT,(X(I),W(I),I=1,NT)

C     USE THE WEIGHTS AND ABSCISSAS TO EVALUATE THE INTEGRAL X**N FOR TESTING
 
      REX=1.0
      ERR=0.0
      DO 2000 J=0,2*NT-1
        RI=0.0
        DO 1000 I=1,NT
          RI=RI+W(I)*X(I)**J
1000    CONTINUE
        ERR=MAX(ERR,ABS(RI-REX)/REX)
    	REX=REX*(J+1)
2000  CONTINUE
 
      WRITE(6,52) NT,ERR
      GO TO 100
      END
 
C     --------------------------------------------------
 
C     To calculate weights and abscissas of a quadrature formula with
C     specified weight function when recurrence relation for orthogonal
C     polynomial is known.
C
C     N : (input) Number of points in the required quadrature formula
C     W : (output) Array of length N, which will contain the weights
C     AB : (output) Array of length N containing the abscissas
C     COF : (input) Array of length 3*N containing the coefficients of
C		the recurrence relation for orthogonal polynomials
C		P_i(x)=(COF(1,i)*x+COF(2,i))P_{i-1}(x) - COF(3,i)*P_{i-2}(x)
C     RI0 : (input) The integral of weight function over the required
C		interval.
C     IER : (output) Error parameter, IER=0 implies successful execution
C		IER=302 implies N.LE.0
C		IER=321 implies that some coefficient becomes imaginary
C			during calculations.
C		In both these cases calculations are abandoned.
C		Other values may be set by TQL2
C     WK : Real array of length N*(N+2) used as scratch space
C
C     Required routines : TQL2
 
      SUBROUTINE GAUSRC(N,W,AB,COF,RI0,IER,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (REPS=1.D-15)
C	For REAL*4 use REPS=6.E-8
C      PARAMETER (REPS=6.E-8)
      DIMENSION WK(N,N+2),COF(3,N),W(N),AB(N)
 
      IER=302
      IF(N.LE.0) RETURN
      LJ=N
 
C     Calculate the coefficients of symmetric tridiagonal matrix
      DO 2000 I=1,N
        WK(I,N+1)=-COF(2,I)/COF(1,I)
        IF(I.LT.N) THEN
          R1=COF(3,I+1)/(COF(1,I)*COF(1,I+1))
          IF(R1.GE.0.0) THEN
            WK(I+1,N+2)=SQRT(R1)
          ELSE
            IER=321
            RETURN
          ENDIF
        ENDIF
        DO 1800 J=1,N
          WK(J,I)=0.0
1800    CONTINUE
        WK(I,I)=1.0
2000  CONTINUE
 
C     Find eigenvalues and eigenvectors of the tridiagonal matrix
      CALL TQL2(WK,N,LJ,WK(1,N+1),WK(1,N+2),REPS,IER)
      IF(IER.GT.0) RETURN
 
C     Calculate the abscissas and weights
      DO 3000 I=1,N
        AB(I)=WK(I,N+1)
        W(I)=WK(1,I)**2*RI0
3000  CONTINUE
 
      END
 
C     --------------------------------------------------
 
C	To find eigenvalues and eigenvectors of Z T Z^T using QL algorithm
C	where T is a symmetric tridiagonal matrix and Z is an orthogonal matrix.
C	If Z is the transformation matrix to reduce original matrix to
C	tridiagonal matrix, it will calculate the eigenvectors of original matrix
C
C	Z : (input/output) Real array of length IZ*N which should contain
C		the transformation matrix required to reduce original real
C		symmetric matrix to tridiagonal form. To find eigenvectors
C		of a symmetric tridiagonal matrix, set Z to unit matrix
C		After execution Z will contain the eigenvector of the original
C		matrix Z T Z^T. Z(i,j) should contain the ith component of
C		jth eigenvector
C	N : (input) Order of matrix
C	IZ : (input) The first dimension of array Z as declared in the
C		calling program. (IZ.GE.N)
C	D : (input/output) Real array of length N, containing the diagonal
C		elements of the tridiagonal matrix, D(i)=T(i,i).
C		After execution it will contain the eigenvalues of the matrix
C	E : (input/output) Real array of length N containing the off-diagonal
C		elements of the tridiagonal matrix, E(i)=T(i,i+1)=T(i+1,i)
C		It is used as scratch space and its contents will be destroyed
C		during execution.
C	REPS : (input) Required tolerance, it should be of order of machine
C		accuracy
C	IER : (output) Error parameter; IER=0 implies successful execution
C		IER=108 implies that N.LE.1 or N>IZ, in which case no
C			calculations are performed
C		IER=143 implies that the QL algorithm failed to converge
C			for some eigenvalue, the calculations are abandoned
C
C	Required routines : None
C	
      SUBROUTINE TQL2(Z,N,IZ,D,E,REPS,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=30)
      DIMENSION Z(IZ,N),D(N),E(N)

      IF(N.LE.0.OR.N.GT.IZ) THEN
        IER=108
        RETURN
      ENDIF
      IER=0
      DO 2000 I=2,N
2000  E(I-1)=E(I)
      E(N)=0
      B=0
      F=0

      DO 4000 L=1,N
        H=REPS*(ABS(D(L))+ABS(E(L)))
        IF(B.LT.H) B=H
C	Look for small off-diagonal elements
        DO 2200 M=L,N
          IF(ABS(E(M)).LE.B) GO TO 2300
2200    CONTINUE
        M=N
2300    IF(M.EQ.L) THEN
C	one eigenvalue is isolated
          D(L)=D(L)+F
          GO TO 4000
        ENDIF

C	Loop for QL transformation 
        DO 3600 IT=1,NIT
C	Find shift
          G=D(L)
          P=(D(L+1)-G)/(2*E(L))
          R=SQRT(P*P+1)
          IF(P.LT.0.0) R=-R
          D(L)=E(L)/(P+R)
          H=G-D(L)
          DO 2600 I=L+1,N
2600      D(I)=D(I)-H
          F=F+H

C	The QL transformation
          P=D(M)
          C=1
          S=0
C	Given's rotations
          DO 3400 I=M-1,L,-1
            G=C*E(I)
            H=C*P
            IF(ABS(P).GE.ABS(E(I))) THEN
              C=E(I)/P
              R=SQRT(C*C+1)
              E(I+1)=S*P*R
              S=C/R
              C=1/R
            ELSE
              C=P/E(I)
              R=SQRT(C*C+1)
              E(I+1)=S*E(I)*R
              S=1./R
              C=C/R
            ENDIF
            P=C*D(I)-S*G
            D(I+1)=H+S*(C*G+S*D(I))

C	Transforming the eigenvectors
            DO 3000 K=1,N
              H=Z(K,I+1)
              Z(K,I+1)=S*Z(K,I)+C*H
              Z(K,I)=C*Z(K,I)-S*H
3000        CONTINUE
3400      CONTINUE

          E(L)=S*P
          D(L)=C*P
          IF(ABS(E(L)).LE.B) THEN
C	One eigenvalue is isolated
            D(L)=D(L)+F
            GO TO 4000
          ENDIF
3600    CONTINUE
C	QL iteration fails to converge
        IER=143
        RETURN
4000  CONTINUE

C	Sort eigenvalues in ascending order by straight selection
      DO 5000 I=1,N-1
        K=I
        P=D(I)
        DO 4200 J=I+1,N
          IF(D(J).LT.P) THEN
            K=J
            P=D(J)
          ENDIF
4200    CONTINUE
        IF(K.NE.I) THEN
C	exchange the eigenvalues and eigenvectors
          D(K)=D(I)
          D(I)=P
          DO 4400 J=1,N
            P=Z(J,I)
            Z(J,I)=Z(J,K)
            Z(J,K)=P
4400      CONTINUE
        ENDIF
5000  CONTINUE
      END
