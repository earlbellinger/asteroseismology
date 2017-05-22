C     PROGRAM TO FIND EIGENVALUES AND EIGENVECTORS OF A HERMITIAN MATRIX
 
      PROGRAM HERMIT
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      DIMENSION ZA(30,30),D(30),E(30),EV(30,30),EI(30),WK(210),ZV(30,30)
C     DATA ((ZA(J,I),I=1,3),J=1,3)/1,(2,-1),(3,-2),(2,1),(4,0),
C     1  (5,-3),(3,2),(5,3),(6,0)/
      DATA ((ZA(J,I),I=1,4),J=1,4)/(7,0),(3,0),(1,2),(-1,2),(3,0),
     1             (7,0),(1,-2),(-1,-2),(1,-2),(1,2),(7,0),(-3,0),
     2             (-1,-2),(-1,2),(-3,0),(7,0)/
 
51    FORMAT(10X,13HTHE MATRIX IS)
52    FORMAT(5(2F8.3,2X))
56    FORMAT(1P6D13.5)
57    FORMAT('    IER =',I4/'   EIGENVALUE',9X,'EIGENVECTOR')
58    FORMAT(1PD14.6,3X,2D14.6,3X,2D14.6/(17X,2D14.6,3X,2D14.6))
 
100   PRINT *,'TYPE N=ORDER,       (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0) STOP
 
      WRITE(6,51)
      PRINT *,'TYPE IN THE COMPLEX MATRIX ROW-WISE'
      DO 500 I=1,N
        PRINT *,I,'TH ROW'
        READ *,(ZA(I,J),J=1,N)
        WRITE(6,52) (ZA(I,J),J=1,N)
500   CONTINUE
 
      IA=30
 
      REPS=1.D-14
 
      CALL HEREVP(ZA,N,IA,EI,ZV,IA,WK,REPS,IER)
      WRITE(6,57) IER
      DO 2000 I=1,N
        WRITE(6,58) EI(I),(ZV(J,I),J=1,N)
2000  CONTINUE
      GO TO 100
      END
 
C     ------------------------------------------------------
 
C	To find eigenvalues and eigenvectors of a complex Hermitian matrix
C
C	ZA : (input) Complex array of length IA*N containing the matrix
C	N : (input) Order of the matrix
C	IA : (input) First dimension of array ZA as declared in the calling program
C	EI : (output) Real array of length N containing the eigenvalues
C	ZV : (output) Complex array of length IZ*N containing the eigenvectors
C		ZV(i,j) is the ith component of Jth eigenvector
C	IZ : (input) First dimension of array ZV as declared in the calling program
C	WK : Real array of length 2N*(2N+2) used as scratch space
C	REPS : (input) Required tolerance, should be of the order of machine
C		accuracy
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=111 implies that N.LE.1 or N>IA or N>IZ
C			in which case no calculations are done
C		Other values may be set by TRED2 or TQL2
C
C	Required routines : TRED2, TQL2
 
      SUBROUTINE HEREVP(ZA,N,IA,EI,ZV,IZ,WK,REPS,IER)
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      DIMENSION ZA(IA,N),WK(2*N,2*N+2),EI(N),ZV(IZ,N)
 
      IF(N.LE.1.OR.N.GT.IA.OR.N.GT.IZ) THEN
        IER=111
        RETURN
      ENDIF

C	Setup the 2N*2N real matrix
      DO 2000 I=1,N
        DO 2000 J=1,N
          WK(I,J)=ZA(I,J)
          WK(I+N,J+N)=WK(I,J)
          WK(I+N,J)=IMAG(ZA(I,J))
          WK(I,J+N)=-WK(I+N,J)
2000  CONTINUE
      N2=2*N
      IW=N2

C	To reduce the 2N*2N matrix to tridiagonal form
      CALL TRED2(WK,N2,IW,WK(1,N2+1),WK(1,N2+2),IER)
      IF(IER.GT.100) RETURN

C	Find eigenvalues and eigenvectors of tridiagonal matrix
      CALL TQL2(WK,N2,IW,WK(1,N2+1),WK(1,N2+2),REPS,IER)
      IF(IER.GT.100) RETURN
 
C	Since all eigenvalues are repeated and sorted in ascending order
C	pick alternate eigenvalues and eigenvectors
C	May create problem for multiple eigenvalues.
      ZI=(0.D0,1.D0)
      DO 3000 I=1,N
        I1=2*I-1
        EI(I)=WK(I1,N2+1)
        DO 2500 J=1,N
          ZV(J,I)=WK(J,I1)+ZI*WK(J+N,I1)
2500    CONTINUE
3000  CONTINUE
 
      END
 
C     -----------------------------------------------------
 
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
 
C     -----------------------------------------------------------
 
C	Reduction of real symmetric matrix to tridiagonal form using
C	Householder's method
C
C	A : (input/output) Real array of length IA*N containing the matrix
C		elements. After execution it will be overwritten by the
C		transformation matrix
C	N : (input) Order of matrix
C	IA : (input) The first dimension of A as specified in the calling program
C	D : (output) Diagonal elements of the transformed tridiagonal matrix
C		D(I) would contain A(I,I)
C	E : (output) Off-diagonal elements of the transformed tridiagonal matrix
C		E(I+1) would contain A(I,I+1)=A(I+1,I)
C	IER : (output) Error parameter, IER=0 implies successful execution
C		IER=107 implies that N.LE.1 or N>IA, in which case no
C			calculations are done
C
C	Required routines : None

      SUBROUTINE TRED2(A,N,IA,D,E,IER)
      IMPLICIT REAL*8(A-H,O-Z)
C	For REAL*4 use REPS=1.E-30
C      PARAMETER(REPS=1.E-30)
      PARAMETER(REPS=1.D-300)
      DIMENSION A(IA,N),D(N),E(N)

      IF(N.LE.1.OR.N.GT.IA) THEN
        IER=107
        RETURN
      ENDIF
      IER=0

      DO 4000 I=N,2,-1
        F=A(I,I-1)
        G=0
        DO 2000 K=1,I-2
2000    G=G+A(I,K)*A(I,K)
        H=G+F*F
        IF(G.LE.REPS) THEN
C	Skip the transformation
          E(I)=F
          H=0
        ELSE

          G=SQRT(H)
          IF(F.GT.0.0) G=-G
          E(I)=G
          H=H-F*G
          A(I,I-1)=F-G
          F=0
          DO 2600 J=1,I-1
C	Elements of u_i/H_i
            A(J,I)=A(I,J)/H
            G=0
C	Form elements of A_iu_i
            DO 2200 K=1,J
2200        G=G+A(J,K)*A(I,K)
            DO 2400 K=J+1,I-1
2400        G=G+A(K,J)*A(I,K)
C	Components of p_i
            E(J)=G/H
            F=F+G*A(J,I)
2600      CONTINUE

C	calculate u_i^Tp_i/2H_i
          HH=0.5*F/H
          DO 3400 J=1,I-1
            F=A(I,J)
            G=E(J)-HH*F
C	Elements of q_i
            E(J)=G
            DO 3000 K=1,J
3000        A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
3400      CONTINUE
        ENDIF
        D(I)=H
4000  CONTINUE

      D(1)=0
      E(1)=0
C	accumulation of transformation matrix Q
      DO 5000 I=1,N
        IF(D(I).NE.0.0) THEN
          DO 4600 J=1,I-1
            G=0
            DO 4200 K=1,I-1
4200        G=G+A(I,K)*A(K,J)
            DO 4400 K=1,I-1
4400        A(K,J)=A(K,J)-G*A(K,I)
4600      CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1
        DO 4800 J=1,I-1
          A(I,J)=0.0
4800    A(J,I)=0.0
5000  CONTINUE
      END
