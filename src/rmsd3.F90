!   COPYRIGHT (C) 2000 MASSIMO MARCHI AND PIERO PROCACCI
!   Full copyright notice at http://www.chim.unifi.it/orac/copyright4.0.html
!   Contact for info M. Marchi, CEA,  Gif Sur Yvette 91191 (FRANCE) 
!   Email:marchi@villon.saclay.cea.fr

!   Rewritten in F90 style by TEODORO LAINO - NEST - INFM 
!   Scuola Normale Superiore di Pisa - t.laino@sns.it
!   CSCS (Manno) - 04/12/02 
!   Copyright Extension 2002 - No modification or partial reproduction of 
!   this document is allowed without a written permission of the authors 
MODULE SUPERPOSITION

  IMPLICIT NONE

CONTAINS

  SUBROUTINE RMSD3( N, R, R0, TRANSL, ROT, ERR, DERR_DR, W, MASSES, IOPT)
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: N, IOPT
    DOUBLE PRECISION, ALLOCATABLE   :: RFIT(:,:)
    DOUBLE PRECISION, INTENT(INOUT) :: R(3,N), R0(3,N), ERR 
    DOUBLE PRECISION, INTENT(IN)    :: W(N), MASSES(N)
    DOUBLE PRECISION  :: MTOT
    DOUBLE PRECISION  :: Q(0:3)
    DOUBLE PRECISION, INTENT(INOUT) :: DERR_DR(3,N), TRANSL(3), ROT(3,3)
    CHARACTER (LEN=80)    :: MSG
    INTEGER :: I, J, K, JOBN, IER, IX, IRET
    DOUBLE PRECISION  :: EPSI, S, WW, DL_AV(3)
    DOUBLE PRECISION  :: M(4,4), LAMBDA(4), Z(4,4), WK(20)
    DOUBLE PRECISION  :: RR0(3), RR(3), DM_R(4,4,3), RP(3,N), R0P(3,N), RRSQ
    DOUBLE PRECISION  :: XX, YY, ZZ
    DATA EPSI / 1.0E-10 /

    ALLOCATE ( RFIT( 3, N ) )

!   ---------------------------------------------------------------------
!   CENTER MOLECULE IN R
!   ---------------------------------------------------------------------
    XX=0.D0
    YY=0.D0
    ZZ=0.D0
    DO I=1,N
       MTOT = MTOT + MASSES(I)
       XX=XX+R(1,I)* MASSES(I)
       YY=YY+R(2,I)* MASSES(I)
       ZZ=ZZ+R(3,I)* MASSES(I)
    ENDDO
    XX=XX/MTOT
    YY=YY/MTOT
    ZZ=ZZ/MTOT
    DO I=1,N
       RP(1,I)=R(1,I)-XX
       RP(2,I)=R(2,I)-YY
       RP(3,I)=R(3,I)-ZZ
    ENDDO
    TRANSL(1)=XX
    TRANSL(2)=YY
    TRANSL(3)=ZZ
!   ---------------------------------------------------------------------
!   CENTER MOLECULE IN R0
!   ---------------------------------------------------------------------
    XX=0.D0
    YY=0.D0
    ZZ=0.D0
    DO I=1,N
       XX=XX+R0(1,I)* MASSES(I)
       YY=YY+R0(2,I)* MASSES(I)
       ZZ=ZZ+R0(3,I)* MASSES(I)
    ENDDO
    XX=XX/MTOT
    YY=YY/MTOT
    ZZ=ZZ/MTOT
    DO I=1,N
       R0P(1,I)=R0(1,I)-XX
       R0P(2,I)=R0(2,I)-YY
       R0P(3,I)=R0(3,I)-ZZ
    ENDDO

    DO I = 1,16
       M(I,1) = 0.D0
    ENDDO
    
    DO I = 1,N
       IF (W(I) .EQ. 0.D0) CYCLE
       RR(1)=RP(1,I)
       RR(2)=RP(2,I)
       RR(3)=RP(3,I)
       RR0(1)=R0P(1,I)
       RR0(2)=R0P(2,I)
       RR0(3)=R0P(3,I)
       RRSQ=W(I)*(RR0(1)**2+RR0(2)**2+RR0(3)**2+RR(1)**2+RR(2)**2+RR(3)**2)
       RR0(1)=W(I)*RR0(1)
       RR0(2)=W(I)*RR0(2)
       RR0(3)=W(I)*RR0(3)
       
       M(1,1)=M(1,1) + RRSQ+2.D0*(-RR0(1)*RR(1)-RR0(2)*RR(2)-RR0(3)*RR(3))
       M(2,2)=M(2,2) + RRSQ+2.D0*(-RR0(1)*RR(1)+RR0(2)*RR(2)+RR0(3)*RR(3))
       M(3,3)=M(3,3) + RRSQ+2.D0*(+RR0(1)*RR(1)-RR0(2)*RR(2)+RR0(3)*RR(3))
       M(4,4)=M(4,4) + RRSQ+2.D0*(+RR0(1)*RR(1)+RR0(2)*RR(2)-RR0(3)*RR(3))
       M(1,2)=M(1,2) +2.D0*(-RR0(2)*RR(3)+RR0(3)*RR(2))
       M(1,3)=M(1,3) +2.D0*(RR0(1)*RR(3)-RR0(3)*RR(1))
       M(1,4)=M(1,4) +2.D0*(-RR0(1)*RR(2)+RR0(2)*RR(1))
       M(2,3)=M(2,3) -2.D0*(RR0(1)*RR(2)+RR0(2)*RR(1))
       M(2,4)=M(2,4) -2.D0*(RR0(1)*RR(3)+RR0(3)*RR(1))
       M(3,4)=M(3,4) -2.D0*(RR0(2)*RR(3)+RR0(3)*RR(2))
    ENDDO
    M(2,1) = M(1,2)
    M(3,1) = M(1,3)
    M(3,2) = M(2,3)
    M(4,1) = M(1,4)
    M(4,2) = M(2,4)
    M(4,3) = M(3,4)
!   ---------------------------------------------------------------------    
!   SOLVE THE EIGENVECTOR PROBLEM FOR M
!   ---------------------------------------------------------------------
    IRET = 0
    JOBN   = 12
    CALL EIGRS(RESHAPE(M(:,:),SHAPE=(/16/)),4,JOBN,LAMBDA,Z,4,WK,IER)
    IF (IER .NE. 0)      STOP ' LSQQTN: FATAL ERROR IN EIGRS'
    IF (WK(1) .GT. 1.D0) IRET = 9
!   ---------------------------------------------------------------------    
!   PICK THE CORRECT EIGENVECTOR(S)
!   ---------------------------------------------------------------------
    S =  1.D0
    IF (Z(1,1) .LT. 0.D0) S = -1.D0
    Q(0) = S*Z(1,1)
    Q(1) = S*Z(2,1)
    Q(2) = S*Z(3,1)
    Q(3) = S*Z(4,1)
    IF (DABS(LAMBDA(1)) .LT. EPSI) THEN
       ERR = 0.D0
    ELSE
       ERR = LAMBDA(1)/DBLE(N)
    ENDIF
    IF (DABS(LAMBDA(1) - LAMBDA(2)) .LT. EPSI)IRET = IRET + 10
    
    IF (IRET .EQ. 0)  MSG = ' LSQQTN: NORMAL EXECUTION, UNIQUE SOLUTION'
    IF (IRET .EQ. 10) MSG = ' LSQQTN: NORMAL EXECUTION, NON-UNIQUE SOLUTION'
    IF (IRET .EQ. 9)  MSG = ' LSQQTN: BAD PERFORM. IN EIGRS, UNIQUE SOLUTION'
    IF (IRET .EQ. 19) MSG = ' LSQQTN: BAD PERFORM. IN EIGRS, NON-UNIQUE SOLUTION'
    IF(IRET.NE.0)WRITE(6,*)MSG
    
    IF (IOPT .EQ. 0) RETURN
!   ---------------------------------------------------------------------    
!   DERIVATIVES OF RMSD WITH RESPECT TO THE POSITIONS
!   ---------------------------------------------------------------------
    DO I = 1,N
       IF (W(I) .EQ. 0.D0) CYCLE
       RR(1)=W(I)*2.*RP(1,I)
       RR(2)=W(I)*2.*RP(2,I)
       RR(3)=W(I)*2.*RP(3,I)
       RR0(1)=W(I)*2.*R0P(1,I)
       RR0(2)=W(I)*2.*R0P(2,I)
       RR0(3)=W(I)*2.*R0P(3,I)
!
       DM_R (1,1,1)=(RR(1)-RR0(1))
       DM_R (1,1,2)=(RR(2)-RR0(2))
       DM_R (1,1,3)=(RR(3)-RR0(3))
! 
       DM_R (1,2,1)=0.D0
       DM_R (1,2,2)= RR0(3)
       DM_R (1,2,3)=-RR0(2)
!
       DM_R (1,3,1)=-RR0(3)
       DM_R (1,3,2)= 0.D0
       DM_R (1,3,3)= RR0(1)
!
       DM_R (1,4,1)= RR0(2)
       DM_R (1,4,2)=-RR0(1)
       DM_R (1,4,3)= 0.D0
!
       DM_R (2,2,1)=(RR(1)-RR0(1))
       DM_R (2,2,2)=(RR(2)+RR0(2))
       DM_R (2,2,3)=(RR(3)+RR0(3))
!
       DM_R (2,3,1)=-RR0(2)
       DM_R (2,3,2)=-RR0(1)
       DM_R (2,3,3)= 0.D0
!
       DM_R (2,4,1)=-RR0(3)
       DM_R (2,4,2)= 0.D0
       DM_R (2,4,3)=-RR0(1)
!
       DM_R (3,3,1)=(RR(1)+RR0(1))
       DM_R (3,3,2)=(RR(2)-RR0(2))
       DM_R (3,3,3)=(RR(3)+RR0(3))
! 
       DM_R (3,4,1)=0.D0
       DM_R (3,4,2)=-RR0(3)
       DM_R (3,4,3)=-RR0(2)
!
       DM_R (4,4,1)=(RR(1)+RR0(1))
       DM_R (4,4,2)=(RR(2)+RR0(2))
       DM_R (4,4,3)=(RR(3)-RR0(3))
!
       DO IX=1,3
          DM_R(2,1,IX)=DM_R(1,2,IX)
          DM_R(3,1,IX)=DM_R(1,3,IX)
          DM_R(4,1,IX)=DM_R(1,4,IX)
          DM_R(3,2,IX)=DM_R(2,3,IX)
          DM_R(4,2,IX)=DM_R(2,4,IX)
          DM_R(4,3,IX)=DM_R(3,4,IX)
       ENDDO
!
       DO IX=1,3
          DERR_DR (IX,I)=0.D0
          DO K=1,4
             DO J=1,4
                DERR_DR (IX,I)=DERR_DR (IX,I)+Q(K-1)*Q(J-1)*DM_R (J,K,IX)
             ENDDO
          ENDDO
          DERR_DR (IX,I)=DERR_DR (IX,I)/DBLE(N)
       ENDDO
    ENDDO
    IF (IOPT .EQ. 1) RETURN
!   ---------------------------------------------------------------------    
!   ROTATION MATRIX IN TERMS OF QUATERNIONS
!   ---------------------------------------------------------------------
    ROT(1,1)=-2.D0*Q(2)**2-2.D0*Q(3)**2+1.D0
    ROT(1,2)=2.D0*(-Q(0)*Q(3)+Q(1)*Q(2))
    ROT(1,3)=2.D0*(Q(0)*Q(2)+Q(1)*Q(3))
    ROT(2,1)=2.D0*(Q(0)*Q(3)+Q(1)*Q(2))
    ROT(2,2)=-2.D0*Q(1)**2-2.D0*Q(3)**2+1.D0
    ROT(2,3)=2.D0*(-Q(0)*Q(1)+Q(2)*Q(3))
    ROT(3,1)=2.D0*(-Q(0)*Q(2)+Q(1)*Q(3))
    ROT(3,2)=2.D0*(Q(0)*Q(1)+Q(2)*Q(3))
    ROT(3,3)=-2.D0*Q(1)**2-2.D0*Q(2)**2+1.D0


!   ---------------------------------------------------------------------    
!   CALCULATE FIT
!   ---------------------------------------------------------------------
    DO I = 1,N
       RFIT(1,I) = ROT(1,1)*R0P(1,I) + ROT(1,2)*R0P(2,I) + ROT(1,3)*R0P(3,I)
       RFIT(2,I) = ROT(2,1)*R0P(1,I) + ROT(2,2)*R0P(2,I) + ROT(2,3)*R0P(3,I)
       RFIT(3,I) = ROT(3,1)*R0P(1,I) + ROT(3,2)*R0P(2,I) + ROT(3,3)*R0P(3,I)
    ENDDO


 
    RETURN
  END SUBROUTINE RMSD3


 
  SUBROUTINE EIGRS  ( A, N, JOBN, D, Z, IZ, WK, IER )
    IMPLICIT NONE
!
!   IMSL ROUTINE NAME   - EIGRS
!
!   --------------------------------------------------------------------
!
!   COMPUTER            - IBM77/DOUBLE
!
!   LATEST REVISION     - JUNE 1, 1980
!
!   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
!                         A REAL SYMMETRIC MATRIX
!
!   USAGE               - CALL EIGRS (A,N,JOBN,D,Z,IZ,WK,IER)
!
!   ARGUMENTS    A      - INPUT REAL SYMMETRIC MATRIX OF ORDER N,
!                         WHOSE EIGENVALUES AND EIGENVECTORS
!                         ARE TO BE COMPUTED. INPUT A IS
!                         DESTROYED IF IJOB IS EQUAL TO 0 OR 1.
!                N      - INPUT ORDER OF THE MATRIX A.
!                JOBN   - INPUT OPTION PARAMETER.  IF JOBN.GE.10
!                         A IS ASSUMED TO BE IN FULL STORAGE MODE
!                         (IN THIS CASE, A MUST BE DIMENSIONED EXACTLY
!                         N BY N IN THE CALLING PROGRAM).
!                         IF JOBN.LT.10 THEN A IS ASSUMED TO BE IN
!                         SYMMETRIC STORAGE MODE.  DEFINE
!                         IJOB=MOD(JOBN,10).  THEN WHEN
!                         IJOB = 0, COMPUTE EIGENVALUES ONLY
!                         IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
!                         VECTORS.
!                         IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
!                         AND PERFORMANCE INDEX.
!                         IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
!                         IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
!                         RETURNED IN WK(1). THE ROUTINES HAVE
!                         PERFORMED (WELL, SATISFACTORILY, POORLY) IF
!                         WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
!                         GREATER THAN 100).
!                D      - OUTPUT VECTOR OF LENGTH N,
!                         CONTAINING THE EIGENVALUES OF A.
!                Z      - OUTPUT N BY N MATRIX CONTAINING
!                         THE EIGENVECTORS OF A.
!                         THE EIGENVECTOR IN COLUMN J OF Z CORRES-
!                         PONDS TO THE EIGENVALUE D(J).
!                         IF IJOB = 0, Z IS NOT USED.
!                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
!                         SPECIFIED IN THE DIMENSION STATEMENT IN THE
!                         CALLING PROGRAM.
!                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
!                         ON THE VALUE OF IJOB, WHEN
!                         IJOB = 0, THE LENGTH OF WK IS AT LEAST N.
!                         IJOB = 1, THE LENGTH OF WK IS AT LEAST N.
!                         IJOB = 2, THE LENGTH OF WK IS AT LEAST
!                         N(N+1)/2+N.
!                         IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
!                IER    - ERROR PARAMETER (OUTPUT)
!                         TERMINAL ERROR
!                         IER = 128+J, INDICATES THAT EQRT2S FAILED
!                         TO CONVERGE ON EIGENVALUE J. EIGENVALUES
!                         AND EIGENVECTORS 1,...,J-1 HAVE BEEN
!                         COMPUTED CORRECTLY, BUT THE EIGENVALUES
!                         ARE UNORDERED. THE PERFORMANCE INDEX
!                         IS SET TO 1000.0
!                         WARNING ERROR (WITH FIX)
!                         IN THE FOLLOWING, IJOB = MOD(JOBN,10).
!                         IER = 66, INDICATES IJOB IS LESS THAN 0 OR
!                         IJOB IS GREATER THAN 3. IJOB SET TO 1.
!                         IER = 67, INDICATES IJOB IS NOT EQUAL TO
!                         ZERO, AND IZ IS LESS THAN THE ORDER OF
!                         MATRIX A. IJOB IS SET TO ZERO.
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - EHOBKS,EHOUSS,EQRT2S,UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                         CONVENTIONS IS AVAILABLE IN THE MANUAL
!                         INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                         APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                         EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!   --------------------------------------------------------------------
!   SPECIFICATIONS FOR ARGUMENTS
!   --------------------------------------------------------------------
    INTEGER            :: N, JOBN, IZ, IER
    DOUBLE PRECISION   :: A(*), D(*), WK(*), Z(IZ,*)
!   --------------------------------------------------------------------
!   SPECIFICATIONS FOR LOCAL VARIABLES
!   --------------------------------------------------------------------
    INTEGER            :: IJOB, IR, JR, IJ, JI, NP1
    INTEGER            :: JER, NA, ND, IIZ, IBEG, IL, KK, LK, I, J, K, L
    DOUBLE PRECISION   :: ANORM, ASUM, PI, SUMZ, SUMR, AN, S
    DOUBLE PRECISION   :: TEN, RDELP, ZERO, ONE, THOUS
    DATA  RDELP / 0.222045D-15 /
    DATA  ZERO, ONE / 0.0D0 , 1.0D0 /, TEN / 10.0D0 /, THOUS / 1000.0D0/
!   --------------------------------------------------------------------
!   INITIALIZE ERROR PARAMETERS
!   FIRST EXECUTABLE STATEMENT
!   --------------------------------------------------------------------
    IER = 0
    JER = 0
    IF (JOBN.LT.10) GO TO 15
!   --------------------------------------------------------------------
!   CONVERT TO SYMMETRIC STORAGE MODE
!   --------------------------------------------------------------------
    K = 1
    JI = N-1
    IJ = 1
    DO 10 J=1,N
       DO 5 I=1,J
          A(K) = A(IJ)
          IJ = IJ+1
          K = K+1
5      CONTINUE
       IJ = IJ + JI
       JI = JI - 1
10  CONTINUE
15  IJOB = MOD(JOBN,10)
    IF (IJOB.GE.0.AND.IJOB.LE.3) GO TO 20
!   --------------------------------------------------------------------
!   WARNING ERROR - IJOB IS NOT IN THE RANGE
!   --------------------------------------------------------------------
    IER = 66
    IJOB = 1
    GO TO 25
20  IF (IJOB.EQ.0) GO TO 35
25  IF (IZ.GE.N)   GO TO 30
!   --------------------------------------------------------------------
!   WARNING ERROR - IZ IS LESS THAN N
!   EIGENVECTORS CAN NOT BE COMPUTED, IJOB SET TO ZERO
!   --------------------------------------------------------------------
    IER = 67
    IJOB = 0
30  IF (IJOB.EQ.3) GO TO 75
35  NA = (N*(N+1))/2
    IF (IJOB.NE.2) GO TO 45
    DO 40 I=1,NA
       WK(I) = A(I)
40  CONTINUE
!   --------------------------------------------------------------------
!   SAVE INPUT A IF IJOB = 2
!   --------------------------------------------------------------------
45  ND = 1
    IF (IJOB.EQ.2) ND = NA+1
!   --------------------------------------------------------------------
!   REDUCE A TO SYMMETRIC TRIDIAGONAL FORM
!   --------------------------------------------------------------------
    CALL EHOUSS (A,N,D,WK(ND),WK(ND))
    IIZ = 1
    IF (IJOB.EQ.0) GO TO 60
    IIZ = IZ
!   --------------------------------------------------------------------
!   SET Z TO THE IDENTITY MATRIX
!   --------------------------------------------------------------------
    DO 55 I=1,N
       DO 50 J=1,N
          Z(I,J) = ZERO
50     CONTINUE
       Z(I,I) = ONE
55  CONTINUE
!   --------------------------------------------------------------------
!   COMPUTE EIGENVALUES AND EIGENVECTORS
!   --------------------------------------------------------------------
60  CALL EQRT2S (D,WK(ND),N,Z,IIZ,JER)
    IF (IJOB.EQ.0)  GO TO 9000
    IF (JER.GT.128) GO TO 65
!   --------------------------------------------------------------------
!   BACK TRANSFORM EIGENVECTORS
!   --------------------------------------------------------------------
    CALL EHOBKS (A,N,1,N,Z,IZ)
65  IF (IJOB.LE.1) GO TO 9000
!   --------------------------------------------------------------------
!   MOVE INPUT MATRIX BACK TO A
!   --------------------------------------------------------------------
    DO 70 I=1,NA
       A(I) = WK(I)
70  CONTINUE
    WK(1) = THOUS
    IF (JER.NE.0) GO TO 9000
!   --------------------------------------------------------------------
!   COMPUTE 1 - NORM OF A
!   --------------------------------------------------------------------
75  ANORM = ZERO
    IBEG = 1
    DO 85 I=1,N
       ASUM = ZERO
       IL = IBEG
       KK = 1
       DO 80 L=1,N
          ASUM = ASUM+DABS(A(IL))
          IF (L.GE.I) KK = L
          IL = IL+KK
80     CONTINUE
       ANORM = DMAX1(ANORM,ASUM)
       IBEG = IBEG+I
85  CONTINUE
    IF (ANORM.EQ.ZERO) ANORM = ONE
!   --------------------------------------------------------------------
!   COMPUTE PERFORMANCE INDEX
!   --------------------------------------------------------------------
    PI = ZERO
    DO 100 I=1,N
       IBEG = 1
       S = ZERO
       SUMZ = ZERO
       DO 95 L=1,N
          LK = IBEG
          KK = 1
          SUMZ = SUMZ+DABS(Z(L,I))
          SUMR = -D(I)*Z(L,I)
          DO 90 K=1,N
             SUMR = SUMR+A(LK)*Z(K,I)
             IF (K.GE.L) KK = K
             LK = LK+KK
90        CONTINUE
          S = S+DABS(SUMR)
          IBEG = IBEG+L
95     CONTINUE
       IF (SUMZ.EQ.ZERO) GO TO 100
       PI = DMAX1(PI,S/SUMZ)
100 CONTINUE
    AN = N
    PI = PI/(ANORM*TEN*AN*RDELP)
    WK(1) = PI
    IF (JOBN.LT.10) GO TO 9000
!   --------------------------------------------------------------------
!   CONVERT BACK TO FULL STORAGE MODE
!   --------------------------------------------------------------------
    NP1 = N+1
    IJ = (N-1)*NP1 + 2
    K = (N*(NP1))/2
    DO 110 JR=1,N
       J = NP1-JR
       DO 105 IR=1,J
          IJ = IJ-1
          A(IJ) = A(K)
          K = K-1
105    CONTINUE
       IJ = IJ-JR
110 CONTINUE
    JI = 0
    K = N-1
    DO 120 I=1,N
       IJ = I-N
       DO 115 J=1,I
          IJ = IJ+N
          JI = JI+1
          A(IJ) = A(JI)
115    CONTINUE
       JI = JI + K
       K = K-1
120 CONTINUE

9000  &
    CONTINUE
    IF (IER.NE.0) CALL UERTST (IER,'EIGRS ')
    IF (JER.EQ.0) GO TO 9005
    IER = JER
    CALL UERTST (IER,'EIGRS ')
9005  & 
    RETURN
  END SUBROUTINE EIGRS


  SUBROUTINE EHOBKS (A, N, M1, M2, Z, IZ) 
    IMPLICIT NONE
!   IMSL ROUTINE NAME   - EHOBKS
!
!   --------------------------------------------------------------------
!
!   COMPUTER            - IBM77/DOUBLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRS
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!   --------------------------------------------------------------------
!
    INTEGER, INTENT(IN) :: IZ, N
    INTEGER :: M1, M2, I, L, IA, J, K
    DOUBLE PRECISION :: A(*),Z(IZ,*),H, S
!   --------------------------------------------------------------------
!   FIRST EXECUTABLE STATEMENT
!   --------------------------------------------------------------------
    IF (N .EQ. 1) GO TO 30
    DO 25 I=2,N
       L = I-1
       IA = (I*L)/2
       H = A(IA+I)
       IF (H.EQ.0.D0) GO TO 25
!   --------------------------------------------------------------------
!   DERIVES EIGENVECTORS M1 TO M2 OF THE ORIGINAL MATRIX FROM 
!   EIGENVECTORS M1 TO M2 OF THE SYMMETRIC TRIDIAGONAL MATRIX
!   --------------------------------------------------------------------
       DO 20 J = M1,M2
          S = 0.0D0
          DO 10 K = 1,L
             S = S+A(IA+K)*Z(K,J)
10        CONTINUE
          S = S/H
          DO 15 K=1,L
             Z(K,J) = Z(K,J)-S*A(IA+K)
15        CONTINUE
20     CONTINUE
25  CONTINUE
30  RETURN
    END SUBROUTINE EHOBKS  
 

    SUBROUTINE EHOUSS (A, N, D, E, E2)
      IMPLICIT NONE
!   IMSL ROUTINE NAME   - EHOUSS
!
!   --------------------------------------------------------------------
!
!   COMPUTER            - IBM77/DOUBLE
!
!   LATEST REVISION     - NOVEMBER 1, 1984
!
!   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRS
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!   --------------------------------------------------------------------
!
      INTEGER, INTENT(IN) :: N
      INTEGER  :: NP1, NN, NBEG, II, I, L, NK, K, JK1, J, IK, JK, JP1
      DOUBLE PRECISION ::  A(*), D(N), E(N), E2(N), ZERO 
      DOUBLE PRECISION ::  H, SCALE, F, G, HH
      DATA  ZERO / 0.0D0 /
!   --------------------------------------------------------------------
!   FIRST EXECUTABLE STATEMENT
!   --------------------------------------------------------------------
      NP1 = N+1
      NN = (N*NP1)/2-1
      NBEG = NN+1-N
      DO 70 II = 1,N
         I = NP1-II
         L = I-1
         H = ZERO
         SCALE = ZERO
         IF (L .LT. 1) GO TO 10
!   --------------------------------------------------------------------
!   SCALE ROW (ALGOL TOL THEN NOT NEEDED)
!   --------------------------------------------------------------------
         NK = NN
         DO 5 K = 1,L
            SCALE = SCALE+DABS(A(NK))
            NK = NK-1
    5    CONTINUE
         IF (SCALE .NE. ZERO) GO TO 15
   10    E(I) = ZERO
         E2(I) = ZERO
         GO TO 65
   15    NK = NN
         DO 20 K = 1,L
            A(NK) = A(NK)/SCALE
            H = H+A(NK)*A(NK)
            NK = NK-1
   20    CONTINUE
         E2(I) = SCALE*SCALE*H
         F = A(NN)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE*G
         H = H-F*G
         A(NN) = F-G
         IF (L .EQ. 1) GO TO 55
         F = ZERO
         JK1 = 1
         DO 40 J = 1,L
            G = ZERO
            IK = NBEG+1
            JK = JK1
!   --------------------------------------------------------------------
!   FORM ELEMENT OF A*U
!   --------------------------------------------------------------------
            DO 25 K = 1,J
               G = G+A(JK)*A(IK)
               JK = JK+1
               IK = IK+1
   25       CONTINUE
            JP1 = J+1
            IF (L .LT. JP1) GO TO 35
            JK = JK+J-1
            DO 30 K = JP1,L
               G = G+A(JK)*A(IK)
               JK = JK+K
               IK = IK+1
   30       CONTINUE
!   --------------------------------------------------------------------
!   FORM ELEMENT OF P
!   --------------------------------------------------------------------
   35       E(J) = G/H
            F = F+E(J)*A(NBEG+J)
            JK1 = JK1+J
   40    CONTINUE
         HH = F/(H+H)
!   --------------------------------------------------------------------
!   FORM REDUCED A
!   --------------------------------------------------------------------
         JK = 1
         DO 50 J = 1,L
            F = A(NBEG+J)
            G = E(J)-HH*F
            E(J) = G
            DO 45 K = 1,J
               A(JK) = A(JK)-F*E(K)-G*A(NBEG+K)
               JK = JK+1
   45       CONTINUE
   50    CONTINUE
   55    DO 60 K = 1,L
            A(NBEG+K) = SCALE*A(NBEG+K)
   60    CONTINUE
   65    D(I) = A(NBEG+I)
         A(NBEG+I) = H*SCALE*SCALE
         NBEG = NBEG-I+1
         NN = NN-I
   70 CONTINUE
   RETURN
   END SUBROUTINE EHOUSS 


   SUBROUTINE EQRT2S ( D, E, N, Z, IZ, IER) 
!   IMSL ROUTINE NAME   - EQRT2S
!
!   --------------------------------------------------------------------
!
!   COMPUTER            - IBM77/DOUBLE
!
!   LATEST REVISION     - NOVEMBER 1, 1984
!
!   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
!                         A SYMMETRIC TRIDIAGONAL MATRIX USING THE
!                         QL METHOD.
!
!   USAGE               - CALL EQRT2S (D,E,N,Z,IZ,IER)
!
!   ARGUMENTS    D      - ON INPUT, THE VECTOR D OF LENGTH N CONTAINS
!                         THE DIAGONAL ELEMENTS OF THE SYMMETRIC
!                         TRIDIAGONAL MATRIX T.
!                         ON OUTPUT, D CONTAINS THE EIGENVALUES OF
!                         T IN ASCENDING ORDER.
!                E      - ON INPUT, THE VECTOR E OF LENGTH N CONTAINS
!                         THE SUB-DIAGONAL ELEMENTS OF T IN POSITION
!                         2,...,N. ON OUTPUT, E IS DESTROYED.
!                N      - ORDER OF TRIDIAGONAL MATRIX T.(INPUT)
!                Z      - ON INPUT, Z CONTAINS THE IDENTITY MATRIX OF
!                         ORDER N.
!                         ON OUTPUT, Z CONTAINS THE EIGENVECTORS
!                         OF T. THE EIGENVECTOR IN COLUMN J OF Z
!                         CORRESPONDS TO THE EIGENVALUE D(J).
!                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
!                         SPECIFIED IN THE DIMENSION STATEMENT IN THE
!                         CALLING PROGRAM. IF IZ IS LESS THAN N, THE
!                         EIGENVECTORS ARE NOT COMPUTED. IN THIS CASE
!                         Z IS NOT USED.
!                IER    - ERROR PARAMETER
!                         TERMINAL ERROR
!                         IER = 128+J, INDICATES THAT EQRT2S FAILED
!                         TO CONVERGE ON EIGENVALUE J. EIGENVALUES
!                         AND EIGENVECTORS 1,...,J-1 HAVE BEEN
!                         COMPUTED CORRECTLY, BUT THE EIGENVALUES
!                         ARE UNORDERED.
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                         CONVENTIONS IS AVAILABLE IN THE MANUAL
!                         INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                         APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                         EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!   --------------------------------------------------------------------
     INTEGER, INTENT(IN) :: N, IZ
     INTEGER :: IER, I, L, J, M, K, L1, MM1, MM1PL, II, IP1
     DOUBLE PRECISION ::  D(*), E(*), Z(IZ,*)
     DOUBLE PRECISION ::  B, C, F, G, H, P, R, S, RDELP, ONE, ZERO
     DATA   RDELP / 0.222045D-15 /
     DATA   ZERO,ONE /0.0D0, 1.0D0/
!   --------------------------------------------------------------------
!   MOVE THE LAST N-1 ELEMENTS OF E INTO THE FIRST N-1 LOCATIONS
!   FIRST EXECUTABLE STATEMENT
!   --------------------------------------------------------------------
     IER  = 0
     IF (N .EQ. 1) GO TO 9005
     DO 5  I=2,N
        E(I-1) = E(I)
5    CONTINUE
     E(N) = ZERO
     B = ZERO
     F = ZERO
     DO  60  L=1,N
         J = 0
         H = RDELP*(DABS(D(L))+DABS(E(L)))
         IF (B.LT.H) B = H
!   --------------------------------------------------------------------
!   LOOK FOR SMALL SUB-DIAGONAL ELEMENT
!   --------------------------------------------------------------------
         DO 10  M=L,N
            K=M
            IF (DABS(E(K)) .LE. B) GO TO 15
   10    CONTINUE
   15    M = K
         IF (M.EQ.L) GO TO 55
   20    IF (J .EQ. 30) GO TO 85
         J = J+1
         L1 = L+1
         G = D(L)
         P = (D(L1)-G)/(E(L)+E(L))
         R = DABS(P)
         IF (RDELP*DABS(P) .LT. 1.0D0) R = DSQRT(P*P+ONE)
         D(L) = E(L)/(P+DSIGN(R,P))
         H = G-D(L)
         DO 25 I = L1,N
            D(I) = D(I)-H
   25    CONTINUE
         F = F+H
!   --------------------------------------------------------------------
!   QL TRANSFORMATION
!   --------------------------------------------------------------------
         P = D(M)
         C = ONE
         S = ZERO
         MM1 = M-1
         MM1PL = MM1+L
         IF (L.GT.MM1) GO TO 50
         DO 45 II=L,MM1
            I = MM1PL-II
            G = C*E(I)
            H = C*P
            IF (DABS(P).LT.DABS(E(I))) GO TO 30
            C = E(I)/P
            R = DSQRT(C*C+ONE)
            E(I+1) = S*P*R
            S = C/R
            C = ONE/R
            GO TO 35
   30       C = P/E(I)
            R = DSQRT(C*C+ONE)
            E(I+1) = S*E(I)*R
            S = ONE/R
            C = C*S
   35       P = C*D(I)-S*G
            D(I+1) = H+S*(C*G+S*D(I))
            IF (IZ .LT. N) GO TO 45
!   --------------------------------------------------------------------
!   FORM VECTOR
!   --------------------------------------------------------------------
            DO 40 K=1,N
               H = Z(K,I+1)
               Z(K,I+1) = S*Z(K,I)+C*H
               Z(K,I) = C*Z(K,I)-S*H
   40       CONTINUE
   45    CONTINUE
   50    E(L) = S*P
         D(L) = C*P
         IF (DABS(E(L)) .GT.B) GO TO 20
   55    D(L) = D(L) + F
60   CONTINUE
!   --------------------------------------------------------------------
!   ORDER EIGENVALUES AND EIGENVECTORS
!   --------------------------------------------------------------------
      DO  80  I=1,N
         K = I
         P = D(I)
         IP1 = I+1
         IF (IP1.GT.N) GO TO 70
         DO 65  J=IP1,N
            IF (D(J) .GE. P) GO TO 65
            K = J
            P = D(J)
   65    CONTINUE
   70    IF (K.EQ.I) GO TO 80
         D(K) = D(I)
         D(I) = P
         IF (IZ .LT. N) GO TO 80
         DO 75 J = 1,N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
   75    CONTINUE
   80 CONTINUE
      GO TO 9005
   85 IER = 128+L
      CALL UERTST(IER,'EQRT2S')
 9005 RETURN
   END SUBROUTINE EQRT2S 

      
   SUBROUTINE UERTST (IER, NAME) 
     IMPLICIT NONE
!   IMSL ROUTINE NAME   - UERTST
!
!   --------------------------------------------------------------------
!
!   COMPUTER            - IBM77/SINGLE
!
!   LATEST REVISION     - MARCH 26, 1982
!
!   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION
!
!   USAGE               - CALL UERTST (IER,NAME)
!
!   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
!                           IER = I+J WHERE
!                             I = 128 IMPLIES TERMINAL ERROR MESSAGE,
!                             I =  64 IMPLIES WARNING WITH FIX MESSAGE,
!                             I =  32 IMPLIES WARNING MESSAGE.
!                             J = ERROR CODE RELEVANT TO CALLING
!                                 ROUTINE.
!                NAME   - A CHARACTER STRING OF LENGTH SIX PROVIDING
!                           THE NAME OF THE CALLING ROUTINE. (INPUT)
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - UGETIO,USPKD
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
!                TO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
!                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
!                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
!                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
!                UGETIO AS FOLLOWS..
!                                NIN = 0
!                                NOUT = NEW OUTPUT UNIT NUMBER
!                                CALL UGETIO(3,NIN,NOUT)
!                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!   --------------------------------------------------------------------
!   SPECIFICATIONS FOR ARGUMENTS
     INTEGER, INTENT(INOUT) :: IER
     CHARACTER (LEN=*)   :: NAME
!   SPECIFICATIONS FOR LOCAL VARIABLES
     INTEGER   ::  I,IEQDF,IOUNIT,LEVEL,LEVOLD,NIN,NMTB
     CHARACTER ::  IEQ,NAMEQ(6),NAMSET(6),NAMUPK(6)
     DATA  NAMSET / 'U', 'E', 'R', 'S', 'E', 'T'/
     DATA  NAMEQ / 6*' ' /
     DATA  LEVEL / 4 /, IEQDF / 0 /, IEQ / '=' /
!   --------------------------------------------------------------------
!   UNPACK NAME INTO NAMUPK FIRST EXECUTABLE STATEMENT
!   --------------------------------------------------------------------
     CALL USPKD (NAME,6,NAMUPK,NMTB)
!   --------------------------------------------------------------------
!   GET OUTPUT UNIT NUMBER
!   --------------------------------------------------------------------
     CALL UGETIO(1,NIN,IOUNIT)
!   --------------------------------------------------------------------
!   CHECK IER
!   --------------------------------------------------------------------
     IF (IER.GT.999) GO TO 25
     IF (IER.LT.-32) GO TO 55
     IF (IER.LE.128) GO TO 5
     IF (LEVEL.LT.1) GO TO 30
!   --------------------------------------------------------------------
!   PRINT TERMINAL MESSAGE
!   --------------------------------------------------------------------
     IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAMUPK
     IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAMUPK
     GO TO 30
5    IF (IER.LE.64) GO TO 10
     IF (LEVEL.LT.2) GO TO 30
!   --------------------------------------------------------------------
!   PRINT WARNING WITH FIX MESSAGE
!   --------------------------------------------------------------------
     IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAMUPK
     IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAMUPK
     GO TO 30
10   IF (IER.LE.32) GO TO 15
!   --------------------------------------------------------------------
!   PRINT WARNING MESSAGE
!   --------------------------------------------------------------------
     IF (LEVEL.LT.3) GO TO 30
     IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAMUPK
     IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAMUPK
     GO TO 30
15   CONTINUE
!   --------------------------------------------------------------------
!   CHECK FOR UERSET CALL
!   --------------------------------------------------------------------
     DO 20 I=1,6
        IF (NAMUPK(I).NE.NAMSET(I)) GO TO 25
20   CONTINUE
     LEVOLD = LEVEL
     LEVEL = IER
     IER = LEVOLD
     IF (LEVEL.LT.0) LEVEL = 4
     IF (LEVEL.GT.4) LEVEL = 4
     GO TO 30
25   CONTINUE
     IF (LEVEL.LT.4) GO TO 30
!   --------------------------------------------------------------------
!   PRINT NON-DEFINED MESSAGE
!   --------------------------------------------------------------------
     IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAMUPK
     IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAMUPK
30   IEQDF = 0
     RETURN
35   FORMAT(" *** TERMINAL ERROR",10X,"(IER = ",I3,") FROM IMSL ROUTINE ",6A1,A1,6A1)
40   FORMAT(" *** WARNING WITH FIX ERROR",2X,"(IER = ",I3,") FROM IMSL ROUTINE ",6A1,A1,6A1)
45   FORMAT(" *** WARNING ERROR",11X,"(IER = ",I3,") FROM IMSL ROUTINE ",6A1,A1,6A1)
50   FORMAT(" *** UNDEFINED ERROR",9X,"(IER = ",I5,") FROM IMSL ROUTINE ",6A1,A1,6A1)
!   --------------------------------------------------------------------
!   SAVE P FOR P = R CASE
!   P IS THE PAGE NAMUPK    -     R IS THE ROUTINE NAMUPK
!   --------------------------------------------------------------------
55   IEQDF = 1
     DO  I=1,6
        NAMEQ(I) = NAMUPK(I)
     END DO
     RETURN
   END SUBROUTINE UERTST 

   SUBROUTINE UGETIO(IOPT, NIN, NOUT) 
     IMPLICIT NONE

!   IMSL ROUTINE NAME   - UGETIO
!
!   --------------------------------------------------------------------
!
!   COMPUTER            - IBM77/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1981
!
!   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
!                           VALUES FOR INPUT AND OUTPUT UNIT
!                           IDENTIFIERS.
!
!   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
!
!   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
!                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT
!                           UNIT IDENTIFIER VALUES ARE RETURNED IN NIN
!                           AND NOUT, RESPECTIVELY.
!                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS
!                           RESET FOR SUBSEQUENT USE.
!                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS
!                           RESET FOR SUBSEQUENT USE.
!                NIN    - INPUT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.
!                NOUT   - OUTPUT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
!                OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT
!                IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR
!                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.
!                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!   --------------------------------------------------------------------
     INTEGER :: IOPT, NIN, NOUT
     INTEGER :: NIND, NOUTD
     DATA   NIND / 5 /, NOUTD / 6 /
     IF (IOPT.EQ.3) GO TO 10
     IF (IOPT.EQ.2) GO TO 5
     IF (IOPT.NE.1) GO TO 9005
     NIN = NIND
     NOUT = NOUTD
     GO TO 9005
5    NIND = NIN
     GO TO 9005
10   NOUTD = NOUT
9005 RETURN
   END SUBROUTINE UGETIO
 
   SUBROUTINE USPKD  (PACKED, NCHARS, UNPAKD, NCHMTB)
     IMPLICIT NONE
!   IMSL ROUTINE NAME   - USPKD
!
!   --------------------------------------------------------------------
!
!   COMPUTER            - IBM77/SINGLE
!
!   LATEST REVISION     - NOVEMBER 1, 1984
!
!   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES THAT HAVE
!                           CHARACTER STRING ARGUMENTS
!
!   USAGE               - CALL USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
!
!   ARGUMENTS    PACKED - CHARACTER STRING TO BE UNPACKED.(INPUT)
!                NCHARS - LENGTH OF PACKED. (INPUT)  SEE REMARKS.
!                UNPAKD - CHARACTER ARRAY TO RECEIVE THE UNPACKED
!                         REPRESENTATION OF THE STRING. (OUTPUT)
!                NCHMTB - NCHARS MINUS TRAILING BLANKS. (OUTPUT)
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - NONE
!
!   REMARKS  1.  USPKD UNPACKS A CHARACTER STRING INTO A CHARACTER ARRAY
!                IN (A1) FORMAT.
!            2.  UP TO 129 CHARACTERS MAY BE USED.  ANY IN EXCESS OF
!                THAT ARE IGNORED.
!
!   COPYRIGHT           - 1984 BY IMSL, INC.  ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE.  NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!   --------------------------------------------------------------------
     INTEGER :: NC, NCHARS, NCHMTB
     INTEGER :: I, N, NN
     CHARACTER   :: UNPAKD(1),IBLANK
     CHARACTER  (LEN=*) :: PACKED
      DATA IBLANK / ' ' /
      NCHMTB = 0
!   --------------------------------------------------------------------
!   RETURN IF NCHARS IS LE ZERO
!   --------------------------------------------------------------------
      IF(NCHARS.LE.0) RETURN
!   --------------------------------------------------------------------
!   SET NC=NUMBER OF CHARS TO BE DECODED
!   --------------------------------------------------------------------
      NC = MIN0 (129,NCHARS)
      READ (PACKED,150) (UNPAKD(I),I=1,NC)
150   FORMAT (129A1)
!   --------------------------------------------------------------------
!   CHECK UNPAKD ARRAY AND SET NCHMTB BASED ON TRAILING BLANKS FOUND
!   --------------------------------------------------------------------
      DO 200 N = 1,NC
         NN = NC - N + 1
         IF(UNPAKD(NN) .NE. IBLANK) GO TO 210
  200 CONTINUE
      NN = 0
210   NCHMTB = NN
      RETURN
   END SUBROUTINE USPKD

END MODULE SUPERPOSITION



MODULE RMSD3_INTERFACE
  USE SUPERPOSITION, ONLY: RMSD3
  USE XYZ, ONLY : POINT
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(3,3) :: ROT
  DOUBLE PRECISION, DIMENSION(3)   :: TRASL
  DOUBLE PRECISION  :: ERR

CONTAINS
  SUBROUTINE OVERLAP_TYPEP ( N, C0, C1, W, MASS, IOPT)
    IMPLICIT NONE
    INTEGER (KIND=8), INTENT(IN)  :: N
    INTEGER, INTENT(IN) :: IOPT
    DOUBLE PRECISION, ALLOCATABLE :: R0 (:,:), R1(:,:), DERR_DR(:,:)
    TYPE(POINT), INTENT(IN) :: C0(N), C1(N)
    DOUBLE PRECISION, INTENT(IN) :: W(N), MASS(N)
    INTEGER :: I, J

    ALLOCATE ( R0(3,N),  R1(3,N), DERR_DR(3,N) )
    DO I=1,3
       DO J=1,3
          ROT(I,J)=0.D0
       END DO
       TRASL(I)=0.D0
    END DO
    ERR = 0.D0

    DO I=1,N
       R0(1,I)=C0(I)%X
       R0(2,I)=C0(I)%Y
       R0(3,I)=C0(I)%Z
       R1(1,I)=C1(I)%X
       R1(2,I)=C1(I)%Y
       R1(3,I)=C1(I)%Z
    END DO

    CALL RMSD3( INT(N,4), R1, R0, TRASL, ROT, ERR, DERR_DR, W, MASS, IOPT)

  END SUBROUTINE OVERLAP_TYPEP

END MODULE RMSD3_INTERFACE
