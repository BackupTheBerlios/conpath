
! Module to perform molecular dynamics in internal coordinates...
! See for theoretical details : 
! P. Pulay, B. Paizs, Chem. Phys. Lett. 353 (2002) 400-406
!
! T. Laino - NEST (INFM) - Scuola Normale Superiore
! ETH Zurich
MODULE DYNA_ZMAT
  USE DYNPREPARE, ONLY:   NA,                  &
                          NB,                  &
                          NC,                  &
                          LOPT,                &
                          GEO,                 &
                          CX,                  &
                          VX,                  &
                          FX,                  &
                          NUM_BONDS,           &
                          NUM_ANGLES,          &
                          NUM_TORSIONS         
  USE START_JOB,   ONLY:  NUMAT,               &
                          ABINIT_INP,          &
                          NSTEP,               &
                          IPRINT,              &
                          TITLE,               &
                          INI_TEMP,            &
                          AREWERESTARTING
  USE CONVFACTORS,  ONLY: TCONV => FSTOAUTIME
  USE UNITS_FILE,   ONLY: FLUSH_FILES
  USE DYNPREPARE,   ONLY: C,                   &                    
                          F,                   &
                          V,                   &
                          M,                   &
                          LAB
  IMPLICIT NONE
CONTAINS

  SUBROUTINE DYNAMICS_ZMAT
     USE START_JOB,   ONLY: DTORIG => TSTEP
     USE FORCE,       ONLY: FORCE_CART
     USE MODKINETIC,  ONLY: KINETIC
     USE CONVFACTORS, ONLY:  ANGTOBOHR
!------------------------------------------------------------------------!
!     FORTRAN PROGRAM TO CONDUCT MOLECULAR DYNAMICS OF ATOMS             !
!     IN INTERNAL COORDINATES                                            !
!                                                                        !
!     REFERENCE:                                                         !
!                                                                        !
!     P. PULAY, B. PAIZS, CHEM. PHYS. LETT., 353 (2002)  400 - 406       !
!------------------------------------------------------------------------!
    IMPLICIT NONE
    ! Local Variables
    DOUBLE PRECISION, ALLOCATABLE  :: QINTN(:), QINTO(:),  QINTV(:)  
    DOUBLE PRECISION, ALLOCATABLE  :: MASS(:,:)
    INTEGER                        :: NUM_VAR
    INTEGER                        :: I, J, L
    INTEGER                        :: RSTEP, STEP
    DOUBLE PRECISION               :: DT, PI
    DOUBLE PRECISION               :: K, E, UPOT, VEXC, VFOND

    DT     = DTORIG * TCONV
!   -------------------------------------------------------------------
!   Zero Accumulators
!   -------------------------------------------------------------------
    STEP = 0
    IF ( IPRINT .LE. 0 ) IPRINT = NSTEP + 1
    NUM_VAR = NUM_BONDS + NUM_ANGLES + NUM_TORSIONS
    ALLOCATE ( QINTN(NUM_VAR), QINTO(NUM_VAR), QINTV(NUM_VAR) )
    ALLOCATE ( MASS(3*NUMAT,3*NUMAT) )
    !
    ! Initialize QINTN with Original values..
    !
    PI = 4.D0 * ATAN(1.D0)
    QINTN  = 0.D0
    QINTO  = 0.D0
    QINTV  = 0.D0
    MASS   = 0.D0
    I = 0
    DO L = 2, NUMAT
       I = I + 1
       GEO(1,L) = GEO(1,L)*ANGTOBOHR
       QINTN(I) = GEO(1,L)
       IF (L.EQ.2) CYCLE
       I = I + 1
       GEO(2,L) = GEO(2,L)*PI/180.D0
       QINTN(I) = GEO(2,L)
       IF (L.EQ.3) CYCLE
       I = I + 1
       GEO(3,L) = GEO(3,L)*PI/180.D0
       QINTN(I) = GEO(3,L)
    END DO
    !
    ! Masses Matrix
    !
    DO I = 1, NUMAT
       MASS((I-1)*3+1,(I-1)*3+1) = 1.D0/M(I)
       MASS((I-1)*3+2,(I-1)*3+2) = 1.D0/M(I)
       MASS((I-1)*3+3,(I-1)*3+3) = 1.D0/M(I)
    END DO
    !
    ! Write Header of Dynamics in Internal Coordinates
    !
    WRITE(6,001)NUM_BONDS, NUM_ANGLES, NUM_TORSIONS
    L = 0
    DO I = 2, NUMAT
       L = L + 1
       WRITE(6,002)L,QINTN(L)/ANGTOBOHR
       IF (I.EQ.2) CYCLE
       L = L + 1
       WRITE(6,003)L,QINTN(L)*180.D0/PI
       IF (I.EQ.3) CYCLE
       L = L + 1
       WRITE(6,004)L,QINTN(L)*180.D0/PI
    END DO
    WRITE(6,005)
    WRITE(6,'(//1X,A)') TITLE
    WRITE(6,'('' NUMBER OF ATOMS  = '',I10  )') NUMAT
    WRITE(6,'('' NUMBER OF STEPS  = '',I10  )') NSTEP
    WRITE(6,'('' OUTPUT FREQUENCY = '',I10  )') IPRINT
    WRITE(6,'('' TIME STEP        = '',F10.6,''  fs'')') DTORIG
    WRITE(6,'('' INITIAL NUCLEAR TEMPERATURE = '',F10.6)')INI_TEMP
    IF (AREWERESTARTING) THEN
       WRITE(6,'('' RESTARTING FROM  FILE! '')')
       STOP ! To implement this part
    ELSE
       RSTEP=0
    END IF
    CALL FLUSH_FILES
    
    CALL KINETIC (  M, V, K, NUMAT )
!   ---------------------------------------------------------------------
!   Calculate initial values 
!   ---------------------------------------------------------------------
    !
    CALL GET_XYZ_GEO
    CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP )
    ! 
    CALL PRINTZMATDYN( UPOT, K, STEP, DTORIG, VFOND, VEXC )

    CALL SINGLE_ZMAT_DYN_STEP(DT, VEXC, VFOND, UPOT, STEP,    & 
                              QINTN, QINTO, QINTV, MASS, NUM_VAR)
!   ---------------------------------------------------------------------
!   Compute Kinetic Energy
!   ---------------------------------------------------------------------
    CALL KINETIC (  M, V, K, NUMAT )  
    WRITE(*,'(//1X,''**** START OF DYNAMICS ****'')')
!   ---------------------------------------------------------------------
    Main: DO  STEP = RSTEP+1, NSTEP
!     -------------------------------------------------------------------
!     Optionally print information
!     -------------------------------------------------------------------
       IF ( MOD( STEP, INT(IPRINT) ) .EQ. 0 ) THEN
          CALL PRINTZMATDYN( UPOT, K, STEP, DTORIG, VFOND, VEXC )
       ENDIF
!     -------------------------------------------------------------------
!     Implementation Algorithm
!     -------------------------------------------------------------------
       CALL  SINGLE_ZMAT_DYN_STEP(DT, VEXC, VFOND, UPOT, STEP,             &
                                  QINTN, QINTO, QINTV, MASS, NUM_VAR)
!     -------------------------------------------------------------------
!     Calculate kinetic energy at current step
!     -------------------------------------------------------------------
       CALL KINETIC ( M, V, K, NUMAT )
    END DO Main
!   ---------------------------------------------------------------------
    RETURN

001 FORMAT("*************************************************************************", &
         /,"*                                                                       *", &
         /,"*  ZMAT DYNAMICAL MODULE:  PERFORMING DYNAMICS IN INTERNAL COORDINATES  *", &
         /,"*                                                                       *", &
         /,"*                                                                       *", &
         /,"*  Ref.:  P. Pulay, B. Paizs, Chem. Phys. Lett., 353 (2002)  400 - 406  *", &
         /,"*                                                                       *", &
         /,"*************************************************************************", &
         /,"*                                                                       *", &
         /,"* NUMBER OF BONDS    :",I6,"                                            *", &
         /,"* NUMBER OF ANGLES   :",I6,"                                            *", &
         /,"* NUMBER OF TORSIONS :",I6,"                                            *", &
         /,"*                                                                       *", &
         /,"*************************************************************************"  &
          )

002 FORMAT("* INTERNAL COORDINATE N.:",I3,"      VALUE:",F12.6," Angstrom * BOND    *"  )
003 FORMAT("* INTERNAL COORDINATE N.:",I3,"      VALUE:",F12.6," Degree   * ANGLE   *"  )
004 FORMAT("* INTERNAL COORDINATE N.:",I3,"      VALUE:",F12.6," Degree   * TORSION *"  )
005 FORMAT("*************************************************************************"  )
  END SUBROUTINE DYNAMICS_ZMAT


  SUBROUTINE PRINTZMATDYN(  UPOT, K, STEP, DTORIG, VFOND, VEXC )
    USE DYNACART,    ONLY:    PUNCHXYZ
    USE CONVFACTORS, ONLY:    FACT => RADTODEGREE, ANGTOBOHR
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: UPOT, K, DTORIG, VFOND, VEXC
    INTEGER, INTENT(IN)          :: STEP
    INTEGER                      :: I

!   ---------------------------------------------------------------------
!   Initial title
!   ---------------------------------------------------------------------
    WRITE(6,100)('*',I=1,80)
    WRITE(6,102)
    WRITE(6,101)STEP,REAL(STEP)*DTORIG
    WRITE(6,102)
    WRITE(6,100)('*',I=1,80)
!   ---------------------------------------------------------------------
!   Coordinates and velocities
!   ---------------------------------------------------------------------
    WRITE(6,102)
    WRITE(6,103)
    WRITE(6,102)      
    WRITE(6,100)('*',I=1,80)
    WRITE(6,102)
    DO I=1,NUMAT
       WRITE(6,106)LAB(I),C(I)%X,C(I)%Y,C(I)%Z,V(I)%X,V(I)%Y,V(I)%Z
    END DO
!   ---------------------------------------------------------------------
!   Energies information
!   ---------------------------------------------------------------------
    WRITE(6,102)
    WRITE(6,100)('*',I=1,80)
    WRITE(6,105)UPOT,K
    WRITE(6,108)UPOT+K,VEXC-VFOND
    WRITE(6,102)
    WRITE(6,107)VFOND,VEXC
!   ---------------------------------------------------------------------
!   Gradients information
!   ---------------------------------------------------------------------
    WRITE(6,100)('*',I=1,80)
    WRITE(6,102)
    WRITE(6,111)
    WRITE(6,102)
    WRITE(6,100)('*',I=1,80)
    WRITE(6,110)
    WRITE(6,100)('*',I=1,80)
    DO I=1,NUMAT
       WRITE(6,104)LAB(I), GEO(1,I)/ANGTOBOHR, NA(I),      &
                           GEO(2,I)*FACT, NB(I),           &
                           GEO(3,I)*FACT, NC(I),           &
                           F(I)%X, F(I)%Y, F(I)%Z
    END DO
    WRITE(6,100)('*',I=1,80)
    CALL PUNCHXYZ(NUMAT,LAB,C)
    CALL FLUSH_FILES
    RETURN

100 FORMAT(80A)
101 FORMAT('*',' STEP NUMBER=',I6,38X,'TIME=',F12.3,' fs',' *')
102 FORMAT('*',78X,'*')
103 FORMAT('*',19X,'GEOMETRY',12X,17X,'VELOCITIES',12X,'*')

104 FORMAT('*',1X,A2,F9.4,I3,F9.4,I3,F9.4,I3,2X,3F12.6,1X,'*')
106 FORMAT('*',1X,A2,3F12.6,2X,3F12.6,1X,'*')
105 FORMAT('*',1X,'POTENTIAL ENERGY=',F12.6,' A.U.',  &
               9X,'KINETIC ENERGY=',F13.6,' A.U.',1X,'*')
107 FORMAT('*',' FUNDAMENTAL STATE ENERGY=',F12.6,6X,     &
               'EXCITED STATE ENERGY=',F12.6,1X,'*')
108 FORMAT('*',' TOTAL ENERGY=',F16.6,1X,'A.U.',9X        &
              ,'( E2 - E1 )=',F16.8,' A.U.',1X,'*')
110 FORMAT('*',11X,'  ZMAT INTERNAL GEOMETRY ',19X, &
                   'XYZ GRADIENTS',10X,'*')
111 FORMAT('*',31X ,'GRADIENTS SECTION',30X,'*')
    
  END SUBROUTINE PRINTZMATDYN


  SUBROUTINE GET_XYZ_GEO !( C,  CX,  GEO,  NUMAT)
    USE  READ_ZMAT,   ONLY:  ZMAT_TO_XYZ
    USE SYSTEM_UTIL,  ONLY:  CPSTOP
    USE READ_ZMAT,    ONLY:  XYZ_TO_INTERNAL
    IMPLICIT NONE
    INTEGER :: I, II, ITERATION
    INTEGER, INTENT(IN) :: NUM_VAR
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NUM_VAR) :: DELTAQ
    DOUBLE PRECISION, INTENT(IN), DIMENSION(NUM_VAR)    :: QINT
    DOUBLE PRECISION, INTENT(IN), DIMENSION(NUM_VAR,3*NUMAT) :: BMAT
    DOUBLE PRECISION, INTENT(IN), DIMENSION(3*NUMAT,3*NUMAT) :: MASS
    DOUBLE PRECISION, ALLOCATABLE :: WORK1(:,:), WORK2(:,:), CXITER(:)
    DOUBLE PRECISION, ALLOCATABLE :: WORKL(:), MATIT(:,:), QINTUP(:), QINTO(:)
    DOUBLE PRECISION  :: CONV
    INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER :: INFO, LWORK

    CALL  ZMAT_TO_XYZ(GEO, CX, NA, NB, NC, INT(NUMAT,KIND=4))
    ! Tranfer values in C
    DO I = 1, NUMAT
       II = (I-1)*3
       C(I)%X = CX(II+1)
       C(I)%Y = CX(II+2)
       C(I)%Z = CX(II+3)
    END DO

    RETURN
    
    ENTRY UPDATE_XYZ ( DELTAQ, NUM_VAR , MASS, BMAT, QINT)
    
    IF (.NOT.ALLOCATED(WORK1))    ALLOCATE(WORK1(3*NUMAT,NUM_VAR))
    IF (.NOT.ALLOCATED(WORK2))    ALLOCATE(WORK2(NUM_VAR,NUM_VAR))
    IF (.NOT.ALLOCATED(MATIT))    ALLOCATE(MATIT(3*NUMAT,NUM_VAR))
    IF (.NOT.ALLOCATED(CXITER))   ALLOCATE(CXITER(3*NUMAT))
    IF (.NOT.ALLOCATED(QINTUP))   ALLOCATE(QINTUP(NUM_VAR))
    WORK1  = 0.D0
    WORK2  = 0.D0
    MATIT  = 0.D0
    CXITER = 0.D0
    QINTUP = 0.D0
    WORK1  = TRANSPOSE(BMAT)
    WORK2  = MATMUL( BMAT, MATMUL(MASS, WORK1))
    ! Invert the matrix
    ! First LU Factorize ....
    IF (.NOT.ALLOCATED(IPIV)) ALLOCATE(IPIV(NUM_VAR))
    !
    CALL DGETRF(NUM_VAR, NUM_VAR, WORK2, NUM_VAR, IPIV, INFO)
    !
    IF (INFO.EQ.0) THEN
       WRITE(6,'(A)')"DGETRF:: FACTORIZATION SUCCESSFULLY!"
    ELSEIF (INFO.LT.0) THEN
       WRITE(6,'(A,I5,A)')"INFO=",INFO,": the i-th argument had an illegal value!"
       CALL CPSTOP("UPDATE_XYZ:: DGETRF")
    ELSE
       WRITE(6,'(A)')"INFO =",INFO,"U(i,i) is exactly zero. The factorization "// &
                 "has been completed, but the factor U is exactly "// &
                 "singular, and division by zero will occur if it is used"// &
                 "to solve a system of equations."
       CALL CPSTOP("UPDATE_XYZ:: DGETRF")
    END IF
    ! Now Invert...
    ! Query the optimal dimension
    ALLOCATE (WORKL(NUM_VAR))
    CALL DGETRI(NUM_VAR, WORK2, NUM_VAR, IPIV, WORKL, -1, INFO)
    LWORK = WORKL(1)
    DEALLOCATE (WORKL)
    ! Invert
    ALLOCATE (WORKL(LWORK))
    CALL DGETRI(NUM_VAR, WORK2, NUM_VAR, IPIV, WORKL, LWORK, INFO)
    DEALLOCATE (WORKL)
    DEALLOCATE (IPIV)
    IF (INFO.EQ.0) THEN
       WRITE(6,'(A)')"DGETRI:: INVERTION SUCCESSFULLY!"
    ELSEIF (INFO.LT.0) THEN
       WRITE(6,'(A,I5,A)')"INFO=",INFO,": the i-th argument had an illegal value!"
       CALL CPSTOP(" UPDATE_XYZ:: DGETRI")
    ELSE
       WRITE(6,'(A)')"INFO =",INFO,"U(i,i) is exactly zero;the matrix is"// &
                     "singular and its inverse could not be computed."
       CALL CPSTOP("UPDATE_XYZ:: DGETRI")
    END IF
    MATIT = MATMUL(MASS, MATMUL(WORK1,WORK2))
    ITERATION = 0
    CXITER = CX 
    CONV = SQRT(DOT_PRODUCT(MATMUL(MATIT, DELTAQ),MATMUL(MATIT, DELTAQ)))
    ! Iterate up to convergence...
    DO WHILE ((CONV.GT.1.D-7).AND.(ITERATION.LT.100))
       CXITER = CXITER - MATMUL(MATIT, DELTAQ)
       ITERATION = ITERATION + 1
       CALL XYZ_TO_INTERNAL(CXITER, QINTUP, INT(NUMAT), NUM_VAR, NA, NB, NC)
       DELTAQ = QINTUP - QINT
       CONV = SQRT(DOT_PRODUCT(MATMUL(MATIT, DELTAQ),MATMUL(MATIT, DELTAQ)))
    END DO
    CXITER = CXITER - MATMUL(MATIT, DELTAQ)
    WRITE(6,'(A,I6,A,ES14.6)')'  NUMBER OF ITERATIONS: ',ITERATION, &
                             '  CONVERGENCE VALUE: ', CONV
    IF (ITERATION.EQ.100) CALL CPSTOP(" BMATRIX ZMAT-> XYZ ITERATION DID NOT CONVERGE!")

    CX = CXITER
    ! Tranfer values in C
    DO I = 1, NUMAT
       II = (I-1)*3
       C(I)%X = CX(II+1)
       C(I)%Y = CX(II+2)
       C(I)%Z = CX(II+3)
    END DO

    DEALLOCATE(QINTUP)
    DEALLOCATE(CXITER)
    DEALLOCATE(MATIT)
    DEALLOCATE(WORK2)
    DEALLOCATE(WORK1)
    RETURN
  END SUBROUTINE GET_XYZ_GEO


  SUBROUTINE SINGLE_ZMAT_DYN_STEP(DT, VEXC, VFOND, UPOT, STEP,             &
                                  QINTN, QINTO, QINTV, MASS, NUM_VAR)
    USE VERLET_ZMAT_INT,  ONLY:   MOVE_ZMAT_A,    &
                                  MOVE_ZMAT_B
    USE FORCE,            ONLY:   FORCE_CART
    USE BUILD_BMAT,       ONLY:   COMPUTE_MATRIX, &
                                  BMAT,           &
                                  CMAT
    IMPLICIT NONE
    ! Arguments
    DOUBLE PRECISION, INTENT(IN)    :: DT
    DOUBLE PRECISION, INTENT(INOUT) :: VEXC, VFOND, UPOT
    INTEGER,          INTENT(IN)    :: STEP, NUM_VAR
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NUM_VAR) :: QINTN, QINTO, QINTV
    DOUBLE PRECISION, INTENT(IN), DIMENSION(3*NUMAT,3*NUMAT) :: MASS
    ! Local Variables
    DOUBLE PRECISION, ALLOCATABLE   :: DELTAQ(:)

    CALL COMPUTE_MATRIX(CX, GEO, NA, NB, NC)

    IF (.NOT.ALLOCATED(DELTAQ)) ALLOCATE(DELTAQ(NUM_VAR))

    QINTO = QINTN
    CALL MOVE_ZMAT_A  (  DT,  MASS, QINTN, QINTV, BMAT, CMAT, NUM_VAR) ! UPDATE QINT

    DELTAQ = QINTN - QINTO
    CALL UPDATE_XYZ ( DELTAQ, NUM_VAR , MASS, BMAT, QINTN) 
    !
    CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP ) 
    !
    CALL MOVE_ZMAT_B  (  DT,  MASS, QINTN, QINTV, BMAT, CMAT, NUM_VAR) ! UPDATE VX
    

    RETURN
  END SUBROUTINE SINGLE_ZMAT_DYN_STEP

END MODULE DYNA_ZMAT

