MODULE OPTIMIZATION_PROCEDURES
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
  USE DYNA_ZMAT,    ONLY: GET_XYZ_GEO
  USE SYSTEM_UTIL,  ONLY: CPSTOP
  IMPLICIT NONE
CONTAINS
  SUBROUTINE INTERNAL_OPTIMIZE_PROC
    USE FORCE,       ONLY: FORCE_CART
    USE CONVFACTORS, ONLY: ANGTOBOHR
    IMPLICIT NONE
    !
    ! THIS PROGRAM PERFORMS A CONICAL INTERSECTION OPTIMIZATION IN INTERNAL COORDINATES
    ! WITH THE POSSIBILITY TO KEEP FIXED SOME INTERNAL COORDINATES..
    ! 
    ! Local Variables
    DOUBLE PRECISION, ALLOCATABLE  :: QINTN(:)
    INTEGER                        :: NUM_VAR
    INTEGER                        :: I, J, L
    INTEGER                        :: RSTEP, STEP
    DOUBLE PRECISION               :: DT, PI
    DOUBLE PRECISION               :: K, E, UPOT, VEXC, VFOND
    DOUBLE PRECISION, ALLOCATABLE  :: LOW(:), UPP(:), NBD(:)
    !
    IF ( IPRINT .LE. 0 ) IPRINT = NSTEP + 1
    NUM_VAR = NUM_BONDS + NUM_ANGLES + NUM_TORSIONS
    ALLOCATE ( QINTN(NUM_VAR))
    PI = 4.D0 * ATAN(1.D0)
    QINTN  = 0.D0
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
    ! Write Header of Optimization in Internal Coordinates
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
    !
    CALL FLUSH_FILES
    CALL driver_Lbfgs ( NUM_VAR, M, QINTN, LOW, UPP, NBD)
    !
    !
    RETURN
001 FORMAT("*************************************************************************", &
         /,"*                                                                       *", &
         /,"* OPTIMIZATION MODULE: PERFORMING OPTIMIZATION IN INTERNAL COORDINATES  *", &
         /,"*                                                                       *", &
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
      ! 
  END SUBROUTINE INTERNAL_OPTIMIZE_PROC

  SUBROUTINE COMP_ENE_FORCES_OPTIMIZE ( XVAL, FUNC, F_INT_OUT, NLENGTH, STEP)
    USE  BUILD_BMAT,   ONLY:  CREATE_BMAT,   &
                              BMAT
    USE FORCE,       ONLY: FORCE_CART
    USE CONVFACTORS, ONLY: ANGTOBOHR
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN)   :: NLENGTH, STEP
    DOUBLE PRECISION, INTENT(INOUT)                     :: FUNC
    DOUBLE PRECISION, DIMENSION(NLENGTH), INTENT(INOUT) :: F_INT_OUT, XVAL
    ! Local Variables
    DOUBLE PRECISION, ALLOCATABLE  :: FINT(:)
    DOUBLE PRECISION :: VFOND, UPOT, PI
    INTEGER :: I, L

    PI = 4.D0 * ATAN(1.D0)
    ALLOCATE ( FINT(NLENGTH) )
    ! Update GEO Vector
    I = 0
    DO L = 2, NUMAT
       I = I + 1
       GEO(1,L) = XVAL(I)
       IF (L.EQ.2) CYCLE
       I = I + 1
       GEO(2,L) = XVAL(I)
       IF (L.EQ.3) CYCLE
       I = I + 1
       GEO(3,L) = XVAL(I)
    END DO
    ! Compute Forces...
    ! Get geometry into XYZ format
    ! Compute forces in Cartesian Coordinates
    CALL GET_XYZ_GEO
    CALL FORCE_CART ( FUNC , VFOND,  UPOT, STEP )
    ! Project out forces into internal coordinates
    CALL CREATE_BMAT(CX, GEO, NA, NB, NC)
    CALL PROJECT_FORCES( BMAT, FX, FINT, NLENGTH )
    ! Nullify gradients according constraints
    CALL IMPOSE_CONSTRAINT( FINT, NLENGTH )
    F_INT_OUT = FINT
    !        
    RETURN
  END SUBROUTINE COMP_ENE_FORCES_OPTIMIZE



  SUBROUTINE PROJECT_FORCES( BMAT, FX, FINT, NUM_VAR)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: NUM_VAR
    DOUBLE PRECISION, DIMENSION(NUM_VAR,3*NUMAT) :: BMAT
    DOUBLE PRECISION, DIMENSION(NUM_VAR) :: FINT
    DOUBLE PRECISION, DIMENSION(3*NUMAT) :: FX
    ! Local Variables
    
    FINT = MATMUL( BMAT, FX )

    RETURN
  END SUBROUTINE PROJECT_FORCES


  SUBROUTINE IMPOSE_CONSTRAINT( FINT, NUM_VAR )
    USE START_JOB, ONLY:   LOPT,         &
                           NINTCONSTR    
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: NUM_VAR
    DOUBLE PRECISION, DIMENSION(NUM_VAR) :: FINT
    ! Local Variables
    INTEGER :: I

    !
    DO I=1, NINTCONSTR
       FINT(LOPT(I)) = 0.D0
    END DO

    RETURN
  END SUBROUTINE IMPOSE_CONSTRAINT

END MODULE OPTIMIZATION_PROCEDURES
