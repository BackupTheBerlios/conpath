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
    ! 
    !
    CALL FLUSH_FILES
    !
    ! Compute Forces...
    !
    ! Sono arrivato qui...
    ! Bisogna:
    ! 1) Calcolare le forze in coordinate cartesiane
    ! 2) Proiettarle nelle coordinate interne
    ! 3) Annullare i gradienti delle coord interne che vuoi congelare
    ! 4) ritrasformale nelle coordinate cartesiane
    ! 5) Chiamare in questa routine il driver per Lbfgs
    ! 6) I punti 1)->4) devono essere fatti da una routine interna a questo modulo...
    !
    CALL GET_XYZ_GEO
    CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP )
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
  END SUBROUTINE INTERNAL_OPTIMIZE_PROC
END MODULE OPTIMIZATION_PROCEDURES
