! Authors: Teodoro Laino && Daniele Passerone      SNS - INFM - ETHZ
!
! Zurich (ETHZ) - Date 07.07.03
!
MODULE BUILD_BMAT
  USE READ_ZMAT, ONLY:    NUMAT,           &
                          NUM_BONDS,       &
                          NUM_ANGLES,      &
                          NUM_TORSIONS
  
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE :: BMAT(:,:), CMAT(:,:,:)
CONTAINS
!    --------------------------------------------------------------------
!    THE B AND C MATRIX ( WILSON - ELYASHEVICH ) ARE DEFINED IN TERMS
!    OF THE EXPANSION OF INTERNAL COORDINATES q_i IN POWERS OF THE 
!    CHANGES OF THE CARTESIAN x_k AT TIME = 0
!
!
!    q_i - q_i^0 = Sum_k ( B_i^k * x_k ) + 
!                  + 1/2 Sum_k Sum_l ( C_i^{kl} x_k x_l ) + ...  
!
!    B CONTAINS THE FIRST DERIVATIVES OF THE INTERNAL COORDINATES WITH 
!      RESPECT TO THE CARTESIAN
!    C IS THE SECOND-ORDER ANALOG, THAT IS IT CONTAINS THE SECOND 
!      DERIVATIVES: C_i^{kl} = d^2 q_i / dx_k dx_l
!
!    --------------------------------------------------------------------
  SUBROUTINE COMPUTE_MATRIX(COORD, GEO, NA, NB, NC)
    USE PRINT_MATRIX, ONLY:  PRINT_MAT2, &
                             PRINT_MAT3
    IMPLICIT NONE
    ! Arguments
    DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: COORD, GEO
    INTEGER, DIMENSION(*), INTENT(IN) :: NA, NB, NC
    ! Local Variables
    INTEGER :: NUM_VAR
    LOGICAL :: DEBUG
    CHARACTER (LEN=80) :: TITLE
    CHARACTER (LEN=3),  ALLOCATABLE :: LABELX(:)
    CHARACTER (LEN=10), ALLOCATABLE :: LABELY(:)
    CHARACTER (LEN=10), ALLOCATABLE :: LABELG(:)

    DEBUG = .FALSE.
    NUM_VAR = NUM_BONDS + NUM_ANGLES + NUM_TORSIONS

    CALL CREATE_BMAT(COORD, GEO, NA, NB, NC)
    CALL CREATE_CMAT(COORD, GEO, NA, NB, NC)

    !
    ! Debug Section
    !
    IF (DEBUG) THEN
       ! DEBUG FIRST DERIVATIVES...
       ALLOCATE ( LABELY(NUM_VAR) )
       ! Debug Headers...
       LABELY( 1) = 'BOND     1'
       LABELY( 2) = 'BOND     2'
       LABELY( 3) = 'ANGLE    2'
       LABELY( 4) = 'BOND     3'
       LABELY( 5) = 'ANGLE    3'
       LABELY( 6) = 'TORSION  3'
       LABELY( 7) = 'BOND     4'
       LABELY( 8) = 'ANGLE    4'
       LABELY( 9) = 'TORSION  4'
       LABELY(10) = 'BOND     5'
       LABELY(11) = 'ANGLE    5'
       LABELY(12) = 'TORSION  5'
       TITLE = "BMATRIX: INTERNAL COORDINATES FIRST DERIVATIVES"

       CALL PRINT_MAT2 ( TITLE=TITLE, LABELY=LABELY, MAT=BMAT,  &
                         IR=NUM_VAR,  IC=3*NUMAT )
       DEALLOCATE ( LABELY )
       ! DEBUG SECOND DERIVATIVES...
        ALLOCATE ( LABELY(3*NUMAT), LABELX(3*NUMAT), LABELG(NUM_VAR) )
        ! Debug Headers
        LABELX( 1) = "DX1" ; LABELY( 1) = "DX1" ;  LABELG( 1) = 'BOND     1'
        LABELX( 2) = "DY1" ; LABELY( 2) = "DY1" ;  LABELG( 2) = 'BOND     2'
        LABELX( 3) = "DZ1" ; LABELY( 3) = "DZ1" ;  LABELG( 3) = 'ANGLE    2'
        LABELX( 4) = "DX2" ; LABELY( 4) = "DX2" ;  LABELG( 4) = 'BOND     3'
        LABELX( 5) = "DY2" ; LABELY( 5) = "DY2" ;  LABELG( 5) = 'ANGLE    3'
        LABELX( 6) = "DZ2" ; LABELY( 6) = "DZ2" ;  LABELG( 6) = 'TORSION  3'
        LABELX( 7) = "DX3" ; LABELY( 7) = "DX3" ;  LABELG( 7) = 'BOND     4'
        LABELX( 8) = "DY3" ; LABELY( 8) = "DY3" ;  LABELG( 8) = 'ANGLE    4'
        LABELX( 9) = "DZ3" ; LABELY( 9) = "DZ3" ;  LABELG( 9) = 'TORSION  4'
        LABELX(10) = "DX4" ; LABELY(10) = "DX4" ;  LABELG(10) = 'BOND     5'
        LABELX(11) = "DY4" ; LABELY(11) = "DY4" ;  LABELG(11) = 'ANGLE    5'
        LABELX(12) = "DZ4" ; LABELY(12) = "DZ4" ;  LABELG(12) = 'TORSION  5'
        LABELX(13) = "DX5" ; LABELY(13) = "DX5" ; 
        LABELX(14) = "DY5" ; LABELY(14) = "DY5" ; 
        LABELX(15) = "DZ5" ; LABELY(15) = "DZ5" ; 
        LABELX(16) = "DX6" ; LABELY(16) = "DX6" ; 
        LABELX(17) = "DY6" ; LABELY(17) = "DY6" ; 
        LABELX(18) = "DZ6" ; LABELY(18) = "DZ6" ; 
        TITLE = "CMATRIX: INTERNAL COORDINATES SECOND DERIVATIVES"

       CALL PRINT_MAT3 ( TITLE=TITLE, LABELX=LABELX, LABELY=LABELY,  &
                         MAT=CMAT, NCOL=18, IR=NUM_VAR,  IC=3*NUMAT, &
                         IC2=3*NUMAT, LABELG=LABELG  )
    ENDIF


    RETURN
  END SUBROUTINE COMPUTE_MATRIX


  SUBROUTINE CREATE_BMAT(COORD,GEO,NA,NB,NC)
    USE DERIVATE_GEO_OBJECTS, ONLY:   DERI_BOND,        &
                                      DERI_ANGLE,       &
                                      DERI_TORSION
    IMPLICIT NONE
    ! Arguments
    DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: COORD, GEO
    INTEGER, DIMENSION(*), INTENT(IN) :: NA, NB, NC
    
    ! Local variables
    INTEGER :: NUM_VAR, INDB, IAT1, IAT2, IAT3, IAT, J, K
    DOUBLE PRECISION, ALLOCATABLE :: F(:)
    DOUBLE PRECISION :: FACT

    ! Number of variables in Bmat Matrix..
    NUM_VAR = NUM_BONDS + NUM_ANGLES + NUM_TORSIONS
    INDB = 0
    IAT  = 0
    IAT1 = 0
    IAT2 = 0
    IAT3 = 0
    FACT = 1.D0

    ! Allocate Bmat
    IF (.NOT.ALLOCATED(F))    ALLOCATE(F(3*NUMAT))
    IF (.NOT.ALLOCATED(BMAT)) ALLOCATE(BMAT(NUM_VAR,NUMAT*3))
    BMAT = 0.D0    
    F    = 0.D0
!------------------------------------------------------------------------
!
!   B_i^k = d q_i / d x_k
!
!   IN MODULE DERIVATE_GEO_OBJECTS :
!     1)   DERI_BOND    ( I, J, C, F, FACT )            
!                                                            ^
!     2)   DERI_ANGLE   ( J, K, I, C, F, FACT )    : Angle <JKI>
!
!     3)   DERI_TORSION ( K, J, M, L, C, F, FACT ) : Dihedral  <KJML>
!
!   B MATRIX DIVIDED IN BLOCKS : 
!     ( BOND1, BOND2, ANGLE2, BOND3, ANGLE3, TORSION3, .., .., TORSIONn)
!     According the order of the GEO  Array (Internal Coordinates Array) 
!------------------------------------------------------------------------
    
    DO J= 2, NUMAT
       ! Bond section
       INDB = INDB +1       
       IAT = NA(J)
       F = 0.D0
       CALL DERI_BOND ( J, IAT, COORD, F, FACT )
       DO K=1,3
          BMAT(INDB, (J-1)*3+K)   = F((J-1)*3+K)
          BMAT(INDB, (IAT-1)*3+K) = F((IAT-1)*3+K) 
       END DO
       IF (J.EQ.2) CYCLE
       ! Angle section
       INDB = INDB +1
       IAT1 = NA(J)
       IAT2 = NB(J)
       F=0.D0
       CALL DERI_ANGLE ( J, IAT1, IAT2, COORD, F, FACT )
       DO K=1,3
          BMAT(INDB,(J-1)*3+K)    = F((J-1)*3+K)
          BMAT(INDB,(IAT1-1)*3+K) = F((IAT1-1)*3+K) 
          BMAT(INDB,(IAT2-1)*3+K) = F((IAT2-1)*3+K) 
       END DO
       If (J.EQ.3) CYCLE
       ! Torsion section
       INDB = INDB +1
       IAT1 = NA(J)
       IAT2 = NB(J)
       IAT3 = NC(J)
       F=0.D0
       CALL DERI_TORSION( J, IAT1, IAT2, IAT3, COORD, F, FACT )
       DO K=1,3
          BMAT(INDB,(J-1)*3+K)    = F((J-1)*3+K)
          BMAT(INDB,(IAT1-1)*3+K) = F((IAT1-1)*3+K) 
          BMAT(INDB,(IAT2-1)*3+K) = F((IAT2-1)*3+K) 
          BMAT(INDB,(IAT3-1)*3+K) = F((IAT3-1)*3+K) 
       END DO
    END DO

    RETURN
  END SUBROUTINE CREATE_BMAT
  
  
  SUBROUTINE CREATE_CMAT(COORD,GEO,NA,NB,NC)
    USE DERIVATE_2_GEO_OBJECTS, ONLY:   DERI_2_BOND,        &
                                        DERI_2_ANGLE,       &
                                        DERI_2_TORSION
    IMPLICIT NONE
    ! Arguments
    DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: COORD, GEO
    INTEGER,          DIMENSION(*), INTENT(IN) :: NA, NB, NC
    ! Local Variables
    DOUBLE PRECISION, ALLOCATABLE :: F(:)
    DOUBLE PRECISION :: FACT
    INTEGER :: NUM_VAR, INDB, INDF
    INTEGER :: IAT, IAT1, IAT2, IAT3
    INTEGER :: I, J, K, II, KK, ATK, ATI

    ! Number of variables in Bmat Matrix..
    NUM_VAR = NUM_BONDS + NUM_ANGLES + NUM_TORSIONS
    INDB = 0
    IAT  = 0
    IAT1 = 0
    IAT2 = 0
    IAT3 = 0
    FACT = 1.D0

    ! Allocate Bmat
    IF (.NOT.ALLOCATED(CMAT)) ALLOCATE(CMAT( NUM_VAR, NUMAT*3, NUMAT*3 ))
    CMAT = 0.D0  

    ! Main Loop for second derivatives computation  
    DO J= 2, NUMAT
       ! Bond section
       IF (.NOT.ALLOCATED(F))    ALLOCATE( F(21) )
       F    = 0.D0
       INDF = 0
       INDB = INDB + 1  
       IAT  = NA(J)
       CALL DERI_2_BOND ( J, IAT, COORD, F, FACT )       
       DO K= 1, 6
          DO I = 1, K
             INDF = INDF + 1
             KK   = MOD (K, 3)
             II   = MOD (I, 3)
             IF (KK.EQ.0) KK = 3
             IF (II.EQ.0) II = 3
             ATK = J
             ATI = J
             IF ( K .GT. 3 ) ATK = IAT
             IF ( I .GT. 3 ) ATI = IAT
             CMAT(INDB, (ATK - 1)*3+KK, (ATI - 1)*3+II )   = F( INDF )
             CMAT(INDB, (ATI - 1)*3+II, (ATK - 1)*3+KK )   = F( INDF ) 
          END DO
       END DO
       DEALLOCATE ( F )

       IF (J.EQ.2) CYCLE
       ! Angle section
       IF (.NOT.ALLOCATED(F))    ALLOCATE( F(45) )
       F    = 0.D0
       INDF = 0
       INDB = INDB + 1
       IAT1 = NA(J)
       IAT2 = NB(J)
       CALL DERI_2_ANGLE ( J, IAT1, IAT2, COORD, F, FACT )
       DO K= 1, 9
          DO I = 1, K
             INDF = INDF + 1
             KK   = MOD (K, 3)
             II   = MOD (I, 3)
             IF (KK.EQ.0) KK = 3
             IF (II.EQ.0) II = 3
             ATK = J
             ATI = J
             IF ( K .GT. 3 ) ATK = IAT1
             IF ( I .GT. 3 ) ATI = IAT1
             IF ( K .GT. 6 ) ATK = IAT2
             IF ( I .GT. 6 ) ATI = IAT2             
             CMAT(INDB, (ATK - 1)*3+KK, (ATI - 1)*3+II )   = F( INDF )
             CMAT(INDB, (ATI - 1)*3+II, (ATK - 1)*3+KK )   = F( INDF )
          END DO
       END DO
       DEALLOCATE ( F )

       IF (J.EQ.3) CYCLE
       ! Torsion section
       IF (.NOT.ALLOCATED(F))    ALLOCATE( F(78) )
       F    = 0.D0
       INDF = 0
       INDB = INDB + 1
       IAT1 = NA(J)
       IAT2 = NB(J)
       IAT3 = NC(J)
       CALL DERI_2_TORSION( J, IAT1, IAT2, IAT3, COORD, F, FACT )
       DO K= 1, 12
          DO I = 1, K
             INDF = INDF + 1
             KK   = MOD (K, 3)
             II   = MOD (I, 3)
             IF (KK.EQ.0) KK = 3
             IF (II.EQ.0) II = 3
             ATK = J
             ATI = J
             IF ( K .GT. 3 ) ATK = IAT1
             IF ( I .GT. 3 ) ATI = IAT1
             IF ( K .GT. 6 ) ATK = IAT2
             IF ( I .GT. 6 ) ATI = IAT2
             IF ( K .GT. 9 ) ATK = IAT3
             IF ( I .GT. 9 ) ATI = IAT3
             CMAT(INDB, (ATK - 1)*3+KK, (ATI - 1)*3+II )   = F( INDF )
             CMAT(INDB, (ATI - 1)*3+II, (ATK - 1)*3+KK )   = F( INDF )
          END DO
       END DO
       DEALLOCATE ( F )

    END DO   
    
    RETURN
  END SUBROUTINE CREATE_CMAT
END MODULE BUILD_BMAT

  
