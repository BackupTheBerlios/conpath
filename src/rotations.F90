MODULE ROTATIONS
  USE XYZ, ONLY:     POINT
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE :: ROT(:,:)

CONTAINS
  SUBROUTINE MATRIX_ROTATION ( ALPHA, BETA, GAMMA)
!   -------------------------------------------------------------------
!   THIS SUBROUTINE BUILD A GENERAL ROTATION MATRIX
!   -------------------------------------------------------------------
!   THE ROTATION GIVEN BY THE EULER ANGLES ALPHA, BETA, AND GAMMA CAN 
!   BE DECOMPOSED INTO A SEQUENCE OF THREE SUCCESSIVE ROTATIONS. 
!   THE FIRST BY ANGLE ALPHA ABOUT THE Z AXIS, THE SECOND BY ANGLE BETA 
!   ABOUT THE X AXIS, AND THE THIRD ABOUT THE Z AXIS (AGAIN) BY ANGLE 
!   GAMMA.  THE ANGLE BETA IS RESTRICTED TO THE RANGE 0 TO PI.   
!   DEFINITION OF EULER ANGLES FOR MAIN ROTATIONS
!
!                                      ALPHA         BETA         GAMMA   
!   ROTATION OF A deg. AROUND X-AXIS :   0            A             0
!   ROTATION OF A deg. AROUND Y-AXIS : -Pi/2          A            Pi/2
!   ROTATION OF A deg. AROUND X-AXIS :   A            0             0
!   -------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: ALPHA, BETA, GAMMA
    
    IF ( ALLOCATED( ROT ) ) DEALLOCATE( ROT )
    ALLOCATE( ROT(3,3) )
    
    ROT=0.D0

    ROT(1,1) =   COS(ALPHA)*COS(GAMMA) - COS(BETA)*SIN(ALPHA)*SIN(GAMMA)
    ROT(1,2) = -(COS(BETA)*COS(GAMMA)*SIN(ALPHA)) - COS(ALPHA)*SIN(GAMMA)  
    ROT(1,3) =   SIN(ALPHA)*SIN(BETA)
 
    ROT(2,1) =   COS(GAMMA)*SIN(ALPHA) + COS(ALPHA)*COS(BETA)*SIN(GAMMA)
    ROT(2,2) =   COS(ALPHA)*COS(BETA)*COS(GAMMA) - SIN(ALPHA)*SIN(GAMMA)
    ROT(2,3) =   -(COS(ALPHA)*SIN(BETA)) 

    ROT(3,1) =   SIN(BETA)*SIN(GAMMA) 
    ROT(3,2) =   COS(GAMMA)*SIN(BETA)
    ROT(3,3) =   COS(BETA)

    RETURN
  END SUBROUTINE MATRIX_ROTATION


  FUNCTION ROTATE3D(V, ALPHA, BETA, GAMMA, VZERO)
!   --------------------------------------------------------------------
!   THIS FUNCTION COMPUTES THE ROTATION VECTOR GIVEN THE THREE EULER
!   ANGLES AROUND POINT V0. 
!   FOR ALPHA, BET AND GAMMA DEFINITION SEE THE PREVIOUS SUBROUTINE
!   --------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(POINT)  :: ROTATE3D
    TYPE(POINT), INTENT(IN)  :: V
    TYPE(POINT), INTENT(IN), OPTIONAL :: VZERO
    TYPE(POINT)  :: V0
    DOUBLE PRECISION, INTENT(IN) :: ALPHA, BETA, GAMMA

    IF (PRESENT(VZERO)) THEN
       V0%X=VZERO%X
       V0%Y=VZERO%Y
       V0%Z=VZERO%Z      
    ELSE
       V0%X=0.D0
       V0%Y=0.D0
       V0%Z=0.D0       
    END IF

    ROTATE3D%X = ( V0%X - V%X)+ (-V0%Z + V%Z)*SIN(BETA)*SIN(GAMMA) +                          & 
                 (-V0%Y + V%Y)*(COS(GAMMA)*SIN(ALPHA) + COS(ALPHA)*COS(BETA)*SIN(GAMMA)) +    &
                 (-V0%X + V%X)*(COS(ALPHA)*COS(GAMMA) - COS(BETA)*SIN(ALPHA)*SIN(GAMMA))
                 
    ROTATE3D%Y = ( V0%Y - V%Y)+ (-V0%Z + V%Z)*COS(GAMMA)*SIN(BETA) +                          &
                 (-V0%X + V%X)*(-(COS(BETA)*COS(GAMMA)*SIN(ALPHA)) - COS(ALPHA)*SIN(GAMMA)) + &
                 (-V0%Y + V%Y)*(COS(ALPHA)*COS(BETA)*COS(GAMMA) - SIN(ALPHA)*SIN(GAMMA))
                 
    ROTATE3D%Z = ( V0%Z - V%Z)+ (-V0%Z + V%Z)*COS(BETA) -   & 
                 (-V0%Y + V%Y)*COS(ALPHA)*SIN(BETA) +       &
                 (-V0%X + V%X)*SIN(ALPHA)*SIN(BETA)

  END FUNCTION ROTATE3D

END MODULE ROTATIONS
