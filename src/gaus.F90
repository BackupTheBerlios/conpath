!-----------------------------------------------------------------------
MODULE CAMPANE
      IMPLICIT NONE
      DOUBLE PRECISION :: VAL1, VAL2
!     ------------------------------------------------------------------
      CONTAINS
!     ==================================================================
      DOUBLE PRECISION FUNCTION GAUSSIAN ( C, C0, G, H, ALPHA,   &
                                SIGMA, NATOM )
!     ==================================================================
      USE XYZ, ONLY: POINT
      IMPLICIT NONE
      TYPE(POINT), INTENT(IN) :: G(*), H(*), C(*), C0(*)
      INTEGER (KIND=8), INTENT(IN) :: NATOM
      INTEGER :: I
      DOUBLE PRECISION :: ALPHA,SIGMA
!     ------------------------------------------------------------------
      VAL1 = 0.D0
      VAL2 = 0.D0
!     ------------------------------------------------------------------
      DO I = 1 , NATOM
         VAL1 = VAL1 + G(I)%X * ( C(I)%X - C0(I)%X ) +    &
                       G(I)%Y * ( C(I)%Y - C0(I)%Y ) +    &
                       G(I)%Z * ( C(I)%Z - C0(I)%Z ) 
         VAL2 = VAL2 + H(I)%X * ( C(I)%X - C0(I)%X ) +    &
                       H(I)%Y * ( C(I)%Y - C0(I)%Y ) +    &
                       H(I)%Z * ( C(I)%Z - C0(I)%Z ) 
      END DO  
      GAUSSIAN=ALPHA*DEXP(-(VAL1*VAL1+VAL2*VAL2)/(2.D0*SIGMA*SIGMA))
      RETURN
!     ------------------------------------------------------------------
      END FUNCTION GAUSSIAN
!     ==================================================================
      SUBROUTINE GAUS_FORCES(C, C0,  G, H, ALPHA, SIGMA, NATOM,    &
                             GAUSSIAN, F)
!     ==================================================================
      USE XYZ, ONLY: POINT
      IMPLICIT NONE
      TYPE(POINT), INTENT(IN) :: G(*), H(*), C(*), C0(*)
      TYPE(POINT), INTENT(INOUT) :: F(*)
      INTEGER (KIND=8), INTENT(IN) :: NATOM
      DOUBLE PRECISION :: ALPHA, SIGMA, GAUSSIAN
      INTEGER :: I
  
      DO I=1,NATOM
         write(6,'(/3f15.9,I5)')F(I)%X,F(I)%Y,F(I)%Z,I
         F(I)%X=F(I)%X-GAUSSIAN*(VAL1+VAL2)/(SIGMA*SIGMA)*(G(I)%X+H(I)%X)  
         F(I)%Y=F(I)%Y-GAUSSIAN*(VAL1+VAL2)/(SIGMA*SIGMA)*(G(I)%Y+H(I)%Y)
         F(I)%Z=F(I)%Z-GAUSSIAN*(VAL1+VAL2)/(SIGMA*SIGMA)*(G(I)%Z+H(I)%Z)
         write(6,'(3f15.9,I5/)')    &   
              GAUSSIAN*(VAL1+VAL2)/(SIGMA*SIGMA)*(G(I)%X+H(I)%X), &
              GAUSSIAN*(VAL1+VAL2)/(SIGMA*SIGMA)*(G(I)%Y+H(I)%Y), &
              GAUSSIAN*(VAL1+VAL2)/(SIGMA*SIGMA)*(G(I)%Z+H(I)%Z), I
      ENDDO
      RETURN
!     ------------------------------------------------------------------
      END SUBROUTINE GAUS_FORCES
!     ==================================================================
END MODULE CAMPANE
