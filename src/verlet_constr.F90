!------------------------------------------------------------------------
MODULE VERLET_CONSTR
      IMPLICIT NONE
!     -------------------------------------------------------------------
      CONTAINS
!     ===================================================================
      SUBROUTINE CONSTRA ( DT, M, V1, V2, C0 )
!     ===================================================================
      USE START_JOB,   ONLY: N => NUMAT
      USE DYNPREPARE,  ONLY: C,               &                        
                             V
      USE XYZ,         ONLY: POINT
!     -------------------------------------------------------------------
!     IT CORRECTS THE POSITIONS AND VELOCITIES WITH THE CONSTRAINT FORCES
!     FIRST PART
!     -------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)                :: DT
      DOUBLE PRECISION                            :: DT2, DTSQ2
      DOUBLE PRECISION, INTENT(IN), DIMENSION(N)  :: M(*)
      DOUBLE PRECISION                            :: A, D, E, F, B
      DOUBLE PRECISION                            :: DEN, LAM1, LAM2
      TYPE(POINT), INTENT(IN)                     :: V1(*), V2(*), C0(*)
      INTEGER                                     :: I
!     -------------------------------------------------------------------
      DT2 = DT / 2.0
      DTSQ2 = DT2 * DT

!
!   IT COMPUTES THE QUANTITIES OF INTEREST
!
      A = 0.d0
      B = 0.d0
      D = 0.d0
      E = 0.d0
      F = 0.d0

      DO I = 1, N
         A = A + (C(I)%X-C0(I)%X)*V1(I)%X
         A = A + (C(I)%Y-C0(I)%Y)*V1(I)%Y
         A = A + (C(I)%Z-C0(I)%Z)*V1(I)%Z
         B = B + (C(I)%X-C0(I)%X)*V2(I)%X
         B = B + (C(I)%Y-C0(I)%Y)*V2(I)%Y
         B = B + (C(I)%Z-C0(I)%Z)*V2(I)%Z
         D = D - V1(I)%X*V2(I)%X/ M(I)
         D = D - V1(I)%Y*V2(I)%Y/ M(I)
         D = D - V1(I)%Z*V2(I)%Z/ M(I)
         E = E - V1(I)%X*V1(I)%X/ M(I)
         E = E - V1(I)%Y*V1(I)%Y/ M(I)
         E = E - V1(I)%Z*V1(I)%Z/ M(I)
         F = F - V2(I)%X*V2(I)%X/ M(I)
         F = F - V2(I)%Y*V2(I)%Y/ M(I)
         F = F - V2(I)%Z*V2(I)%Z/ M(I)
      ENDDO

       
      DEN = DTSQ2 * (D * D - F * E)
      LAM1 = (F * A - D * B) / DEN
      LAM2 = (E * B - D * A) / DEN

! 
!  CORRECTION OF THE POSITIONS AND VELOCITIES
! 


      DO I = 1 , N
         C(I)%X =  C(I)%X - DTSQ2 * (V1(I)%X * LAM1 + V2(I)%X * LAM2) / M(I)  
         C(I)%Y =  C(I)%Y - DTSQ2 * (V1(I)%Y * LAM1 + V2(I)%Y * LAM2) / M(I)  
         C(I)%Z =  C(I)%Z - DTSQ2 * (V1(I)%Z * LAM1 + V2(I)%Z * LAM2) / M(I)  
         V(I)%X =  V(I)%X - DTSQ2 * (V1(I)%X * LAM1 + V2(I)%X * LAM2) / M(I)  
         V(I)%Y =  V(I)%Y - DTSQ2 * (V1(I)%Y * LAM1 + V2(I)%Y * LAM2) / M(I)  
         V(I)%Z =  V(I)%Z - DTSQ2 * (V1(I)%Z * LAM1 + V2(I)%Z * LAM2) / M(I)  
      END DO
 
!
!  CHECK 
!
      A = 0.d0
      B = 0.D0
      DO I = 1, N
         A = A + (C(I)%X-C0(I)%X)*V1(I)%X
         A = A + (C(I)%Y-C0(I)%Y)*V1(I)%Y
         A = A + (C(I)%Z-C0(I)%Z)*V1(I)%Z
         B = B + (C(I)%X-C0(I)%X)*V2(I)%X
         B = B + (C(I)%Y-C0(I)%Y)*V2(I)%Y
         B = B + (C(I)%Z-C0(I)%Z)*V2(I)%Z
      END DO

      WRITE (6,*)'Constraints: ',A,B
      RETURN
!     -------------------------------------------------------------------
      END SUBROUTINE CONSTRA
!     ===================================================================
      SUBROUTINE CONSTRB ( DT, M, V1, V2, C0 )
!     ===================================================================
      USE START_JOB,   ONLY: N => NUMAT
      USE DYNPREPARE,  ONLY: V
      USE XYZ,         ONLY: POINT
!     -------------------------------------------------------------------
!     IT CORRECTS THE VELOCITIES WITH THE CONSTRAINT FORCES
!     SECOND PART
!     -------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)                :: DT
      DOUBLE PRECISION                            :: DT2, DTSQ2
      DOUBLE PRECISION, INTENT(IN), DIMENSION(N)  :: M(*)
      DOUBLE PRECISION                            :: A, D, E, F, B
      DOUBLE PRECISION                            :: DEN, LAM1, LAM2
      TYPE(POINT), INTENT(IN)                     :: V1(*), V2(*), C0(*)
      INTEGER                                     :: I
!     -------------------------------------------------------------------
      DT2 = DT / 2.0
      DTSQ2 = DT2 * DT

!
!   IT COMPUTES THE QUANTITIES OF INTEREST
!
      A = 0.d0
      B = 0.d0
      D = 0.d0
      E = 0.d0
      F = 0.d0

      DO I = 1, N
         A = A + V(I)%X*V1(I)%X
         A = A + V(I)%Y*V1(I)%Y
         A = A + V(I)%Z*V1(I)%Z
         B = B + V(I)%X*V2(I)%X
         B = B + V(I)%Y*V2(I)%Y
         B = B + V(I)%Z*V2(I)%Z
         D = D - V1(I)%X*V2(I)%X/ M(I)
         D = D - V1(I)%Y*V2(I)%Y/ M(I)
         D = D - V1(I)%Z*V2(I)%Z/ M(I)
         E = E - V1(I)%X*V1(I)%X/ M(I)
         E = E - V1(I)%Y*V1(I)%Y/ M(I)
         E = E - V1(I)%Z*V1(I)%Z/ M(I)
         F = F - V2(I)%X*V2(I)%X/ M(I)
         F = F - V2(I)%Y*V2(I)%Y/ M(I)
         F = F - V2(I)%Z*V2(I)%Z/ M(I)
      ENDDO

       
      DEN = DTSQ2 * (D * D - F * E)
      LAM1 = (F * A - D * B) / DEN
      LAM2 = (E * B - D * A) / DEN

! 
!  CORRECTION OF VELOCITIES
! 


      DO I = 1 , N
         V(I)%X =  V(I)%X - DTSQ2 * (V1(I)%X * LAM1 + V2(I)%X * LAM2) / M(I)  
         V(I)%Y =  V(I)%Y - DTSQ2 * (V1(I)%Y * LAM1 + V2(I)%Y * LAM2) / M(I)  
         V(I)%Z =  V(I)%Z - DTSQ2 * (V1(I)%Z * LAM1 + V2(I)%Z * LAM2) / M(I)  
      END DO
      RETURN
!     -------------------------------------------------------------------
      END SUBROUTINE CONSTRB
!     ===================================================================    

END MODULE VERLET_CONSTR
