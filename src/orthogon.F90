!------------------------------------------------------------------------
MODULE VECTORALGEBRA                                                    !
      IMPLICIT NONE                                                     !
!     ------------------------------------------------------------------!
      CONTAINS                                                          !
!     ===================================================================
      SUBROUTINE GRAMSCHMIDT(A,DIMVEC,NVEC)
!     ===================================================================
      IMPLICIT NONE                                                     !
      DOUBLE PRECISION, PARAMETER :: THRS=1.D-15                        !
      INTEGER, INTENT(IN) :: DIMVEC, NVEC                               !
      DOUBLE PRECISION, INTENT(INOUT) :: A(DIMVEC,*)                    !
      DOUBLE PRECISION :: SUM, SUM2, NORM                               !
      INTEGER          :: N, M, I, J, K                                 !
                                                                        !
      DO K=1,NVEC-1                                                     !
         SUM = 0.D0                                                     !
         DO I=1,DIMVEC                                                  !
            SUM = SUM + A(I, K)*A(I, K)                                 !
         END DO                                                         !
         NORM=SQRT(SUM)                                                 !
         IF (NORM.LT.THRS) NORM=1.D0                                    !
         DO I=1,DIMVEC                                                  !
            A(I, K) = A(I, K) / NORM                                    !
         END DO                                                         !
         DO J=K+1,NVEC                                                  !
            SUM=0.D0                                                    !
            SUM2=0.D0                                                   !
            DO I=1,DIMVEC                                               !
               SUM = SUM + A(I, K) * A(I, J)                            !
            END DO                                                      !
            DO I=1,DIMVEC                                               !
               A(I, J) = A(I, J) - A(I, K) * SUM                        !
               SUM2 = SUM2 + A(I, J) * A(I, J)                          !
            END DO                                                      !
            NORM=SQRT(SUM2)                                             !
            IF (NORM.LT.THRS) NORM=1.D0                                 !
            DO I=1,DIMVEC                                               !
               A(I, J) = A(I, J) / NORM                                 !
            END DO                                                      !
         END DO                                                         !
      END DO                                                            !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE GRAMSCHMIDT                                        !
!     ===================================================================
      SUBROUTINE PROJECTOUT(V1,A,DIMVEC)      
!     ===================================================================
      IMPLICIT NONE                                                     !
      INTEGER, INTENT(IN) :: DIMVEC                                     !
      DOUBLE PRECISION, INTENT(INOUT) :: A(DIMVEC,*), V1(*)             !
      DOUBLE PRECISION :: SUMG, SUMH                                    !
      INTEGER :: I                                                      !
                                                                        !
      SUMG = 0.D0                                                       !
      SUMH = 0.D0                                                       !
      DO I=1,DIMVEC                                                     !
         SUMG = SUMG + V1(I) * A(I,1)                                   !
         SUMH = SUMH + V1(I) * A(I,2)                                   !
      END DO                                                            !
                                                                        !
      DO I=1, DIMVEC                                                    !
         V1(I) = V1(I) - SUMG * A(I,1) - SUMH * A(I,2)                  !
      END DO                                                            !
                                                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE PROJECTOUT                                         !
!     ==================================================================!
      
END MODULE VECTORALGEBRA
