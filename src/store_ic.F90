!------------------------------------------------------------------------
MODULE CONICAL                                                          !
      USE XYZ,          ONLY: POINT                                     !
      USE START_JOB,    ONLY: NUMAT                                     !
      USE SYSTEM_UTIL,  ONLY: CPSTOP                                    !
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      TYPE(POINT), ALLOCATABLE :: V1(:,:), V2(:,:),CSAVE(:,:)           !
      TYPE(POINT), ALLOCATABLE :: V1_TMP(:,:), V2_TMP(:,:),C_TMP(:,:)   !
!     ------------------------------------------------------------------!
      CONTAINS                                                          !
!     ===================================================================
      SUBROUTINE STORE_IC ( NIC )
!     ===================================================================
      USE DYNPREPARE,  ONLY: FSAVE,               &                     !
                             A,                   &                     !
                             C                                          !
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      INTEGER :: NDIM, I, J                                             !
      INTEGER, INTENT(IN) :: NIC                                        !
                                                                        !
      NDIM = NUMAT                                                      !
                                                                        !
      IF (NIC.EQ.1) THEN                                                !
         ALLOCATE (  V1(NDIM, NIC), &                                   !
                     V2(NDIM, NIC), &
                  CSAVE(NDIM, NIC)     )
                                                                        !
         DO J=1,NUMAT                                                   !
            V1(J,1)%X = FSAVE(J)%X                                      !
            V1(J,1)%Y = FSAVE(J)%Y                                      !
            V1(J,1)%Z = FSAVE(J)%Z                                      !
                                                                        !
            V2(J,1)%X = A(J)%X                                          !
            V2(J,1)%Y = A(J)%Y                                          !
            V2(J,1)%Z = A(J)%Z                                          !

            CSAVE(J,1)%X = C(J)%X                         
            CSAVE(J,1)%Y = C(J)%Y                         
            CSAVE(J,1)%Z = C(J)%Z                                 

         END DO                                                         !
                                                                        !
      ELSEIF (NIC.GT.1) THEN                                            !
         ALLOCATE (  V1_TMP(NDIM, NIC-1), &                             !
              V2_TMP(NDIM, NIC-1),        &
              C_TMP (NDIM, NIC-1)   )                                   !
                                                                        !
         DO I=1,NIC-1                                                   !
            DO J=1,NUMAT                                                !
               V1_TMP(J,I)%X=V1(J,I)%X                                  !
               V1_TMP(J,I)%Y=V1(J,I)%Y                                  !
               V1_TMP(J,I)%Z=V1(J,I)%Z                                  !
                                                                        !
               V2_TMP(J,I)%X=V2(J,I)%X                                  !
               V2_TMP(J,I)%Y=V2(J,I)%Y                                  !
               V2_TMP(J,I)%Z=V2(J,I)%Z                                  !

               C_TMP(J,I)%X=CSAVE(J,I)%X 
               C_TMP(J,I)%Y=CSAVE(J,I)%Y 
               C_TMP(J,I)%Z=CSAVE(J,I)%Z 

            END DO                                                      !
         END DO                                                         !
                                                                        !
         DEALLOCATE ( V1, V2, CSAVE )  
         ALLOCATE   ( V1(NDIM, NIC), &                                  !
              V2(NDIM, NIC),         &
              CSAVE(NDIM, NIC) )                                        !
                                                                        !
         DO I=1,NIC-1                                                   !
            DO J=1,NUMAT                                                !
               V1(J,I)%X=V1_TMP(J,I)%X                                 !
               V1(J,I)%Y=V1_TMP(J,I)%Y                                 !
               V1(J,I)%Z=V1_TMP(J,I)%Z                                 !
                                                                        !
               V2(J,I)%X=V2_TMP(J,I)%X                                  !
               V2(J,I)%Y=V2_TMP(J,I)%Y                                  !
               V2(J,I)%Z=V2_TMP(J,I)%Z                                  !

               CSAVE(J,I)%X=C_TMP(J,I)%X
               CSAVE(J,I)%Y=C_TMP(J,I)%Y
               CSAVE(J,I)%Z=C_TMP(J,I)%Z
 
            END DO                                                      !
         END DO                                                         !
                                                                        !
         DEALLOCATE ( V1_TMP, V2_TMP, C_TMP )  
                                                                        !
         DO J=1,NUMAT                                                   !
            V1(J,NIC)%X = FSAVE(J)%X                                    !
            V1(J,NIC)%Y = FSAVE(J)%Y                                    !
            V1(J,NIC)%Z = FSAVE(J)%Z                                    !
                                                                        !
            V2(J,NIC)%X = A(J)%X                                        !
            V2(J,NIC)%Y = A(J)%Y                                        !
            V2(J,NIC)%Z = A(J)%Z                                        !

            CSAVE(J,NIC)%X=C(J)%X
            CSAVE(J,NIC)%Y=C(J)%Y
            CSAVE(J,NIC)%Z=C(J)%Z             

         END DO                                                         !
                                                                        !
      ELSE                                                              !
         WRITE(6,'(A)')'ERROR IN NIC: NUMBER OF CONICAL INTERSECTIONS!' !
         CALL CPSTOP('STORE_IC')                                        !
      END IF                                                            !
      WRITE(6,10)NIC                                                    !
      RETURN                                                            !
10    FORMAT(1X,'NUMBER OF CONICAL INTERSECTION FOUND = ',I5)           !
!     ------------------------------------------------------------------!
      END SUBROUTINE STORE_IC                                           !
!     ===================================================================
END MODULE CONICAL
