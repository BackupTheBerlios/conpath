!-----------------------------------------------------------------------!
MODULE GEOMODULE                                                        !
                                                                        !
      IMPLICIT NONE                                                     !
!     ------------------------------------------------------------------!     
      CONTAINS                                                          !
!     ===================================================================
      SUBROUTINE CREATE_GEOFILE ( GEOFILE )
!     ===================================================================
      USE DYNPREPARE,  ONLY: C,               &                         !
                             LAB                                        !
      USE START_JOB,   ONLY: NUMAT,           &                         !
                             ABINIPROG                                  !
      USE SYSTEM_UTIL, ONLY: CPSTOP,          &                         !
                             EATW,            &                         !
                             LFBL,            &                         !
                             CP_OPEN,         &                         !
                             CP_CLOSE,        &                         !
                             PSYST,           &                         !
                             COMMAND                                    !
      USE XYZ,         ONLY: POINT                                      !
      USE PARSING,     ONLY: CLEANSTRING,     &                         !  
                             FREESPACE                                  !
      USE CONVFACTORS, ONLY: CV => BOHRTOANG                            !
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      INTEGER             :: ICHANNEL,  I                               !
      CHARACTER (LEN=*)   :: GEOFILE                                    !
      CHARACTER (LEN=256) :: LFILE,LFILE2                               !
                                                                        !
      CALL CLEANSTRING(LFILE,LEN(LFILE))                                !
      CALL CLEANSTRING(LFILE2,LEN(LFILE2))                              !
      SELECT CASE (ABINIPROG)                                           !
      CASE ('MOLPRO    ')                                               !
         COMMAND='rm -f '//GEOFILE                                      !
         WRITE(6,'(A)')' DELETING OLD GEOFILE AND CREATING THE NEW ONE!'! 
         CALL PSYST                                                     !
         ICHANNEL=14                                                    !
         CALL CP_OPEN(ICHANNEL,GEOFILE,'UNKNOWN','FORMATTED')           !
         WRITE(ICHANNEL,'(A)')'geometry={nosym;noorient;'               !
                                                                        !
         WRITE(LFILE,'(I5,A2)')NUMAT,';'                                !
         CALL FREESPACE(LFILE,LFILE2,LEN(LFILE),';')                    !
         WRITE(ICHANNEL,'(A)')LFILE2(1:LFBL(LFILE2,LEN(LFILE2)))        !
         WRITE(ICHANNEL,'(A)')   &                                      !
                      ' XYZ Geometry - ConPath program - Angs unit ;'   !
         DO I=1,NUMAT                                                   !
            WRITE(LFILE,10)LAB(I),C(I)%X*CV,C(I)%Y*CV,C(I)%Z*CV         !
            CALL CLEANSTRING(LFILE2,LEN(LFILE2))                        !
            CALL FREESPACE(LFILE,LFILE2,LEN(LFILE),';')                 !
            WRITE(ICHANNEL,'(A)')LFILE2(1:LFBL(LFILE2,LEN(LFILE2)))     ! 
         END DO                                                         !
         WRITE(ICHANNEL,'(A)')'}'                                       !
         CALL CP_CLOSE(ICHANNEL, GEOFILE)                               !
      CASE DEFAULT                                                      !
         WRITE(6,'(A)')'ABINITIO PROGRAM NOT IMPLEMENTED ON CONPATH!'   !
         CALL CPSTOP('GEOMODULE')                                       !
      END SELECT                                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 10   FORMAT(A2,', ',F15.9,', ',F15.9,', ',F15.9,'; ')                  !
!     ------------------------------------------------------------------!
      END SUBROUTINE CREATE_GEOFILE                                     !
!     ===================================================================

END MODULE GEOMODULE
