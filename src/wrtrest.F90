!------------------------------------------------------------------------
MODULE RESTART                                                          !
      USE UNITS_FILE, ONLY: RESFILE                                     !
      USE MODPROPERTIES, ONLY: COMPUTE_PROPERTIES                       !
      USE DYNPREPARE,    ONLY: M,   &                                   !
                               ASSIGNMASS                               !
      USE START_JOB,     ONLY: IPRINT,      &                           !
                               PROPERTIES,  &                           !
                               PROP_AND_STOP                            !
      IMPLICIT NONE                                                     !
      INTEGER, PARAMETER :: IBUFFERS=4                                  ! 
      INTEGER  :: RSTEP,RN_IC                                           ! 
!     ------------------------------------------------------------------!
      CONTAINS                                                          !
!     ===================================================================
      SUBROUTINE WRITE_RESTART ( STEP , N_IC )
!     ===================================================================
      USE DYNPREPARE, ONLY:    C,       &                               !
                               V,       &                               !
                               LAB                                      !
      USE START_JOB,  ONLY:    NUMAT                                    !
      USE CONICAL,    ONLY:    V1, V2, CSAVE                            !
!     ------------------------------------------------------------------!
!     Structure of Restart File :                                       !
!     When you modify this subroutine remember to increase the Ibuffers !
!     variable....                                                      !
!                                                                       !
!     Ibuffer = 1   : 1, Numat, Step                                    !
!     Ibuffer = 2   : 2, Labels                                         !
!     Ibuffer = 3   : 3, Coordinates                                    !
!     Ibuffer = 4   : 4, Velocities                                     !
!     Ibuffer = 5   : 5, Number of conical intersections found          !
!     Ibuffer = 6   : 6, V1                                             !
!     Ibuffer = 7   : 7, V2                                             !
!     Ibuffer = 8   : 8, CSAVE                                          !
!     EndofBuffer   : 999999999                                         !
!                                                                       !
!     ------------------------------------------------------------------!
      INTEGER, INTENT(IN) :: STEP, N_IC                                 ! 
      INTEGER :: I ,J                                                   ! 
                                                                        !
      WRITE(RESFILE)1,NUMAT,STEP                                        !
      WRITE(RESFILE)2,(LAB(I),I=1,NUMAT)                                !
      WRITE(RESFILE)3,(C(I)%X,C(I)%Y,C(I)%Z,I=1,NUMAT)                  !
      WRITE(RESFILE)4,(V(I)%X,V(I)%Y,V(I)%Z,I=1,NUMAT)                  !
      WRITE(RESFILE)5,N_IC                                              !
      IF (N_IC.NE.0)THEN                                                !
         WRITE(RESFILE)6,         &                                     !
             ((V1(I,J)%X,V1(I,J)%Y,V1(I,J)%Z,I=1,NUMAT),J=1,N_IC)       !
         WRITE(RESFILE)7,         &                                     !
             ((V2(I,J)%X,V2(I,J)%Y,V2(I,J)%Z,I=1,NUMAT),J=1,N_IC)       !
         WRITE(RESFILE)8,         &                                     !
          ((CSAVE(I,J)%X,CSAVE(I,J)%Y,CSAVE(I,J)%Z,I=1,NUMAT),J=1,N_IC) !
      ENDIF                                                             !
!     ------------------------------------------------------------------!
!     End of buffer                                                     !
!     ------------------------------------------------------------------!
      WRITE(RESFILE)999999999                                           !
!     ------------------------------------------------------------------!
!     Compute Properties                                                !
!     ------------------------------------------------------------------!
      IF (PROPERTIES.NE.0)  &                                           !
           CALL COMPUTE_PROPERTIES(C,V,NUMAT,INT(IPRINT,4),STEP)        !
!     ------------------------------------------------------------------!
      CALL FLUSH(RESFILE)
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE WRITE_RESTART                                      !
!     ===================================================================
      SUBROUTINE READ_RESTART 
!     ===================================================================
      USE DYNPREPARE,  ONLY:    C,       &                              !
                                V,       &                              !
                                LAB                                     !
      USE START_JOB,   ONLY:    N => NUMAT                              !
      USE SYSTEM_UTIL, ONLY:    CPSTOP                                  !
      USE CONICAL, ONLY:        V1, V2, CSAVE                           !
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      INTEGER (KIND=8) :: RNATOM                                        ! 
      INTEGER :: IBUF, ENDBUFF,I,J,NTIMES,IFIRST                        ! 
      NTIMES=0                                                          !
      IFIRST=0                                                          !
                                                                        !
                                                                        !
      DO WHILE (.TRUE.)   ! Not looping.. don't worry... ;)             !
         READ(RESFILE,END=99)IBUF,RNATOM,RSTEP                          !
         IF (RNATOM.NE.N) GO TO 999                                     !
         READ(RESFILE,END=999)IBUF,(LAB(I),I=1,N)                       !
         IF (IFIRST.EQ.0) THEN                                          !
            CALL ASSIGNMASS(LAB, M, N)                                  !
            IFIRST=1                                                    !
         END IF                                                         !
         READ(RESFILE,END=999)IBUF,(C(I)%X,C(I)%Y,C(I)%Z,I=1,N)         !
         READ(RESFILE,END=999)IBUF,(V(I)%X,V(I)%Y,V(I)%Z,I=1,N)         !
         READ(RESFILE,END=999)IBUF,RN_IC                                !
         IF (RN_IC.NE.0)THEN                                            !
              IF (NTIMES.EQ.0)THEN                                      !
                 ALLOCATE    (V1(RNATOM, RN_IC), &                      !
                 V2(RNATOM, RN_IC))                                     !
                 ALLOCATE  (CSAVE (RNATOM, RN_IC))                      !
              ELSE                                                      !
                 DEALLOCATE  (V1, V2, CSAVE)                            !
                 ALLOCATE   (V1(RNATOM, RN_IC), &                       !
                 V2(RNATOM, RN_IC))                                     !
                 ALLOCATE  (CSAVE (RNATOM, RN_IC))                      !
              ENDIF                                                     !
              NTIMES=1                                                  !
            READ(RESFILE,END=999)IBUF,           &                      !
            ((V1(I,J)%X,V1(I,J)%Y,V1(I,J)%Z,I=1,RNATOM),J=1,RN_IC)      !
            READ(RESFILE,END=999)IBUF,           &                      !
            ((V2(I,J)%X,V2(I,J)%Y,V2(I,J)%Z,I=1,RNATOM),J=1,RN_IC)      !
            READ(RESFILE,END=999)IBUF,           &                      !
      ((CSAVE(I,J)%X,CSAVE(I,J)%Y,CSAVE(I,J)%Z,I=1,RNATOM),J=1,RN_IC)   !
         ENDIF                                                          !
         READ(RESFILE,END=999)ENDBUFF                                   !
!     ------------------------------------------------------------------!
!     Compute Properties                                                !
!     ------------------------------------------------------------------!
         IF (PROPERTIES.NE.0)  &                                        !
              CALL COMPUTE_PROPERTIES(C,V,N,INT(IPRINT,4),RSTEP)        !
         IF (ENDBUFF.NE.999999999) GO TO 999                            !
      END DO                                                            !
                                                                        !
 99   WRITE(6,'(1X,A//)')   &                                           !
      'RESTART FILE READ! GOING ON WITH THE ORIGINAL JOB!'              !
      IF (PROP_AND_STOP) THEN                                           !
         WRITE(6,'(A)')'SIMPLE PROPERTIES CALCULATION ON RESTART FILE!' !
         WRITE(6,'(A)')'NORMAL TERMINATION OF CONPATH PROGRAM...'       !
         STOP 0                                                         !
      END IF                                                            !
      WRITE(6,'(80A)')('*',I=1,80)                                      !
      RETURN                                                            !
!     ------------------------------------------------------------------!
999   WRITE(6,'(A)')'RESTART FILE CORRUPTED!'                           !
      CALL CPSTOP('READ_RESTART')                                       !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE READ_RESTART                                       !
!     ===================================================================

END MODULE RESTART
