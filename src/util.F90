!------------------------------------------------------------------------
MODULE EMPTYLINE                                                        !
      IMPLICIT NONE                                                     !
      PUBLIC                                                            !
      CONTAINS                                                          !
!     ===================================================================
      LOGICAL FUNCTION EMPTY(STRING,LENS)                               
!     =================================================================== 
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=*), INTENT(IN) :: STRING                           !
      INTEGER, INTENT(IN) :: LENS                                       !
      INTEGER :: I                                                      !
                                                                        !
      EMPTY=.TRUE.                                                      !
      DO I=1,LENS                                                       !
         IF (STRING(I:I).NE.' ') EMPTY=.FALSE.                          !
      END DO                                                            !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END FUNCTION EMPTY                                                !
!     =================================================================== 
END MODULE EMPTYLINE

!------------------------------------------------------------------------
MODULE SYSTEM_UTIL                                                      !
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=256) :: COMMAND                                    !
!     ------------------------------------------------------------------!
      CONTAINS
!     ===================================================================
      SUBROUTINE PSYST                                                  
!     ===================================================================
      IMPLICIT NONE                                                     !
      CALL SYSTEM (COMMAND)                                             !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE PSYST                                              !
!     ===================================================================
      SUBROUTINE GOABINITIO(ABINIEXE,ABINIT_INP)
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=256), INTENT(IN):: ABINIEXE, ABINIT_INP            !
      INTEGER :: I,K                                                    !
                                                                        !
      I=INDEX(ABINIEXE,' ')  -1                                         !
      K=INDEX(ABINIT_INP,' ')-1                                         !
                                                                        !
      COMMAND=ABINIEXE(1:I)//' '//ABINIT_INP(1:K)                       ! 
      CALL PSYST                                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE GOABINITIO                                         !
!     ===================================================================
      SUBROUTINE CP_OPEN_SYS(CHANNEL,FILE,STATUS,FORM,POSIT)
!     ===================================================================
      CHARACTER (LEN=*), INTENT(IN) :: FILE, STATUS, FORM, POSIT        !
      INTEGER, INTENT(IN) :: CHANNEL                                    !
      INTEGER :: I                                                      !
                                                                        !
      I=INDEX(FILE,' ')-1                                               !
      WRITE(6,'(3A)')' OPENING FILE : ',FILE(1:I),' .'                  !
      OPEN(CHANNEL,FILE=FILE(1:I),STATUS=STATUS,FORM=FORM,   &          !
           POSITION=POSIT,ERR=99)                                       !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 99   WRITE(6,'(3A)')' ERROR IN OPENING FILE: ',FILE(1:I),' !'          !
      CALL CPSTOP('CP_OPEN_SYS')                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE CP_OPEN_SYS                                        !
!     =================================================================== 
      SUBROUTINE CP_OPEN(CHANNEL,FILE,STATUS,FORM)
!     ===================================================================
      CHARACTER (LEN=*), INTENT(IN) :: FILE, STATUS, FORM               !
      INTEGER, INTENT(IN) :: CHANNEL                                    !
      INTEGER :: I                                                      !
                                                                        !
      I=INDEX(FILE,' ')-1                                               !
      WRITE(6,'(3A)')' OPENING FILE : ',FILE(1:I),' .'                  !
      OPEN(CHANNEL,FILE=FILE(1:I),STATUS=STATUS,FORM=FORM,ERR=99)       !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 99   WRITE(6,'(3A)')' ERROR IN OPENING FILE: ',FILE(1:I),' !'          !
      CALL CPSTOP('CP_OPEN')                                            !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE CP_OPEN                                            !
!     =================================================================== 
      SUBROUTINE CP_CLOSE(CHANNEL,FILE)
!     ===================================================================
      CHARACTER (LEN=*), INTENT(IN) :: FILE                             !
      INTEGER, INTENT(IN) :: CHANNEL                                    ! 
      INTEGER :: I                                                      !
                                                                        !
      I=INDEX(FILE,' ')-1                                               !
      WRITE(6,'(3A)')' CLOSING FILE : ',FILE(1:I),' .'                  !
      CLOSE(CHANNEL)                                                    !
!     ------------------------------------------------------------------!
      END SUBROUTINE CP_CLOSE                                           !
!     =================================================================== 
      SUBROUTINE CPSTOP(KEY)
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=*), INTENT(IN) :: KEY                              !
                                                                        !
      WRITE(6,'(2A)')'EXECUTION ABORTED IN SUBROUTINE: ',KEY            !
      WRITE(6,'(A)')'STOPPING NOW!'                                     !
      STOP 999                                                         !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE CPSTOP                                             !
!     ===================================================================  
      SUBROUTINE SEARCH(TARGET,ICHANNEL,REWIND,DIRECTION,CONFERMA,LINE)
!     ===================================================================
!     TARGET :: STRING TO SEARCH IN FILE CORRESPONDING TO ICHANNEL      !
!     REWIND :: .TRUE. OR .FALSE. WHETHER TO REWIND OR NOT THE FILE     !
!     DIRECTION :: "AHEA" or "BACK"                                     !
!     CONFERMA :: .TRUE. IF SEARCH FOUND THE TARGET                     !
!     ------------------------------------------------------------------!
      USE EMPTYLINE                                                     !
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=20), INTENT(IN)  :: TARGET                         !
      CHARACTER (LEN=256),INTENT(OUT) :: LINE                           !
      CHARACTER (LEN=256) :: LINEOLD                                    !
      INTEGER, INTENT(IN) :: ICHANNEL                                   !
      LOGICAL, INTENT(IN) :: REWIND                                     !
      CHARACTER (LEN=4), INTENT(IN) :: DIRECTION                        !
      LOGICAL, INTENT(INOUT) :: CONFERMA                                !
      INTEGER :: IIND, I                                                !
                                                                        !
      IF (REWIND) REWIND ICHANNEL                                       !
                                                                        !
      SELECT CASE (DIRECTION)                                           !
      CASE ('AHEA')                                                     !
         DO WHILE (.NOT.CONFERMA)                                       !
            READ(ICHANNEL,'(A)',END=99,ERR=999)LINE                     !
            IIND=INDEX(LINE,TARGET(1:LFBL(TARGET,LEN(TARGET))))         !
            IF (IIND.NE.0) CONFERMA=.TRUE.                              !
         END DO                                                         !
      CASE ('BACK')                                                     !
         DO WHILE (.NOT.CONFERMA)                                       !
            DO I=1,2                                                    !
               BACKSPACE (ICHANNEL,ERR=999)                             !
            END DO                                                      !
            LINEOLD=LINE                                                !
            READ(ICHANNEL,'(A)')LINE                                    !
            IF ((LINEOLD.EQ.LINE).AND.(.NOT.EMPTY(LINE,256)))  GO TO 99 !
            IIND=INDEX(LINE,TARGET(1:LFBL(TARGET,LEN(TARGET))))         !
            IF (IIND.NE.0) CONFERMA=.TRUE.                              !
         END DO                                                         !
      CASE DEFAULT                                                      !
         WRITE(6,'(A)')'SEARCH IS POSSIBLE ONLY BACKWARD AND FORWARD!'  !
         CALL CPSTOP('SEARCH')                                          !
      END SELECT                                                        !
                                                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 99   CONTINUE                                                          !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 999  WRITE(6,'(A)')'AN UNEXPECTED ERROR OCCURED WHILE SEARCHING ...'   !
      CALL CPSTOP('SEARCH')                                             !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE SEARCH                                             !
!     =================================================================== 
!     ===================================================================  
      SUBROUTINE SEARCHT(TARGET1,TARGET2,TARGET3,ICHANNEL,REWIND, &
                         DIRECTION,CONFERMA,LINE)
!     ===================================================================
!     TARGET :: STRING TO SEARCH IN FILE CORRESPONDING TO ICHANNEL      !
!     REWIND :: .TRUE. OR .FALSE. WHETHER TO REWIND OR NOT THE FILE     !
!     DIRECTION :: "AHEA" or "BACK"                                     !
!     CONFERMA :: .TRUE. IF SEARCH FOUND THE TARGET                     !
!     ------------------------------------------------------------------!
      USE EMPTYLINE                                                     !
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=10), INTENT(IN) :: TARGET1,TARGET2,TARGET3         !
      CHARACTER (LEN=256),INTENT(OUT) :: LINE                           !
      CHARACTER (LEN=256) :: LINEOLD                                    !
      INTEGER, INTENT(IN) :: ICHANNEL                                   !
      LOGICAL, INTENT(IN) :: REWIND                                     !
      CHARACTER (LEN=4), INTENT(IN) :: DIRECTION                        !
      LOGICAL, INTENT(INOUT) :: CONFERMA                                !
      INTEGER :: IIND, I                                                !
                                                                        !
      IF (REWIND) REWIND ICHANNEL                                       !
                                                                        !
      SELECT CASE (DIRECTION)                                           !
      CASE ('AHEA')                                                     !
         DO WHILE (.NOT.CONFERMA)                                       !
            READ(ICHANNEL,'(A)',END=99,ERR=999)LINE                     !
            IIND=INDEX(LINE,TARGET1(1:LFBL(TARGET1,10)))                !
            IF (IIND.NE.0) THEN                                         !
               IIND=INDEX(LINE,TARGET2(1:LFBL(TARGET2,10)))             !
               IF (IIND.NE.0) THEN                                      !
                  IIND=INDEX(LINE,TARGET3(1:LFBL(TARGET3,10)))          !
                  IF (IIND.NE.0) CONFERMA=.TRUE.                        !
               END IF                                                   !
            END IF                                                      !
         END DO                                                         !
      CASE ('BACK')                                                     !
         DO WHILE (.NOT.CONFERMA)                                       !
            DO I=1,2                                                    !
               BACKSPACE (ICHANNEL,ERR=999)                             !
            END DO                                                      !
            LINEOLD=LINE                                                !
            READ(ICHANNEL,'(A)')LINE                                    !
            IF ((LINEOLD.EQ.LINE).AND.(.NOT.EMPTY(LINE,256)))  GO TO 99 !
            IIND=INDEX(LINE,TARGET1(1:LFBL(TARGET1,10)))                !
            IF (IIND.NE.0) THEN                                         !
               IIND=INDEX(LINE,TARGET2(1:LFBL(TARGET2,10)))             !
               IF (IIND.NE.0) THEN                                      !
                  IIND=INDEX(LINE,TARGET3(1:LFBL(TARGET3,10)))          !
                  IF (IIND.NE.0) CONFERMA=.TRUE.                        !
               END IF                                                   !
            END IF                                                      !            
         END DO                                                         !
      CASE DEFAULT                                                      !
         WRITE(6,'(A)')'SEARCH IS POSSIBLE ONLY BACKWARD AND FORWARD!'  !
         CALL CPSTOP('SEARCHT')                                         !
      END SELECT                                                        !
                                                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 99   CONTINUE                                                          !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 999  WRITE(6,'(A)')'AN UNEXPECTED ERROR OCCURED WHILE SEARCHING ...'   !
      CALL CPSTOP('SEARCH')                                             !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE SEARCHT                                            !
!     ===================================================================
      CHARACTER (LEN=10) FUNCTION GETLABEL(LINE,CHAR,LENS)                        
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=*) :: LINE                                         !
      CHARACTER (LEN=1) :: CHAR                                         !
      INTEGER :: IND                                                    !
      INTEGER, INTENT(IN) :: LENS                                       !
                                                                        !
      CALL EATW(LINE,LENS)                                              !
      IND=INDEX(LINE,CHAR)                                              !
      IF (IND.EQ.0) THEN                                                !
         GETLABEL='NULL'                                                !
      ELSE                                                              !
         GETLABEL=LINE(1:IND-1)                                         !
      END IF                                                            !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END FUNCTION GETLABEL                                             !
!     ===================================================================
      SUBROUTINE EATW(LINE,LENS)
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=*) :: LINE                                         !
      INTEGER :: I                                                      !
      INTEGER,  INTENT(IN) :: LENS                                      !
                                                                        !
      DO I=1,LENS                                                       !
         IF (LINE(I:I).NE.' ') THEN                                     !
            LINE(1:)=LINE(I:)                                           !
            EXIT                                                        !
         END IF                                                         !
      END DO                                                            !
      IF (I.EQ.LENS+1)  LINE='EMPTY'                                    !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE EATW                                               !
!     =================================================================== 
      SUBROUTINE EATB(LINE,LENS)
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=*) :: LINE                                         !
      INTEGER :: I                                                      !
      INTEGER,  INTENT(IN) :: LENS                                      !
                                                                        !
      DO I=1,LENS                                                       !
         IF (LINE(I:I).EQ.' ') THEN                                     !
            LINE(1:)=LINE(I:)                                           !
            EXIT                                                        !
         END IF                                                         !
      END DO                                                            !
      IF (I.EQ.LENS+1)  LINE='EMPTY'                                    !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE EATB                                               !
!     =================================================================== 
      INTEGER FUNCTION  LFBL(LINE,LENS)
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=*) :: LINE                                         !
      INTEGER :: I                                                      !
      INTEGER,  INTENT(IN) :: LENS                                      !
                                                                        !
      DO I=LENS,1,-1                                                    !
         IF (LINE(I:I).NE.' ') THEN                                     !
            EXIT                                                        !
         END IF                                                         !
      END DO                                                            !
      LFBL=I                                                            !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END FUNCTION  LFBL                                                !
!     =================================================================== 
      SUBROUTINE JLINE(ICHANNEL, NLINE, CONFERMA)
!     ===================================================================
      IMPLICIT NONE                                                     ! 
      INTEGER, INTENT(IN) :: ICHANNEL, NLINE                            !
      INTEGER             :: I                                          !
      CHARACTER (LEN=256) :: LINE                                       !
      LOGICAL             :: CONFERMA                                   !
                                                                        !
      DO I=1,NLINE                                                      !
         READ(ICHANNEL,'(A)',END=99,ERR=999)LINE                        !
      END DO                                                            !
                                                                        !
      CONFERMA=.TRUE.                                                   !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 99   WRITE(6,'(A,I5)')' REACHED END OF FILE RELATED TO CHANNEL:', &    !
                       ICHANNEL                                         !
      CONFERMA=.FALSE.                                                  !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 999  WRITE(6,'(A,I5)')' AN ERROR OCCURRED READING FILE ON CHANNEL:', & !
                       ICHANNEL                                         !
      CONFERMA=.FALSE.                                                  !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE JLINE                                              !
!     ===================================================================
END MODULE SYSTEM_UTIL



!------------------------------------------------------------------------
MODULE XYZ                                                              !
      IMPLICIT NONE                                                     !
!     ------------------------------------------------------------------!
      TYPE POINT                                                        !
        REAL (KIND=8) :: X,Y,Z                                          !
      END TYPE POINT                                                    !
!     ------------------------------------------------------------------!
      CONTAINS                                                          !
!     ===================================================================
      DOUBLE PRECISION FUNCTION PROD2(A,B,DIM)
!     ===================================================================
      IMPLICIT NONE                                                     !
      TYPE(POINT), INTENT(IN) :: A(*) , B(*)                            !
      INTEGER, INTENT(IN)     :: DIM                                    !
      INTEGER                 :: I                                      !
                                                                        !
      PROD2=0.D0                                                        !
      DO I=1,DIM                                                        !
         PROD2 = PROD2 + A(I)%X * B(I)%X                                !
         PROD2 = PROD2 + A(I)%Y * B(I)%Y                                !
         PROD2 = PROD2 + A(I)%Z * B(I)%Z                                !
      END DO                                                            !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END FUNCTION PROD2                                                !
!     ===================================================================
      DOUBLE PRECISION FUNCTION PROD3(A,B,C,DIM)
!     ===================================================================
      IMPLICIT NONE                                                     !
      TYPE(POINT), INTENT(IN) :: B(*),  C(*)                            !
      DOUBLE PRECISION, INTENT(IN) :: A(*)                              !
      INTEGER, INTENT(IN)     :: DIM                                    !
      INTEGER                 :: I                                      !
                                                                        !
      PROD3=0.D0                                                        !
      DO I=1,DIM                                                        !
         PROD3 = PROD3 + A(I) * B(I)%X * C(I)%X                         !
         PROD3 = PROD3 + A(I) * B(I)%Y * C(I)%Y                         !
         PROD3 = PROD3 + A(I) * B(I)%Z * C(I)%Z                         !
      END DO                                                            !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END FUNCTION PROD3                                                !
!     ===================================================================

      SUBROUTINE COPY_VEC ( N , V1 , V2 )
!     ===================================================================
      IMPLICIT NONE
      INTEGER (KIND=8)                   :: N, I
      TYPE (POINT), INTENT (IN)          :: V1 (*)
      TYPE (POINT), INTENT (OUT)         :: V2 (*)
 
      DO I=1,N
        V2 (I)%X = V1 (I)%X
        V2 (I)%Y = V1 (I)%Y
        V2 (I)%Z = V1 (I)%Z
      END DO
      RETURN
      END SUBROUTINE COPY_VEC
!     ===================================================================
END MODULE XYZ



!------------------------------------------------------------------------
MODULE PARSING                                                          !
      IMPLICIT NONE                                                     !
!     ------------------------------------------------------------------!
      CONTAINS                                                          !
!     ===================================================================       
      SUBROUTINE GATHEREXP(ICHANNEL,TAR1,TAR2,LINE,CONFERMA)
!     =================================================================== 
      USE SYSTEM_UTIL, ONLY: SEARCH,       &                            !
                             CPSTOP,       &                            !
                             LFBL,         &                            !
                             EATW                                       !
      IMPLICIT NONE                                                     !
      INTEGER, INTENT(IN)               :: ICHANNEL                     !
      CHARACTER (LEN=1), INTENT(IN)     :: TAR1, TAR2                   !
      CHARACTER (LEN=1024), INTENT(OUT) :: LINE                         !
      CHARACTER (LEN=256)               :: CLINE                        !
      CHARACTER (LEN=20)                :: SLINE                        !
      LOGICAL :: CONFERMA, CONF2, FOUND                                 !
      INTEGER :: I, IENT, INDT1, INDT2, IL1, IL2, LF, ILTEMP            !
                                                                        !
      DO I=1,LEN(LINE)                                                  !
          LINE(I:I)=' '    ! Write blank on line                        !
      END DO                                                            !
      DO I=1,LEN(CLINE)                                                 !
         CLINE(I:I)=' '    ! Write blank on cline                       !
      END DO                                                            !
      CONF2=.FALSE.                                                     !
      SLINE='                    '                                      !
      SLINE=TAR1                                                        !
      CALL SEARCH(SLINE,ICHANNEL,.FALSE.,'AHEA',CONF2,CLINE)            !
      IF (.NOT.CONF2) CALL CPSTOP('GATHEREXP')                          !
      INDT1=INDEX(CLINE,TAR1)                                           !
      INDT2=INDEX(CLINE(INDT1+1:),TAR2)                                 !
      FOUND=.FALSE.                                                     !
      IF (INDT2.NE.0) THEN                                              !
         LINE(1:)=CLINE(INDT1:INDT2)                                    !
      ELSE                                                              !
         LF=LFBL(CLINE,256)                                             !
         IL2=LF-INDT1+1                                                 !
         IF (CLINE(LF:LF).EQ.';') THEN                                  !
            LINE(1:IL2)=CLINE(INDT1:LF)                                 !
         ELSE                                                           !
            IL2=IL2+1                                                   !
            LINE(1:IL2)=CLINE(INDT1:LF)//';'                            !
         ENDIF                                                          !
cycle1:  DO WHILE (.NOT.FOUND)                                          !
            IL2=IL2+1                                                   !
            READ(ICHANNEL,'(A)')CLINE                                   !
            CALL EATW(CLINE,256)                                        !
            IF (INDEX(CLINE,'EMPTY').NE.0) CYCLE cycle1                 !
            INDT2=INDEX(CLINE,TAR2)                                     !
            IF (INDT2.NE.0) THEN                                        !
               FOUND=.TRUE.                                             !
               ILTEMP=INDT2                                             !
            ELSE                                                        !
               ILTEMP=LFBL(CLINE,256)                                   !
            ENDIF                                                       !
            IL1=IL2                                                     !
            IF ((CLINE(ILTEMP:ILTEMP).EQ.';').OR.FOUND) THEN            !
               IL2=ILTEMP+IL1-1                                         !
               LINE(IL1:IL2)=CLINE(1:ILTEMP)                            !
            ELSE                                                        !
               IL2=ILTEMP+IL1                                           !
               LINE(IL1:IL2)=CLINE(1:ILTEMP)//';'                       !
            END IF                                                      !
         END DO  cycle1                                                 !
      ENDIF                                                             !
      IF (IL2.GT.LEN(LINE)) THEN                                        !
         WRITE(6,'(A)')'RECORD LONGER THAN VARIABLE IN MEMORY!'         !
         RETURN                                                         !
      ENDIF                                                             !
      CONFERMA=.TRUE.                                                   !
      RETURN                                                            !
!     -------------------------------------------------------------------
      END SUBROUTINE GATHEREXP                                          !
!     ===================================================================       
      SUBROUTINE GETFIELD(FIELD,LFIELD,TAR1,TAR2,LINE,POS)
!     =================================================================== 
      USE        SYSTEM_UTIL, ONLY: EATW,      &                        !
                                    EATB,      &                        !
                                    CPSTOP                              !
      USE        EMPTYLINE,   ONLY: EMPTY                               !
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=*), INTENT(IN)    :: LINE                          !
      CHARACTER (LEN=1024)             :: CLINE                         !
      CHARACTER (LEN=*), INTENT(OUT)   :: FIELD                         !
      INTEGER :: LFIELD, ICHANNEL, ITAR1, ITAR2, ILONG, I               !
      CHARACTER (LEN=1), INTENT(IN)   :: TAR1,TAR2                      !
      INTEGER, INTENT(IN) :: POS                                        !
                                                                        !
      CLINE=LINE                                                        !
      CALL EATW(CLINE,1024)                                             !
      ITAR1=0                                                           !
      ITAR2=0                                                           !
      DO I=1,POS                                                        !
         ITAR1=INDEX(LINE(ITAR1+1:),TAR1)+ITAR1+1                       !
         IF (ITAR1.EQ.0) THEN                                           !
            WRITE(6,'(3A)')'CANNOT FIND DELIMITER "',TAR1,'" ON LINE:'  !
            WRITE(6,'(3A)')LINE                                         !
            CALL CPSTOP('GETFIELD')                                     !
         ENDIF                                                          !
         ITAR2=INDEX(LINE(ITAR1+1:),TAR2)                               !
         IF (ITAR2.EQ.0) THEN                                           !
            WRITE(6,'(3A)')'CANNOT FIND DELIMITER "',TAR2,'" ON LINE:'  !
            WRITE(6,'(3A)')LINE(ITAR1+1:)                               !
            CALL CPSTOP('GETFIELD')                                     !
         ENDIF                                                          !
         ILONG=ITAR2-1                                                  !
         IF (ILONG.GT.LFIELD) THEN                                      !
            WRITE(6,'(A)')'FIELD BEETWEN "',TAR1,'" AND "',TAR2, &      !
            ' IS LONGER THAN DEFAULT VARIABLE!'                         !
            CALL CPSTOP('GETFIELD')                                     !
         ENDIF                                                          !
      END DO                                                            !
                                                                        !
      FIELD=CLINE(ITAR1:ITAR1+ITAR2-1)                                  !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE GETFIELD                                           !
!     ===================================================================
      SUBROUTINE CLEANSTRING(STRING,LSTR)
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER  (LEN=*), INTENT(INOUT)  :: STRING                      !
      INTEGER, INTENT(IN)                :: LSTR                        !
      INTEGER                            :: I                           !
                                                                        !
      DO I=1,LSTR                                                       !
         STRING(I:I)=' '                                                !
      END DO                                                            !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE CLEANSTRING                                        !
!     ===================================================================
      SUBROUTINE FREESPACE(STRING,STRING2,LSTR,TEND)
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=*), INTENT(IN)    :: STRING                        !
      CHARACTER (LEN=*), INTENT(OUT)   :: STRING2                       !
      CHARACTER (LEN=*), INTENT(IN)    :: TEND                          !
      INTEGER, INTENT(IN)              :: LSTR                          !
      INTEGER                          :: I, J, IFIN                    !
                                                                        !
      IFIN=LSTR                                                         !
      IF (TEND.NE.'$') IFIN=INDEX(STRING,TEND)                          !
                                                                        !
      J=0                                                               !
      DO I=1,IFIN                                                       !
         IF (STRING(I:I).NE.' ') THEN                                   !
            J=J+1                                                       !
            STRING2(J:J)=STRING(I:I)                                    !
         END IF                                                         !
      END DO                                                            !
      STRING2(J+1:)=STRING(I:)                                          !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE FREESPACE                                          !
!     ===================================================================

END MODULE PARSING

!------------------------------------------------------------------------
MODULE CONVFACTORS                                                      !
      IMPLICIT NONE                                                     ! 
                                                                        !
      DOUBLE PRECISION, PARAMETER :: BOHRTOANG=0.529177249D0            !
      DOUBLE PRECISION, PARAMETER :: ANGTOBOHR=1.D0 / BOHRTOANG         !
      DOUBLE PRECISION, PARAMETER :: HARTREETOEV=27.2113962D0           !
      DOUBLE PRECISION, PARAMETER :: EVTOKCAL=23.060542301389D0         !
      DOUBLE PRECISION, PARAMETER :: FSTOAUTIME=41.34137221718D0        !
      DOUBLE PRECISION, PARAMETER :: AMUTOAU=1822.888515D0              !
      DOUBLE PRECISION, PARAMETER :: KBOLTZ=1.380658D-23  ! (J/K)       !
      DOUBLE PRECISION, PARAMETER :: NAVOGADRO=6.0221D-23               !
      DOUBLE PRECISION, PARAMETER :: EVTOKELVIN=1.1604D4                !
      DOUBLE PRECISION, PARAMETER :: CALTOJOULE=4.1868D0                !
      DOUBLE PRECISION, PARAMETER :: AUTOKELVIN=HARTREETOEV*EVTOKELVIN  !
      DOUBLE PRECISION, PARAMETER :: RADTODEGREE=57.295779513082D0      !
END MODULE CONVFACTORS                                                  !


 
!------------------------------------------------------------------------
MODULE PERIODICTABLE                                                    !
      IMPLICIT NONE                                                     !
      CHARACTER  (LEN=2), ALLOCATABLE  ::  ELEMENTS(:)                  !
      DOUBLE PRECISION,   ALLOCATABLE  ::  MASSES(:)                    !
      INTEGER, PARAMETER :: NSTORED=109                                 !
!-----------------------------------------------------------------------!
      CONTAINS                                                          !
!     ==================================================================!
      SUBROUTINE SETPERIODICTABLE                                       !
      IMPLICIT NONE                                                     !
                                                                        !
      ALLOCATE  ( ELEMENTS(NSTORED) ,   &                               !
                   MASSES (NSTORED)        )                            !
                                                                        !
      CALL SETELEMENT                                                   !
      CALL SETMASSES                                                    !
                                                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE SETPERIODICTABLE                                   !
!     ===================================================================
      SUBROUTINE SETELEMENT 
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=2), DIMENSION(NSTORED) :: TMP                      !
      DATA TMP     /' H','He',                                         &! 
      'Li','Be',' B',' C',' N',' O',' F','Ne',                         &!
      'Na','Mg','Al','Si',' P',' S','Cl','Ar',                         &!
      ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu',          &!
      'Zn','Ga','Ge','As','Se','Br','Kr',                              &!
      'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',          &!
      'Cd','In','Sn','Sb','Te',' I','Xe',                              &!
      'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',     &!
      'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir','Pt',     &!
      'Au','Hg','Tl','Pb','Bi','Po','At','Rn',                         &!
      'Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es',&!
      'Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt'/                !
                                                                        !
      ELEMENTS=TMP                                                      !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE SETELEMENT                                         !
!     ===================================================================
      SUBROUTINE SETMASSES 
!     ===================================================================
      IMPLICIT NONE                                                     !
      DOUBLE PRECISION, DIMENSION(NSTORED) :: TMP                       !
      DATA TMP    / 1.00790D0,  4.00260D0,  6.94000D0,  9.01218D0,     &!
       10.81000D0, 12.01100D0, 14.00670D0, 15.99940D0, 18.99840D0,     &!
       20.17900D0, 22.98977D0, 24.30500D0, 26.98154D0, 28.08550D0,     &!
       30.97376D0, 32.06000D0, 35.45300D0, 39.94800D0, 39.09830D0,     &!
       40.08000D0, 44.95590D0, 47.90000D0, 50.94150D0, 51.99600D0,     &!
       54.93800D0, 55.84700D0, 58.93320D0, 58.71000D0, 63.54600D0,     &!
       65.38000D0, 69.73500D0, 72.59000D0, 74.92160D0, 78.96000D0,     &!
       79.90400D0, 83.80000D0, 85.46780D0, 87.62000D0, 88.90590D0,     &!
       91.22000D0, 92.90640D0, 95.94000D0, 98.90620D0, 101.0700D0,     &!
       102.9055D0, 106.4000D0, 107.8680D0, 112.4100D0, 114.8200D0,     &!
       118.6900D0, 121.7500D0, 127.6000D0, 126.9045D0, 131.3000D0,     &!
       132.9054D0, 137.3300D0, 15*0.000D0, 178.4900D0, 180.9479D0,     &!
       183.8500D0, 186.2070D0, 190.2000D0, 192.2200D0, 195.0900D0,     &!
       196.9665D0, 200.5900D0, 204.3700D0, 207.2000D0, 208.9804D0,     &! 
       26*0.000D0/                                                      !         
       MASSES=TMP                                                       !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE SETMASSES                                          !
!     ===================================================================
END MODULE PERIODICTABLE      
      
