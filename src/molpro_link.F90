MODULE MOLPRO_LINK
      IMPLICIT NONE
      LOGICAL :: COMP_CHARGE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CHARGE
!     -------------------------------------------------------------------
      CONTAINS
!     ===================================================================
      SUBROUTINE GET_VECTOR_MOLPRO ( ICHANNEL, G,  NUMAT,  LABEL ) 
!     ===================================================================
      USE  SYSTEM_UTIL, ONLY:  SEARCH,           &
                               CPSTOP,           &
                               JLINE
      USE  XYZ,         ONLY:  POINT
!     -------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN)                    :: ICHANNEL
      INTEGER (KIND=8), INTENT(IN)           :: NUMAT
      TYPE(POINT), INTENT(INOUT)             :: G(*)
      CHARACTER  (LEN=20),  INTENT(IN)       :: LABEL
      CHARACTER  (LEN=256)                   :: LINE
      LOGICAL                                :: CONFERMA
      INTEGER                                :: I, ATDUM
      
      CONFERMA=.FALSE.
      CALL SEARCH(LABEL,ICHANNEL,.TRUE.,'AHEA',CONFERMA,LINE)
      IF (.NOT.CONFERMA) THEN
         WRITE(6,100)'CANNOT FIND LABEL :',LABEL,' ON CHANNEL FILE:', &
                     ICHANNEL
         CALL CPSTOP('GET_VECTOR_MOLPRO')
      END IF
      CONFERMA=.FALSE.
      CALL JLINE(ICHANNEL,3,CONFERMA)
      IF (.NOT.CONFERMA) CALL CPSTOP('JLINE by GET_GRAD_DIFF_MOLPRO')
      DO I=1,NUMAT
         READ(ICHANNEL,*,ERR=99)ATDUM,G(I)%X,G(I)%Y,G(I)%Z
      END DO
      RETURN
!     -------------------------------------------------------------------
 99   WRITE(6,101)'AN ERROR OCCURED WHILE READING FILE ON CHANNEL :',&
                  ICHANNEL,' FOR :',LABEL
      CALL CPSTOP('GET_VECTOR_MOLPRO')
      RETURN
!     -------------------------------------------------------------------
 100  FORMAT (A,/A/,A,I5)
 101  FORMAT (A,I5,/2A)
!     -------------------------------------------------------------------
      END SUBROUTINE GET_VECTOR_MOLPRO
!     ===================================================================
      SUBROUTINE GET_ENERGIES_MOLPRO ( ICHANNEL, VEXC, VFOND )
!     ===================================================================
      USE  SYSTEM_UTIL, ONLY:  SEARCHT,           &
                               CPSTOP,            &
                               SEARCH
      USE  XYZ,         ONLY:  POINT
      USE  PARSING,     ONLY:  CLEANSTRING 
!     -------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER  (LEN=10)                    :: LABEL1,LABEL2,LABEL3
      CHARACTER  (LEN=256)                   :: LINE
      LOGICAL                                :: CONFERMA
      INTEGER                                :: IND1
      INTEGER,  INTENT(IN)                   :: ICHANNEL
      DOUBLE PRECISION,  INTENT(OUT)         :: VEXC, VFOND

      LABEL1='MC'
      LABEL2='1.1'
      LABEL3='ENERGY'
      CALL CLEANSTRING(LINE,LEN(LINE))
      CONFERMA=.FALSE.
      CALL SEARCHT(LABEL1,LABEL2,LABEL3,ICHANNEL,.TRUE.,'AHEA', & 
                   CONFERMA,LINE)
      IF (.NOT.CONFERMA) THEN
         WRITE(6,100)LABEL1,LABEL2,LABEL3,ICHANNEL
         CALL CPSTOP('GET_ENERGIES_MOLPRO')
      END IF      
      IND1=INDEX(LINE,'ENERGY')+10
      READ(LINE(IND1:),*)VFOND

      LABEL1='MC'
      LABEL2='2.1'
      LABEL3='ENERGY'
      CALL CLEANSTRING(LINE,LEN(LINE))
      CONFERMA=.FALSE.
      CALL SEARCHT(LABEL1,LABEL2,LABEL3,ICHANNEL,.TRUE.,'AHEA', & 
                   CONFERMA,LINE)
      IF (.NOT.CONFERMA) THEN
         WRITE(6,100)LABEL1,LABEL2,LABEL3,ICHANNEL
         CALL CPSTOP('GET_ENERGIES_MOLPRO')
      END IF      
      IND1=INDEX(LINE,'ENERGY')+10
      READ(LINE(IND1:),*)VEXC

      RETURN
!     -------------------------------------------------------------------
 100  FORMAT ('CANNOT FIND LABELS :',3(/A/),' ON CHANNEL FILE:',I5)
!     -------------------------------------------------------------------
      END SUBROUTINE GET_ENERGIES_MOLPRO
!     ===================================================================
      SUBROUTINE GET_ENERGY_MOLPRO ( ICHANNEL, VEXC, VFOND )
!     ===================================================================
      USE  SYSTEM_UTIL, ONLY:  SEARCHT,           &
                               CPSTOP,            &
                               SEARCH
      USE  XYZ,         ONLY:  POINT
      USE  PARSING,     ONLY:  CLEANSTRING 
!     -------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER  (LEN=10)                    :: LABEL1,LABEL2,LABEL3
      CHARACTER  (LEN=256)                   :: LINE
      LOGICAL                                :: CONFERMA
      INTEGER                                :: IND1
      INTEGER,  INTENT(IN)                   :: ICHANNEL
      DOUBLE PRECISION,  INTENT(OUT)         :: VEXC, VFOND

      LABEL1='RHF'
      LABEL2='1.1'
      LABEL3='ENERGY'
      CALL CLEANSTRING(LINE,LEN(LINE))
      CONFERMA=.FALSE.
      CALL SEARCHT(LABEL1,LABEL2,LABEL3,ICHANNEL,.TRUE.,'AHEA', & 
                   CONFERMA,LINE)
      IF (.NOT.CONFERMA) THEN
         WRITE(6,100)LABEL1,LABEL2,LABEL3,ICHANNEL
         CALL CPSTOP('GET_ENERGIES_MOLPRO')
      END IF      
      IND1=INDEX(LINE,'ENERGY')+10
      READ(LINE(IND1:),*)VFOND

      VEXC = 0.D0
      RETURN
!     -------------------------------------------------------------------
 100  FORMAT ('CANNOT FIND LABELS :',3(/A/),' ON CHANNEL FILE:',I5)
!     -------------------------------------------------------------------
       END SUBROUTINE GET_ENERGY_MOLPRO
!     ===================================================================
      SUBROUTINE GET_CHARGES ( ICHANNEL, NUMAT )
!     ===================================================================
      USE  SYSTEM_UTIL, ONLY:  SEARCHT,           &
                               CPSTOP,            &
                               SEARCH
      USE  XYZ,         ONLY:  POINT
      USE  PARSING,     ONLY:  CLEANSTRING 
!     -------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER  (LEN=3)                     :: SIGN
      CHARACTER  (LEN=10)                    :: LABEL1,LABEL2,LABEL3
      CHARACTER  (LEN=256)                   :: LINE
      LOGICAL                                :: CONFERMA
      INTEGER                                :: IND1, I
      INTEGER,  INTENT(IN)                   :: ICHANNEL
      INTEGER (KIND=8), INTENT(IN)           :: NUMAT

      LABEL1='Unique'
      LABEL2='Total'
      LABEL3='Charge'
      COMP_CHARGE = .FALSE.
      CALL CLEANSTRING(LINE,LEN(LINE))
      CONFERMA=.FALSE.
      CALL SEARCHT(LABEL1,LABEL2,LABEL3,ICHANNEL,.TRUE.,'AHEA', & 
                   CONFERMA,LINE)
      IF (.NOT.CONFERMA) THEN
         WRITE(6,'(A)')'NO CHARGE AVAILABLE'
         RETURN
      END IF 
      COMP_CHARGE = .TRUE.
      IF (.NOT.ALLOCATED(CHARGE)) ALLOCATE(CHARGE(NUMAT))
      CHARGE = 0.D0
      IND1=INDEX(LINE,'Total')+5
      DO I = 1, NUMAT
         READ(ICHANNEL,'(A)')LINE
         READ(LINE(IND1:IND1+3),'(A3)')SIGN
         READ(LINE(IND1+3:),*)CHARGE(I)
         IF (INDEX(SIGN,'-').NE.0) CHARGE(I) = CHARGE(I) * (-1.D0)
      END DO

      REWIND(ICHANNEL)
      RETURN
!     -------------------------------------------------------------------
      END SUBROUTINE GET_CHARGES
!     ===================================================================
END MODULE MOLPRO_LINK
