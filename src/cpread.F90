!------------------------------------------------------------------------
MODULE UNITS_FILE                                                       !
      IMPLICIT NONE                                                     !
!     ------------------------------------------------------------------!
      INTEGER :: XYZFILE, GHPLANE, RESFILE, PRPFILE
      CONTAINS                                                          !
!     ===================================================================
      SUBROUTINE SETFILES(INPUTFILE, RESTART )
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=256), INTENT(IN) :: INPUTFILE                      !      
      LOGICAL, INTENT(IN) :: RESTART                                    ! 
                                                                        !
      XYZFILE = 20                                                      !
      GHPLANE = 21                                                      !
      RESFILE = 22                                                      !
      PRPFILE = 23
      IF (RESTART) THEN                                                 !
         CALL OPENFILES(INPUTFILE,'xyz',XYZFILE,'UNKNOWN','FORMATTED'  &! 
              ,'APPEND')                                                !      
         CALL OPENFILES(INPUTFILE,'ghp',GHPLANE,'UNKNOWN','FORMATTED'  &! 
              ,'APPEND')                                                !
         CALL OPENFILES(INPUTFILE,'prp',PRPFILE,'UNKNOWN','FORMATTED'  &! 
              ,'APPEND')                                                !
      ELSE                                                              !
         CALL OPENFILES(INPUTFILE,'xyz',XYZFILE,'UNKNOWN','FORMATTED'  &!
              ,'REWIND')                                                !
         CALL OPENFILES(INPUTFILE,'ghp',GHPLANE,'UNKNOWN','FORMATTED'  &!
              ,'REWIND')                                                !
         CALL OPENFILES(INPUTFILE,'prp',PRPFILE,'UNKNOWN','FORMATTED'  &!
              ,'REWIND')                                                !
      END IF                                                            !
      CALL OPENFILES(INPUTFILE,'res',RESFILE,'UNKNOWN','UNFORMATTED'&   !
                              ,'REWIND')                                !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE SETFILES                                           !
!     ===================================================================
      SUBROUTINE OPENFILES(INPUTFILE,TARG,ICHANNEL,STATUS,FORM,ACCESS)         
!     ===================================================================
      USE SYSTEM_UTIL, ONLY : CP_OPEN_SYS                               !
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=256), INTENT(IN) :: INPUTFILE                      ! 
      CHARACTER (LEN=3), INTENT(IN) :: TARG                             ! 
      CHARACTER (LEN=*), INTENT(IN) :: STATUS, FORM, ACCESS             ! 
      CHARACTER (LEN=256) :: FILE                                       ! 
      INTEGER,  INTENT(IN) :: ICHANNEL                                  ! 
      INTEGER :: I                                                      ! 
                                                                        !
      I=INDEX(INPUTFILE,'.')                                            !
      IF (I.EQ.0) I=INDEX(INPUTFILE,' ')                                !
      FILE=INPUTFILE(1:I-1)//'.'//TARG                                  !
                                                                        !
      CALL CP_OPEN_SYS(ICHANNEL,FILE,STATUS,FORM,ACCESS)                !
                                                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE OPENFILES                                          !
!     ===================================================================
      SUBROUTINE FLUSH_FILES
!     ===================================================================
      IMPLICIT NONE                                                     !
                                                                        !
      CALL FLUSH(6)                                                     !
      CALL FLUSH(XYZFILE)                                               !
      CALL FLUSH(GHPLANE)                                               !
      CALL FLUSH(RESFILE)                                               !
      CALL FLUSH(PRPFILE)                                               !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE FLUSH_FILES                                        !
!     ===================================================================    
END MODULE


!------------------------------------------------------------------------
MODULE START_JOB                                                        !
      IMPLICIT NONE                                                     !
!     -------------------------------------------------------------------
      TYPE ZMATYPE                                                      !
        SEQUENCE                                                        ! 
        CHARACTER (LEN=2)  :: LABEL                                     ! 
        INTEGER   (KIND=2) :: AT1                                       ! 
        INTEGER   (KIND=2) :: AT2                                       ! 
        INTEGER   (KIND=2) :: AT3                                       ! 
        REAL      (KIND=8) :: BOND                                      ! 
        REAL      (KIND=8) :: ANGLE                                     ! 
        REAL      (KIND=8) :: TORSION                                   ! 
      END TYPE ZMATYPE                                                  ! 
!     ------------------------ MAIN VAR --------------------------------- 
      CHARACTER (LEN=256) :: INPUTFILE                                  ! 
      INTEGER   (KIND=8)  :: NUMAT  , IRESVEL                           ! 
      CHARACTER (LEN=256) :: GEOFILE                                    ! 
      CHARACTER (LEN=10)  :: ABINIPROG                                  ! 
      CHARACTER (LEN=256) :: ABINIT_INP,ABINIEXE,ABINIT_OUT             ! 
      REAL      (KIND=8)  :: TSTEP,ENECONV                              ! 
      CHARACTER (LEN=20)  :: METHOD                                     ! 
      CHARACTER (LEN=10)  :: DYNTYPE                                    ! 
      CHARACTER (LEN=9 )  :: COORDINATES                                !
      CHARACTER (LEN=9 )  :: TYPEZMAT                                   !
      INTEGER   (KIND=8)  :: NSTEP,IPRINT                               ! 
      INTEGER             :: CLX,CLY,CLZ,CPX,CPY,CPZ,CVG,CVH,KK,NCONSTR !
      CHARACTER (LEN=256) :: TITLE                                      ! 
      LOGICAL             :: USEGRADEX,USESQUARENE                      ! 
      LOGICAL             :: AREWERESTARTING                            ! 
      LOGICAL             :: PROP_AND_STOP                              ! 
      LOGICAL             :: NOROT,NOELLE,NOPI,NOGH                     !
      LOGICAL             :: KUTTEH                                     ! 
      LOGICAL             :: V1V2UPDATE                                 !
      INTEGER             :: PROPERTIES                                 ! 
      INTEGER             :: IMOM, IBOND, IANGLE, ITORSION, IPVEC       ! 
      CHARACTER (LEN=10), ALLOCATABLE :: PROPERTY(:)                    ! 
      INTEGER, ALLOCATABLE :: BOND(:,:), ANGLE(:,:), TORSION(:,:)       ! 
      INTEGER             :: UPV12                                      !
      DOUBLE PRECISION    :: INI_TEMP                                   !
      DOUBLE PRECISION    :: ALPHA,SIGMA                                !
      DOUBLE PRECISION    :: TEMP_SCALE                                 !
      DOUBLE PRECISION    :: DELTAT                                     !
      DOUBLE PRECISION    :: RANVEC                                     !
      DOUBLE PRECISION    :: TEMP_MAX                                   !
      DOUBLE PRECISION    :: LFACT                                      !
      DOUBLE PRECISION    :: DAMP                                       !
      INTEGER             :: NREP1, NREP2                               !
      LOGICAL             :: ADDREPULSIVE                               !
      DOUBLE PRECISION    :: AREP, SIGREP                               !
!     -------------------------------------------------------------------
      CONTAINS
!     ===================================================================
      SUBROUTINE CONPATH_START
!     ===================================================================
      USE SYSTEM_UTIL, ONLY: CP_OPEN,      &                            !  
                             CP_CLOSE                                   !
      USE UNITS_FILE,  ONLY: SETFILES                                   !
      IMPLICIT NONE                                                     !
      INTEGER :: ICHANNEL                                               ! 
      INTEGER, EXTERNAL :: IARGC
!     ------------------------------------------------------------------!
!     Set Default values                                                !
!     ------------------------------------------------------------------!
      NUMAT=0                                                           !
      GEOFILE='geo.inp'                                                 !
      TSTEP=0.1                 ! femtoseconds                          !
      ENECONV=0.0001            ! A.U.                                  !
      METHOD='DYNAMIC'                                                  !
      ABINIT_INP='molpro.com'                                           !
      ABINIT_OUT='molpro.out'                                           !
      COORDINATES='CARTESIAN'                                           !
      NSTEP=100                                                         !
      IPRINT=1                                                          !
      TITLE='NO TITLE'                                                  !
      ABINIEXE='/apps/molpro2000.1/bin/molpro'                          !
      ABINIPROG='MOLPRO'                                                !
      DAMP=1.D0                                                         !
      DYNTYPE='CONICAL   '                                              !
      USEGRADEX=.FALSE.                                                 !
      USESQUARENE=.FALSE.                                               !
      AREWERESTARTING=.FALSE.                                           !
      PROP_AND_STOP=.FALSE.                                             !
      LFACT=0.D0                                                        !
      TEMP_SCALE=1.D0                                                   !
      DELTAT=0.D0                                                       !
      RANVEC=0.D0                                                       !
      TEMP_MAX=5000.D0                                                  !
      UPV12=1                                                           !
      ALPHA=0.D0                                                        !
      IRESVEL=0                                                         !
      PROPERTIES=0                                                      !
      SIGMA=1.D0                                                        !
      IMOM = 0                                                          !
      IPVEC = 0                                                         !
      IBOND = 0                                                         !
      IANGLE = 0                                                        !
      ITORSION = 0                                                      !
      NOROT=.FALSE.                                                     !
      NOPI=.TRUE.                                                       !
      NOELLE=.TRUE.                                                     !
      NOGH=.TRUE.                                                       !
      KUTTEH=.FALSE.                                                    !
      V1V2UPDATE=.FALSE.                                                !
      ADDREPULSIVE=.FALSE.                                              !
      TYPEZMAT="MOLPRO"                                                 !
!     ------------------------------------------------------------------!
!     Get Input file                                                    !
!     ------------------------------------------------------------------!
#ifdef HPUX                                          
      CALL GETARG(2,INPUTFILE)                                          !
      IF (INDEX(INPUTFILE,'restart').NE.0) AREWERESTARTING=.TRUE.       !
      CALL GETARG(1,INPUTFILE)                                          ! 
#elif IFC                                            
      IF (IARGC().GT.1) THEN                                            !
         CALL GETARG(2,INPUTFILE)                                       !
         IF (INDEX(INPUTFILE,'restart').NE.0) AREWERESTARTING=.TRUE.    !
      ENDIF                                                             !
      CALL GETARG(1,INPUTFILE)                                          !
#else                                                  
      CALL GETARG(2,INPUTFILE)                                          !
      IF (INDEX(INPUTFILE,'restart').NE.0) AREWERESTARTING=.TRUE.       !
      CALL GETARG(1,INPUTFILE)                                          !
#endif                                                 
      ICHANNEL=10                                                       !
!     ------------------------------------------------------------------!
!     Open Input file and read all the stuff you need                   !
!     ------------------------------------------------------------------!
      CALL   CP_OPEN( ICHANNEL, INPUTFILE, 'OLD', 'FORMATTED' )         !
      CALL READ_DATA( ICHANNEL, 'MAIN' )                                !
      CALL READ_DATA( ICHANNEL, 'ABINITIO')                             !
      CALL  CP_CLOSE( ICHANNEL, INPUTFILE )                             !
      CALL SETFILES ( INPUTFILE, AREWERESTARTING )                      !
!     ------------------------------------------------------------------!
      END SUBROUTINE CONPATH_START                                      !
!     ===================================================================
      SUBROUTINE READ_DATA ( ICHANNEL, TARGET )
!     ===================================================================
      USE SYSTEM_UTIL, ONLY: SEARCH,       &                            !
                             GETLABEL,     &                            !
                             CPSTOP,       &                            !
                             CP_OPEN,      &                            !
                             CP_CLOSE,     &                            !
                             LFBL                                       !
      USE PARSING,     ONLY: GETFIELD                                   !
      USE EMPTYLINE                                                     !
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=*), INTENT(IN) :: TARGET                           ! 
      CHARACTER (LEN=256) :: LINE                                       ! 
      INTEGER, INTENT(IN) :: ICHANNEL                                   ! 
      INTEGER :: INDU, ICHANNEL2, IND1, IND2, IPROP, J, PRT, ICL        ! 
      LOGICAL :: CONFERMA                                               ! 
      CHARACTER (LEN=10) :: LABEL                                       ! 
      CHARACTER (LEN=20) :: SLINE, SCRLIN                               ! 
      INTEGER, ALLOCATABLE :: BONDS(:,:), ANGLES(:,:), TORSIONS(:,:)    ! 
                                                                        !
      IPROP = 0                                                         !
      CONFERMA = .FALSE.                                                !
      SLINE='                    '                                      !
      SLINE='&'//TARGET                                                 !
      CALL SEARCH(SLINE,ICHANNEL,.TRUE.,'AHEA',CONFERMA,LINE)           !
      IF (.NOT.CONFERMA) GO TO 9                                        !
                                                                        !
      READ(ICHANNEL,'(A)')LINE                                          !
      SELECT CASE (TARGET)                                              !
      CASE ('MAIN')                                                     !
         DO WHILE (INDEX(LINE,'&END').EQ.0)                             !
            IF (.NOT.EMPTY(LINE,256)) THEN                              !
               LABEL=GETLABEL(LINE,'=',256)                             !
               SELECT CASE (LABEL)                                      !
               CASE ('PROPERTIES')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)PROPERTIES                     !
               CASE ('PROPERTY  ')                                      !
                  IF (PROPERTIES.EQ.0) THEN                             !
                     WRITE(6,'(A)')  &                                  !
                     ' NUMBER OF PROPERTIES REQUESTED IS ZERO!'         !
                     WRITE(6,'(A)')'CHECK YOUR INPUT AND RERUN THE JOB!'!
                     CALL CPSTOP('READ_DATA')                           !
                  ELSE                                                  !
                     IPROP = IPROP + 1                                  !
                     IF (IPROP.EQ.1) ALLOCATE ( PROPERTY(PROPERTIES) )  !
                  END IF                                                !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,'(')-1                                !
                  IF (IND2.EQ.-1) IND2=INDEX(LINE,',')-1                !
                  READ(LINE(IND1:IND2),'(A)')PROPERTY(IPROP)            !
                  DO ICL=IND2-IND1+2,10
                     PROPERTY(IPROP)(ICL:ICL)=' '
                  END DO 
                  SELECT CASE(PROPERTY(IPROP))                          !
                  CASE ('BOND      ')                                   !
                     IBOND = IBOND + 1                                  !
                     IF (IBOND.EQ.1) THEN                               !
                        ALLOCATE ( BOND(1,2) )                          !
                     ELSE                                               !
                        ALLOCATE ( BONDS(IBOND-1,2) )                   !
                        BONDS=BOND                                      !
                        DEALLOCATE ( BOND )                             !
                        ALLOCATE ( BOND(IBOND,2) )                      !
                        DO J=1,IBOND-1                                  !
                           BOND(J,1)=BONDS(J,1)                         !
                           BOND(J,2)=BONDS(J,2)                         !
                        END DO                                          !
                        DEALLOCATE (BONDS)                              !
                     END IF                                             !
                     IND1=IND2+2                                        !
                     IND2=INDEX(LINE,')')-1                             !
                     CALL GETFIELD(SCRLIN,LEN(SCRLIN),',',',', &        !
                          ','//LINE(IND1:IND2)//',',1)                  !
                     READ(SCRLIN,*)BOND(IBOND,1)                        !
                     CALL GETFIELD(SCRLIN,LEN(SCRLIN),',',',', &        !
                          ','//LINE(IND1:IND2)//',',2)                  !
                     READ(SCRLIN,*)BOND(IBOND,2)                        !
                                                                        !
                  CASE ('ANGLE     ')                                   !
                     IANGLE = IANGLE + 1                                !
                     IF (IANGLE.EQ.1) THEN                              !
                        ALLOCATE ( ANGLE(1,3) )                         !
                     ELSE                                               !
                        ALLOCATE ( ANGLES(IANGLE-1,3) )                 !
                        ANGLES=ANGLE                                    !
                        DEALLOCATE ( ANGLE )                            !
                        ALLOCATE ( ANGLE(IANGLE,3) )                    !
                        DO J=1,IANGLE-1                                 !
                           ANGLE(J,1)=ANGLES(J,1)                       !
                           ANGLE(J,2)=ANGLES(J,2)                       !
                           ANGLE(J,3)=ANGLES(J,3)                       !
                        END DO                                          !
                        DEALLOCATE (ANGLES)                             !
                     END IF                                             !
                     IND1=IND2+2                                        !
                     IND2=INDEX(LINE,')')-1                             !
                     CALL GETFIELD(SCRLIN,LEN(SCRLIN),',',',', &        !
                          ','//LINE(IND1:IND2)//',',1)                  !
                     READ(SCRLIN,*)ANGLE(IANGLE,1)                      !
                     CALL GETFIELD(SCRLIN,LEN(SCRLIN),',',',', &        !
                          ','//LINE(IND1:IND2)//',',2)                  !
                     READ(SCRLIN,*)ANGLE(IANGLE,2)                      !
                     CALL GETFIELD(SCRLIN,LEN(SCRLIN),',',',', &        !
                          ','//LINE(IND1:IND2)//',',3)                  !
                     READ(SCRLIN,*)ANGLE(IANGLE,3)                      !
                                                                        !
                  CASE ('TORSION   ')                                   !
                     ITORSION = ITORSION + 1                            !
                     IF (ITORSION.EQ.1) THEN                            !
                        ALLOCATE ( TORSION(1,4) )                       !
                     ELSE                                               !
                        ALLOCATE ( TORSIONS(ITORSION-1,4) )             !
                        TORSIONS = TORSION                              !
                        DEALLOCATE ( TORSION )                          !
                        ALLOCATE ( TORSION(ITORSION,4) )                !
                        DO J=1,ITORSION-1                               !
                           TORSION(J,1)=TORSIONS(J,1)                   !
                           TORSION(J,2)=TORSIONS(J,2)                   !
                           TORSION(J,3)=TORSIONS(J,3)                   !
                           TORSION(J,4)=TORSIONS(J,4)                   !
                        END DO                                          !
                        DEALLOCATE (TORSIONS)                           !
                     END IF                                             !
                     IND1=IND2+2                                        !
                     IND2=INDEX(LINE,')')-1                             !
                     CALL GETFIELD(SCRLIN,LEN(SCRLIN),',',',', &        !
                          ','//LINE(IND1:IND2)//',',1)                  !
                     READ(SCRLIN,*)TORSION(ITORSION,1)                  !
                     CALL GETFIELD(SCRLIN,LEN(SCRLIN),',',',', &        !
                          ','//LINE(IND1:IND2)//',',2)                  !
                     READ(SCRLIN,*)TORSION(ITORSION,2)                  !
                     CALL GETFIELD(SCRLIN,LEN(SCRLIN),',',',', &        !
                          ','//LINE(IND1:IND2)//',',3)                  !
                     READ(SCRLIN,*)TORSION(ITORSION,3)                  !
                     CALL GETFIELD(SCRLIN,LEN(SCRLIN),',',',', &        !
                          ','//LINE(IND1:IND2)//',',4)                  !
                     READ(SCRLIN,*)TORSION(ITORSION,4)                  !
                                                                        !
                  CASE ('ANGULARMOM')                                   !
                                                                        !
                     IMOM = IMOM + 1                                    !
                                                                        !
                  CASE ('MOMENTUM  ')                                   !
                                                                        !
                     IPVEC = IPVEC + 1                                  !
                                                                        !
                  CASE DEFAULT                                          !
                     WRITE(6,'(2A)')'PROPERTY NOT YET IMPLEMENTED!', &  !
                          PROPERTY(IPROP)                               !
                     CALL CPSTOP('READ_DATA')                           !
                  END SELECT                                            !
                                                                        !
               CASE ('NOROTATION')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)NOROT                          !
               CASE ('LCONSTR   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)NOELLE
               CASE ('PCONSTR   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)NOPI
               CASE ('GHCONSTR  ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)NOGH
               CASE ('KUTTEH    ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)KUTTEH                         !
               CASE ('V1V2UPDATE    ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)V1V2UPDATE
               CASE ('VELREST   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)IRESVEL                        !
               CASE ('JUSTPROP  ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)PROP_AND_STOP                  !
               CASE ('ADDREPULS ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)ADDREPULSIVE                   !
               CASE ('NREP1     ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)NREP1
               CASE ('NREP2     ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)NREP2
               CASE ('AREP      ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)AREP
               CASE ('SIGREP    ')                                    !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)SIGREP
               CASE ('ALPHA     ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)ALPHA                          ! 
               CASE ('SIGMA     ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)SIGMA                          ! 
               CASE ('INITEMP   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)INI_TEMP                       !
               CASE ('SCALETEMP   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)TEMP_SCALE
               CASE ('DELTATEMP ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)DELTAT
               CASE ('RANDVEC   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)RANVEC
               CASE ('MAXTEMP   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)TEMP_MAX
               CASE ('STEPV1V2 ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)UPV12
               CASE ('LFACT     ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)LFACT                          !
               CASE ('USESQUAREN')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)USESQUARENE                    !
               CASE ('USEGRADEX ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)USEGRADEX                      !
               CASE ('DYNTYPE   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')DYNTYPE                    !
               CASE ('DAMP      ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)DAMP                           !
               CASE ('ABINIPROG ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')ABINIPROG                  !
               CASE ('ABINIEXE  ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')ABINIEXE                   !
               CASE ('NSTEP     ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)NSTEP                          !
               CASE ('IPRINT    ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)IPRINT                         !
               CASE ('COORDINATE')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')COORDINATES                !
               CASE ('TYPEZMAT')                                        !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')TYPEZMAT                   !
               CASE ('TITLE     ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')TITLE                      !
               CASE ('NUMAT     ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)NUMAT                          !
               CASE ('GEOFILE   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')GEOFILE                    !
               CASE ('TSTEP     ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)TSTEP                          !
               CASE ('ENECONV   ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),*)ENECONV                        !
               CASE ('METHOD    ')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')METHOD                     !
               CASE ('ABINIT_INP')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')ABINIT_INP                 !
               CASE ('ABINIT_OUT')                                      !
                  IND1=INDEX(LINE,'=')+1                                !
                  IND2=INDEX(LINE,',')-1                                !
                  READ(LINE(IND1:IND2),'(A)')ABINIT_OUT                 !
               CASE DEFAULT                                             !
                  WRITE(6,'(A)')'ERROR IN ConPath INPUT FILE.'          !
                  WRITE(6,'(A)') &                                      !
                  'UNKNOWN KEYWORD OR SYNTAX ERROR ON LINE :'           !
                  WRITE(6,'(A)') LINE                                   !
                  CALL CPSTOP('READ_DATA')                              !
               END SELECT                                               !
            END IF                                                      !
            READ(ICHANNEL,'(A)')LINE                                    !
         END DO                                                         !
      CASE ('ABINITIO')                                                 !
         ICHANNEL2=11                                                   !
         CALL CP_OPEN(ICHANNEL2,ABINIT_INP,'UNKNOWN','FORMATTED')       !
         DO WHILE (INDEX(LINE,'&END').EQ.0)                             !
            WRITE(ICHANNEL2,'(A)')LINE(1:LFBL(LINE,LEN(LINE)))          !
            READ(ICHANNEL,'(A)')LINE                                    !
         END DO                                                         !
         CALL CP_CLOSE(ICHANNEL2,ABINIT_INP)                            !
      CASE DEFAULT                                                      !
         WRITE(6,'(A)')'AVAILABLE ENTRIES : MAIN, ABINITIO !'           !
         CALL CPSTOP('READ_DATA')                                       !
      END SELECT                                                        !
      IF (KUTTEH) THEN                                                  !
         KK=1                                                           !
         IF (NOGH)THEN                                                  !
           WRITE (6,*)'Constraining G'                                  !
           CVG=KK                                                       !
           KK=KK+1                                                      !
         ENDIF                                                          !
         IF (NOGH) THEN                                                 !
           WRITE (6,*)'Constraining H'                                  !
           CVH=KK                                                       !
           KK=KK+1                                                      !
         ENDIF                                                          !
         IF (NOELLE) THEN                                               !
           WRITE (6,*)'Constraining Lx'                                 !
           CLX=KK                                                       !
           KK=KK+1                                                      !
         ENDIF                                                          !
         IF (NOELLE) THEN                                               !
           WRITE (6,*)'Constraining Ly'                                 !
           CLY=KK                                                       !
           KK=KK+1                                                      !
         ENDIF                                                          !
         IF (NOELLE) THEN                                               !
           WRITE (6,*)'Constraining Lz'                                 !
           CLZ=KK                                                       !
           KK=KK+1                                                      !
         ENDIF                                                          !
         IF (NOPI) THEN                                                 !
           WRITE (6,*)'Constraining Px'                                 !
           CPX=KK                                                       !
           KK=KK+1                                                      !
         ENDIF                                                          !
         IF (NOPI) THEN                                                 !
           WRITE (6,*)'Constraining Py'                                 !
           CPY=KK                                                       !
           KK=KK+1                                                      !
         ENDIF                                                          !
         IF (NOPI) THEN                                                 !
           WRITE (6,*)'Constraining Pz'                                 !
           CPZ=KK                                                       !
           KK=KK+1                                                      !
         ENDIF                                                          !
         KK=KK-1                                                        !
         NCONSTR=KK                                                     !
      ENDIF                                                             !
                                                                        !
      PRT = IMOM+IBOND+IANGLE+ITORSION+IPVEC                            !
      IF (KUTTEH) NOROT=.FALSE.                                         !
      IF (PRT.NE.PROPERTIES) THEN                                       !
         WRITE(6,'(A,I5,2A,/,A,I5,A)') &                                !
              'YOU HAVE REQUESTED ',PROPERTIES,' PROPERTIES ',&         !
              'CALCULATIONS',' BUT YOU HAVE SPECIFIED ONLY ',PRT, &     !
              ' PROPERTIES.'                                            !
         CALL CPSTOP('READ_DATA')                                       !
      END IF                                                            !
                                                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 9    WRITE(6,'(4A)')'TARGET : ',TARGET,' NOT FOUND IN INPUT FILE :', & !
                    INPUTFILE                                           !
      CALL CPSTOP('READ_DATA')                                          !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE READ_DATA                                          !
!     ===================================================================
      SUBROUTINE CONPATH_LOGO
!     ===================================================================
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=10)  :: VERSION                                    ! 
      CHARACTER (LEN=20)  :: AUTHOR1, AUTHOR2                           ! 
      CHARACTER (LEN=50)  :: INSTITUTION1, INSTITUTION2, INSTITUTION3   ! 
      CHARACTER (LEN=20)  :: LOGIN                                      ! 
                                                                        !
      VERSION='alpha 0.00'                                              !
      AUTHOR1='Teodoro Laino'                                           !
      AUTHOR2='Daniele Passerone'                                       !
      INSTITUTION1='Scuola Normale Superiore di Pisa'                   !
      INSTITUTION2='NEST - INFM'                                        !
      INSTITUTION3='CSCS (Manno Switzerland)'                           !
      WRITE(6,999)VERSION,AUTHOR1,INSTITUTION1,INSTITUTION2,  &         !
                  AUTHOR2,INSTITUTION3                                  !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 999  FORMAT(9X,'Program to perform search of conical Intersections',&  !
     //,9X,' #####                  ######                         ',&  !   
      /,9X,'#     #   ####   #    # #     #    ##     #####  #    #',&  !   
      /,9X,'#        #    #  ##   # #     #   #  #      #    #    #',&  !   
      /,9X,'#        #    #  # #  # ######   #    #     #    ######',&  !   
      /,9X,'#        #    #  #  # # #        ######     #    #    #',&  !   
      /,9X,'#     #  #    #  #   ## #        #    #     #    #    #',&  !   
      /,9X,' #####    ####   #    # #        #    #     #    #    #',&  !
     //,' Version :',A/,' Authors :',A,10X,A/,61X,A/,'          ',A ,&  !
        18X,A/ )                                                        !
!     ------------------------------------------------------------------!
      END SUBROUTINE CONPATH_LOGO                                       !
!     ===================================================================

END MODULE START_JOB
