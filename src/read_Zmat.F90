MODULE READ_ZMAT
  USE PERIODICTABLE, ONLY: AMS   => MASSES,         &
                           ELEMENTS,                &
                           NSTORED
  USE DYNPREPARE,    ONLY: NUM_BONDS,               &
                           NUM_ANGLES,              &
                           NUM_TORSIONS
  IMPLICIT NONE
  INTEGER :: NUMAT
  CHARACTER (LEN=16), ALLOCATABLE :: TXTATM(:)
  CHARACTER (LEN=12), ALLOCATABLE :: TGEO(:,:)
  INTEGER,  ALLOCATABLE           :: LGEO(:,:), LOC(:,:)
  CHARACTER (LEN=10), ALLOCATABLE :: SIMBOL(:)
  
CONTAINS

  SUBROUTINE READ_ZMAT_INPUT
    USE DYNPREPARE,   ONLY:  LABELS => LAB,         &
                             NA,                    &
                             NB,                    &
                             NC,                    &
                             ATMASS => M,           &
                             LOPT,                  &
                             GEO
    USE START_JOB,    ONLY:  NUMATM => NUMAT,       &
                             GEOFILE
    USE SYSTEM_UTIL,  ONLY:  CP_OPEN,               &
                             CP_CLOSE,              &
                             CPSTOP
    IMPLICIT NONE 
    !Local variables
    CHARACTER (LEN=2), ALLOCATABLE :: ELEMNT(:)
    CHARACTER (LEN=80)             :: LINE, STRING
    DOUBLE PRECISION               :: SUM
    INTEGER, DIMENSION(20)         :: ISTART
    INTEGER                        :: NERR, NA1, MAXTXT, ILINE, NATOMS, NVALUE
    INTEGER                        :: IVAR, NVAR, KERR, MERR
    INTEGER                        :: K, J, L, N, LERR
    INTEGER                        :: I, IREAD
    LOGICAL                        :: LEADSP, UNRECOG, LINE_IS_BLANK

    IREAD=12                                       
    CALL CP_OPEN(IREAD,GEOFILE,'OLD','FORMATTED')  
    REWIND IREAD                                      

    !Allocate Arrays
    IF (.NOT.ALLOCATED(TGEO))         ALLOCATE(TGEO(3,NUMATM))
    IF (.NOT.ALLOCATED(TXTATM))       ALLOCATE(TXTATM(NUMATM))
    IF (.NOT.ALLOCATED(LGEO))         ALLOCATE(LGEO(3,NUMATM))
    IF (.NOT.ALLOCATED(ELEMNT))       ALLOCATE(ELEMNT(NSTORED))

    !Set ELEMNT vector
    !------------------------------------!
    !                                    !
    !   GET ALL CHARACTERS IN UPPER CASE !
    !                                    !
    !------------------------------------!
    ELEMNT = ELEMENTS
    DO J = 1, NSTORED
       DO I = 1, 2
          ILINE=ICHAR(ELEMNT(J)(I:I))
          IF(ILINE.GE.ICHAR('a').AND.ILINE.LE.ICHAR('z')) THEN
             ELEMNT(J)(I:I)=CHAR(ILINE+ICHAR('A')-ICHAR('a'))
          ENDIF
       END DO
    END DO

    ! Begin Instructions
    NERR=0
    NUMAT=0
    NA1=0
    NA=0
    NB=0
    NC=0
    GEO=0.D0
    MAXTXT=0
    ATMASS=0.D0
    LOPT=0

    NUM_BONDS    = 0
    NUM_ANGLES   = 0
    NUM_TORSIONS = 0
    NUMAT        = 0

    WRITE(6,'(A/)')'ZMATRIX READ IN INPUT FILE:'

    NATOMS = 0
    LINE_IS_BLANK = .FALSE.
    Read_zmat_info:  DO WHILE (.TRUE.)
       READ(IREAD,'(A)',END=70,ERR=100) LINE
       WRITE(6,'(A)')LINE
       LINE_IS_BLANK = (LINE.EQ.' ')
       IF(LINE_IS_BLANK) EXIT
       NATOMS = NATOMS + 1
!------------------------------------------------------------------------
!
!   SEE IF TEXT IS ASSOCIATED WITH THIS ELEMENT
!
!------------------------------------------------------------------------
       I=INDEX(LINE,'(')
       IF(I.NE.0)THEN
!------------------------------------------------------------------------
!
!   YES, ELEMENT IS LABELLED.
!
!------------------------------------------------------------------------
          K=INDEX(LINE,')')
          TXTATM(NATOMS)=LINE(I:K)
          MAXTXT=MIN(13,MAX(MAXTXT,K-I+1))
          STRING=LINE(1:I-1)//LINE(K+1:)
          LINE=STRING
       ELSE
          TXTATM(NATOMS)=' '
       ENDIF
!------------------------------------------------------------------------
!
!   GET ALL CHARACTERS IN UPPER CASE
!
!------------------------------------------------------------------------
       DO I = 1, 80
          ILINE=ICHAR(LINE(I:I))
          IF(ILINE.GE.ICHAR('a').AND.ILINE.LE.ICHAR('z')) THEN
             LINE(I:I)=CHAR(ILINE+ICHAR('A')-ICHAR('a'))
          ENDIF
       END DO
!------------------------------------------------------------------------
!
!   GET ELEMENT INFORMATION (ELEMENT LABEL AND MASSES)   
!
!------------------------------------------------------------------------
       NVALUE=0
       LEADSP=.TRUE.
       DO I = 1, 80
          IF (LEADSP.AND.LINE(I:I).NE.' ') THEN
             NVALUE=NVALUE+1
             ISTART(NVALUE)=I
          END IF
          LEADSP=(LINE(I:I).EQ.' ')
       END DO
       
       UNRECOG = .TRUE.
       DO J=1,100
          IF(INDEX(' '//LINE(ISTART(1):ISTART(1)+2),ELEMNT(J)//' ').NE.0) THEN
             UNRECOG = .FALSE.
             EXIT
          END IF
       END DO
       IF(INDEX(' '//LINE(ISTART(1):ISTART(1)+2),' X').NE.0) THEN
          J=99
          UNRECOG = .FALSE.
       ENDIF
       IF (UNRECOG) THEN
          WRITE(6,'(2A)')' ELEMENT NOT RECOGNIZED: ',LINE(ISTART(1):ISTART(1)+2)
          NERR=NERR+1
       ELSE
          LABELS(NATOMS)=ELEMNT(J)
       END IF

       IF(J.NE.99)THEN
          NUMAT=NUMAT+1
          ATMASS(NUMAT)=READA(LINE(1:MAX(ISTART(2)-1,1)),ISTART(1))
          IF(ATMASS(NUMAT).GT.1.D-15)THEN
             WRITE(6,'('' FOR ATOM'',I4,''  ISOTOPIC MASS:'',F15.5)')NATOMS, ATMASS(NUMAT)
          ELSE
             ATMASS(NUMAT)=AMS(J)
          ENDIF
       ENDIF

       TGEO(1,NATOMS)=' '
       TGEO(2,NATOMS)=' '
       TGEO(3,NATOMS)=' '
       IF(NATOMS.EQ.1) CYCLE
       NA(NATOMS)=NINT(READA(LINE,ISTART(8)))
       NUM_BONDS = NUM_BONDS + 1
       CALL GETVAL(LINE(ISTART(2):),GEO(1,NATOMS),TGEO(1,NATOMS))
       LOPT((NATOMS-1)*3+1)=NINT(READA(LINE,ISTART(3)))
       IF(NATOMS.EQ.2) CYCLE
       NB(NATOMS)=NINT(READA(LINE,ISTART(9)))
       NUM_ANGLES = NUM_ANGLES + 1
       CALL GETVAL(LINE(ISTART(4):),GEO(2,NATOMS),TGEO(2,NATOMS))
       LOPT((NATOMS-1)*3+2)=NINT(READA(LINE,ISTART(5)))
       IF(NATOMS.EQ.3) CYCLE
       NC(NATOMS)=NINT(READA(LINE,ISTART(10)))
       NUM_TORSIONS = NUM_TORSIONS + 1
       CALL GETVAL(LINE(ISTART(6):),GEO(3,NATOMS),TGEO(3,NATOMS))
       LOPT((NATOMS-1)*3+3)=NINT(READA(LINE,ISTART(7)))

    END DO Read_zmat_info

70  CONTINUE ! THE FILE IS ENDED... CHECK ON ATOMS NUMBER AND CONTINUE... 
    IF (NUMATM.NE.NATOMS)   GO TO 100
    IF (NERR.NE.0)          GO TO 100  ! Exit whethear there are errors in execution time..

    DO I=1,NATOMS
       DO  J=1,3
          LGEO(J,I)=-1
       END DO
    END DO

    !At this point we check for variable symbol definition..

    IVAR=-1
    NVAR=0
    KERR=0

    !Check how many symbols have been defined

    MERR=0
    DO I=1,NATOMS
       DO  J=1,3
          IF(GEO(J,I).LT.-998) MERR=MERR+1
       END DO
    END DO
!------------------------------------------------------------------------
!
!  IF ALL SYMBOLS ARE DEFINED, THEN DO NOT READ 'FIXED' SYMBOLS
!  MERR : number of symbols in ZMAT definition
!
!------------------------------------------------------------------------
    IF(MERR.EQ.0) GO TO 80
    IVAR=NVAR
    IF (.NOT.ALLOCATED(LOC))          ALLOCATE(LOC(2,MERR))
    IF (.NOT.ALLOCATED(SIMBOL))       ALLOCATE(SIMBOL(MERR))

90  READ(IREAD,'(A)',END=80,ERR=100)LINE
    IF(LINE.EQ.' ')GOTO 90
    
!------------------------------------------------------------------------
!
!  UPPER CASE CONVERSION..
!
!------------------------------------------------------------------------
    DO I=1,80
       ILINE=ICHAR(LINE(I:I))
       IF(ILINE.GE.ICHAR('a').AND.ILINE.LE.ICHAR('z')) THEN
          LINE(I:I)=CHAR(ILINE+ICHAR('A')-ICHAR('a'))
       ENDIF
    END DO
    

    DO I=1,80
       IF(LINE(I:I).NE.' ') EXIT
    END DO
    DO  L=I,I+12
       IF(LINE(L:L).EQ.' ') EXIT
    END DO
    SUM=READA(LINE,L)
    N=0
    LERR=0
    
    DO  J=1,NATOMS
       DO  K=1,3
          IF(TGEO(K,J).EQ.LINE(I:L) .OR. TGEO(K,J)(2:).EQ.LINE(I:L).AND.TGEO(K,J)(1:1).EQ.'-')THEN
             IF(LGEO(K,J).NE.-1)LERR=1
             LGEO(K,J)=LGEO(K,J)+1
             N=N+1
             GEO(K,J)=SUM

             NVAR=NVAR+1
             LOC(1,NVAR)=J
             LOC(2,NVAR)=K
             SIMBOL(NVAR)=TGEO(K,J)(:10)
             IF(SIMBOL(NVAR)(1:1).EQ.'-')THEN
                GEO(K,J)=-SUM
             ENDIF
          ENDIF
       END DO
    END DO

    KERR=KERR+LERR
    IF(LERR.EQ.1)THEN
       WRITE(6,'(2A)')' THE FOLLOWING SYMBOL HAS BEEN DEFINED MORE'//' THAN ONCE:',LINE(I:L)
       NERR=NERR+1
       GO TO 100
    ENDIF
    IF(N.EQ.0)THEN
       WRITE(6,'(2A)')' THE FOLLOWING SYMBOLIC WAS NOT USED:',LINE(I:L)
       NERR=NERR+1
       GO TO 100
    ENDIF
    GOTO 90

    MERR=0
    DO  I=1,NATOMS
       DO  J=1,3
          IF(GEO(J,I).LT.-998)MERR=MERR+1
       END DO
    END DO
!    
    IF(MERR.NE.0)WRITE(6,'(I4,A)')MERR,' GEOMETRY VARIABLES WERE NOT'//' DEFINED'
    IF(MERR+KERR+NERR.NE.0)THEN
       WRITE(6,'(A,I3,A)')' THE GEOMETRY DATA-SET CONTAINED',MERR+KERR+NERR,' ERRORS'
       CALL CPSTOP('ERROR IN ZMAT_READ')
       RETURN
    ENDIF

80  WRITE(6,'(A)')'ZMATRIX READ! CONTINUE CALCULATION'

    !DeAllocate Arrays
    IF (ALLOCATED(TGEO))         DEALLOCATE(TGEO)
    IF (ALLOCATED(TXTATM))       DEALLOCATE(TXTATM)
    IF (ALLOCATED(LGEO))         DEALLOCATE(LGEO)

    DO I=1, NATOMS
       GEO(2,I) = GEO(2,I) * 4.D0 * ATAN(1.D0) / 180.D0
       GEO(3,I) = GEO(3,I) * 4.D0 * ATAN(1.D0) / 180.D0
    END DO

    IF (NUMAT.NE.NUMATM) THEN
       CALL CPSTOP('ZMAT_READ| NUMBER OF ATOMIC SPECIES DOESN''T MATCH!')
    END IF

    RETURN
100 WRITE(6,'(A)')'AN ERROR OCCURED IN READ_ZMAT_INPUT! ZMAT FORMAT ERROR!'
    CALL CPSTOP('ERROR IN ZMAT_READ')
    RETURN
  END SUBROUTINE READ_ZMAT_INPUT



  SUBROUTINE GETVAL(LINE,X,T)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: LINE
    CHARACTER (LEN=12), INTENT(INOUT) :: T
    CHARACTER (LEN=1) :: CH1, CH2
    DOUBLE PRECISION, INTENT(INOUT) :: X
!    INTEGER, EXTERNAL :: ICHAR
    INTEGER :: I

    CH1=LINE(1:1)
    CH2=LINE(2:2)
    IF((ICHAR(CH1).LT.ICHAR('A').OR.ICHAR(CH1).GT.ICHAR('Z')) .AND.   &
       (ICHAR(CH2).LT.ICHAR('A').OR.ICHAR(CH2).GT.ICHAR('Z')))THEN
!------------------------------------------------------------------------
!
!   IS A NUMBER
!
!------------------------------------------------------------------------
       X=READA(LINE,1)
       T=' '
    ELSE
!------------------------------------------------------------------------
!
!   IS A STRING. WE SET X TO -999.DO TO TRACE BACK LATER THE REAL VALUE!
!
!------------------------------------------------------------------------
       I=INDEX(LINE,' ')
       T=LINE(:I)
       X=-999.D0
    ENDIF
    RETURN
  END SUBROUTINE GETVAL


  DOUBLE PRECISION FUNCTION READA(STRING,ISTART)
!------------------------------------------------------------------------
!
!     FORTRAN FUNCTION TO EXTRACT NUMBER FROM STRING
!
!------------------------------------------------------------------------
    IMPLICIT NONE
    ! Arguments
    CHARACTER (LEN=*), INTENT(IN) :: STRING
    INTEGER, INTENT(IN) :: ISTART
    ! Local Variables
    LOGICAL          :: EXPNNT
    INTEGER :: I0, I9, IDOT, INEG, IPOS, ICAPD, ICAPE, ISMLD, ISMLE
    INTEGER :: L, I, J, IADD, N
!------------------------------------------------------------------------
!
!     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
!
!------------------------------------------------------------------------
    I0=ICHAR('0')
    I9=ICHAR('9')
    IDOT=ICHAR('.')
    INEG=ICHAR('-')
    IPOS=ICHAR('+')
    ICAPD=ICHAR('D')
    ICAPE=ICHAR('E')
    ISMLD=ICHAR('d')
    ISMLE=ICHAR('e')
    L=LEN(STRING)
!------------------------------------------------------------------------
!
!     FIND THE START OF THE NUMERIC FIELD
!
!------------------------------------------------------------------------
    DO I=ISTART,L
       IADD=0
       N=ICHAR(STRING(I:I))
!------------------------------------------------------------------------
!
!       SIGNAL START OF NUMERIC FIELD IF DIGIT FOUND
!
!------------------------------------------------------------------------
       IF (N.GE.I0.AND.N.LE.I9) GOTO 20
!------------------------------------------------------------------------
!
!       ACCOUNT FOR CONSECUTIVE SIGNS [- AND(OR) +]
!
!------------------------------------------------------------------------
       IF (N.EQ.INEG.OR.N.EQ.IPOS) THEN
          IADD=IADD+1
          IF(I+IADD.GT.L)GOTO 50
          N=ICHAR(STRING(I+IADD:I+IADD))
          IF(N.GE.I0.AND.N.LE.I9)GOTO 20
       ENDIF
!------------------------------------------------------------------------
!
!       ACCOUNT FOR CONSECUTIVE DECIMAL POINTS (.)
!
!------------------------------------------------------------------------
       IF(N.EQ.IDOT)THEN
          IADD=IADD+1
          IF(I+IADD.GT.L)GOTO 50
          N=ICHAR(STRING(I+IADD:I+IADD))
          IF(N.GE.I0.AND.N.LE.I9)GOTO 20
       ENDIF
    END DO
    GOTO 50
!------------------------------------------------------------------------
!
!     FIND THE END OF THE NUMERIC FIELD
!
!------------------------------------------------------------------------
20  EXPNNT=.FALSE.
    DO J=I+1,L
       IADD=0
       N=ICHAR(STRING(J:J))
!------------------------------------------------------------------------
!
!       CONTINUE SEARCH FOR END IF DIGIT FOUND
!
!------------------------------------------------------------------------
       IF(N.GE.I0.AND.N.LE.I9) CYCLE
!------------------------------------------------------------------------
!
!       CONTINUE SEARCH FOR END IF SIGN FOUND AND EXPNNT TRUE
!
!------------------------------------------------------------------------
       IF(N.EQ.INEG.OR.N.EQ.IPOS)THEN
          IF(.NOT.EXPNNT)GOTO 40
          IADD=IADD+1
          IF(J+IADD.GT.L)GOTO 40
          N=ICHAR(STRING(J+IADD:J+IADD))
          IF(N.GE.I0.AND.N.LE.I9) CYCLE
       ENDIF
       IF(N.EQ.IDOT)THEN
          IADD=IADD+1
          IF(J+IADD.GT.L)GOTO 40
          N=ICHAR(STRING(J+IADD:J+IADD))
          IF(N.GE.I0.AND.N.LE.I9) CYCLE
          IF(N.EQ.ICAPE.OR.N.EQ.ISMLE.OR.N.EQ.ICAPD.OR.N.EQ.ISMLD) CYCLE
       ENDIF
       IF(N.EQ.ICAPE.OR.N.EQ.ISMLE.OR.N.EQ.ICAPD.OR.N.EQ.ISMLD)THEN
          IF(EXPNNT)GOTO 40
          EXPNNT=.TRUE.
          CYCLE
       ENDIF
       GOTO 40
    END DO
    J=L+1
40  N=ICHAR(STRING(J-1:J-1))
    IF(N.EQ.ICAPE.OR.N.EQ.ISMLE.OR.N.EQ.ICAPD.OR.N.EQ.ISMLD)J=J-1
!------------------------------------------------------------------------
!
!     FOUND THE END OF THE NUMERIC FIELD (IT RUNS 'I' THRU 'J-1')
!
!------------------------------------------------------------------------
    N=0
    N=N+INDEX(STRING(I:J-1),'e')
    N=N+INDEX(STRING(I:J-1),'E')
    N=N+INDEX(STRING(I:J-1),'d')
    N=N+INDEX(STRING(I:J-1),'D')
    IF(N.EQ.0)THEN
       READA=DIGIT(STRING(I:J-1),1)
    ELSE
       READA=DIGIT(STRING(:I+N-2),I)*1.D1**DIGIT(STRING(:J-1),I+N)
    ENDIF
    RETURN

!
!     DEFAULT VALUE RETURNED BECAUSE NO NUMERIC FIELD FOUND
!
!------------------------------------------------------------------------
50  READA=0.D0
    RETURN

END FUNCTION READA


DOUBLE PRECISION FUNCTION DIGIT(STRING,ISTART)
!------------------------------------------------------------------------
!
!     FORTRAN FUNCTION TO CONVERT NUMERIC FIELD TO DOUBLE PRECISION
!     NUMBER.  THE STRING IS ASSUMED TO BE CLEAN (NO INVALID DIGIT
!     OR CHARACTER COMBINATIONS FROM ISTART TO THE FIRST NONSPACE,
!     NONDIGIT, NONSIGN, AND NONDECIMAL POINT CHARACTER).
!
!------------------------------------------------------------------------
  CHARACTER (LEN=*) :: STRING
  DOUBLE PRECISION  :: C1, C2, DECIML
  LOGICAL           :: SIGN
  INTEGER           :: ISTART
  INTEGER           :: I0, I9, INEG, IPOS, IDOT, ISPC, L, IDIG, N, I
  INTEGER           :: J
!------------------------------------------------------------------------
!
!     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
!
!------------------------------------------------------------------------
  I0=ICHAR('0')
  I9=ICHAR('9')
  INEG=ICHAR('-')
  IPOS=ICHAR('+')
  IDOT=ICHAR('.')
  ISPC=ICHAR(' ')
  C1=0.D0
  C2=0.D0
  SIGN=.TRUE.
  L=LEN(STRING)
!------------------------------------------------------------------------
!
!     DETERMINE THE CONTRIBUTION TO THE NUMBER GREATER THAN ONE
!
!------------------------------------------------------------------------
  IDIG=0
  DO I=ISTART,L
     N=ICHAR(STRING(I:I))
     IF(N.GE.I0.AND.N.LE.I9)THEN
        IDIG=IDIG+1
        C1=C1*1.D1+N-I0
     ELSEIF(N.EQ.INEG.OR.N.EQ.IPOS.OR.N.EQ.ISPC)THEN
        IF(N.EQ.INEG)SIGN=.FALSE.
     ELSEIF(N.EQ.IDOT)THEN
        GOTO 20
     ELSE
        GOTO 40
     ENDIF
  END DO
!------------------------------------------------------------------------
!
!     DETERMINE THE CONTRIBUTION TO THE NUMBER LESS THAN THAN ONE
!
!------------------------------------------------------------------------
20 DECIML=1.D0
  DO  J=I+1,L
     N=ICHAR(STRING(J:J))
     IF(N.GE.I0.AND.N.LE.I9)THEN
        DECIML=DECIML/1.D1
        C2=C2+(N-I0)*DECIML
     ELSEIF(N.NE.ISPC)THEN
        GOTO 40
     ENDIF
  END DO
!------------------------------------------------------------------------
!
!     PUT THE PIECES TOGETHER
!
!------------------------------------------------------------------------
40 DIGIT=C1+C2
  IF(.NOT.SIGN)DIGIT=-DIGIT
  RETURN
END FUNCTION DIGIT





SUBROUTINE ZMAT_TO_XYZ(GEO,COORD,NA,NB,NC,NATOMS)
  IMPLICIT NONE
  !ARGUMENTS
  INTEGER,   INTENT(IN)  :: NATOMS
  INTEGER,   INTENT(INOUT), DIMENSION(NATOMS) ::  NA, NB, NC
  DOUBLE PRECISION, INTENT(IN), DIMENSION(3,NATOMS)  :: GEO
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,NATOMS) :: COORD
  DOUBLE PRECISION, DIMENSION(3,3) :: TVEC
  !LOCAL VARIABLES
  INTEGER :: I, J, K, MA, MB, MC
  DOUBLE PRECISION :: CCOS, COSA, XB, YB, ZB, XA, YA, ZA, RBC, XYB
  DOUBLE PRECISION :: XPA, YPA, XPB, COSTH, SINTH, COSKH, SINKH 
  DOUBLE PRECISION :: ZQA, YZA, SINA, SIND, COSD, XD, YD, ZD
  DOUBLE PRECISION :: XPQ, YPD, ZPD, XQD, YQD, ZQD, XRD
  DOUBLE PRECISION :: SINPH, COSPH, XPD
!------------------------------------------------------------------------
!
!    ZMAT_TO_XYZ  COMPUTES COORDINATES FROM BOND-ANGLES AND LENGTHS.
!    IT IS ADAPTED FROM THE PROGRAM WRITTEN BY M.J.S. DEWAR.
!
!        NORMAL CONVERSION FROM INTERNAL TO CARTESIAN COORDINATES
!        IS DONE.
! 
!  ON INPUT:
!         GEO    = ARRAY OF INTERNAL COORDINATES.
!         NATOMS = NUMBER OF ATOMS, INCLUDING DUMMIES.
!         NA     = ARRAY OF ATOM LABELS FOR BOND LENGTHS.
!                  
!
!
!
!  ON OUTPUT:
!         COORD  = ARRAY OF CARTESIAN COORDINATES
!
!------------------------------------------------------------------------

! Copy atomic position of first atom in COORD array
  DO I=1,MIN(3,NATOMS)
     DO J=1,3
        ! First Atom..
        COORD(J,I)=GEO(J,I)
     END DO
  END DO
  
  IF(NATOMS.GT.1)THEN
     IF(NA(2).EQ.1)THEN
! Move second atom along X direction
        ! Second Atom...
        COORD(1,2)=COORD(1,1)+GEO(1,2)
        COORD(2,2)=COORD(2,1)
        COORD(3,2)=COORD(3,1)
     ENDIF
  ELSE
! Finish the job
     GOTO 120
  ENDIF

  IF(NATOMS.EQ.2) GOTO 120
  CCOS=COS(GEO(2,3))
  ! Third Atom...
  IF(NA(3).EQ.1)THEN
     COORD(1,3)=COORD(1,1)+GEO(1,3)*CCOS
     COORD(2,3)=COORD(2,1)+GEO(1,3)*SIN(GEO(2,3))
     COORD(3,3)=COORD(3,1)
  ELSEIF(NA(3).EQ.2)THEN
     COORD(1,3)=COORD(1,2)-GEO(1,3)*CCOS
     COORD(2,3)=COORD(2,2)+GEO(1,3)*SIN(GEO(2,3))
     COORD(3,3)=COORD(3,2)
  ELSE
     COORD(1,3)=GEO(1,3)
     COORD(2,3)=GEO(2,3)
     COORD(3,3)=GEO(3,3)
  ENDIF

  DO I=4,NATOMS
     MC=NA(I)
     IF(MC.NE.0)THEN
        COSA=COS(GEO(2,I))
        MB=NB(I)
        XB=COORD(1,MB)-COORD(1,MC)
        YB=COORD(2,MB)-COORD(2,MC)
        ZB=COORD(3,MB)-COORD(3,MC)
        RBC=XB*XB+YB*YB+ZB*ZB
        IF(RBC.LT.1.D-16)THEN
!------------------------------------------------------------------------
!
!     TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.
!
!------------------------------------------------------------------------
           WRITE(6,'(A,I4,A,I4,A)')' ATOMS',MB,' AND',MC,' ARE COINCIDENT'
           WRITE(6,'(A)')' THIS IS A FATAL ERROR, RUN STOPPED IN ZMAT_TO_XYZ!'
           STOP
        ELSE
           RBC=1.0D00/DSQRT(RBC)
        ENDIF
50      MA=NC(I)
        XA=COORD(1,MA)-COORD(1,MC)
        YA=COORD(2,MA)-COORD(2,MC)
        ZA=COORD(3,MA)-COORD(3,MC)
!------------------------------------------------------------------------
!
!     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
!     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
!
!------------------------------------------------------------------------
        XYB=SQRT(XB*XB+YB*YB)
        K=-1
        IF (XYB.GT.0.1D00) GO TO 60
        XPA=ZA
        ZA=-XA
        XA=XPA
        XPB=ZB
        ZB=-XB
        XB=XPB
        XYB=SQRT(XB*XB+YB*YB)
        K=+1
!------------------------------------------------------------------------
!
!     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
!
!------------------------------------------------------------------------
60      COSTH=XB/XYB
        SINTH=YB/XYB
        XPA=XA*COSTH+YA*SINTH
        YPA=YA*COSTH-XA*SINTH
        SINPH=ZB*RBC
        COSPH=SQRT(ABS(1.D00-SINPH*SINPH))
        ZQA=ZA*COSPH-XPA*SINPH
!------------------------------------------------------------------------
!
!     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
!
!------------------------------------------------------------------------
        YZA=SQRT(YPA**2+ZQA**2)
        IF(YZA.LT.1.D-4)GOTO 80
        IF(ABS(COSA).LT.0.9998D0)THEN
           IF(YZA.LT.2.D-2)THEN
!------------------------------------------------------------------------
!
!   FAULTY ATOM: I; BOND LENGTH TO I: NA(I); ANGLE TO I: NB(I)
!
!   SET A NEW NC(I) AND RE-TRY
!
!------------------------------------------------------------------------
              J=NC(I)
              CALL RENUM(COORD,NA,NB,NC,I,NATOMS)
              IF (NC(I).EQ.0.OR.NC(I).EQ.J) THEN
                 WRITE(6,'(//20X,'' CALCULATION ABANDONED AT THIS POINT'')')
                 IF (NC(I).EQ.0) THEN
                    WRITE(6,'(//10X,'' THREE ATOMS BEING USED TO DEFINE THE'',' &
                         //'/10X,'' COORDINATES OF A FOURTH ATOM, WHOSE BOND-ANGLE IS'')')
                    WRITE(6,'(10X,'' NOT ZERO OR 180 DEGREES, ARE IN AN ALMOST STRAIGHT'')')
                    WRITE(6,'(10X,'' LINE.  THERE IS A HIGH PROBABILITY THAT THE'', ' & 
                         //'/10X,'' COORDINATES OF THE ATOM WILL BE INCORRECT.'')')
                    WRITE(6,'(//20X,''THE FAULTY ATOM IS ATOM NUMBER'',I6)')I
                 ELSE
                    WRITE(6,'(A)')' THE CARTESIAN GEOMETRY HAS BECOME UNREASONABLE !'
                 ENDIF
                 IF(NC(I).EQ.0)THEN
                    WRITE(6,'(//20X,''CARTESIAN COORDINATES UP TO FAULTY ATOM'')')
                    WRITE(6,'(//5X,''I'',12X,''X'',12X,''Y'',12X,''Z'')')
                    DO J=1,I
                       WRITE(6,'(I6,F16.5,2F13.5)')J,(COORD(K,J),K=1,3)
                    END DO
                    WRITE(6,'(//6X,'' ATOMS'',I3,'','',I3,'', AND'',I3,'' ARE WITHIN'','  &
                         //'F7.4,'' ANGSTROMS OF A STRAIGHT LINE'')')MC,MB,MA,YZA
                 ENDIF
                 STOP
                 RETURN
              ENDIF
              GOTO 50
           ENDIF
        ENDIF
        COSKH=YPA/YZA
        SINKH=ZQA/YZA
        GOTO 90
80      CONTINUE
!------------------------------------------------------------------------
!
!   ANGLE TOO SMALL TO BE IMPORTANT
!
!------------------------------------------------------------------------
        COSKH=1.D0
        SINKH=0.D0
90      CONTINUE
!------------------------------------------------------------------------
!
!     COORDINATES :-   A=(???,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
!     NONE ARE NEGATIVE.
!     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
!
!------------------------------------------------------------------------
        SINA=SIN(GEO(2,I))
        SIND=-SIN(GEO(3,I))
        COSD=COS(GEO(3,I))
        XD=GEO(1,I)*COSA
        YD=GEO(1,I)*SINA*COSD
        ZD=GEO(1,I)*SINA*SIND
!------------------------------------------------------------------------
!
!     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
!
!------------------------------------------------------------------------
        YPD=YD*COSKH-ZD*SINKH
        ZPD=ZD*COSKH+YD*SINKH
        XPD=XD*COSPH-ZPD*SINPH
        ZQD=ZPD*COSPH+XD*SINPH
        XQD=XPD*COSTH-YPD*SINTH
        YQD=YPD*COSTH+XPD*SINTH
        IF (K.LT.1) GO TO 100
        XRD=-ZQD
        ZQD=XQD
        XQD=XRD
100     COORD(1,I)=XQD+COORD(1,MC)
        COORD(2,I)=YQD+COORD(2,MC)
        COORD(3,I)=ZQD+COORD(3,MC)
     ELSE
        COORD(1,I)=GEO(1,I)
        COORD(2,I)=GEO(2,I)
        COORD(3,I)=GEO(3,I)
     ENDIF
  END DO

120 CONTINUE !The End..
  RETURN
END SUBROUTINE ZMAT_TO_XYZ

SUBROUTINE RENUM(COORD,NA,NB,NC,II,NATOMS)
  IMPLICIT NONE
  !ARGUMENTS
  INTEGER, INTENT(IN),    DIMENSION(NATOMS) ::  NA, NB
  INTEGER, DIMENSION(NATOMS), INTENT(INOUT) ::  NC
  INTEGER, INTENT(IN) :: II, NATOMS
  DOUBLE PRECISION, INTENT(IN), DIMENSION(3,NATOMS) :: COORD
  !LOCAL VARIABLES
  INTEGER :: NAI, NBI, JJ, I
  DOUBLE PRECISION :: THETA, RMIN, ANGLE, RAB
!------------------------------------------------------------------------
!
!  RENUMBER THE NC OF ATOM II.  ON INPUT, THE ANGLE NA(II)-NB(II)-NC(II)
!  IS TOO NEAR TO 0 OR 180 DEGREES.  FIND A NEW ATOM FOR NC(II), SO THAT
!  THE ANGLE WILL BE ACCEPTABLE (AS LARGE AS POSSIBLE)
!
!------------------------------------------------------------------------

  NAI=NA(II)
  NBI=NB(II)
!------------------------------------------------------------------------
!
!   THETA = 45 DEGREES
!
!------------------------------------------------------------------------
  THETA=0.7853D0
  JJ=0
  RMIN=1.D10
10 DO 20 I=1,II-1
     IF (I.EQ.NAI.OR.I.EQ.NBI) GOTO 20
     CALL BANGLE(COORD,NAI,NBI,I,ANGLE)
     IF (ANGLE.GT.1.5707963D0) ANGLE=2.D0*ASIN(1.D0)-ANGLE
     IF (ANGLE.LT.THETA) GOTO 20
!------------------------------------------------------------------------
!
!   ANGLE IS OK.  NOW FIND ATOM OF LOWEST DISTANCE
!
!------------------------------------------------------------------------
     RAB=(COORD(1,NBI)-COORD(1,I))**2+ &
          (COORD(2,NBI)-COORD(2,I))**2+ &
          (COORD(3,NBI)-COORD(3,I))**2
     IF(RAB.LT.RMIN)THEN
        JJ=I
        RMIN=RAB
     ENDIF
20 END DO
  IF(JJ.EQ.0)THEN
!------------------------------------------------------------------------
!
!   NO ATOM INSIDE THE ALLOWED ANGLE - REDUCE THE ANGLE
!
!------------------------------------------------------------------------
     THETA=THETA*0.5D0
     IF (THETA.LT.0.0174533D0) THEN
        NC(II)=0
        RETURN
     ENDIF
     GOTO 10
  ENDIF
!------------------------------------------------------------------------
!
!  BEST NC IS JJ; BEST ANGLE IS THMIN
!
!------------------------------------------------------------------------
  NC(II)=JJ
  RETURN
END SUBROUTINE RENUM

SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)
  IMPLICIT NONE
  ! ARGUMENTS
  DOUBLE PRECISION, INTENT(IN), DIMENSION(3,*) :: XYZ
  INTEGER, INTENT(IN) :: I, J, K
  DOUBLE PRECISION, INTENT(OUT) :: ANGLE
  ! LOCAL VARIABLES
  DOUBLE PRECISION :: D2IJ, D2JK, D2IK, XY, TEMP
!------------------------------------------------------------------------
!
! BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
!        CARTESIAN COORDINATES ARE IN XYZ.
!
!------------------------------------------------------------------------
  D2IJ = (XYZ(1,I)-XYZ(1,J))**2+    &
         (XYZ(2,I)-XYZ(2,J))**2+    & 
         (XYZ(3,I)-XYZ(3,J))**2
  D2JK = (XYZ(1,J)-XYZ(1,K))**2+    &
         (XYZ(2,J)-XYZ(2,K))**2+    &
         (XYZ(3,J)-XYZ(3,K))**2
  D2IK = (XYZ(1,I)-XYZ(1,K))**2+    &
         (XYZ(2,I)-XYZ(2,K))**2+    &
         (XYZ(3,I)-XYZ(3,K))**2
  XY = SQRT(D2IJ*D2JK)
  TEMP = 0.5D0 * (D2IJ+D2JK-D2IK) / XY
  IF (TEMP .LT. -1.D0+1.D-12)TEMP=-1.0D0
  IF (TEMP .GT.  1.D0-1.D-12)TEMP=1.D0
  ANGLE = ACOS( TEMP )
  RETURN
END SUBROUTINE BANGLE

SUBROUTINE XYZ_TO_INTERNAL(COORD, QINT, NUMAT, NUM_VAR, NA, NB, NC)
  USE GEO_OBJECTS, ONLY: BOND_VAL,     &
                         ANGLE_VAL,    &
                         TORS_VAL
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(IN) :: NUMAT, NUM_VAR
  INTEGER, INTENT(IN), DIMENSION(NUMAT) :: NA, NB, NC
  DOUBLE PRECISION, INTENT(IN), DIMENSION(3*NUMAT)    :: COORD
  DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NUM_VAR) :: QINT 
  ! Local Variables
  INTEGER :: I, IND

  IND = 0
  DO I = 2, NUMAT
     IND=IND+1
     QINT(IND) = BOND_VAL(I, NA(I), COORD)
     IF (I.EQ.2) CYCLE
     IND=IND+1
     QINT(IND) = ANGLE_VAL(I, NA(I), NB(I), COORD)
     IF (I.EQ.3) CYCLE
     IND=IND+1
     QINT(IND) = TORS_VAL(I, NA(I), NB(I), NC(I), COORD)
     IF (I.EQ.4) CYCLE
  END DO

  RETURN
END SUBROUTINE XYZ_TO_INTERNAL

END MODULE READ_ZMAT
