!------------------------------------------------------------------------
MODULE DYNPREPARE                                                       !
      USE START_JOB,     ONLY: ZMATYPE                                  !
      USE XYZ,           ONLY: POINT                                    !
      IMPLICIT NONE                                                     !
!     ------------------------------------------------------------------!
!     Define Coordinates vector, Force, Velocities, Acceleration, Masses!
!     ------------------------------------------------------------------! 
      TYPE(POINT), ALLOCATABLE       :: C(:),             &! Cart. Coord!
                                        F(:),             &! Cart. Forcs!
                                        V(:),             &! Cart. Veloc!
                                        A(:),             &! Cart. Accel!
                                        F2(:),            &! Exc.St.Grad!
                                        FSAVE(:)           ! Old.Val.Fs.!
      TYPE(POINT), ALLOCATABLE       :: OLDF(:),          &! Old Forces !
                                        OLDA(:)            ! Old Acceler!
      DOUBLE PRECISION, ALLOCATABLE  :: M(:),             &! Masses Vect!
                                        GHMATRIX(:,:)      ! New Matr GH!
      DOUBLE PRECISION, ALLOCATABLE  :: GHMATRIXOLD(:,:)   ! Old Matr GH!
      CHARACTER (LEN=2), ALLOCATABLE :: LAB(:)             ! Label Vectr!
      ! Internal Coordinates                               
      INTEGER, ALLOCATABLE           :: NA(:),            &! Bond  Flags!
                                        NB(:),            &! Angle Flags!
                                        NC(:)              ! Torsion Fl.!
      INTEGER, ALLOCATABLE           :: LOPT(:)            ! Opt. Flags !     
      DOUBLE PRECISION, ALLOCATABLE  :: GEO(:,:)           ! Int. Coords!
      DOUBLE PRECISION, ALLOCATABLE  :: CX(:)              ! Cart.Coords!
      DOUBLE PRECISION, ALLOCATABLE  :: VX(:)              ! Cart.Veloc.!
      DOUBLE PRECISION, ALLOCATABLE  :: FX(:)              ! Cart.Forces!
      INTEGER                        :: NUM_BONDS,        &             !
                                        NUM_ANGLES,       &             !
                                        NUM_TORSIONS                    !
!     ------------------------------------------------------------------!
      CONTAINS                                                          !
!     ===================================================================
      SUBROUTINE DYN_ALLOCATE                                           
!     ===================================================================
      USE START_JOB,     ONLY: NUMAT,       &                           !
                               GEOFILE,     &                           !
                               COORDINATES                              !
      IMPLICIT NONE                                                     !
      INTEGER :: I                                                      !
                                                                        !
      ALLOCATE  (  LAB(NUMAT),      &                                   !
                     M(NUMAT)           )                               !
!     ------------------------------------------------------------------!
!     Initialize with zero values                                       ! 
!     ------------------------------------------------------------------! 
      DO I=1,NUMAT                                                      !
         LAB(I)='  '                                                    !
         M(I)=0.D0                                                      !
      END DO                                                            !
!     ------------------------------------------------------------------!
                                                                        !
      SELECT CASE (COORDINATES)                                         !
      CASE ('CARTESIAN')                                                !
                                                                        !
         ALLOCATE (    C(NUMAT),    &                                   !
                       F(NUMAT),    &                                   !
                      F2(NUMAT),    &                                   !
                   FSAVE(NUMAT),    &                                   !
                       V(NUMAT),    &                                   !
                       A(NUMAT),    &                                   !
                    OLDF(NUMAT),    &                                   !
                    OLDA(NUMAT),    &                                   !
         GHMATRIXOLD(2,3*NUMAT),    &                                   !
            GHMATRIX(3*NUMAT,2)          )                              !
!     ------------------------------------------------------------------!
!     Initialize with zero values                                       ! 
!     ------------------------------------------------------------------! 
         DO I=1,NUMAT                                                   !
            C(I)%X=0.D0                                                 !
            C(I)%Y=0.D0                                                 !
            C(I)%Z=0.D0                                                 !
                                                                        !
            F(I)%X=0.D0                                                 !
            F(I)%Y=0.D0                                                 !
            F(I)%Z=0.D0                                                 !
                                                                        !
            F2(I)%X=0.D0                                                !
            F2(I)%Y=0.D0                                                !
            F2(I)%Z=0.D0                                                !
                                                                        !
            V(I)%X=0.D0                                                 !
            V(I)%Y=0.D0                                                 !
            V(I)%Z=0.D0                                                 !
                                                                        !
            A(I)%X=0.D0                                                 !
            A(I)%Y=0.D0                                                 !
            A(I)%Z=0.D0                                                 !
                                                                        !
            OLDF(I)%X=0.D0                                              !
            OLDF(I)%Y=0.D0                                              !
            OLDF(I)%Z=0.D0                                              !
                                                                        !
            OLDA(I)%X=0.D0                                              !
            OLDA(I)%Y=0.D0                                              !
            OLDA(I)%Z=0.D0                                              !
         END DO                                                         !
         GHMATRIX=0.D0                                                  !
         GHMATRIXOLD=0.D0                                               !
                                                                        !
      CASE ('ZMATCOORD')                                                !
                                                                        !
         ALLOCATE (  NA(NUMAT),     &                                   !
                     NB(NUMAT),     &                                   !
                     NC(NUMAT),     &                                   !
                     LOPT(3*NUMAT), &                                   !
                     GEO(3,NUMAT),  &                                   !
                     CX(3*NUMAT),   &                                   !
                     VX(3*NUMAT),   &                                   !
                     FX(3*NUMAT)       )                                !
                                                                        !
         ALLOCATE (    C(NUMAT),    &                                   !
                       F(NUMAT),    &                                   !
                       V(NUMAT)        )                                !
                                                                        !
         DO I=1,NUMAT                                                   !
            C(I)%X=0.D0                                                 !
            C(I)%Y=0.D0                                                 !
            C(I)%Z=0.D0                                                 !
                                                                        !
            F(I)%X=0.D0                                                 !
            F(I)%Y=0.D0                                                 !
            F(I)%Z=0.D0                                                 !
                                                                        !
            V(I)%X=0.D0                                                 !
            V(I)%Y=0.D0                                                 !
            V(I)%Z=0.D0                                                 !
         END DO                                                         !
                                                                        !
         NA     = 0                                                     !
         NB     = 0                                                     !
         NC     = 0                                                     !
         LOPT   = 0                                                     !
         GEO    = 0.D0                                                  !
         CX     = 0.D0                                                  !
         VX     = 0.D0                                                  !
         FX     = 0.D0                                                  !
         NUM_BONDS    = 0                                               !
         NUM_ANGLES   = 0                                               !
         NUM_TORSIONS = 0                                               !
                                                                        !
                                                                        !
      END SELECT                                                        !
!     ===================================================================
      END SUBROUTINE DYN_ALLOCATE                                       !
!     ===================================================================
      SUBROUTINE DYN_RDGEO                                              
!     ===================================================================      
      USE SYSTEM_UTIL,  ONLY:  CPSTOP,      &                           !
                               CP_OPEN,     &                           !
                               CP_CLOSE,    &                           !
                               EATW,        &                           !
                               LFBL                                     !
      USE START_JOB,    ONLY:  GEOFILE,     &                           !
                               NUMAT,       &                           !
                               COORDINATES, &                           !
                               TYPEZMAT                                 !
      USE PARSING,      ONLY:  GATHEREXP,   &                           !
                               GETFIELD                                 !
      USE CONVFACTORS,  ONLY:  ANGTOBOHR                                !
      IMPLICIT NONE                                                     !
      INTEGER :: ICHANNEL, INDT, IAT, IFIE, I, K                        !
      LOGICAL :: CONFERMA                                               !
      CHARACTER (LEN=512) :: LINE                                       !
      CHARACTER (LEN=256) :: FIELD,FIELD2                               !
                                                                        !
      ICHANNEL=12                                                       !
      CONFERMA=.FALSE.                                                  !
      CALL CP_OPEN(ICHANNEL,GEOFILE,'OLD','FORMATTED')                  !
      REWIND ICHANNEL                                                   !
      CALL GATHEREXP(ICHANNEL,'{','}',LINE,CONFERMA)                    !
      IF (.NOT.CONFERMA) THEN                                           !
         SELECT CASE (TYPEZMAT)                                         !
         CASE ('MOLPRO')                                                !
            WRITE(6,'(A)')'UNKNOWN ZMAT FORMAT.. MOLPRO TYPE EXPECTED!' !
            CALL CPSTOP('DYN_RDZMAT')                                   !
         CASE DEFAULT                                                   !
            ! Read all the stuff in Main directory                      !
            RETURN                                                      !
         END SELECT                                                     !
      ENDIF                                                             !
!     ------------------------------------------------------------------!
!     Sligtly modification of line record to improve the reading        !
!     ------------------------------------------------------------------!
      INDT=INDEX(LINE,'}')                                              !
      LINE(INDT:INDT)=';'                                               !
!     ------------------------------------------------------------------!
      CALL GETFIELD(FIELD,LEN(FIELD),';',';',LINE,1)                    !
      SELECT CASE (COORDINATES)                                         !
      CASE ('CARTESIAN')                                                !
         CALL GETFIELD(FIELD,LEN(FIELD),';',';',LINE,1)                 !
      END SELECT                                                        !
      WRITE(6,'(72A)')('-',I=1,72)                                      !
      WRITE(6,'(/15X,A/)')'STARTING GEOMETRY READ FROM ABINITIO INPUT'  !
      DO IAT=1, NUMAT                                                   !
         SELECT CASE(COORDINATES)                                       !
         CASE ('CARTESIAN')                                             !
            CALL GETFIELD(FIELD,LEN(FIELD),';',';',LINE,IAT+2)          !
         CASE ('ZMATCOORD')                                             !
            CALL GETFIELD(FIELD,LEN(FIELD),';',';',LINE,IAT  )          !
         END SELECT                                                     !
!     ------------------------------------------------------------------!
!     Sligtly modification of field record to improve the next reading  !
!     ------------------------------------------------------------------!
         CALL EATW(FIELD,LEN(FIELD))                                    !
         IFIE=LFBL(FIELD,LEN(FIELD))                                    !
         FIELD(1:)=','//FIELD(1:IFIE)//','                              !
!     ------------------------------------------------------------------!
         SELECT CASE (COORDINATES)                                      !
         CASE ('CARTESIAN')                                             !
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,1)           !
            CALL EATW(FIELD2,LEN(FIELD2))                               !
            LAB(IAT)=FIELD2(1:2)                                        !
            IF (LFBL(FIELD2,LEN(FIELD2)).EQ.1) LAB(IAT)=' '//FIELD2(1:1)!
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,2)           !
            READ(FIELD2,*)C(IAT)%X                                      !
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,3)           !
            READ(FIELD2,*)C(IAT)%Y                                      !
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,4)           !
            READ(FIELD2,*)C(IAT)%Z                                      !
            C(IAT)%X=C(IAT)%X*ANGTOBOHR                                 !
            C(IAT)%Y=C(IAT)%Y*ANGTOBOHR                                 !
            C(IAT)%Z=C(IAT)%Z*ANGTOBOHR                                 !
            WRITE(6,'(A2,3F15.9)')LAB(IAT),C(IAT)%X,C(IAT)%Y,C(IAT)%Z   !
         CASE ('ZMATCOORD')                                             !
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,1)           !
            CALL EATW(FIELD2,LEN(FIELD2))                               !
            LAB(IAT)=FIELD2(1:2)                                        !
            IF (LFBL(FIELD2,LEN(FIELD2)).EQ.1) LAB(IAT)=' '//FIELD2(1:1)!
            IF (IAT.EQ.1) GO TO 50                                      !
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,2)           !
            READ(FIELD2,*)NA(IAT)
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,3) 
            READ(FIELD2,*)GEO(1,IAT)                                    !
            NUM_BONDS = NUM_BONDS + 1
            IF (IAT.EQ.2) GO TO 50                                      !
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,4)           !
            READ(FIELD2,*)NB(IAT)
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,5)           !
            READ(FIELD2,*)GEO(2,IAT)                                    !
            NUM_ANGLES = NUM_ANGLES + 1
            IF (IAT.EQ.3) GO TO 50                                      !
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,6)           !
            READ(FIELD2,*)NC(IAT)
            CALL GETFIELD(FIELD2,LEN(FIELD2),',',',',FIELD,7)           !
            READ(FIELD2,*)GEO(3,IAT)                                    !
            NUM_TORSIONS = NUM_TORSIONS + 1
50          WRITE(6,100)LAB(IAT),GEO(1,IAT),NA(IAT),    &               !
                                 GEO(2,IAT),NB(IAT),    &               !
                                 GEO(3,IAT),NC(IAT)                     !
100         FORMAT(A2,5X,F12.6,I5,F12.6,I5,F12.6,I5)                    !
         END SELECT                                                     !
      END DO                                                            !
      WRITE(6,'(72A)')('-',I=1,72)                                      !
!     ------------------------------------------------------------------!
!     Get information on atomic masses                                  !
!     ------------------------------------------------------------------!
      CALL ASSIGNMASS(LAB,M,NUMAT)                                      !
!     ------------------------------------------------------------------!
      CALL CP_CLOSE(ICHANNEL,GEOFILE)                                   ! 
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE DYN_RDGEO                                          !
!     ===================================================================
      SUBROUTINE ASSIGNMASS(LAB,MPRV,NUMAT)
!     ===================================================================
      USE SYSTEM_UTIL,     ONLY  :  CPSTOP
      USE PERIODICTABLE,   ONLY  :  MASSES,    &                        !
                                    ELEMENTS,  &                        !
                                    NSTORED                             !
      USE CONVFACTORS,     ONLY  :  AMUTOAU                             ! 
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      CHARACTER (LEN=2), INTENT(IN)   :: LAB(*)                         !
      DOUBLE PRECISION, INTENT(OUT)   :: MPRV(*)                        !
      INTEGER (KIND=8), INTENT(IN)    :: NUMAT                          !
      INTEGER :: I,J                                                    !
      LOGICAL :: FOUND                                                  !
                                                                        !
      DO I=1,NUMAT                                                      !
         FOUND=.FALSE.                                                  !
         DO J=1,NSTORED                                                 !
            IF (LAB(I).EQ.ELEMENTS(J)) EXIT                             !
         END DO                                                         !
         IF (J.EQ.NSTORED+1) GO TO 99                                   !
         MPRV(I)=MASSES(J)*AMUTOAU                                      !
      END DO                                                            !
                                                                        !
      RETURN                                                            !
!     ------------------------------------------------------------------!
 99   WRITE(*,'(A)')'UNABLE TO FIND ATOM''S LABEL:',LAB(I)              !
      CALL CPSTOP('ASSIGNMASS')                                         !
!     ------------------------------------------------------------------!
      END SUBROUTINE ASSIGNMASS                                         !
!     ===================================================================
END MODULE DYNPREPARE
