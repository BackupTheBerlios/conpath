!------------------------------------------------------------------------
MODULE MODKINETIC                                                       !
      IMPLICIT NONE                                                     !
!     ------------------------------------------------------------------!
      CONTAINS                                                          !
!     ===================================================================
      SUBROUTINE KINETIC ( M, V, K, N )                                 
!     ===================================================================
      USE START_JOB, ONLY: IRESVEL, &                                   !
                           INI_TEMP, &                                     !
                           TEMP_SCALE, &
                           TEMP_MAX
      USE  XYZ,   ONLY:   POINT,    &                                   !
                          PROD3                                         !
      USE DYNLIB, ONLY:   VEL_RESCALE                                   !
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      INTEGER (KIND=8), INTENT(IN)     :: N                             !
      DOUBLE PRECISION, INTENT(IN)     :: M(*)                          !
      TYPE(POINT),      INTENT(INOUT)  :: V(*)                          !
      REAL    (KIND=8)                 :: K                             !
      INTEGER                          :: I                             !
       
                                                                        !
      IF (INI_TEMP.LT.TEMP_MAX) INI_TEMP=INI_TEMP*TEMP_SCALE
      IF (IRESVEL.GT.0) CALL VEL_RESCALE(INI_TEMP, V, N)                !
      K=0.D0                                                            !
      K=PROD3(M,V,V,INT(N,4))                                           !
      K=0.5D0*K                                                         !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE KINETIC                                            !
!     ===================================================================      
END MODULE MODKINETIC

!------------------------------------------------------------------------
MODULE FORCE                                                            !
      IMPLICIT NONE                                                     !
!     ------------------------------------------------------------------!     
      CONTAINS                                                          !
!     ===================================================================
      SUBROUTINE FORCE_CART( VEXC , VFOND,  UPOT, STEP )
!     ===================================================================
      USE DYNPREPARE,  ONLY: C,               &                         !
                             F,               &    ! Diff.Energy Grad.  !
                             F2,              &    ! Excited State Grad.!
                             M,               &                         !
                             A,               &    ! NAC Gradient       !
                             LAB,             &                         !
                             FSAVE,           &                         !
                             OLDF,            &                         !
                             OLDA,            &                         !
                             FX                                         !
      USE SYSTEM_UTIL, ONLY: PSYST,           &                         !
                             COMMAND,         &                         !
                             CPSTOP,          &                         !
                             EATW,            &                         !
                             LFBL,            &                         !
                             CP_OPEN,         &                         !
                             CP_CLOSE,        &                         !
                             SEARCH                                     !
      USE START_JOB,   ONLY: NUMAT,           &                         !
                             GEOFILE,         &                         !
                             ABINIEXE,        &                         !
                             ABINIT_INP,      &                         !
                             ABINIT_OUT,      &                         !
                             ABINIPROG,       &                         !
                             IPRINT,          &                         !
                             DYNTYPE,         &                         !
                             USEGRADEX,       &                         !
                             USESQUARENE,     &                         !
                             LFACT,           &                         !
                             ADDREPULSIVE                               !
      USE MOLPRO_LINK, ONLY: GET_VECTOR_MOLPRO,   &                     !
                             GET_ENERGIES_MOLPRO,   &                   !
                             GET_ENERGY_MOLPRO,   &                     !
                             GET_CHARGES                                !
      USE GEOMODULE,   ONLY: CREATE_GEOFILE                             !
!     ------------------------------------------------------------------!
!     ROUTINE TO COMPUTE FORCES WITH AN EXTERNAL ABINITIO PROGRAM.      !
!     THE ROUTINE IS COMBINED WITH THE LEAPFROG ALGORITHM               !
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      DOUBLE PRECISION      ::   VEXC, VFOND, UPOT                      !
      INTEGER, INTENT(IN)   ::   STEP                                   !
      INTEGER               ::   I, J, ICHANNEL, K                      !
      CHARACTER  (LEN=256)  ::   OUTPUTFILE                             !
      CHARACTER  (LEN=20)   ::   LABEL                                  !
      LOGICAL               ::   CONFERMA                               !
      CHARACTER  (LEN=256)  ::   LINE                                   !
      CHARACTER  (LEN=20)   ::   TARGET                                 !
      CHARACTER  (LEN=256)  ::   WFUDIR                                 !
      INTEGER, SAVE :: IFAIL                                            !
      DATA IFAIL / 0 /                                                  !
!     Get enviromental variable                                         !
      CALL GETENV('wfudir',WFUDIR)                                      !
      K=LFBL(WFUDIR,LEN(WFUDIR))                                        !
      WFUDIR(1:)='-W '//WFUDIR(1:K)                                     !
      K=K+3                                                             !  
                                                                        !
      CALL CREATE_GEOFILE ( GEOFILE )                                   !
      I=LFBL(ABINIEXE,LEN(ABINIEXE))                                    !
      J=LFBL(ABINIT_INP,LEN(ABINIT_INP))                                !
      COMMAND=ABINIEXE(1:I)//' -s '//WFUDIR(1:K)//'  '//ABINIT_INP(1:J) !
      write (6,*)Command
      CALL PSYST                                                        !
!     Sometime MCSCF doesn't converge... That's an intrinsic problem    !
!     related to the space defined for MCSCF correlation...             !
!     When this appen I choose to :                                     !
!       (1) -  Print on screen a warning message ...                    !
!       (2) -  Keep the gradient of the previous step hoping in good    !
!              convergence at the next step ...                         !
!       (3) -  After 5 failed attempts kill the job..                   !
      CONFERMA=.FALSE.                                                  !
      ICHANNEL=13                                                       !
      CALL CP_OPEN(ICHANNEL, ABINIT_OUT, 'OLD', 'FORMATTED')            !
                                                                        !
      TARGET='** NO CONVERGENCE **'                                     !
      CALL SEARCH(TARGET,ICHANNEL,.TRUE.,'AHEA',CONFERMA,LINE)          !
      IF (CONFERMA)  THEN                                               !
         IFAIL = IFAIL + 1                                              !
         WRITE(6,'(A)')'>>>WARNING IN MCSCF PROCEDURE<<<'               !
         WRITE(6,'(2A)')' MCSCF DIDN''T CONVERGE... CONTINUING WITH', & !
                       ' OLD GRADIENTS!!!'                              !
         WRITE(6,'(A,I5/)')' FAILURE NUMBER : ',IFAIL                   !
         WRITE(6,'(A)')' DELETING OLD RESTART FILE...'
         COMMAND='rm -f  '//WFUDIR(3:K)//'/'//ABINIT_INP(1:J-3)//'res ' !
         write (6,*)Command
         CALL PSYST                                                     !
         IF (IFAIL.GT.5) THEN                                           !
            WRITE(6,'(A)')' ERROR IN MCSCF CALCULATION!!!'              !
            WRITE(6,'(A)')' PROBABILY YOUR ACTIVE SPACE IS ILL-DEFINED!'!
            WRITE(6,'(A,I5,2A)')' THERE HAVE BEEN :',IFAIL,   &         !
                 ' CONSECUTIVE FAILED ATTEMPTS TO ', &                  !
                 ' COMPUTE ELECTRONIC WAVEFUNCTION ...'                 !
            CALL CPSTOP('FORCE_CART')                                   !
         END IF                                                         !
         CALL CP_CLOSE(ICHANNEL, ABINIT_OUT)                            !
         RETURN                                                         !
      END IF                                                            !
!     ------------------------------------------------------------------!
!     At this point gradients and potential have been computed...       !
!     We gather them from ABINIT_OUT file                               !
!     ------------------------------------------------------------------!
      IFAIL = 0                                                         !
      REWIND ICHANNEL                                                   !
                                                                        !
      SELECT CASE (ABINIPROG)                                           !
      CASE ('MOLPRO    ')                                               !
         SELECT CASE (DYNTYPE)                                          !
         CASE ('CONICAL   ')                                            !
            CALL STORE_OLD_VALUE(FSAVE,A,OLDF,OLDA,NUMAT)               !
            CALL GET_ENERGIES_MOLPRO  ( ICHANNEL, VEXC, VFOND )         !
            LABEL='SA-MC DIFF. GRADIENT'                                !
            CALL GET_VECTOR_MOLPRO ( ICHANNEL, F,  NUMAT,  LABEL )      !
            LABEL='SA-MC NACME FOR STAT'                                !
            CALL GET_VECTOR_MOLPRO ( ICHANNEL, A,  NUMAT,  LABEL )      !
            UPOT=VEXC-VFOND                                             !
                                                                        !
            IF (USEGRADEX) THEN                                         !
!     ------------------------------------------------------------------!
!     Use gradient of Excited state                                     !
!     ------------------------------------------------------------------!
               LABEL='GRADIENT FOR STATE 2'                             !
               CALL GET_VECTOR_MOLPRO ( ICHANNEL, F2, NUMAT,  LABEL )   !
!     ------------------------------------------------------------------!
!     Orthonormalize the two vectors G & H and  project out of it  the  !
!     gradient of the excited state                                     !
!     ------------------------------------------------------------------!
               CALL CREATEXCVEC( F2, F, A, OLDF, OLDA, NUMAT, STEP )    ! 
            ENDIF                                                       !
                                                                        !
            IF (USESQUARENE) THEN                                       !
!     ------------------------------------------------------------------!
!     We get as potential the Form (E2-E1)^2 . Modify Gradients and U   !
!     ------------------------------------------------------------------!
               DO I=1,NUMAT                                             !
                  F(I)%X=F(I)%X*UPOT*2.D0                               !
                  F(I)%Y=F(I)%Y*UPOT*2.D0                               !
                  F(I)%Z=F(I)%Z*UPOT*2.D0                               !
               END DO                                                   !
               UPOT=UPOT*UPOT                                           !
            ENDIF                                                       !
                                                                        !
            DO I=1,NUMAT                                                !
               FSAVE(I)%X=F(I)%X                                        !
               FSAVE(I)%Y=F(I)%Y                                        !
               FSAVE(I)%Z=F(I)%Z                                        !
            END DO                                                      !
                                                                        !
            IF (USEGRADEX) THEN                                         !
!     ------------------------------------------------------------------!
!     Sum up the gradients of Excited state projected out of g-h plane  !
!     ------------------------------------------------------------------!
               DO I=1,NUMAT                                             !
                  F(I)%X = ( 1.D0 - LFACT ) *F(I)%X + LFACT * F2(I)%X   !
                  F(I)%Y = ( 1.D0 - LFACT ) *F(I)%Y + LFACT * F2(I)%Y   !
                  F(I)%Z = ( 1.D0 - LFACT ) *F(I)%Z + LFACT * F2(I)%Z   !
               END DO                                                   !
            END IF                                                      !
         CASE ('CONSTRAIN ')                                            !
            CALL STORE_OLD_VALUE(FSAVE,A,OLDF,OLDA,NUMAT)               !
            CALL GET_CHARGES ( ICHANNEL, NUMAT )                        !
            CALL GET_ENERGIES_MOLPRO  ( ICHANNEL, VEXC, VFOND )         !
            LABEL='SA-MC DIFF. GRADIENT'                                !
            CALL GET_VECTOR_MOLPRO ( ICHANNEL, F,  NUMAT,  LABEL )      !
            LABEL='SA-MC NACME FOR STAT'                                !
            CALL GET_VECTOR_MOLPRO ( ICHANNEL, A,  NUMAT,  LABEL )      !
            UPOT=VEXC-VFOND                                             !
                                                                        !
!     ------------------------------------------------------------------!
!     Use gradient of Excited state                                     !
!     ------------------------------------------------------------------!
            LABEL='GRADIENT FOR STATE 2'                                !
            CALL GET_VECTOR_MOLPRO ( ICHANNEL, F2, NUMAT,  LABEL )      !
!     ------------------------------------------------------------------!
!     Adding a weak repulsive force on the bond                         !
!     ------------------------------------------------------------------!
                                                                        !
            IF (ADDREPULSIVE) CALL ADDREPULSIVE_ROUTINE( F, NUMAT)  
            CALL CREATEXCVEC( F2, F, A, OLDF, OLDA, NUMAT, STEP )       !
                                                                        !
            DO I=1,NUMAT                                                !
               FSAVE(I)%X=F(I)%X                                        !
               FSAVE(I)%Y=F(I)%Y                                        !
               FSAVE(I)%Z=F(I)%Z                                        !
            END DO                                                      !
                                                                        !
!     ------------------------------------------------------------------!
!     Sum up the gradients of Excited state projected out of g-h plane  !
!     ------------------------------------------------------------------!
            DO I=1,NUMAT                                                !
               F(I)%X =  F2(I)%X                                        !
               F(I)%Y =  F2(I)%Y                                        !
               F(I)%Z =  F2(I)%Z                                        !
            END DO                                                      !
         CASE ('ONESTATE  ')                                            !
            CALL GET_ENERGY_MOLPRO  ( ICHANNEL, VEXC, VFOND )           !
            IF (USEGRADEX) THEN                                         !
!     ------------------------------------------------------------------!
!     Use gradient of Excited state                                     !
!     ------------------------------------------------------------------!
               LABEL='GRADIENT FOR STATE 2'                             !
               CALL GET_VECTOR_MOLPRO ( ICHANNEL, F , NUMAT,  LABEL )   !
               UPOT = VEXC                                              !
            ELSE                                                        !
!     ------------------------------------------------------------------!
!     Use gradient of Foundamental state                                !
!     ------------------------------------------------------------------!
               LABEL='GRADIENT FOR STATE 1'                             !
               CALL GET_VECTOR_MOLPRO ( ICHANNEL, F , NUMAT,  LABEL )   !
               UPOT = VFOND                                             !
            ENDIF                                                       !
            DO  I = 1, NUMAT                                            !
               FX((I-1)*3+1) = F(I)%X                                   !
               FX((I-1)*3+2) = F(I)%Y                                   !
               FX((I-1)*3+3) = F(I)%Z                                   !
            END DO                                                      !
         CASE DEFAULT                                                   !
            WRITE(6,'(3A)')'CONICAL, CONSTRAIN AND  ONESTATE TYPE OF', &!
                          ' DYNAMICS ARE',                             &!
                          ' POSSIBLE IN THIS VERSION OF CONPATH!'       !
            CALL CPSTOP('FORCE_CART')                                   !
         END SELECT                                                     !
      CASE DEFAULT                                                      !
         WRITE(6,'(A)')'ABINITIO PROGRAM NOT IMPLEMENTED ON CONPATH!'   !
         CALL CPSTOP('FORCE_CART')                                      !
      END SELECT                                                        !
                                                                        !
      CALL CP_CLOSE(ICHANNEL, ABINIT_OUT)                               !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE FORCE_CART                                         !
!     ===================================================================
      SUBROUTINE CREATEXCVEC( F2, F, A, OLDF, OLDA,  NUMAT, STEP ) 
!     ===================================================================
      USE VECTORALGEBRA, ONLY    :  GRAMSCHMIDT,        &               !
                                    PROJECTOUT                          !
      USE XYZ,           ONLY    :  POINT                               !
      USE DYNPREPARE,    ONLY    :  MAT => GHMATRIX,    &               !
                                    MATOLD => GHMATRIXOLD               !
      USE UNITS_FILE,    ONLY    :  GHPLANE                             !
      USE START_JOB,     ONLY    :  DYNTYPE                             !
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      TYPE(POINT), INTENT(IN) ::  F(*), A(*), OLDF(*), OLDA(*)          ! 
      TYPE(POINT), INTENT(INOUT) :: F2(*)                               ! 
      INTEGER (KIND=8), INTENT(IN) :: NUMAT                             ! 
      INTEGER, INTENT(IN) :: STEP                                       ! 
      DOUBLE PRECISION, ALLOCATABLE :: V1(:)                            !
      INTEGER  :: I, VECDIM, J                                          ! 
      DOUBLE PRECISION :: SUP1, SUP2, SUP12                             !
                                                                        !
!     ------------------------------------------------------------------!
!     Allocate Vectors to use them instead of Point type                !
!     ------------------------------------------------------------------!
      VECDIM = NUMAT * 3                                                !
      ALLOCATE   (  V1(VECDIM)  )                                       !
                                                                        !
      DO I=1,NUMAT                                                      !
         DO J=1,2                                                       !
            MATOLD(J,(I-1)*3+1)=MAT((I-1)*3+1,J)                        !
            MATOLD(J,(I-1)*3+2)=MAT((I-1)*3+2,J)                        !
            MATOLD(J,(I-1)*3+3)=MAT((I-1)*3+3,J)                        !
         END DO                                                         !
         MAT((I-1)*3+1,1) = F(I)%X                                      !
         MAT((I-1)*3+2,1) = F(I)%Y                                      !
         MAT((I-1)*3+3,1) = F(I)%Z                                      !
         MAT((I-1)*3+1,2) = A(I)%X                                      !
         MAT((I-1)*3+2,2) = A(I)%Y                                      !
         MAT((I-1)*3+3,2) = A(I)%Z                                      !
         V1((I-1)*3+1) =  F2(I)%X                                       !
         V1((I-1)*3+2) =  F2(I)%Y                                       !
         V1((I-1)*3+3) =  F2(I)%Z                                       !
      END DO                                                            !
                                                                        !
      CALL GRAMSCHMIDT( MAT, VECDIM, 2  )                               !
!     ------------------------------------------------------------------!
!     Writes information on gh-plane on file <nomefile>.ghp             !
!     ------------------------------------------------------------------!
      WRITE(GHPLANE,'(A,I6,A)')'STEP NUMBER =',STEP                     !
      WRITE(GHPLANE,'(9X,A,14X,A,13X,A)')'G','H','DE2'                  !
      WRITE(GHPLANE,'(3F15.9)')(MAT(I,1),MAT(I,2),V1(I),I=1,VECDIM)     !
      WRITE(GHPLANE,'(A,2F15.9)')'PROJECTION OF DE2 ON TWO VECTORS:',&  !
           MATMUL(V1,MAT)                                               !
!     ------------------------------------------------------------------!
      CALL PROJECTOUT ( V1, MAT, VECDIM )                               !
!     ------------------------------------------------------------------!
      WRITE(GHPLANE,'(A,2F15.9)')'RESIDUAL PROJECTION OF DE2 :',&       !
           MATMUL(V1,MAT)                                               !
!     ------------------------------------------------------------------!
!     Compute the superposition of vectors G, H of time step i and  i-1 !
!     ------------------------------------------------------------------!
      SUP1=0.D0                                                         !
      SUP2=0.D0                                                         !
      SUP12=0.D0                                                        !
      DO I=1,NUMAT                                                      !
         SUP1  = SUP1  + OLDF(I)%X*F(I)%X     &                         !
                       + OLDF(I)%Y*F(I)%Y     &                         !
                       + OLDF(I)%Z*F(I)%Z                               !
         SUP2  = SUP2  + OLDA(I)%X*A(I)%X     &                         !
                       + OLDA(I)%Y*A(I)%Y     &                         !
                       + OLDA(I)%Z*A(I)%Z                               ! 
         SUP12 = SUP12 + F(I)%X*A(I)%X        &                         !
                       + F(I)%Y*A(I)%Y        &                         !
                       + F(I)%Z*A(I)%Z                                  !    
      END DO                                                            !
      WRITE(GHPLANE,'(A,F15.9)')'PROJECTION OF G_{I} ON G_{I-1} =',SUP1 !
      WRITE(GHPLANE,'(A,F15.9)')'PROJECTION OF H_{I} ON H_{I-1} =',SUP2 !
      WRITE(GHPLANE,'(A,F15.9)')'PROJECTION OF G_{I} ON H_{I}   =',SUP12!
      WRITE(GHPLANE,'(A/)')'SUPERPOSITION OF GH NORMALIZED MATRIX'      !
      WRITE(GHPLANE,'(/22X,A,14X,A)')'G','H'                            !
      WRITE(GHPLANE,'(4X,''G'',5X,2F15.6,/4X,''H'',5X,2F15.6)') &       !
                                               MATMUL(MATOLD,MAT)       !
!     ------------------------------------------------------------------!
      WRITE(GHPLANE,'(A)')'---'                                         !
!     ------------------------------------------------------------------!
      IF (INDEX(DYNTYPE,'CONICAL').NE.0)  THEN                          !
         DO I=1,NUMAT                                                   !
            F2(I)%X = V1((I-1)*3+1)                                     !
            F2(I)%Y = V1((I-1)*3+2)                                     !
            F2(I)%Z = V1((I-1)*3+3)                                     !
         END DO                                                         !
      END IF                                                            !
                                                                        !
      DEALLOCATE ( V1 )                                                 !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE  CREATEXCVEC                                       !
!     ===================================================================
      SUBROUTINE STORE_OLD_VALUE(F,A,OLDF,OLDA,NUMAT)
!     ===================================================================
      USE XYZ,   ONLY:  POINT                                           !
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      TYPE(POINT), INTENT(IN)    :: F(*),A(*)                           !
      TYPE(POINT), INTENT(INOUT) :: OLDF(*),OLDA(*)                     !
      INTEGER (KIND=8), INTENT(IN) :: NUMAT                             !
      INTEGER :: I                                                      !
                                                                        !
      DO I=1,NUMAT                                                      !
         OLDF(I)%X=F(I)%X                                               !
         OLDF(I)%Y=F(I)%Y                                               !
         OLDF(I)%Z=F(I)%Z                                               !
                                                                        !
         OLDA(I)%X=A(I)%X                                               !
         OLDA(I)%Y=A(I)%Y                                               !
         OLDA(I)%Z=A(I)%Z                                               !
      END DO                                                            !
      RETURN                                                            !
!     ------------------------------------------------------------------!
      END SUBROUTINE STORE_OLD_VALUE                                    !
!     ===================================================================
      
!     ===================================================================
      SUBROUTINE ADDREPULSIVE_ROUTINE(F,NUMAT)
!     ===================================================================
      USE XYZ,   ONLY:  POINT                                           !
      USE START_JOB,   ONLY: NREP1, NREP2, AREP, SIGREP             !
      USE DYNPREPARE,  ONLY: C
!     ------------------------------------------------------------------!
      IMPLICIT NONE                                                     !
      TYPE(POINT), INTENT(INOUT)    :: F(*)                             !
      INTEGER (KIND=8), INTENT(IN) :: NUMAT                             !
!     ------------------------------------------------------------------!
      !Local Variables 
      DOUBLE PRECISION :: XVAL, YVAL, ZVAL, RVAL, MODREP, FMOD
      
      XVAL   = C(NREP1)%X - C(NREP2)%X
      YVAL   = C(NREP1)%Y - C(NREP2)%Y
      ZVAL   = C(NREP1)%Z - C(NREP2)%Z
      RVAL   = SQRT( XVAL * XVAL + YVAL * YVAL + ZVAL * ZVAL ) 
      MODREP = AREP * EXP(-(RVAL/SIGREP)**2)
      write (6,*)'Modrep = ',modrep

      FMOD   = - ( 2.D0 / SIGREP**2 ) * MODREP 
      ! Add Forces to F
      ! NREP1
      write (6,*)'ADDREPULSIVE: ',NREP1,NREP2,AREP,SIGREP
      F(NREP1)%X  =  F(NREP1)%X + FMOD *  XVAL
      F(NREP1)%Y  =  F(NREP1)%Y + FMOD *  YVAL
      F(NREP1)%Z  =  F(NREP1)%Z + FMOD *  ZVAL
      ! NREP2
      F(NREP2)%X  =  F(NREP2)%X - FMOD *  XVAL
      F(NREP2)%Y  =  F(NREP2)%Y - FMOD *  YVAL
      F(NREP2)%Z  =  F(NREP2)%Z - FMOD *  ZVAL
      ! Done...
      RETURN
      END SUBROUTINE ADDREPULSIVE_ROUTINE
!     ===================================================================
END MODULE FORCE
