!------------------------------------------------------------------------
MODULE DYNACART
      IMPLICIT NONE
!     -------------------------------------------------------------------
      CONTAINS
!     ===================================================================
      SUBROUTINE DYNAMICS_CART
!     ===================================================================
      USE START_JOB,   ONLY: N => NUMAT,      &                        
                             ABINIT_INP,      &                        
                             ENECONV,         &                        
                             DT => TSTEP,     &
                             NSTEP,           &
                             IPRINT,          &
                             TITLE,           &
                             LFACT,           &
                             INI_TEMP,        &
                             ALPHA,           &
                             SIGMA,           &
                             AREWERESTARTING
      USE DYNPREPARE,  ONLY: C,               &                        
                             F,               &      
                             FSAVE,           &
                             F2,              &
                             M,               &                        
                             V,               &                        
                             A,               &                        
                             LAB
      USE VERLET_INT,  ONLY: MOVEA,           &
                             MOVEB
      USE FORCE,       ONLY: FORCE_CART
      USE MODKINETIC,  ONLY: KINETIC
      USE CONVFACTORS, ONLY: TCONV => FSTOAUTIME
      USE DYNLIB,      ONLY: BOLTZMANN_DIST
      USE UNITS_FILE,  ONLY: FLUSH_FILES
      USE RESTART,     ONLY: WRITE_RESTART,   &
                             RSTEP,RN_IC
      USE CONICAL,     ONLY: STORE_IC, V1, V2, CSAVE
      USE CAMPANE,     ONLY: GAUSSIAN, GAUS_FORCES
      IMPLICIT NONE
!     -------------------------------------------------------------------
!     FORTRAN PROGRAM TO CONDUCT MOLECULAR DYNAMICS OF ATOMS        
!     IN CARTESIAN COORDINATES                                      
!                                                                   
!     REFERENCE:                                                    
!                                                                   
!     FINCHAM AND HEYES, CCP5 QUARTERLY, 6, 4, 1982.                
!                                                                   
!     ROUTINES REFERENCED:                                          
!                                                                   
!     SUBROUTINE FORCE ( DT, F, A, STEP )                                 
!        CALCULATES THE ACCELERATIONS                               
!     SUBROUTINE MOVEA ( DT, M )                                    
!        ADVANCES POSITIONS FROM T TO T + DT. AND VELOCITIES FROM 
!        T + DT/2
!     SUBROUTINE MOVEB ( DT, M )
!        ADVANCES VELOCITIES FROM T + DT/2 TO T + DT  
!     SUBROUTINE KINETIC ( M, V, K, N )                                     
!        CALCULATES KINETIC ENERGY.     
!                                                                   
!     PRINCIPAL VARIABLES:                                          
!                                                                   
!     INTEGER N                 NUMBER OF ATOMS                     
!     REAL    DT                TIMESTEP                            
!     REAL    C(N)X,C(N)Y,C(N)Z ATOMIC POSITIONS                    
!     REAL    V(N)X,V(N)Y,V(N)Z ATOMIC VELOCITIES                   
!     REAL    F(N)X,F(N)Y,F(N)Z ATOMIC FORCES                       
!                                                                   
!     USAGE:                                                        
!                                                                   
!     THE LEAPFROG ALGORITHM IS OF THE FORM                         
!     VX(T + DT/2) = VX(T - DT/2) + DT * AX(T)  (SIMILARLY FOR Y,Z) 
!     RX(T + DT)   = RX(T) + DT * VX(T + DT/2)  (SIMILARLY FOR Y,Z) 
!     -------------------------------------------------------------------
      INTEGER            :: NORM, STEP, I, N_IC
      INTEGER, PARAMETER :: FREE=3
      LOGICAL            :: FLAG_ADDED
      DOUBLE PRECISION   :: K, E, UPOT, VEXC, VFOND, UGAUS
      DOUBLE PRECISION   :: VN, KN, EN, TEMP
      DOUBLE PRECISION   :: ACV, ACK, ACE, ACT
      DOUBLE PRECISION   :: AVV, AVK, AVE, AVT, LFACT0
      DOUBLE PRECISION   :: ACVSQ, ACKSQ, ACESQ, ACTSQ
      DOUBLE PRECISION   :: FLV, FLK, FLE, FLT, DTORIG
!     -------------------------------------------------------------------
!     Convert Time to Atomic Units
!     -------------------------------------------------------------------
      DTORIG=DT
      DT=DT*TCONV
!     -------------------------------------------------------------------
!     Zero Accumulators
!     -------------------------------------------------------------------
      ACV  = 0.0
      ACK  = 0.0
      ACE  = 0.0
      ACT  = 0.0
      ACVSQ  = 0.0
      ACKSQ  = 0.0
      ACESQ  = 0.0
      ACTSQ  = 0.0
      FLV  = 0.0
      FLK  = 0.0
      FLE  = 0.0
      FLT  = 0.0
      STEP = 0
      IF (.NOT.AREWERESTARTING) N_IC = 0
      FLAG_ADDED=.FALSE.
      IF ( IPRINT .LE. 0 ) IPRINT = NSTEP + 1
!     -------------------------------------------------------------------
!     Print some fancy stuff on stdout
!     -------------------------------------------------------------------
      WRITE(6,'(27X,'' **** PROGRAM FROGGY **** '',//)')
      WRITE(6,'(20X,'' MOLECULAR DYNAMICS     CARTESIAN MODULE   '')')
      WRITE(6,'(20X,'' LEAPFROG ALGORITHM WITH MINIMAL STORAGE   '')')
      
            
      WRITE(6,'(//1X,A)') TITLE
      WRITE(6,'('' NUMBER OF ATOMS  = '',I10  )') N
      WRITE(6,'('' NUMBER OF STEPS  = '',I10  )') NSTEP
      WRITE(6,'('' OUTPUT FREQUENCY = '',I10  )') IPRINT
      WRITE(6,'('' TIME STEP        = '',F10.6,''  fs'')') DTORIG 
      WRITE(6,'('' MIXING GRADIENTS FACTOR = '',F10.6)')LFACT
      WRITE(6,'('' INITIAL NUCLEAR TEMPERATURE = '',F10.6)')INI_TEMP
      IF (AREWERESTARTING) THEN
         WRITE(6,'('' RESTARTING FROM  FILE! '')')
         N_IC=RN_IC
         WRITE(6,'('' COMPUTING GRADIENTS IN THE STARTING POINT '')')
         CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP )         
      ELSE
         RSTEP=0
      END IF
      LFACT0=LFACT
      UGAUS=0.d0
      CALL FLUSH_FILES
!     -------------------------------------------------------------------
!     Restarting way...
!     -------------------------------------------------------------------
      IF (AREWERESTARTING) GOTO 11
      CALL BOLTZMANN_DIST(INI_TEMP,V,N)
      CALL KINETIC (  M, V, K, N )
!     -------------------------------------------------------------------
!     Set Nuclei velocities to the Boltzmann distribution of INI_TEMP
!     -------------------------------------------------------------------
!     -------------------------------------------------------------------
!     Calculate initial values if NORESTART available!
!     -------------------------------------------------------------------
      CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP )
      CALL PRINTGEOVEL( FSAVE,A,C,V,F,FSAVE,F2,A,N,UPOT,K,STEP,DTORIG, &
                        VFOND,VEXC )
      CALL MOVEA (  DT,  M )
      CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP ) 
      CALL MOVEB (  DT,  M )
!     -------------------------------------------------------------------
!     Compute Kinetic Energy
!     -------------------------------------------------------------------
11    CALL KINETIC (  M, V, K, N )  
      WRITE(*,'(//1X,''**** START OF DYNAMICS ****'')')
!     -------------------------------------------------------------------
!     Main Loop begins     
!     -------------------------------------------------------------------
Main: DO  STEP = RSTEP+1, NSTEP
!     -------------------------------------------------------------------
!     Save Restart Information
!     -------------------------------------------------------------------
      CALL WRITE_RESTART ( STEP, N_IC )
!     -------------------------------------------------------------------
!     Checking convergence (E2-E1)
!     -------------------------------------------------------------------

      IF ((VEXC-VFOND).GE.ENECONV)FLAG_ADDED=.FALSE.
      IF ((VEXC-VFOND).LT.ENECONV.AND..NOT.FLAG_ADDED) THEN
         FLAG_ADDED=.TRUE.
         N_IC=N_IC+1
         CALL STORE_IC ( N_IC )
      END IF
!     -------------------------------------------------------------------
!     Optionally print information
!     -------------------------------------------------------------------
         IF ( MOD( STEP, INT(IPRINT) ) .EQ. 0 ) THEN
            CALL PRINTGEOVEL( FSAVE,A,C,V,F,FSAVE,F2,A,N,UPOT,K,STEP,DTORIG, & 
                              VFOND,VEXC )
         ENDIF
!     -------------------------------------------------------------------
!     Implementation Algorithm
!     -------------------------------------------------------------------
         CALL MOVEA ( DT, M )
 
!
!  modify lfact if ugaus is <> 0
!
         IF (UGAUS.GT.0.D0)THEN
            IF (UGAUS.LE.ALPHA)THEN
               LFACT=LFACT0+(1.d0-LFACT0)*(UGAUS/ALPHA)**2
            ELSE
               LFACT=1.d0
            ENDIF
            WRITE (6,*)'lfact, lfact0 = ',lfact,lfact0
         ENDIF
         CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP )
!     -------------------------------------------------------------------
!     Locate gaussian repulsion potential at every conical point found..
!     -------------------------------------------------------------------
         DO I=1,N_IC            
            UGAUS=GAUSSIAN(C,      &
            RESHAPE(CSAVE(1:N,I:I),SHAPE=(/N/)),&     
            RESHAPE(V1(1:N,I:I),SHAPE=(/N/)),&
            RESHAPE(V2(1:N,I:I),SHAPE=(/N/)),ALPHA, SIGMA,N)
            UPOT=UPOT+UGAUS

            WRITE(6,*)'ugaus, i = ',ugaus,i
            CALL GAUS_FORCES(C,   &
            RESHAPE(CSAVE(1:N,I:I),SHAPE=(/N/)),&      
            RESHAPE(V1(1:N,I:I),SHAPE=(/N/)),&
            RESHAPE(V2(1:N,I:I),SHAPE=(/N/)),ALPHA, SIGMA, N, UGAUS, F)
         END DO
         CALL MOVEB ( DT, M )
!     -------------------------------------------------------------------
!     Calculate kinetic energy at current step
!     -------------------------------------------------------------------
         CALL KINETIC ( M, V, K, N )
!     -------------------------------------------------------------------
!     Calculate instantaneous values
!     -------------------------------------------------------------------
         E    = K + UPOT
         VN   = UPOT  / REAL ( N )
         KN   = K  / REAL ( N )
         EN   = E  / REAL ( N )
         TEMP = 2.0 * KN / FREE
!     -------------------------------------------------------------------
!     Increment accumulators
!     -------------------------------------------------------------------
         ACE    = ACE  + EN
         ACK    = ACK  + KN
         ACV    = ACV  + VN
         ACESQ  = ACESQ  + EN  ** 2
         ACKSQ  = ACKSQ  + KN  ** 2
         ACVSQ  = ACVSQ  + VN  ** 2
      END DO Main
!     -------------------------------------------------------------------
!     Main Loop Ends.. Print Some averages..
!     -------------------------------------------------------------------
      WRITE(*,'(/1X,''**** END OF DYNAMICS **** ''//)')
      NORM   = REAL ( STEP - 1 )
      AVE  = ACE  / NORM
      AVK  = ACK  / NORM
      AVV  = ACV  / NORM
      ACESQ  = ( ACESQ  / NORM ) - AVE  ** 2
      ACKSQ  = ( ACKSQ  / NORM ) - AVK  ** 2
      ACVSQ  = ( ACVSQ  / NORM ) - AVV  ** 2
      IF ( ACESQ  .GT. 0.0 ) FLE  = SQRT ( ACESQ  )
      IF ( ACKSQ  .GT. 0.0 ) FLK  = SQRT ( ACKSQ  )
      IF ( ACVSQ  .GT. 0.0 ) FLV  = SQRT ( ACVSQ  )
      AVT = AVK * 2.0 / FREE
      FLT = FLK * 2.0 / FREE
      WRITE(*,'('' AVERAGES'',6(2X,F10.5))')AVE, AVK, AVV, AVT
      WRITE(*,'('' FLUCTS  '',6(2X,F10.5))')FLE, FLK, FLV, FLT
      RETURN
!     -------------------------------------------------------------------
      END SUBROUTINE DYNAMICS_CART                       
!     ===================================================================
      SUBROUTINE PRINTGEOVEL(VV1,VV2,C,V,F,FSAVE,F2,A,N,UPOT,K,STEP,DT, &
                             VFOND,VEXC)
!     ===================================================================
      USE  XYZ,        ONLY:  POINT
      USE  DYNPREPARE, ONLY: LAB
      USE  START_JOB,  ONLY: USESQUARENE,            &
                             USEGRADEX
      USE UNITS_FILE,  ONLY: FLUSH_FILES
      USE CONICAL,     ONLY: STORE_IC, V1, V2, CSAVE
      USE MOLPRO_LINK, ONLY: COMP_CHARGE, CHARGE
!     -------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(POINT), INTENT(IN)       :: C(*), V(*), F(*), A(*),    &
                                       FSAVE(*),F2(*), VV1(*), VV2(*)
      INTEGER (KIND=8), INTENT(IN)  :: N
      INTEGER, INTENT(IN)           :: STEP
      DOUBLE PRECISION, INTENT(IN)  :: UPOT, K
      DOUBLE PRECISION, INTENT(IN)  :: DT,VFOND,VEXC
      INTEGER :: I
!     -------------------------------------------------------------------
!     Initial title
!     -------------------------------------------------------------------
      WRITE(6,100)('*',I=1,80)
      WRITE(6,102)
      WRITE(6,101)STEP,REAL(STEP)*DT
      WRITE(6,102)
      WRITE(6,100)('*',I=1,80)
!     -------------------------------------------------------------------
!     Coordinates and velocities
!     -------------------------------------------------------------------
      WRITE(6,102)
      WRITE(6,103)
      WRITE(6,102)      
      WRITE(6,100)('*',I=1,80)
      WRITE(6,102)
      DO I=1,N
         WRITE(6,104)LAB(I),C(I)%X,C(I)%Y,C(I)%Z,V(I)%X,V(I)%Y,V(I)%Z
      END DO
!     -------------------------------------------------------------------
!     Energies information
!     -------------------------------------------------------------------
      WRITE(6,102)
      WRITE(6,100)('*',I=1,80)
!     -------------------------------------------------------------------
!     If charges are available, let's print them
!     -------------------------------------------------------------------
      IF (COMP_CHARGE) THEN
         WRITE(6,102)
         WRITE(6,112)
         WRITE(6,102)
         WRITE(6,100)('*',I=1,80)
         WRITE(6,102)
         DO I=1,N
            WRITE(6,113)LAB(I),CHARGE(I)
         END DO
         WRITE(6,102)
         WRITE(6,100)('*',I=1,80)         
      END IF
      IF (USESQUARENE) THEN
         WRITE(6,109)UPOT,K
      ELSE
         WRITE(6,105)UPOT,K
      END IF
      WRITE(6,108)UPOT+K,VEXC-VFOND
      WRITE(6,102)
      WRITE(6,107)VFOND,VEXC
!     -------------------------------------------------------------------
!     Gradients information
!     -------------------------------------------------------------------
      WRITE(6,100)('*',I=1,80)
      WRITE(6,102)
      WRITE(6,111)
      WRITE(6,102)
      WRITE(6,100)('*',I=1,80)
      WRITE(6,106)
      WRITE(6,100)('*',I=1,80)
      DO I=1,N
         WRITE(6,104)LAB(I),VV1(I)%X,VV1(I)%Y,VV1(I)%Z,   &
                            VV2(I)%X,    VV2(I)%Y,    VV2(I)%Z
      END DO      
      WRITE(6,100)('*',I=1,80)
      IF (USEGRADEX) THEN 
         WRITE(6,110)
         WRITE(6,100)('*',I=1,80)
         DO I=1,N
            WRITE(6,104)LAB(I),F2(I)%X,F2(I)%Y,F2(I)%Z,   &
                                F(I)%X, F(I)%Y, F(I)%Z
         END DO
         WRITE(6,100)('*',I=1,80)
      END IF

      CALL PUNCHXYZ(N,LAB,C)
      CALL FLUSH_FILES
      RETURN
!     -------------------------------------------------------------------
 100  FORMAT(80A)
 101  FORMAT('*',' STEP NUMBER=',I6,38X,'TIME=',F12.3,' fs',' *')
 102  FORMAT('*',78X,'*')
 103  FORMAT('*',19X,'GEOMETRY',12X,17X,'VELOCITIES',12X,'*')
 104  FORMAT('*',1X,A2,3F12.6,2X,3F12.6,1X,'*')
 105  FORMAT('*',1X,'POTENTIAL ENERGY=',F12.6,' A.U.',  &
                9X,'KINETIC ENERGY=',F13.6,' A.U.',1X,'*')
 106  FORMAT('*',8X,'GRADIENT OF (E2-E1)^2 (A.U.)',2X,9X  &
                ,'NONADIABATIC COUPLING (A.U.)',3X,'*')
 107  FORMAT('*',' FUNDAMENTAL STATE ENERGY=',F12.6,6X,     &
                 'EXCITED STATE ENERGY=',F12.6,1X,'*')
 108  FORMAT('*',' TOTAL ENERGY=',F16.6,1X,'A.U.',        &
                9X,'( E2 - E1 )=',F16.8,' A.U.',1X,'*')
 109  FORMAT('*',1X,'POTENTIAL ENERGY=',F12.6,' A.U.^2',  &
                7X,'KINETIC ENERGY=',F13.6,' A.U.',1X,'*')
 110  FORMAT('*',11X,'EXCITED   STATE  GRADIENT',19X, &
                    'FULL GRADIENT',10X,'*',/, &
             '*',11X,'PROJECTED OUT OF GH-PLANE',42X,'*')
 111  FORMAT('*',31X ,'GRADIENTS SECTION',30X,'*')
      !
 112  FORMAT('*',6X,"CHARGE",66X,'*')
 113  FORMAT('*',1X,A2,F12.6,63X,'*')
!     -------------------------------------------------------------------
      END SUBROUTINE  PRINTGEOVEL                     
!     ===================================================================
      SUBROUTINE PUNCHXYZ(NUMAT, LAB, C)
!     ===================================================================
      USE  UNITS_FILE,   ONLY:      XYZFILE
      USE  XYZ,          ONLY:      POINT
      USE  CONVFACTORS,  ONLY:      CNV => BOHRTOANG 
!     -------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER (KIND=8), INTENT(IN) :: NUMAT
      CHARACTER (LEN=2), INTENT(IN) :: LAB(*)
      TYPE(POINT), INTENT(IN) :: C(*)
      INTEGER :: I


      WRITE(XYZFILE,'(I6/)') NUMAT
      
      DO I=1,NUMAT
         WRITE(XYZFILE,'(A2,5X,3F15.9)')LAB(I),                     &
                                        C(I)%X*CNV,C(I)%Y*CNV,C(I)%Z*CNV
      END DO

      RETURN
!     -------------------------------------------------------------------
      END SUBROUTINE PUNCHXYZ                    
!     ===================================================================

END MODULE DYNACART
