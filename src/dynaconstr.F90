!------------------------------------------------------------------------
MODULE DYNACONSTRAINT
      USE XYZ,         ONLY: POINT,    &
                             COPY_VEC
      IMPLICIT NONE
      TYPE(POINT), ALLOCATABLE :: V1(:), V2(:), C0(:)
!     -------------------------------------------------------------------
      CONTAINS
!     ===================================================================
      SUBROUTINE DYNAMICS_CONSTR
!     ===================================================================
      USE START_JOB,     ONLY: N => NUMAT,      &                        
                               ABINIT_INP,      &                        
                               ENECONV,         &                        
                               DTORIG => TSTEP, &
                               NSTEP,           &
                               IPRINT,          &
                               TITLE,           &
                               LFACT,           &
                               INI_TEMP,        &
                               UPV12,           &
                               RANVEC,          &
                               AREWERESTARTING, &
                               NOROT,           &
                               KUTTEH,          &
                               V1V2UPDATE
      USE DYNPREPARE,    ONLY: C,               &                        
                               F,               &      
                               FSAVE,           &
                               F2,              &
                               M,               &                        
                               V,               &                        
                               A,               &                        
                               LAB
      USE VERLET_INT,    ONLY: MOVEA,           &
                               MOVEB
      USE FORCE,         ONLY: FORCE_CART
      USE MODKINETIC,    ONLY: KINETIC
      USE CONVFACTORS,   ONLY: TCONV => FSTOAUTIME
      USE DYNLIB,        ONLY: BOLTZMANN_DIST
      USE MARSAGLIAS,    ONLY: AMRSET, AMRAND
      USE UNITS_FILE,    ONLY: FLUSH_FILES
      USE RESTART,       ONLY: WRITE_RESTART,   &
                               RSTEP
      USE DYNACART,      ONLY: PRINTGEOVEL
      USE VERLET_CONSTR, ONLY: CONSTRA,       &
                               CONSTRB
      USE CONSTRAINTS,   ONLY: MAKEROT,       &
                               CONSTRUN,      &
                               MCONS,         &
                               MCONSTR

      USE CONSTRAINTS_KUTTEH, ONLY: INIT_VECT,  &
                                    R0V0INIT,   &
                                    CONSTKUTTEH,&
                                    MCONS_KUTTEH => MCONS
      
      IMPLICIT NONE
!     -------------------------------------------------------------------
!     FORTRAN PROGRAM TO CONDUCT MOLECULAR DYNAMICS OF ATOMS        
!     IN CARTESIAN COORDINATES WITH CONSTRAINT ON (E2-E1) VALUE 
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
      INTEGER            :: NORM, STEP, I
      INTEGER, PARAMETER :: FREE=3
      DOUBLE PRECISION   :: K, E, UPOT, VEXC, VFOND
      DOUBLE PRECISION   :: VN, KN, EN, TEMP
      DOUBLE PRECISION   :: ACV, ACK, ACE, ACT
      DOUBLE PRECISION   :: AVV, AVK, AVE, AVT
      DOUBLE PRECISION   :: ACVSQ, ACKSQ, ACESQ, ACTSQ
      DOUBLE PRECISION   :: FLV, FLK, FLE, FLT, DT
!     -------------------------------------------------------------------
!     Convert Time to Atomic Units
!     -------------------------------------------------------------------
      DT=DTORIG
      DT=DT*TCONV
!     -------------------------------------------------------------------
!     Zero Accumulators
!     -------------------------------------------------------------------
      ACV  = 0.D0
      ACK  = 0.D0
      ACE  = 0.D0
      ACT  = 0.D0
      ACVSQ  = 0.D0
      ACKSQ  = 0.D0
      ACESQ  = 0.D0
      ACTSQ  = 0.D0
      FLV  = 0.D0
      FLK  = 0.D0
      FLE  = 0.D0
      FLT  = 0.D0
      STEP = 0
      K = 0.D0
!     -------------------------------------------------------------------
!     Allocate Vectors V1 and V2 ...
!     -------------------------------------------------------------------
      ALLOCATE  (  V1(N) , &
                   V2(N)    )
      ALLOCATE  (  C0(N)    )
!     -------------------------------------------------------------------
      IF ( IPRINT .LE. 0 ) IPRINT = NSTEP + 1
!     -------------------------------------------------------------------
!     Print some fancy stuff on stdout
!     -------------------------------------------------------------------
      WRITE(6,'(27X,'' *** CONSTRAINED  DYNAMICS *** '',//)')
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
         WRITE(6,'('' COMPUTING GRADIENTS IN THE STARTING POINT '')')
         CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP )         
         CALL COPY_VEC ( N , FSAVE , V1 )
         CALL COPY_VEC ( N , A , V2 )
!     -------------------------------------------------------------------
!     Create a set of vector to describe rotations and traslations
!     -------------------------------------------------------------------
         IF ( NOROT  )  CALL MAKEROT ( V1, V2, C, N, M )
         IF ( KUTTEH )  CALL INIT_VECT( V1, V2, N) 
      ELSE
         RSTEP=0
      END IF
!     -------------------------------------------------------------------
!     Storing the original positions into C0 vector
!     -------------------------------------------------------------------
      CALL COPY_VEC ( N , C , C0 )
      CALL FLUSH_FILES
!     -------------------------------------------------------------------
!     Restarting way...
!     -------------------------------------------------------------------
      IF (AREWERESTARTING) GO TO 11
!     -------------------------------------------------------------------
!     Calculate initial values if NORESTART available!
!     -------------------------------------------------------------------
      CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP )
      CALL COPY_VEC ( N , FSAVE , V1 )
      CALL COPY_VEC ( N , A , V2 )
      CALL PRINTGEOVEL( V1,V2,C,V,F,FSAVE,F2,A,N,UPOT,K,STEP,DTORIG, &
                        VFOND,VEXC )
!     -------------------------------------------------------------------
!     Create a set of vector to describe rotations and traslations
!     -------------------------------------------------------------------
      IF ( NOROT )  CALL MAKEROT ( V1, V2, C, N, M )
      IF ( KUTTEH)  CALL INIT_VECT( V1, V2, N ) 
!     -------------------------------------------------------------------
!     Set Nuclei velocities to the Boltzmann distribution of INI_TEMP
!     -------------------------------------------------------------------
      CALL BOLTZMANN_DIST(INI_TEMP, V, N)
      CALL KINETIC (  M, V, K, N )
!     -------------------------------------------------------------------
      IF ( KUTTEH ) CALL R0V0INIT ( C, V )
      CALL MOVEA   (  DT,  M )
!     -------------------------------------------------------------------
!     SET CONSTRAINTS - PART A - CONSTRAINTS ON COORDINATES
!     -------------------------------------------------------------------
      IF (.NOT.KUTTEH) THEN
         IF ( NOROT ) THEN
            CALL CONSTRUN( DT, MCONS, MCONSTR, C, V, C0)
         ELSE
            CALL CONSTRA ( DT, M, V1, V2, C0 )
         END IF
         CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP ) 
      END IF
!     -------------------------------------------------------------------
      CALL MOVEB   (  DT,  M )
!     -------------------------------------------------------------------
!     SET CONSTRAINTS - PART B - CONSTRAINTS ON VELOCITIES
!     -------------------------------------------------------------------
      IF (.NOT.KUTTEH) THEN
         IF ( NOROT ) THEN
            CALL CONSTRUN( DT, MCONS, MCONSTR, C, V )
         ELSE
            CALL CONSTRB (  DT,  M, V1, V2, C0 )      
         END IF
      ELSE
!     -------------------------------------------------------------------
!     Randomize iniziatization
!     -------------------------------------------------------------------
         CALL AMRSET (1230)
         CALL CONSTKUTTEH ( DT, MCONS_KUTTEH,  C, V, C0 )
         CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP ) 
      END IF
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
      CALL WRITE_RESTART ( STEP, 0 )
!     -------------------------------------------------------------------
!     Optionally print information
!     -------------------------------------------------------------------
      IF ( MOD( STEP, INT(IPRINT) ) .EQ. 0 ) THEN
         CALL PRINTGEOVEL( V1,V2,C,V,F,FSAVE,F2,A,N,UPOT,K,STEP,DTORIG, & 
              VFOND,VEXC )
      ENDIF
!     -------------------------------------------------------------------
!     Implementation Algorithm
!     -------------------------------------------------------------------
      IF ( KUTTEH ) CALL R0V0INIT ( C, V )
      CALL MOVEA ( DT, M )
!     -------------------------------------------------------------------
!     SET CONSTRAINTS - PART A - CONSTRAINTS ON COORDINATES
!     -------------------------------------------------------------------
      IF (.NOT.KUTTEH) THEN
         IF ( NOROT ) THEN
            CALL CONSTRUN( DT, MCONS, MCONSTR, C, V, C0)
         ELSE
            CALL CONSTRA ( DT, M, V1, V2, C0 )
         END IF
         CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP )
      END IF
!     -------------------------------------------------------------------
      CALL MOVEB ( DT, M )
!     -------------------------------------------------------------------
!     SET CONSTRAINTS - PART B - CONSTRAINTS ON VELOCITIES
!     -------------------------------------------------------------------
      IF (KUTTEH) THEN
         DO I=1,N
            MCONS_KUTTEH%X=MCONS_KUTTEH%X+AMRAND()*RANVEC
            MCONS_KUTTEH%Y=MCONS_KUTTEH%Y+AMRAND()*RANVEC
            MCONS_KUTTEH%Z=MCONS_KUTTEH%Z+AMRAND()*RANVEC
         ENDDO
         CALL CONSTKUTTEH ( DT, MCONS_KUTTEH,  C, V, C0 )
         CALL FORCE_CART ( VEXC , VFOND,  UPOT, STEP ) 
      ELSE
         IF ( NOROT ) THEN
            CALL CONSTRUN( DT, MCONS, MCONSTR, C, V )
         ELSE
            CALL CONSTRB (  DT,  M, V1, V2, C0 )      
         END IF
      END IF


      IF (V1V2UPDATE.AND.MOD(STEP,UPV12).EQ.0)THEN
         
         WRITE (*,*)'Updating V1 and V2'
         CALL COPY_VEC ( N , FSAVE , V1 )
         CALL COPY_VEC ( N , A , V2 )
         IF ( NOROT )  CALL MAKEROT ( V1, V2, C, N, M )
         IF ( KUTTEH)  CALL INIT_VECT( V1, V2, N ) 
         CALL COPY_VEC ( N , C , C0 )
      ENDIF



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
      DEALLOCATE ( V1, V2 )
      RETURN
!     -------------------------------------------------------------------
      END SUBROUTINE DYNAMICS_CONSTR                      
!     ===================================================================

END MODULE DYNACONSTRAINT
