MODULE VERLET_INT
      IMPLICIT NONE
!     -------------------------------------------------------------------
!                   VELOCITY VERSION OF VERLET ALGORITHM              
!     -------------------------------------------------------------------
!     TWO ROUTINES THAT TOGETHER IMPLEMENT VELOCITY VERLET METHOD. 
!                                                                  
!     REFERENCE:                                                   
!                                                                  
!     SWOPE ET AL., J. CHEM. PHYS. 76, 637, 1982.                  
!                                                                  
!     ROUTINES SUPPLIED:                                           
!                                                                  
!     SUBROUTINE MOVEA ( DT, M )                                   
!        MOVES POSITIONS AND PARTIALLY UPDATES VELOCITIES.         
!     SUBROUTINE MOVEB ( DT, M, K )                                
!        COMPLETES VELOCITY MOVE AND CALCULATES KINETIC ENERGY.    
!                                                                  
!     PRINCIPAL VARIABLES:                                         
!                                                                  
!     INTEGER N                   NUMBER OF MOLECULES              
!     REAL    DT                  TIMESTEP                         
!     REAL    M                   ATOMIC MASS                      
!     REAL    RX(N),RY(N),RZ(N)   POSITIONS                        
!     REAL    VX(N),VY(N),VZ(N)   VELOCITIES                       
!     REAL    FX(N),FY(N),FZ(N)   FORCES                           
!                                                                  
!     USAGE:                                                       
!                                                                  
!     AT THE START OF A TIMESTEP, MOVEA IS CALLED TO ADVANCE THE   
!     POSITIONS AND 'HALF-ADVANCE' THE VELOCITIES.  THEN THE FORCE 
!     ROUTINE IS CALLED, AND THIS IS FOLLOWED BY MOVEB WHICH       
!     COMPLETES THE ADVANCEMENT OF VELOCITIES.                     
!     -------------------------------------------------------------------
      CONTAINS
!     ===================================================================
      SUBROUTINE MOVEA ( DT, M )
!     ===================================================================
      USE START_JOB,   ONLY: N => NUMAT
      USE DYNPREPARE,  ONLY: C,               &                        
                             F,               &                        
                             V
      USE XYZ,         ONLY: POINT,           &
                             PROD2,           &
                             PROD3
!     -------------------------------------------------------------------
!     FIRST PART OF VELOCITY VERLET ALGORITHM                      
!                                                                  
!     USAGE:                                                       
!                                                                  
!     THE FIRST PART OF THE ALGORITHM IS A TAYLOR SERIES WHICH     
!     ADVANCES POSITIONS FROM T TO T + DT AND VELOCITIES FROM      
!     T TO T + DT/2.  AFTER THIS, THE FORCE ROUTINE IS CALLED.     
!     -------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)           :: DT
      DOUBLE PRECISION                       :: DT2, DTSQ2
      DOUBLE PRECISION, INTENT(IN), DIMENSION(N)  :: M(*)
      INTEGER                                :: I
!     -------------------------------------------------------------------
      DT2   = DT / 2.0
      DTSQ2 = DT * DT2

      DO  I = 1, N
         C(I)%X =  C(I)%X + DT  * V(I)%X - DTSQ2 *  F(I)%X / M(I)
         C(I)%Y =  C(I)%Y + DT  * V(I)%Y - DTSQ2 *  F(I)%Y / M(I)
         C(I)%Z =  C(I)%Z + DT  * V(I)%Z - DTSQ2 *  F(I)%Z / M(I)
         V(I)%X =  V(I)%X - DT2 * F(I)%X / M(I)
         V(I)%Y =  V(I)%Y - DT2 * F(I)%Y / M(I)
         V(I)%Z =  V(I)%Z - DT2 * F(I)%Z / M(I)
      END DO
      
      RETURN
!     -------------------------------------------------------------------
      END SUBROUTINE MOVEA
!     ===================================================================
      SUBROUTINE MOVEB ( DT, M )
!     ===================================================================
      USE START_JOB,   ONLY: N => NUMAT,      &
                             DAMP
      USE DYNPREPARE,  ONLY: C,               &                        
                             F,               &                        
                             V
      USE XYZ,         ONLY: POINT,           &
                             PROD2,           &
                             PROD3      
!     -------------------------------------------------------------------
!     SECOND PART OF VELOCITY VERLET ALGORITHM                    
!                                                                 
!     USAGE:                                                      
!                                                                 
!     THE SECOND PART OF THE ALGORITHM ADVANCES VELOCITIES FROM   
!     T + DT/2 TO T + DT. THIS ASSUMES THAT FORCES HAVE BEEN      
!     COMPUTED IN THE FORCE ROUTINE AND STORED IN FX, FY, FZ.     
!     -------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)                :: DT
      DOUBLE PRECISION                            :: DT2, DTSQ2
      DOUBLE PRECISION, INTENT(IN), DIMENSION(N)  :: M(*)
      INTEGER                                     :: I
!     -------------------------------------------------------------------
      DT2 = DT / 2.0

      DO I = 1, N
         V(I)%X = ( V(I)%X - DT2 * F(I)%X / M(I) ) * DAMP
         V(I)%Y = ( V(I)%Y - DT2 * F(I)%Y / M(I) ) * DAMP
         V(I)%Z = ( V(I)%Z - DT2 * F(I)%Z / M(I) ) * DAMP
      END DO
      RETURN
!     -------------------------------------------------------------------
      END SUBROUTINE MOVEB
!     ===================================================================

END MODULE VERLET_INT
