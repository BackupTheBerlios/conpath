
MODULE DYNLIB
      IMPLICIT NONE
      INTEGER (KIND=4) :: DUMMY
!     -------------------------------------------------------------------     
      CONTAINS
!     ===================================================================
      SUBROUTINE BOLTZMANN_DIST ( TEMP, V, N )
!     ===================================================================
      USE XYZ,         ONLY:  POINT,PROD3
      USE DYNPREPARE,  ONLY:  M
      USE CONVFACTORS, ONLY:  TFAC => AUTOKELVIN,    &
                              AMUTOAU
      USE MARSAGLIAS,  ONLY:  GAUSS
!     -------------------------------------------------------------------
!     TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  
!                                                                   
!     THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (UNIT) MASS.
!     THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    
!     MOLECULES, AND NON-LINEAR MOLECULES.                          
!                                                                   
!     ROUTINE REFERENCED:                                           
!                                                                   
!     REAL FUNCTION GAUSS ( DUMMY )                                 
!        RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             
!        DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             
!     -------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(POINT), INTENT(INOUT)   :: V(*)
      DOUBLE PRECISION, INTENT(IN) :: TEMP
      INTEGER (KIND=8), INTENT(IN) :: N
      INTEGER :: I
      DOUBLE PRECISION :: EKIN, TNOW, RTEMP, SUMX, SUMY, SUMZ, SUMM, SD
      
      RTEMP = TEMP / ( TFAC )
      
      DO  I = 1, N
         SD = SQRT ( RTEMP / M(I) )
         V(I)%X = GAUSS (0.D0, SD)
         V(I)%Y = GAUSS (0.D0, SD)
         V(I)%Z = GAUSS (0.D0, SD)
      END DO
!     -------------------------------------------------------------------
!     Remove Net Momentum
!     -------------------------------------------------------------------
      SUMX = 0.0
      SUMY = 0.0
      SUMZ = 0.0
      SUMM = 0.0
      DO I = 1, N
         SUMM = SUMM + M(I)/ AMUTOAU
         SUMX = SUMX + V(I)%X * M(I)/ AMUTOAU
         SUMY = SUMY + V(I)%Y * M(I)/ AMUTOAU
         SUMZ = SUMZ + V(I)%Z * M(I)/ AMUTOAU
      END DO
      SUMX = SUMX / ( SUMM  )
      SUMY = SUMY / ( SUMM  )
      SUMZ = SUMZ / ( SUMM  )
      EKIN=PROD3 (M,V,V,INT(N,4))
      TNOW= EKIN * TFAC / ( 3.D0 * N )
      WRITE(6,10)('*',I=1,80)
      WRITE(6,11) TNOW
!      DO I = 1, N
!         V(I)%X = V(I)%X - SUMX
!         V(I)%Y = V(I)%Y - SUMY
!         V(I)%Z = V(I)%Z - SUMZ
!      END DO
      EKIN = 0.0
      EKIN=PROD3 (M,V,V,INT(N,4))
      TNOW= EKIN * TFAC / ( 3.D0 * N )
      WRITE(6,12) TNOW
      WRITE(6,10)('*',I=1,80)
      RETURN
10    FORMAT(80A)
11    FORMAT('*',' INITAL TEMPERATURE FROM BOLTZMANN DISTRIBUTION:', &
                 4X,F12.6,14X,'*')
12    FORMAT('*',' INITAL TEMPERATURE AFTER REMOVE OF NET MOMENTUM:', &
                 3X,F12.6,14X,'*')
!     -------------------------------------------------------------------
      END SUBROUTINE BOLTZMANN_DIST
!     ===================================================================
      SUBROUTINE VEL_RESCALE ( TEMP, V, N )
!     ===================================================================
      USE XYZ,         ONLY:  POINT, PROD3
      USE DYNPREPARE,  ONLY:  M
      USE START_JOB,    ONLY:  DELTAT
      USE CONVFACTORS, ONLY:  TFAC => AUTOKELVIN,    &
                              AMUTOAU
!     -------------------------------------------------------------------
!     RESCALE THE VELOCITIES ACCORDING TO TEMP
!                                                                   
!     -------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(POINT), INTENT(INOUT)   :: V(*)
      DOUBLE PRECISION, INTENT(IN) :: TEMP
      INTEGER (KIND=8), INTENT(IN) :: N
      INTEGER :: I
      DOUBLE PRECISION :: TNOW, EKIN, RTEMP, FACT
      DOUBLE PRECISION :: SUMM, SUMX, SUMY, SUMZ
!     -------------------------------------------------------------------
!     DO NOT Remove Net Momentum
!     -------------------------------------------------------------------
 
!      SUMX = 0.0
!      SUMY = 0.0
!      SUMZ = 0.0
!      SUMM = 0.0
!      DO I = 1, N
!         SUMM = SUMM + M(I)/ AMUTOAU
!         SUMX = SUMX + V(I)%X * M(I)/ AMUTOAU
!         SUMY = SUMY + V(I)%Y * M(I)/ AMUTOAU
!         SUMZ = SUMZ + V(I)%Z * M(I)/ AMUTOAU
!      END DO
!      SUMX = SUMX / (REAL ( N ) * SUMM  )
!      SUMY = SUMY / (REAL ( N ) * SUMM  )
!      SUMZ = SUMZ / (REAL ( N ) * SUMM  )
!      DO I = 1, N
!         V(I)%X = V(I)%X - SUMX
!         V(I)%Y = V(I)%Y - SUMY
!         V(I)%Z = V(I)%Z - SUMZ
!      END DO
      EKIN = 0.0
      EKIN=PROD3 (M,V,V,INT(N,4))
      TNOW=EKIN* TFAC/ ( 3.D0 * N )
      WRITE(6,10)('*',I=1,80)
      WRITE(6,11) TNOW
      write (6,*)'I want ',temp,deltat
!     -------------------------------------------------------------------
!     Calculate kinetic energy
!     -------------------------------------------------------------------
 

      IF (DABS(TEMP-TNOW).GT.DELTAT)THEN
         FACT=DSQRT(TEMP/TNOW)
         DO I = 1, N
            V(I)%X =  FACT * V(I)%X
            V(I)%Y =  FACT * V(I)%Y
            V(I)%Z =  FACT * V(I)%Z
         END DO

         EKIN = 0.0
         EKIN=PROD3 (M,V,V,INT(N,4))
         TNOW=EKIN*TFAC/ ( 3.D0 * N )
         WRITE(6,12) TNOW
         WRITE(6,10)('*',I=1,80)
      ENDIF
!     
      RETURN
10    FORMAT(80A)
11    FORMAT('*',' TEMPERATURE BEFORE VELOCITIES RESCALING:', &
                 4X,F12.6,21X,'*')
12    FORMAT('*',' TEMPERATURE AFTER  VELOCITIES RESCALING:', &
                 4X,F12.6,21X,'*')
!     -------------------------------------------------------------------
      END SUBROUTINE VEL_RESCALE
!
!     ===================================================================

END MODULE DYNLIB

