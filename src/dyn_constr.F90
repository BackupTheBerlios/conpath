MODULE CONSTRAINTS

      USE SYSTEM_UTIL, ONLY: CPSTOP
      USE XYZ,         ONLY: POINT
      USE ROTATIONS,   ONLY: ROTATE3D
      IMPLICIT NONE

      INTERFACE
         SUBROUTINE DGESV ( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
           INTEGER, INTENT(IN) :: INFO, LDA, LDB, N, NRHS
           INTEGER, INTENT(OUT):: IPIV( * )
           DOUBLE PRECISION, INTENT(INOUT) ::  A( LDA, *), B( LDB, *)
         END SUBROUTINE DGESV
      END INTERFACE

      INTERFACE OPERATOR (*)
         MODULE PROCEDURE SCALARPOINT
      END INTERFACE
      
      DOUBLE PRECISION, ALLOCATABLE :: VV(:,:), XV(:,:)
      TYPE(POINT), ALLOCATABLE      :: MCONS(:,:)
      INTEGER                       :: MCONSTR
CONTAINS
!     ===================================================================
      SUBROUTINE CONSTRUN ( DT, VCONSTR, NCONSTR, C, V, C0 )
!     ===================================================================
      USE START_JOB,   ONLY: N => NUMAT
      USE DYNPREPARE,  ONLY: MASS => M
      USE CONVFACTORS,     ONLY  :  AMUTOAU
!     -------------------------------------------------------------------
!     IT CORRECTS THE POSITIONS AND VELOCITIES WITH THE CONSTRAINT FORCES
!     FIRST PART
!     -------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)                :: DT
      DOUBLE PRECISION                            :: DT2, DTSQ2
      TYPE(POINT), INTENT(IN), OPTIONAL           :: C0(*)
      TYPE(POINT), INTENT(INOUT)                  :: C(*), V(*)
      TYPE(POINT), INTENT(IN)                     :: VCONSTR(N,NCONSTR)
      INTEGER, INTENT(IN)                         :: NCONSTR
!     -------------------------------------------------------------------
      INTEGER, ALLOCATABLE :: IPIV(:)
      TYPE(POINT), ALLOCATABLE :: CLANG(:), CZERO(:), W1(:), W2(:)
      INTEGER :: INFO, I, J, K
      DOUBLE PRECISION :: RESIDUE
!     -------------------------------------------------------------------
      IF (ALLOCATED ( VV )) DEALLOCATE ( VV )
      IF (ALLOCATED ( XV )) DEALLOCATE ( XV )

      ALLOCATE ( VV( NCONSTR, NCONSTR),   &
                 XV( NCONSTR,       1),   &
                         IPIV(NCONSTR)       )

      ALLOCATE ( CZERO(N),   &
                    W1(N),   &
                    W2(N)       )

      IF (PRESENT(C0)) THEN
         DO I=1,N
            CZERO(I)%X = ( C(I)%X - C0(I)%X ) 
            CZERO(I)%Y = ( C(I)%Y - C0(I)%Y ) 
            CZERO(I)%Z = ( C(I)%Z - C0(I)%Z ) 
         END DO
      ELSE
          DO I=1,N
            CZERO(I)%X = V(I)%X  
            CZERO(I)%Y = V(I)%Y  
            CZERO(I)%Z = V(I)%Z  
         END DO        
      END IF
!     -------------------------------------------------------------------
!     TO FIND LAGRANGE MULTIPLIER I HAVE TO SOLVE THE LINEAR
!     SYSTEM :
!                   VV (:,:) * LAMBDA(:) = XV(:,:)
!     WHERE :
!                - VV(I,J)   = <V_I | V_J>
!                - LAMBDA(I)  = Ith LAGRANGE MULTIPLIER
!                - XV(I,1)     = < Xnow | V_I >
!     -------------------------------------------------------------------
      DO I=1,NCONSTR
         DO K=1,N
            W1(K)%X = VCONSTR(K,I)%X
            W1(K)%Y = VCONSTR(K,I)%Y
            W1(K)%Z = VCONSTR(K,I)%Z
         END DO
         DO J=I,NCONSTR
            DO K=1,N
               W2(K)%X = VCONSTR(K,J)%X / MASS(K) 
               W2(K)%Y = VCONSTR(K,J)%Y / MASS(K) 
               W2(K)%Z = VCONSTR(K,J)%Z / MASS(K) 
            END DO
!     -------------------------------------------------------------------
!     BUILDS THE MATRIX VV - MATRIX OF LINEAR SYSTEM 
!     -------------------------------------------------------------------
            VV(I,J) = W1 * W2
            VV(J,I) = VV(I,J)
         END DO
!     -------------------------------------------------------------------
!     BUILDS THE VECTOR XV
!     -------------------------------------------------------------------
         XV(I,1) = W1 * CZERO
      END DO
      IPIV = 0
!     -------------------------------------------------------------------
!     LAPACK ROUTINE TO SOLVE THE LINEAR SYSTEM
!     -------------------------------------------------------------------
      CALL DGESV ( NCONSTR, 1 , VV, NCONSTR, IPIV, XV, NCONSTR, INFO )
!     -------------------------------------------------------------------
!     ANALYZE RESULTS FROM LAPACK LINEAR SOLVER ROUTINE
!     -------------------------------------------------------------------
      IF (INFO.NE.0) THEN
         WRITE(6,'(A,I5)')' LINEAR SOLVER LAPACK3.0 ROUTINE.INFO VALUE =',&
              INFO
         WRITE(6,'(A)')'IF INFO :',                                       &
         '< 0: IF INFO = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE    ', &
         '> 0: IF INFO = K, U(K,K) IS EXACTLY ZERO.  THE FACTORIZATION ', &
         '     HAS BEEN COMPLETED, BUT THE FACTOR U IS EXACTLY         ', &
         '     SINGULAR, SO THE SOLUTION COULD NOT BE COMPUTED.        '
      END IF
      DEALLOCATE ( IPIV )
!     -------------------------------------------------------------------
!     WE HAVE THE SOLUTION OF THE LINEAR SYSTEM .. NOW WE CAN CORRECT THE
!     CONSTRAINTS
!     -------------------------------------------------------------------      
      ALLOCATE ( CLANG(N) )

      DO I=1, N
         CLANG(I)%X = 0.D0
         CLANG(I)%Y = 0.D0
         CLANG(I)%Z = 0.D0
      END DO

      DO J=1,NCONSTR
         DO I=1,N
            CLANG(I)%X = CLANG(I)%X + XV(J,1) * VCONSTR(I,J)%X / MASS(I)
            CLANG(I)%Y = CLANG(I)%Y + XV(J,1) * VCONSTR(I,J)%Y / MASS(I)             
            CLANG(I)%Z = CLANG(I)%Z + XV(J,1) * VCONSTR(I,J)%Z / MASS(I) 
         END DO
      END DO

      DO I=1,N
         IF (PRESENT(C0)) THEN
            C(I)%X =  C(I)%X - CLANG(I)%X
            C(I)%Y =  C(I)%Y - CLANG(I)%Y
            C(I)%Z =  C(I)%Z - CLANG(I)%Z
         END IF
         V(I)%X =  V(I)%X - CLANG(I)%X
         V(I)%Y =  V(I)%Y - CLANG(I)%Y
         V(I)%Z =  V(I)%Z - CLANG(I)%Z
      END DO

      DEALLOCATE (CLANG)
      DEALLOCATE (CZERO)
!     -------------------------------------------------------------------
!     VERIFY THE CONSTRAINT EFFICACY
!     -------------------------------------------------------------------
      DO I=1,N
         W1(I)%X = 0.D0
         W1(I)%Y = 0.D0
         W1(I)%Z = 0.D0
      END DO
      write(*,*)'teo ================>',NCONSTR
      DO I=1,NCONSTR
         DO J=1,N
            W1(J)%X = W1(J)%X +VCONSTR(J,I)%X 
            W1(J)%Y = W1(J)%Y +VCONSTR(J,I)%Y 
            W1(J)%Z = W1(J)%Z +VCONSTR(J,I)%Z            
         END DO
      END DO
      RESIDUE = 0.D0
      IF (PRESENT(C0)) THEN
         DO I=1,N
            RESIDUE = RESIDUE +  ( C(I)%X - C0(I)%X ) * W1(I)%X
            RESIDUE = RESIDUE +  ( C(I)%Y - C0(I)%Y ) * W1(I)%Y
            RESIDUE = RESIDUE +  ( C(I)%Z - C0(I)%Z ) * W1(I)%Z
         END DO
         WRITE(6,'(A,F15.9)')'CONSTRAINED VALUE ON COORDINATES= ',RESIDUE
      ELSE
          DO I=1,N
            RESIDUE = RESIDUE + V(I)%X * W1(I)%X
            RESIDUE = RESIDUE + V(I)%Y * W1(I)%Y
            RESIDUE = RESIDUE + V(I)%Z * W1(I)%Z
         END DO
         WRITE(6,'(A,F15.9)')'CONSTRAINED VALUE ON  VELOCITIES= ',RESIDUE
      END IF

      DEALLOCATE ( W1, W2)
      RETURN
      END SUBROUTINE CONSTRUN


      DOUBLE PRECISION FUNCTION SCALARPOINT ( V1, V2 )
      USE START_JOB,   ONLY: N => NUMAT

      IMPLICIT NONE
      TYPE(POINT), INTENT(IN) :: V1(N), V2(N)
      INTEGER :: I, DIM1
      DOUBLE PRECISION :: PARTIAL

      DIM1 = SIZE(V1)
      IF (DIM1.NE.SIZE(V2)) THEN
         WRITE(6,'(A)')' ERROR IN VECTORS DIMENSION!'
         CALL CPSTOP('SCALARPOINT')
      END IF

      PARTIAL = 0.D0
      DO I=1,DIM1
         PARTIAL = PARTIAL + V1(I)%X * V2(I)%X
         PARTIAL = PARTIAL + V1(I)%Y * V2(I)%Y
         PARTIAL = PARTIAL + V1(I)%Z * V2(I)%Z
      END DO
      
      SCALARPOINT = PARTIAL
      RETURN
      END FUNCTION SCALARPOINT

      
      SUBROUTINE MAKEROT ( V1, V2, C, N, M )
        USE ROTATIONS,   ONLY:  ROTATE3D
        USE DYNPREPARE,  ONLY: MASS => M
        IMPLICIT NONE
        TYPE(POINT), INTENT(IN) :: V1(*), V2(*), C(*)
        INTEGER (KIND=8), INTENT(IN) :: N
        DOUBLE PRECISION, INTENT(IN) :: M(*)
        TYPE(POINT) :: WRK, VZERO
        INTEGER :: I, J
        DOUBLE PRECISION :: DELTA,  PI2, MTOT, NORM

        DELTA = 0.001
        PI2 = ATAN(1.D0) * 2.D0
        MCONSTR = 2  + 3 + 3   ! V1, V2 - TRASLATIONS - ROTATIONS
        IF (.NOT.ALLOCATED (MCONS) ) ALLOCATE ( MCONS(N,MCONSTR) )
!     -------------------------------------------------------------------
!     COMPUTE CENTRE OF MASS OF THE SYSTEM
!     -------------------------------------------------------------------
        VZERO%X = 0.D0
        VZERO%Y = 0.D0
        VZERO%Z = 0.D0       
        MTOT = 0.D0
        DO I=1,N
           MTOT = MTOT + M(I)
           VZERO%X = VZERO%X + C(I)%X * M(I) 
           VZERO%Y = VZERO%Y + C(I)%Y * M(I) 
           VZERO%Z = VZERO%Z + C(I)%Z * M(I) 
        END DO
        VZERO%X = VZERO%X/MTOT 
        VZERO%Y = VZERO%Y/MTOT  
        VZERO%Z = VZERO%Z/MTOT                
!     -------------------------------------------------------------------
!     STORE V1 & V2 IN MCONS
!     -------------------------------------------------------------------
        DO I=1,N
           MCONS(I,1)%X=V1(I)%X
           MCONS(I,1)%Y=V1(I)%Y
           MCONS(I,1)%Z=V1(I)%Z

           MCONS(I,2)%X=V2(I)%X
           MCONS(I,2)%Y=V2(I)%Y
           MCONS(I,2)%Z=V2(I)%Z
        END DO
!      -------------------------------------------------------------------
!      MAKE TRASLATIONS
!      -------------------------------------------------------------------
        DO I=1,N
          MCONS(I,3)%X = 1.D0 * MASS(I)
          MCONS(I,3)%Y = 0.D0
          MCONS(I,3)%Z = 0.D0

          MCONS(I,4)%X = 0.D0
          MCONS(I,4)%Y = 1.D0 * MASS(I)
          MCONS(I,4)%Z = 0.D0

          MCONS(I,5)%X = 0.D0
          MCONS(I,5)%Y = 0.D0
          MCONS(I,5)%Z = 1.D0 * MASS(I)
        END DO
!      -------------------------------------------------------------------
!      MAKE ROTATION ALONG  X AXIS
!      -------------------------------------------------------------------
        DO I=1,N
           WRK%X = C(I)%X
           WRK%Y = C(I)%Y
           WRK%Z = C(I)%Z
           
           MCONS(I,6) = ROTATE3D( WRK, 0.D0 , DELTA, 0.D0, VZERO )
        END DO
!      -------------------------------------------------------------------
!      MAKE ROTATION ALONG  Y AXIS
!      -------------------------------------------------------------------
        DO I=1,N
           WRK%X = C(I)%X
           WRK%Y = C(I)%Y
           WRK%Z = C(I)%Z
           
           MCONS(I,7) = ROTATE3D( WRK, -PI2 , DELTA, PI2, VZERO )
        END DO
!      -------------------------------------------------------------------
!      MAKE ROTATION ALONG  Z AXIS
!      -------------------------------------------------------------------
        DO I=1,N
           WRK%X = C(I)%X
           WRK%Y = C(I)%Y
           WRK%Z = C(I)%Z
           
           MCONS(I,8) = ROTATE3D( WRK, DELTA, 0.D0, 0.D0, VZERO )
        END DO
!      -------------------------------------------------------------------
!      NORMALIZE THE VECTORS
!      -------------------------------------------------------------------
        ! NORMALIZE CONSTRAINTS VECTORS

        DO I=1,N
           MCONS(I,6)%X = MCONS(I,6)%X   * MASS(I)
           MCONS(I,6)%Y = MCONS(I,6)%Y   * MASS(I)
           MCONS(I,6)%Z = MCONS(I,6)%Z   * MASS(I)

           MCONS(I,7)%X = MCONS(I,7)%X   * MASS(I)
           MCONS(I,7)%Y = MCONS(I,7)%Y   * MASS(I)
           MCONS(I,7)%Z = MCONS(I,7)%Z   * MASS(I)

           MCONS(I,8)%X = MCONS(I,8)%X   * MASS(I)
           MCONS(I,8)%Y = MCONS(I,8)%Y   * MASS(I)
           MCONS(I,8)%Z = MCONS(I,8)%Z   * MASS(I)
        END DO

NORM_cycle:   DO J=1,MCONSTR
                 NORM = 0.D0
                 DO I=1,N
                    NORM = NORM + MCONS(I,J)%X * MCONS(I,J)%X + &
                                  MCONS(I,J)%Y * MCONS(I,J)%Y + &
                                  MCONS(I,J)%Z * MCONS(I,J)%Z
                 END DO
                 NORM = SQRT ( NORM )
                 DO I=1,N
                    MCONS(I,J)%X = MCONS(I,J)%X / NORM
                    MCONS(I,J)%Y = MCONS(I,J)%Y / NORM
                    MCONS(I,J)%Z = MCONS(I,J)%Z / NORM      
                 END DO          
              END DO NORM_cycle

        RETURN
      END SUBROUTINE MAKEROT

END MODULE CONSTRAINTS
