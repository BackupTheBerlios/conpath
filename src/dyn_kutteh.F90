MODULE CONSTRAINTS_KUTTEH

      USE SYSTEM_UTIL,    ONLY: CPSTOP
      USE XYZ,            ONLY: POINT
      USE ROTATIONS,      ONLY: ROTATE3D
      USE MODPROPERTIES,  ONLY: ANGULARMOM,     &
                                LMOM,           &
                                COMPMOMENTUM

      IMPLICIT NONE

      INTERFACE
         SUBROUTINE DGESV ( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
           INTEGER, INTENT(IN) :: INFO, LDA, LDB, N, NRHS
           INTEGER, INTENT(OUT):: IPIV( * )
           DOUBLE PRECISION, INTENT(INOUT) ::  A( LDA, *), B( LDB, *)
         END SUBROUTINE DGESV
      END INTERFACE
  
      DOUBLE PRECISION, ALLOCATABLE :: R0(:,:), V0(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: R(:,:), V(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: VV(:,:), XV(:,:), XI(:,:), XVO(:,:), VVO(:,:)
      TYPE(POINT), ALLOCATABLE      :: MCONS(:,:)
      INTEGER                       :: MCONSTR
      DOUBLE PRECISION, DIMENSION(3):: LL0, PP0
      DOUBLE PRECISION              :: LL2, PP2
CONTAINS

  SUBROUTINE INIT_VECT (V1, V2, N)
    IMPLICIT NONE
    TYPE(POINT), INTENT(IN) :: V1(*), V2(*)
    INTEGER (KIND=8), INTENT(IN) :: N
    INTEGER :: I,J,K
    DOUBLE PRECISION :: NORM
    
    MCONSTR  = 2
    IF (.NOT.ALLOCATED(MCONS)) ALLOCATE (MCONS(N,MCONSTR)) 
    DO I=1,N
       MCONS(I,1)%X=V1(I)%X
       MCONS(I,1)%Y=V1(I)%Y
       MCONS(I,1)%Z=V1(I)%Z
       
       MCONS(I,2)%X=V2(I)%X
       MCONS(I,2)%Y=V2(I)%Y
       MCONS(I,2)%Z=V2(I)%Z
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
  END SUBROUTINE INIT_VECT

! ===================================================================
  SUBROUTINE R0V0INIT (  C, CP )
! ===================================================================
    USE START_JOB,   ONLY: N => NUMAT
    USE DYNPREPARE,  ONLY: MASS => M
    IMPLICIT NONE
    TYPE(POINT), INTENT(IN)   ::  C(*), CP(*)
    INTEGER                   ::  I
    
    IF (ALLOCATED ( R0 )) DEALLOCATE ( R0 )
    IF (ALLOCATED ( V0 )) DEALLOCATE ( V0 )
    
    ALLOCATE ( R0( 3, N), V0( 3, N))
    
    DO I=1,N
       R0(1, I)=C(I)%X
       R0(2, I)=C(I)%Y
       R0(3, I)=C(I)%Z
       V0(1, I)=CP(I)%X
       V0(2, I)=CP(I)%Y
       V0(3, I)=CP(I)%Z
    ENDDO
    
    CALL ANGULARMOM(R0, V0, N)
    LL0=LMOM
    WRITE (6,*)'ANGMOM BEFORE - R0V0INIT',LL0(1),LL0(2),LL0(3)
    CALL COMPMOMENTUM(R0, V0, N, 0)
    WRITE (6,*)'LINMOM BEFORE - R0V0INIT',LMOM(1),LMOM(2),LMOM(3)
    PP0=LMOM
    
    RETURN
  END SUBROUTINE R0V0INIT
  
      

!     ===================================================================
  SUBROUTINE CONSTKUTTEH ( DT, VCONSTR,  C, CP, C0 )
!     ===================================================================
    USE START_JOB,   ONLY: N => NUMAT,  &
                           NOPI, NOELLE, NOGH, &
                           CLX,CLY,CLZ,CPX,CPY,CPZ,CVG,CVH,NCONSTR
    USE DYNPREPARE,  ONLY: MASS => M
    USE CONVFACTORS, ONLY: AMUTOAU
!     -------------------------------------------------------------------
!     IT CORRECTS THE POSITIONS AND VELOCITIES WITH THE CONSTRAINT FORCES
!     FIRST PART
!     -------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)                :: DT
    DOUBLE PRECISION                            :: DT2
    TYPE(POINT), INTENT(IN)                     :: C0(*)
    TYPE(POINT), INTENT(INOUT)                  :: C(*), CP(*)
    TYPE(POINT), INTENT(IN)                     :: VCONSTR(N,2)
!     -------------------------------------------------------------------
    INTEGER, ALLOCATABLE :: IPIV(:)
    DOUBLE PRECISION, ALLOCATABLE :: CCC(:), QQ(:), XT(:)
    TYPE(POINT), ALLOCATABLE :: CLANG(:), W1(:), W2(:)
    INTEGER :: INFO, I, J, K, II, K1, K2, L
    DOUBLE PRECISION :: RESIDUE, TEMP, PPX,  TTOT, T2, T3, TERM
    DOUBLE PRECISION, DIMENSION (3) :: LL, PP
    DOUBLE PRECISION :: L2, P2
    LOGICAL :: KVEC
    PARAMETER (KVEC=.TRUE.)
    INTEGER, DIMENSION(3,3) :: I3, I4 
    INTEGER, DIMENSION(3) :: I1, I2
    INTEGER :: K3, K4, KK
    DATA I1(1)/2/I2(1)/3/I1(2)/3/I2(2)/1/I1(3)/1/I2(3)/2/
!
!
!  NUMBER OF CONSTRAINTS
!
    WRITE (6,*)'CHECK ',CPX,CPY,CPZ
    IF (ALLOCATED ( VV )) DEALLOCATE ( VV )
    IF (ALLOCATED ( CCC )) DEALLOCATE ( CCC )
    IF (ALLOCATED ( QQ )) DEALLOCATE ( QQ )
    IF (ALLOCATED ( VVO )) DEALLOCATE ( VVO )
    IF (ALLOCATED ( XI )) DEALLOCATE ( XI )
    IF (ALLOCATED ( XV )) DEALLOCATE ( XV )
    IF (ALLOCATED ( XVO )) DEALLOCATE ( XVO )

    ALLOCATE ( VV( NCONSTR, NCONSTR),   &
         VVO( NCONSTR, NCONSTR),   &
         XV( NCONSTR,       1),   &
         XVO( NCONSTR,       1),   &
         IPIV(NCONSTR)      ,  &
         XI(NCONSTR,NCONSTR ), CCC(NCONSTR), QQ(NCONSTR), XT(NCONSTR))

    IF (ALLOCATED ( R )) DEALLOCATE ( R )
    IF (ALLOCATED ( V )) DEALLOCATE ( V )

    ALLOCATE ( R( 3, N), V( 3, N))

    ALLOCATE (    W1(N),   &
                  W2(N)       )

    I3=0
    I4=0
    I3(1,2)=1
    I3(2,1)=-1
    I3(1,3)=-1
    I3(3,1)=1
    I3(2,3)=1
    I3(3,2)=-1
    I4(1,2)=3
    I4(2,1)=3
    I4(1,3)=2
    I4(3,1)=2
    I4(2,3)=1
    I4(3,2)=1

    DT2= DT * DT

       DO I=1,N
          R( 1, I)=C(I)%X
          R( 2, I)=C(I)%Y
          R( 3, I)=C(I)%Z

          V( 1, I)=CP(I)%X
          V( 2, I)=CP(I)%Y
          V( 3, I)=CP(I)%Z
       END DO

!     -------------------------------------------------------------------
!     TO FIND LAGRANGE MULTIPLIER I HAVE TO SOLVE THE LINEAR
!     SYSTEM :
!                   VV (:,:) * LAMBDA(:) = XV(:,:)
!     WHERE :
!                - VV(I,J)   = (\partial Sigma / \partial LAMBDA )_{I,J}
!                - LAMBDA(I)  = Ith LAGRANGE MULTIPLIER
!                - XV(I,1)     = Sigma_I (R',V')
!     -------------------------------------------------------------------

    ALLOCATE ( CLANG(N) )
!     -------------------------------------------------------------------
!     COMPUTES THE ANGULAR MOMENTUM AND LINEAR MOMENTUM
!     -------------------------------------------------------------------     
    DO II=1,1
        XV=0.D0
       LMOM=0.0D0
       CALL ANGULARMOM(R, V, N)
       WRITE (6,*)'ANG MOM BEFORE ',LMOM(1),LMOM(2),LMOM(3)
       LL=LMOM
       LL2=LL(1)*LL(1)+LL(2)*LL(2)+LL(3)*LL(3)
       IF (CLX.NE.0) XV( CLX, 1)=-LL(1)
       IF (CLY.NE.0) XV( CLY, 1)=-LL(2)
       IF (CLZ.NE.0) XV( CLZ, 1)=-LL(3)
       CALL COMPMOMENTUM(R, V, N, 0)
       PP=LMOM
       WRITE (6,*)'LINMOM BEFORE ',LMOM(1),LMOM(2),LMOM(3)
       PP2=PP(1)*PP(1)+PP(2)*PP(2)+PP(3)*PP(3)
          
       IF (CPX.NE.0) XV( CPX, 1)=-PP(1)
       IF (CPY.NE.0) XV( CPY, 1)=-PP(2)
       IF (CPZ.NE.0) XV( CPZ, 1)=-PP(3)
!     -------------------------------------------------------------------
!     COMPUTES THE CONSTRAINTS AT R' AND V'
!     -------------------------------------------------------------------

       DO I=1,N
          IF (.NOT.KVEC)THEN
             IF (CVG.NE.0)THEN
             XV( 1, 1)=XV( 1, 1)-(R( 1, I)-C0(I)%X)*VCONSTR( I, 1)%X 
             XV( 1, 1)=XV( 1, 1)-(R( 2, I)-C0(I)%Y)*VCONSTR( I, 1)%Y 
             XV( 1, 1)=XV( 1, 1)-(R( 3, I)-C0(I)%Z)*VCONSTR( I, 1)%Z 
             ENDIF

             IF (CVH.NE.0)THEN
             XV( 2, 1)=XV( 2, 1)-(R( 1, I)-C0(I)%X)*VCONSTR( I, 2)%X 
             XV( 2, 1)=XV( 2, 1)-(R( 2, I)-C0(I)%Y)*VCONSTR( I, 2)%Y 
             XV( 2, 1)=XV( 2, 1)-(R( 3, I)-C0(I)%Z)*VCONSTR( I, 2)%Z 
             ENDIF
          ELSE
             IF (CVG.NE.0)THEN
             XV( 1, 1)=XV( 1, 1)-(V( 1, I))*VCONSTR( I, 1)%X 
             XV( 1, 1)=XV( 1, 1)-(V( 2, I))*VCONSTR( I, 1)%Y 
             XV( 1, 1)=XV( 1, 1)-(V( 3, I))*VCONSTR( I, 1)%Z 

             ENDIF

             IF (CVH.NE.0)THEN
             XV( 2, 1)=XV( 2, 1)-(V( 1, I))*VCONSTR( I, 2)%X 
             XV( 2, 1)=XV( 2, 1)-(V( 2, I))*VCONSTR( I, 2)%Y 
             XV( 2, 1)=XV( 2, 1)-(V( 3, I))*VCONSTR( I, 2)%Z 
             ENDIF
          ENDIF
             
       ENDDO


!     -------------------------------------------------------------------
!     BUILDS THE MATRIX VV - MATRIX OF LINEAR SYSTEM 
!     -------------------------------------------------------------------
       VV=0.D0

!     -------------------------------------------------------------------
!     FIRST AND SECOND ROW: X1 AND X2 CONSTRAINT
!     VV(X1-X2,X1-X2)
!     -------------------------------------------------------------------
       IF (CVG.NE.0)THEN
          DO K=1,2
             DO J=1,2
               TEMP=0.D0 
               DO I=1,N
                  TEMP=TEMP+VCONSTR( I, K)%X*VCONSTR( I, J)%X/MASS(I) 
                  TEMP=TEMP+VCONSTR( I, K)%Y*VCONSTR( I, J)%Y/MASS(I) 
                  TEMP=TEMP+VCONSTR( I, K)%Z*VCONSTR( I, J)%Z/MASS(I) 
               ENDDO
               IF (KVEC)THEN
                  VV( K, J)=-2.D0*TEMP*DT
               ELSE
                  VV( K, J)=-TEMP*DT2
               ENDIF
            ENDDO
!        -------------------------------------------------------------------
!        VV(X1-X2,Lx)  
!        -------------------------------------------------------------------
            IF (CLX.NE.0)THEN 
               TEMP=0.D0
               DO I=1,N
                  TEMP=TEMP+(-R0(3,I))*VCONSTR( I, K)%Y 
                  TEMP=TEMP+(R0(2,I))*VCONSTR( I, K)%Z 
               ENDDO
               IF (KVEC)THEN
                 VV( K, CLX)=-2.D0*TEMP*DT
               ELSE
                 VV( K, CLX)=-TEMP*DT2
               ENDIF
            ENDIF
!        -------------------------------------------------------------------
!        VV(X1-X2,Ly) 
!        -------------------------------------------------------------------

            IF (CLY.NE.0)THEN 
              TEMP=0.D0
              DO I=1,N
                 TEMP=TEMP+(R0(3,I))*VCONSTR( I, K)%X 
                 TEMP=TEMP+(-R0(1,I))*VCONSTR( I, K)%Z 
              ENDDO
              IF (KVEC)THEN
                VV( K, CLY)=-2.D0*TEMP*DT
              ELSE
                VV( K, CLY)=-TEMP*DT2
              ENDIF
            ENDIF

!        -------------------------------------------------------------------
!        VV(X1-X2,Lz) 
!        -------------------------------------------------------------------
            IF (CLZ.NE.0)THEN 
               TEMP=0.D0
               DO I=1,N
                  TEMP=TEMP+(-R0(2,I))*VCONSTR( I, K)%X 
                  TEMP=TEMP+(R0(1,I))*VCONSTR( I, K)%Y 
               ENDDO
               IF (KVEC)THEN
                 VV( K, CLZ)=-2.D0*TEMP*DT
               ELSE
                 VV( K, CLZ)=-TEMP*DT2
               ENDIF
            ENDIF
 
 
!        -------------------------------------------------------------------
!        VV(X1-X2,Px) 
!        -------------------------------------------------------------------
            IF (CPX.NE.0)THEN 
               TEMP=0.D0
               DO I=1,N
                  TEMP=TEMP+VCONSTR( I, K)%X 
               ENDDO
               IF (KVEC)THEN
                 VV( K, CPX)=-2.D0*TEMP*DT
               ELSE
                 VV( K, CPX)=-TEMP*DT2
               ENDIF
            ENDIF
 
           
!        -------------------------------------------------------------------
!        VV(X1-X2,Py) 
!        -------------------------------------------------------------------
            IF (CPY.NE.0)THEN 
               TEMP=0.D0
               DO I=1,N
                  TEMP=TEMP+VCONSTR( I, K)%Y 
               ENDDO
               IF (KVEC)THEN
                 VV( K, CPY)=-2.D0*TEMP*DT
               ELSE
                 VV( K, CPY)=-TEMP*DT2
               ENDIF
            ENDIF

!        -------------------------------------------------------------------
!        VV(X1-X2,Pz) 
!        -------------------------------------------------------------------
            IF (CPZ.NE.0)THEN 
               TEMP=0.D0
               DO I=1,N
                  TEMP=TEMP+VCONSTR( I, K)%Z 
               ENDDO
               IF (KVEC)THEN
                 VV( K, CPZ)=-2.D0*TEMP*DT
               ELSE
                 VV( K, CPZ)=-TEMP*DT2
               ENDIF
            ENDIF
         ENDDO
      ENDIF

!     -------------------------------------------------------------------
!     ROW 3-5: ANGULAR MOMENTUM (X,Y,Z COMPONENT)
!     -------------------------------------------------------------------

       IF (CLX.NE.0)THEN
          IF (CVG.NE.0)THEN
            DO L=1,3
               DO K=1,2
!         -------------------------------------------------------------------
!         VV(L[x-z],X1-X2) 
!         -------------------------------------------------------------------
                  TEMP=0.D0
                  DO I=1,N
                     IF (I3(1,L).NE.0) TEMP=TEMP-I3(1,L)*V(I4(1,L),I)*VCONSTR( I, K)%X 
                     IF (I3(2,L).NE.0) TEMP=TEMP-I3(2,L)*V(I4(2,L),I)*VCONSTR( I, K)%Y 
                     IF (I3(3,L).NE.0) TEMP=TEMP-I3(3,L)*V(I4(3,L),I)*VCONSTR( I, K)%Z 
                  ENDDO
                  DO I=1,N
                     IF (I3(1,L).NE.0) TEMP=TEMP+2.D0*I3(1,L)*R(I4(1,L),I)*VCONSTR( I, K)%X/DT 
                     IF (I3(2,L).NE.0) TEMP=TEMP+2.D0*I3(2,L)*R(I4(2,L),I)*VCONSTR( I, K)%Y/DT 
                     IF (I3(3,L).NE.0) TEMP=TEMP+2.D0*I3(3,L)*R(I4(3,L),I)*VCONSTR( I, K)%Z/DT 
                  ENDDO
                  VV( L+CLX-1, K)=-TEMP*DT2
               ENDDO
            ENDDO
          ENDIF

!        -------------------------------------------------------------------
!        VV(L[x-z],L[x-z]) 
!        -------------------------------------------------------------------
          DO J=1,3
             DO L=1,3
                TEMP=0.D0
                DO I=1,N
                   DO K=1,3
                      IF (I3(K,J)*I3(K,L).NE.0)THEN
                         TEMP=TEMP+I3(K,J)*I3(K,L)*(-V(I4(K,J),I))*(R0(I4(K,L),I))*MASS(I)
                      ENDIF
                   ENDDO
                ENDDO
                DO I=1,N
                   DO K=1,3
                      IF (I3(K,J)*I3(K,L).NE.0)THEN
                         TEMP=TEMP+2.D0*I3(K,J)*I3(K,L)*(R(I4(K,J),I))*(R0(I4(K,L),I))/DT*MASS(I)
                      ENDIF
                   ENDDO
                ENDDO
                VV( J+CLX-1, L+CLX-1)=-TEMP*DT2
             ENDDO
          ENDDO

!        -------------------------------------------------------------------
!        VV(L[x-z],P[x-z]) 
!        -------------------------------------------------------------------

          IF (CPX.NE.0) THEN
             DO L=1,3
                DO K=1,3
                   TEMP=0.D0
                   DO I=1,N
                      IF (I3(K,L).NE.0)TEMP=TEMP+MASS(I)*I3(K,L)*(-V(I4(K,L),I))
                   ENDDO
                   DO I=1,N
                      IF (I3(K,L).NE.0)TEMP=TEMP+2.D0*MASS(I)*I3(K,L)*(R(I4(K,L),I))/DT
                   ENDDO
                   VV(L+CLX-1,K+CPX-1)=-TEMP*DT2
                ENDDO
             ENDDO
          ENDIF
       ENDIF
      IF (CPX.NE.0)THEN
!        -------------------------------------------------------------------
!        ROW  6-8: LINEAR MOMENTUM
!        -------------------------------------------------------------------
         IF (CVG.NE.0)THEN
            DO K=1,2
!           -------------------------------------------------------------------
!           VV(Px,X1-X2) 
!           -------------------------------------------------------------------
               TEMP=0.D0
               DO I=1,N
                  TEMP=TEMP+2.D0*VCONSTR( I, K)%X 
!                  TEMP=TEMP+VCONSTR( I, K)%X 
               ENDDO
               VV( CPX, K)=-TEMP*DT
!               VV( CPX, K)=-TEMP*DT2
!           -------------------------------------------------------------------
!           VV(Py,X1-X2) 
!           -------------------------------------------------------------------
               TEMP=0.D0
               DO I=1,N
                  TEMP=TEMP+2.D0*VCONSTR( I, K)%Y 
!                  TEMP=TEMP+VCONSTR( I, K)%Y 
               ENDDO
               VV( CPY, K)=-TEMP*DT
!               VV( CPY, K)=-TEMP*DT2
!           -------------------------------------------------------------------
!           VV(Pz,X1-X2) 
!           -------------------------------------------------------------------
               TEMP=0.D0
               DO I=1,N
                  TEMP=TEMP+2.D0*VCONSTR( I, K)%Z 
!                  TEMP=TEMP+VCONSTR( I, K)%Z 
               ENDDO
               VV( CPZ, K)=-TEMP*DT
!               VV( CPZ, K)=-TEMP*DT2
            ENDDO
         ENDIF

!        -------------------------------------------------------------------
!        VV(Px-z,Lx-z) 
!        -------------------------------------------------------------------
         IF (CLX.NE.0)THEN
            DO L=1,3
               DO K=1,3
                  TEMP=0.D0
                  DO I=1,N
                      IF (I3(K,L).NE.0)TEMP=TEMP+2.D0*MASS(I)*I3(K,L)*(R0(I4(K,L),I))/DT
                  ENDDO
                  VV(K+CPX-1,L+CLX-1)=-TEMP*DT2
!                  VV(K+CPX-1,L+CLX-1)=-TEMP*DT2
               ENDDO
            ENDDO
         ENDIF


!        -------------------------------------------------------------------
!        VV(P[x-z],P[x-z]) 
!        -------------------------------------------------------------------
          
          TEMP=0.D0
          DO I=1,N
             TEMP=TEMP+2.D0*MASS(I)
!             TEMP=TEMP+MASS(I)
          ENDDO
          VV( CPX, CPX)=-TEMP*DT
          VV( CPY, CPY)=-TEMP*DT
          VV( CPZ, CPZ)=-TEMP*DT
!          VV( CPX, CPX)=-TEMP*DT2
!          VV( CPY, CPY)=-TEMP*DT2
!          VV( CPZ, CPZ)=-TEMP*DT2
        ENDIF






!       write(*,*)' VVVVVVVV '
!       write(*,'(8G18.9)')VV
!       write(*,'(8G18.9)')XV

!
!
!     -------------------------------------------------------------------
!     LAPACK ROUTINE TO SOLVE THE LINEAR SYSTEM
!     -------------------------------------------------------------------

       VVO=VV
       XVO=XV
       
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
!     -------------------------------------------------------------------
!     TRY TO SEE IF IT WORKED
!     -------------------------------------------------------------------
      DO I=1,NCONSTR
         XT(I)=0
         DO J=1,NCONSTR
            XT(I)=XT(I)+VVO(I,J)*XV(J,1)
         ENDDO
      ENDDO

        
         
       
!     -------------------------------------------------------------------
!     INVERT VV
!     -------------------------------------------------------------------

       XI=0.D0
        DO I=1,NCONSTR
           DO J=1,NCONSTR
              IF (I.EQ.J)XI(I,J)=1.D0
           ENDDO
        ENDDO

        CALL DGESV (NCONSTR, NCONSTR, VVO, NCONSTR, IPIV, XI, NCONSTR, INFO)


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


 
       QQ=0.D0
!     -------------------------------------------------------------------
!     Q CORRECTION
!     -------------------------------------------------------------------
!     TERMS PROPORTIONAL TO DT^2: NO CONTRIBUTION
!     -------------------------------------------------------------------
                        
!     -------------------------------------------------------------------
!     TERMS PROPORTIONAL TO DT^3: THERE IS A TERM IN THE ANGULAR MOMENTUM
!     -------------------------------------------------------------------
!     ANGULAR MOMENTUM PART
      IF (CLX.NE.0)THEN
         DO L=1,3
            TTOT=0.D0
            DO I=1,N
               DO K1=1,3
                  DO J=1,N
                     DO K2=1,3
                        IF (K1.NE.K2.AND.K1.NE.L.AND.K2.NE.L.AND.I.EQ.J)THEN
!           -------------------------------------------------------------------
!           I HAVE A TERM
!           -------------------------------------------------------------------
                           TERM=I3(L,K1)*MASS(I)
                        ELSE
                           TERM=0.D0
                        ENDIF
       
                        TERM=TERM*2.D0*DT*DT2/MASS(I)/MASS(J)

                        DO K3=1,NCONSTR
                           T2=0.D0
                           IF (K3.LT.CLX)THEN
                              IF (K1.EQ.1)T2=XV(K3,1)*VCONSTR(I,K3)%X
                              IF (K1.EQ.2)T2=XV(K3,1)*VCONSTR(I,K3)%Y
                              IF (K1.EQ.3)T2=XV(K3,1)*VCONSTR(I,K3)%Z
                           ELSEIF (K3.GE.CLX.AND.K3.LE.CLZ)THEN
                              IF(I3(K1,K3-CLX+1).NE.0)T2=XV(K3,1)*I3(K1,K3-CLX+1)*MASS(I)*R0(I4(K1,K3-CLX+1),I)
                           ELSEIF (K3.GE.CPX.AND.K3.LE.CPZ)THEN
                              T2=XV(K3,1)*MASS(I)
                           ENDIF

                           DO K4=1,NCONSTR
                              T3=0.D0
                              IF (K4.LT.CLX)THEN
                                 IF (K2.EQ.1)T3=XV(K4,1)*VCONSTR(J,K4)%X
                                 IF (K2.EQ.2)T3=XV(K4,1)*VCONSTR(J,K4)%Y
                                 IF (K2.EQ.3)T3=XV(K4,1)*VCONSTR(J,K4)%Z
                              ELSEIF (K4.GE.CLX.AND.K4.LE.CLZ)THEN
                                 IF(I3(K2,K4-CLX+1).NE.0)T3=XV(K4,1)*I3(K2,K4-CLX+1)*MASS(J)*R0(I4(K2,K4-CLX+1),J)
                              ELSEIF (K4.GE.CPX.AND.K4.LE.CPZ)THEN
                                 T3=XV(K4,1)*MASS(J)
                              ENDIF
                              TTOT=TTOT+T3*T2*TERM
                           ENDDO
                        ENDDO   
                     ENDDO   
                  ENDDO
               ENDDO   
            ENDDO
         QQ(L+CLX-1)=QQ(L+CLX-1)+TTOT
         ENDDO
      ENDIF
!     -------------------------------------------------------------------
!     TERMS PROPORTIONAL TO DT^4
!     NO TERM
!     -------------------------------------------------------------------



      DO I=1,NCONSTR
          CCC(I)=0.D0
          DO J=1,NCONSTR
             CCC(I)=CCC(I)-XI(I,J)*(-XVO(J,1)+QQ(J))
          ENDDO
          XV(I,1)=CCC(I)
       ENDDO
       DO I=1,NCONSTR
!          WRITE (6,*)'CCC ',I,CCC(I),XV(I,1)
!          WRITE (6,*)'XT, Q ',I,XT(I),XVO(I,1),QQ(I)
       ENDDO

!     -------------------------------------------------------------------
!     CYCLE ON THE SECOND ORDER CORRECTIONS
!     -------------------------------------------------------------------
          

         
                  

!     -------------------------------------------------------------------
!     WE HAVE THE SOLUTION OF THE LINEAR SYSTEM .. NOW WE CAN CORRECT THE
!     CONSTRAINTS
!     -------------------------------------------------------------------      

!      ---------------------------------------------------------------------------
!  COMMENTS
!      ---------------------------------------------------------------------------
       DO I=1, N
          CLANG(I)%X = 0.D0
          CLANG(I)%Y = 0.D0
          CLANG(I)%Z = 0.D0
       END DO

       IF (CVG.NE.0)THEN
         DO J=1,2
            DO I=1,N
               CLANG(I)%X = CLANG(I)%X + CCC(J) * VCONSTR(I,J)%X / MASS(I) 
               CLANG(I)%Y = CLANG(I)%Y + CCC(J) * VCONSTR(I,J)%Y / MASS(I) 
               CLANG(I)%Z = CLANG(I)%Z + CCC(J) * VCONSTR(I,J)%Z / MASS(I) 
            END DO
         END DO
       ENDIF
       IF (CLX.NE.0)THEN
          DO J=1,3
             DO I=1,N
                IF(I3(1,J).NE.0)CLANG(I)%X = CLANG(I)%X + CCC(CLX-1+J) *  I3(1,J)*( R0(I4(1,J),I))
                IF(I3(2,J).NE.0)CLANG(I)%Y = CLANG(I)%Y + CCC(CLX-1+J) *  I3(2,J)*( R0(I4(2,J),I))
                IF(I3(3,J).NE.0)CLANG(I)%Z = CLANG(I)%Z + CCC(CLX-1+J) *  I3(3,J)*( R0(I4(3,J),I))
             END DO
          ENDDO
       ENDIF
       IF (CPX.NE.0)THEN
          DO J=CPX,CPZ
             DO I=1,N
                IF (J.EQ.CPX)CLANG(I)%X = CLANG(I)%X + CCC(J) 
                IF (J.EQ.CPY)CLANG(I)%Y = CLANG(I)%Y + CCC(J) 
                IF (J.EQ.CPZ)CLANG(I)%Z = CLANG(I)%Z + CCC(J) 
             END DO
          ENDDO
       ENDIF

       DO I=1,N
          C(I)%X =  C(I)%X - CLANG(I)%X * DT2 
          C(I)%Y =  C(I)%Y - CLANG(I)%Y * DT2
          C(I)%Z =  C(I)%Z - CLANG(I)%Z * DT2
          CP(I)%X = CP(I)%X - CLANG(I)%X * 2.D0 * DT 
          CP(I)%Y = CP(I)%Y - CLANG(I)%Y * 2.D0 * DT
          CP(I)%Z = CP(I)%Z - CLANG(I)%Z * 2.D0 * DT
       END DO
!   --------------------------------------------------------------------   
!   CONTROLS
       DO I=1, N
          R( 1, I)=C(I)%X
          R( 2, I)=C(I)%Y
          R( 3, I)=C(I)%Z

          V( 1, I)=CP(I)%X
          V( 2, I)=CP(I)%Y
          V( 3, I)=CP(I)%Z
       END DO
       LMOM=0.0D0
       CALL ANGULARMOM(R, V, N)
       WRITE (6,*)'ANG MOM AFTER ',LMOM(1),LMOM(2),LMOM(3)
       CALL COMPMOMENTUM (R ,V  ,N, 0 )
       WRITE (6,*)'LINMOM AFTER ',LMOM(1),LMOM(2),LMOM(3)

       ENDDO
!   --------------------------------------------------------------------   
       DEALLOCATE (CLANG)
!     -------------------------------------------------------------------
!     VERIFY THE CONSTRAINT EFFICACY
!     -------------------------------------------------------------------
       DO I=1,N
          W1(I)%X = 0.D0
          W1(I)%Y = 0.D0
          W1(I)%Z = 0.D0
       END DO
       DO I=1,2
          DO J=1,N
             W1(J)%X = W1(J)%X +VCONSTR(J,I)%X 
             W1(J)%Y = W1(J)%Y +VCONSTR(J,I)%Y 
             W1(J)%Z = W1(J)%Z +VCONSTR(J,I)%Z            
          END DO
       END DO
       RESIDUE = 0.D0
       DO I=1,N
          RESIDUE = RESIDUE +  ( C(I)%X - C0(I)%X ) * W1(I)%X
          RESIDUE = RESIDUE +  ( C(I)%Y - C0(I)%Y ) * W1(I)%Y
          RESIDUE = RESIDUE +  ( C(I)%Z - C0(I)%Z ) * W1(I)%Z
       END DO
       WRITE(6,'(A,F21.15)')'CONSTRAINED VALUE ON COORDINATES= ',RESIDUE
       RESIDUE = 0.D0
       DO I=1,N
          RESIDUE = RESIDUE + CP(I)%X * W1(I)%X
          RESIDUE = RESIDUE + CP(I)%Y * W1(I)%Y
          RESIDUE = RESIDUE + CP(I)%Z * W1(I)%Z
       END DO
       WRITE(6,'(A,F21.15)')'CONSTRAINED VALUE ON  VELOCITIES= ',RESIDUE
    DEALLOCATE ( IPIV )
    DEALLOCATE ( W1, W2)
    RETURN
  END SUBROUTINE CONSTKUTTEH
  
  
  
END MODULE CONSTRAINTS_KUTTEH
