MODULE PRINT_MATRIX
  USE READ_ZMAT, ONLY:    NUMAT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE PRINT_MAT2(TITLE, LABELX, LABELY, NCOL, MAT, IR, IC)
    IMPLICIT NONE
    ! Arguments
    INTEGER,           INTENT(IN)                          :: IR, IC
    INTEGER,           INTENT(IN),                OPTIONAL :: NCOL
    DOUBLE PRECISION,  INTENT(IN)                          :: MAT(IR, IC)
    CHARACTER (LEN=*), INTENT(IN),                OPTIONAL :: TITLE
    CHARACTER (LEN=3), INTENT(IN), DIMENSION(IC), OPTIONAL :: LABELX
    CHARACTER (LEN=10),INTENT(IN), DIMENSION(IR), OPTIONAL :: LABELY
    ! Local Variables
    INTEGER :: I, J, TIMES, T, K, ADJ, ISTRT, IEND
    INTEGER :: NUMBER_OF_COLUMN_TO_PRINT, NCP
    INTEGER :: NUMBER_OF_ATOMS_TO_PRINT, NAP
    DOUBLE PRECISION :: HERM
    
    IF (PRESENT(NCOL)) THEN
       IF (NCOL.GT.30) THEN
          WRITE(6, '(A)')'PRINT_MAT2: MAXIMUM NUMBER OF COLUMNS SET TO 30!'
          STOP
       ENDIF
       NUMBER_OF_ATOMS_TO_PRINT  = NCOL / 3  ! We print columns in group of 3..
    ELSE
       NUMBER_OF_ATOMS_TO_PRINT  = 4
    ENDIF
    NUMBER_OF_COLUMN_TO_PRINT    = NUMBER_OF_ATOMS_TO_PRINT * 3
    NAP   = NUMBER_OF_ATOMS_TO_PRINT
    NCP   = NUMBER_OF_COLUMN_TO_PRINT
    ADJ   = MOD(IC, NCP)
    IF (ADJ.NE.0) ADJ=1
    TIMES = IC / NCP + ADJ
    ISTRT = 1
    IEND  = NAP

    ! Print Header...
    IF ( PRESENT( TITLE ) ) WRITE(6,'(/,/A/)')TITLE
    IF (.NOT. PRESENT (LABELX)) THEN
       IF ( PRESENT(LABELY) ) THEN
          WRITE(6,150)('X',K,'Y',K,'Z',K, K = ISTRT, IEND)
       ELSE
          WRITE(6,100)('X',K,'Y',K,'Z',K, K = ISTRT, IEND)
       ENDIF
    ELSE
       WRITE(6,200)(LABELX((K-1)*3+1),LABELX((K-1)*3+2),LABELX((K-1)*3+3), K = ISTRT, IEND)
    ENDIF
    ! Printing Main Loop... 
    DO T = 1, TIMES
       ! Print NCOL Columns
       DO I = 1, IR
          IF (.NOT.PRESENT(LABELY)) THEN
             WRITE(6,'(30F10.6)')(MAT( I, J), J = (T-1)*NCP+1, MIN(T*NCP,IC))
          ELSE
             WRITE(6,'(A10,30F10.6)')LABELY(I),(MAT( I, J), J = (T-1)*NCP+1, MIN(T*NCP,IC))
          ENDIF
       END DO
       ! Print Next Header
       ISTRT = ISTRT + NAP
       IEND  = IEND  + NAP
       IF (T .NE.TIMES) THEN
          IF (.NOT. PRESENT (LABELX)) THEN
             IF ( PRESENT(LABELY) ) THEN
                WRITE(6,150)('X',K,'Y',K,'Z',K, K = ISTRT, MIN(IEND,NUMAT))
             ELSE
                WRITE(6,100)('X',K,'Y',K,'Z',K, K = ISTRT, MIN(IEND,NUMAT))
             ENDIF
          ELSE
             WRITE(6,200)(LABELX((K-1)*3+1),LABELX((K-1)*3+2),LABELX((K-1)*3+3), K = ISTRT, IEND)
          ENDIF
       ENDIF
       ! Cycle ...
    END DO

    IF ( IR .EQ. IC ) THEN
    ! Hermiticity Test
       HERM = 0.D0
       DO I = 1, IR
          DO J = 1, I - 1
             HERM = HERM + MAT(I,J)-MAT(J,I)
          END DO
       END DO
       WRITE(6,'(/A)')'HERMITICITY TEST!'
       WRITE(6,'(A)' )'THIS IS A SQUARE MATRIX.. AND THE VALUE OF MAT(I,J)-MAT(J,I) '
       WRITE(6,'(A,f15.9/)')'SUMMED ALL OVER THE TERMS BUT NOT DIAGONALS IS :: ', HERM
    END IF
    RETURN


    100 FORMAT(/,5X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                 4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                 4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                 4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                 4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                 4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                 4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                 4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                 4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                 4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,/ )

    150 FORMAT(/,15X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                  4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                  4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                  4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                  4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                  4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                  4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                  4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                  4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,  &
                  4X,A1,I1,8X,A1,I1,8X,A1,I1,4X,/ )

    200 FORMAT(/,15X,A3,7X,A3,7X,A3,3X,  &
                  4X,A3,7X,A3,7X,A3,3X,  &
                  4X,A3,7X,A3,7X,A3,3X,  &
                  4X,A3,7X,A3,7X,A3,3X,  &
                  4X,A3,7X,A3,7X,A3,3X,  &
                  4X,A3,7X,A3,7X,A3,3X,  &
                  4X,A3,7X,A3,7X,A3,3X,  &
                  4X,A3,7X,A3,7X,A3,3X,  &
                  4X,A3,7X,A3,7X,A3,3X,  &
                  4X,A3,7X,A3,7X,A3,3X,/ )

  END SUBROUTINE PRINT_MAT2

  SUBROUTINE PRINT_MAT3(TITLE, LABELG, LABELX, LABELY, NCOL, MAT, IR, IC, IC2)
    IMPLICIT NONE
    ! Arguments
    INTEGER,           INTENT(IN)                           :: IR, IC, IC2
    INTEGER,           INTENT(IN),                 OPTIONAL :: NCOL
    DOUBLE PRECISION,  INTENT(IN)                           :: MAT(IR, IC, IC2)
    CHARACTER (LEN=*), INTENT(IN),                 OPTIONAL :: TITLE
    CHARACTER (LEN=10),INTENT(IN), DIMENSION(IR ), OPTIONAL :: LABELG
    CHARACTER (LEN=10),INTENT(IN), DIMENSION(IC2), OPTIONAL :: LABELY
    CHARACTER (LEN=3), INTENT(IN), DIMENSION(IC ), OPTIONAL :: LABELX


    ! Local Variables
    INTEGER :: IVAR


    ! Loop on all internal variables...
    DO IVAR = 1, IR
       WRITE(6,100) LABELG(IVAR)
       CALL PRINT_MAT2 ( TITLE=TITLE, LABELY=LABELY, &
                         NCOL=NCOL, IR=IC,  IC=IC2 , &
                         LABELX=LABELX,              & 
                         MAT=RESHAPE(MAT(IVAR:IVAR,1:IC,1:IC2),SHAPE=(/IC,IC2/)))
    END DO

    RETURN

    100 FORMAT("DEBUGGING INTERNAL COORDINATE::",A)
  END SUBROUTINE PRINT_MAT3

END MODULE PRINT_MATRIX
