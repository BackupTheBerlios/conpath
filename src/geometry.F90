!
! Author : Teodoro Laino             Scuola Normale Superiore -  NEST  - INFM
!
! A - No part of this software package  can  be reproduced  or transmitted in 
!     any form or by any means: electronic, mechanical, or otherwise, without 
!     prior written permission from the author therein.
! B - No part of this source code can be altered unless permission is granted 
!     by the author itself.
! C - This included program can not be used in any  unlawful  manner  or in a 
!     way that would be construed as unlawful by the government or equivalent. 
! D - No part of the program can be copied and re-sold under another name. 
! E - Title, ownership rights and intellectual property rights  in and to the 
!     Software and Documentation shall remain to the author.
!     This Software is protected by the copyright laws of Italian  government
!     and international copyright treaties.
!
! Copyright 2002-2004 - Date 06.12.02 
! History:
! CSCS (Manno- Switzerland) - Lugano
! ETHZ (Zurich) - Date 04.11.02 - Second derivatives Implemented
! ETHZ (Zurich) - Date 29.11.02 - Fixed Second Derivatives for small angles 
!                                 or \pi
MODULE GEO_OBJECTS
  IMPLICIT NONE
CONTAINS
  DOUBLE PRECISION FUNCTION BOND_VAL ( I, J, C )
! COMPUTES THE BOND VALUE BEETWEN ATOMS I AND J
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I, J
    DOUBLE PRECISION, INTENT(IN) :: C(*)  
    INTEGER :: II, JJ
    DOUBLE PRECISION :: XDIST, YDIST, ZDIST

    II=(I-1)*3
    JJ=(J-1)*3
    XDIST=C(II+1)-C(JJ+1)
    YDIST=C(II+2)-C(JJ+2)
    ZDIST=C(II+3)-C(JJ+3)
    BOND_VAL=SQRT(XDIST*XDIST+YDIST*YDIST+ZDIST*ZDIST)

    RETURN
  END FUNCTION BOND_VAL

  DOUBLE PRECISION FUNCTION ANGLE_VAL ( J, K, I, C )
! COMPUTES THE ANGLE VALUE BEETWEN ATOMS J, K, I  (IN THIS ORDER!!!)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I, J, K
    DOUBLE PRECISION, INTENT(IN) :: C(*)
    INTEGER :: II, JJ, KK
    DOUBLE PRECISION :: X1DIST, Y1DIST, Z1DIST, R1IN
    DOUBLE PRECISION :: X2DIST, Y2DIST, Z2DIST, R2IN

    II=(I-1)*3
    JJ=(J-1)*3
    KK=(K-1)*3

    X1DIST=C(KK+1)-C(JJ+1)
    Y1DIST=C(KK+2)-C(JJ+2)
    Z1DIST=C(KK+3)-C(JJ+3)
    R1IN=SQRT(X1DIST*X1DIST+Y1DIST*Y1DIST+Z1DIST*Z1DIST)

    X2DIST=C(KK+1)-C(II+1)
    Y2DIST=C(KK+2)-C(II+2)
    Z2DIST=C(KK+3)-C(II+3)
    R2IN=SQRT(X2DIST*X2DIST+Y2DIST*Y2DIST+Z2DIST*Z2DIST)

    ANGLE_VAL = ACOS((X1DIST*X2DIST+Y1DIST*Y2DIST+Z1DIST*Z2DIST)/(R1IN*R2IN))
    RETURN
  END FUNCTION ANGLE_VAL

  DOUBLE PRECISION FUNCTION TORS_VAL ( K, J, M, L, C )
! COMPUTES THE TORSION ANGLE VALUE BEETWEN ATOMS K, J, M, L (ORDER!!!)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: K, J, M, L
    DOUBLE PRECISION, INTENT(IN) :: C(*)
    INTEGER :: MM, JJ, KK, LL
    DOUBLE PRECISION, DIMENSION(3) :: X, Y, Z, V1, V2
    DOUBLE PRECISION :: V1NORM, V2NORM, TNORM, PRDT, COSPHI
    DOUBLE PRECISION :: PRDT_SIGN, PI

    PI = ATAN(1.D0)*4.D0
    KK=(K-1)*3
    JJ=(J-1)*3
    MM=(M-1)*3
    LL=(L-1)*3
    X(1)=C(JJ+1)-C(KK+1)
    X(2)=C(JJ+2)-C(KK+2)    
    X(3)=C(JJ+3)-C(KK+3)    
    
    Y(1)=C(LL+1)-C(MM+1)
    Y(2)=C(LL+2)-C(MM+2)    
    Y(3)=C(LL+3)-C(MM+3)

    Z(1)=C(MM+1)-C(JJ+1)
    Z(2)=C(MM+2)-C(JJ+2)    
    Z(3)=C(MM+3)-C(JJ+3)    
    
    CALL CROSSPRODUCT(X,Z,V1)
    CALL CROSSPRODUCT(Z,Y,V2)
    V1NORM = SQRT( DDOT(V1,V1,3) )
    V2NORM = SQRT( DDOT(V2,V2,3) )
    TNORM = V1NORM * V2NORM
    PRDT = DDOT(V1, V2, 3)
    PRDT_SIGN = DDOT(X,V2,3)
    COSPHI = PRDT / (V1NORM*V2NORM)
    IF (ABS(COSPHI).GT.1) COSPHI=1.D0*COSPHI/ABS(COSPHI)
    TORS_VAL = SIGN(1.D0,PRDT_SIGN) * ACOS ( COSPHI )
    IF (TORS_VAL.LT.0) TORS_VAL = 2.D0*PI + TORS_VAL

    RETURN
  END FUNCTION TORS_VAL


  SUBROUTINE CROSSPRODUCT(V1, V2, F)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: V1, V2
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: F

    F(1)=V1(2)*V2(3)-V1(3)*V2(2)
    F(2)=V1(3)*V2(1)-V1(1)*V2(3)
    F(3)=V1(1)*V2(2)-V1(2)*V2(1)
    RETURN
  END SUBROUTINE CROSSPRODUCT

  DOUBLE PRECISION FUNCTION DDOT(X,Y,N)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN), DIMENSION(N) :: X, Y
    INTEGER, INTENT(IN) :: N
    INTEGER :: I

    DDOT=0.D0
    DO I=1,N
       DDOT = DDOT + X(I) * Y(I)
    END DO
    RETURN
  END FUNCTION DDOT

END MODULE GEO_OBJECTS



MODULE DERIVATE_GEO_OBJECTS
  IMPLICIT NONE
CONTAINS

  SUBROUTINE DERI_BOND ( I, J, C, F, FACT )
!   CALCOLA LE DERIVATE PER UN BOND
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I, J
    DOUBLE PRECISION, INTENT(IN) :: C(*), FACT
    DOUBLE PRECISION, INTENT(INOUT) :: F(*)
    INTEGER :: II, JJ
    DOUBLE PRECISION :: XDIST, YDIST, ZDIST, B, FMOD

    II=(I-1)*3
    JJ=(J-1)*3
    XDIST=C(II+1)-C(JJ+1)
    YDIST=C(II+2)-C(JJ+2)
    ZDIST=C(II+3)-C(JJ+3)

    B = SQRT(XDIST*XDIST+YDIST*YDIST+ZDIST*ZDIST)
    
    FMOD = 1.D0 / B

!   STORE THE FORCES IN F
    F(II+1) =   XDIST * FMOD * FACT + F(II+1)
    F(II+2) =   YDIST * FMOD * FACT + F(II+2)
    F(II+3) =   ZDIST * FMOD * FACT + F(II+3)

    F(JJ+1) = - XDIST * FMOD * FACT + F(JJ+1)
    F(JJ+2) = - YDIST * FMOD * FACT + F(JJ+2)
    F(JJ+3) = - ZDIST * FMOD * FACT + F(JJ+3)

    RETURN
  END SUBROUTINE DERI_BOND

    
  SUBROUTINE DERI_ANGLE ( J, K, I, C, F, FACT )
! CALCOLA LE DERIVATE PER UN ANGLE IN ( J, K, I ) ORDER!!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I, J, K
    DOUBLE PRECISION, INTENT(IN) :: C(*), FACT
    DOUBLE PRECISION, INTENT(INOUT) :: F(*)
    INTEGER :: II, JJ, KK
    DOUBLE PRECISION :: X1DIST, Y1DIST, Z1DIST
    DOUBLE PRECISION :: X2DIST, Y2DIST, Z2DIST
    DOUBLE PRECISION :: R1IN, R2IN, ANG, FMOD, PRODSC, SQ, F1, F2
    DOUBLE PRECISION, DIMENSION(9) :: DER
    DOUBLE PRECISION :: PI
    
    PI=4.D0*ATAN(1.D0)
    II=(I-1)*3
    JJ=(J-1)*3
    KK=(K-1)*3
    X1DIST=C(KK+1)-C(JJ+1)
    Y1DIST=C(KK+2)-C(JJ+2)
    Z1DIST=C(KK+3)-C(JJ+3)
    R1IN=SQRT(X1DIST*X1DIST+Y1DIST*Y1DIST+Z1DIST*Z1DIST)

    X2DIST=C(KK+1)-C(II+1)
    Y2DIST=C(KK+2)-C(II+2)
    Z2DIST=C(KK+3)-C(II+3)
    R2IN=SQRT(X2DIST*X2DIST+Y2DIST*Y2DIST+Z2DIST*Z2DIST)

    PRODSC = X1DIST*X2DIST+Y1DIST*Y2DIST+Z1DIST*Z2DIST
    SQ = 1.D0 / (R1IN*R2IN)
    ANG = ACOS( PRODSC * SQ )

    F1 = PRODSC / ( R2IN * R2IN )
    F2 = PRODSC / ( R1IN * R1IN )

!   DERIVATIVES OF COS(ANG) RESPECT TO ATOM K
    DER(1)=( (-F1+1.D0) * X2DIST + (-F2+1.D0) * X1DIST ) * SQ
    DER(2)=( (-F1+1.D0) * Y2DIST + (-F2+1.D0) * Y1DIST ) * SQ
    DER(3)=( (-F1+1.D0) * Z2DIST + (-F2+1.D0) * Z1DIST ) * SQ    
!   DERIVATIVES OF COS(ANG) RESPECT TO ATOM J
    DER(4)=( F2 * X1DIST - X2DIST ) * SQ
    DER(5)=( F2 * Y1DIST - Y2DIST ) * SQ
    DER(6)=( F2 * Z1DIST - Z2DIST ) * SQ
!   DERIVATIVES OF COS(ANG) RESPECT TO ATOM I
    DER(7)=( F1 * X2DIST - X1DIST ) * SQ
    DER(8)=( F1 * Y2DIST - Y1DIST ) * SQ
    DER(9)=( F1 * Z2DIST - Z1DIST ) * SQ    
!   CHECK FOR NUMERICAL INCONSISTENCIES
    IF ((ABS(ANG).LT.0.00001D0).OR.(ABS(ANG-PI).LT.0.00001D0)) THEN
       FMOD=   0.D0
    ELSE
       FMOD= - 1.D0 / SIN(ANG)
    ENDIF

    F(KK+1) =  DER(1) * FMOD * FACT + F(KK+1)
    F(KK+2) =  DER(2) * FMOD * FACT + F(KK+2)
    F(KK+3) =  DER(3) * FMOD * FACT + F(KK+3)
                             
    F(JJ+1) =  DER(4) * FMOD * FACT + F(JJ+1)
    F(JJ+2) =  DER(5) * FMOD * FACT + F(JJ+2)
    F(JJ+3) =  DER(6) * FMOD * FACT + F(JJ+3)
                         
    F(II+1) =  DER(7) * FMOD * FACT + F(II+1)
    F(II+2) =  DER(8) * FMOD * FACT + F(II+2)
    F(II+3) =  DER(9) * FMOD * FACT + F(II+3)

    RETURN
  END SUBROUTINE DERI_ANGLE

  SUBROUTINE DERI_TORSION ( K, J, M, L, C, F, FACT )
    USE GEO_OBJECTS, ONLY :  CROSSPRODUCT, &
                             DDOT
! CALCOLA LE DERIVATE PER UN TORSION IN ( K, J, M, L) ORDER!!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: K, J, M, L
    DOUBLE PRECISION, INTENT(IN) :: C(*), FACT
    DOUBLE PRECISION, INTENT(INOUT) :: F(*)
    INTEGER :: MM, JJ, KK, LL, IK
    DOUBLE PRECISION, DIMENSION(3) :: X, Y, Z, V1, V2
    DOUBLE PRECISION :: V1NORM, V2NORM, TNORM, PRDT, COSPHI, TOR
    DOUBLE PRECISION, DIMENSION(9)  :: DDOTD, DNORMD, DCOSD
    DOUBLE PRECISION, DIMENSION(12) :: DER
    DOUBLE PRECISION :: FD1, FD2, XNORMQ, FMOD
    DOUBLE PRECISION :: PI
    DOUBLE PRECISION :: PRDT_SIGN
    
    PI=4.D0*ATAN(1.D0)
    KK=(K-1)*3
    JJ=(J-1)*3
    MM=(M-1)*3
    LL=(L-1)*3
    X(1)=C(JJ+1)-C(KK+1)
    X(2)=C(JJ+2)-C(KK+2)    
    X(3)=C(JJ+3)-C(KK+3)    
    
    Y(1)=C(LL+1)-C(MM+1)
    Y(2)=C(LL+2)-C(MM+2)    
    Y(3)=C(LL+3)-C(MM+3)

    Z(1)=C(MM+1)-C(JJ+1)
    Z(2)=C(MM+2)-C(JJ+2)    
    Z(3)=C(MM+3)-C(JJ+3)    
    
    CALL CROSSPRODUCT(X,Z,V1)
    CALL CROSSPRODUCT(Z,Y,V2)
    V1NORM    = SQRT( DDOT(V1,V1,3) )
    V2NORM    = SQRT( DDOT(V2,V2,3) )
    TNORM     = V1NORM * V2NORM
    PRDT      = DDOT(V1, V2, 3)
    PRDT_SIGN = DDOT(X,V2,3)
    COSPHI    = PRDT / TNORM
    IF (ABS(COSPHI).GT.1) COSPHI=COSPHI/ABS(COSPHI)
    TOR = SIGN(1.D0,PRDT_SIGN) * ACOS ( COSPHI )    

!   CALCOLO DERIVATE DELL'ANGOLO DI TORSIONE <TOR>
!   DDOT PART DERIVATIVES
    DDOTD(1)=   Z(2)*V2(3) - Z(3)*V2(2)                             ! DX1
    DDOTD(2)= - Z(1)*V2(3) + Z(3)*V2(1)                             ! DX2
    DDOTD(3)=   Z(1)*V2(2) - Z(2)*V2(1)                             ! DX3
    DDOTD(4)= - Z(2)*V1(3) + Z(3)*V1(2)                             ! DY1
    DDOTD(5)=   Z(1)*V1(3) - Z(3)*V1(1)                             ! DY2
    DDOTD(6)= - Z(1)*V1(2) + Z(2)*V1(1)                             ! DY3
    DDOTD(7)=   Y(2)*V1(3) - X(2)*V2(3) - Y(3)*V1(2) + X(3)*V2(2)   ! DZ1
    DDOTD(8)= - Y(1)*V1(3) + X(1)*V2(3) + Y(3)*V1(1) - X(3)*V2(1)   ! DZ2
    DDOTD(9)=   Y(1)*V1(2) - X(1)*V2(2) - Y(2)*V1(1) + X(2)*V2(1)   ! DZ3    
!   NORM PART DERIVATIVES
    FD1=V2NORM/V1NORM
    DNORMD(1)=( Z(2)*V1(3) - Z(3)*V1(2) )* FD1
    DNORMD(2)=(-Z(1)*V1(3) + Z(3)*V1(1) )* FD1
    DNORMD(3)=( Z(1)*V1(2) - Z(2)*V1(1) )* FD1
    FD2=V1NORM/V2NORM
    DNORMD(4)=(-Z(2)*V2(3) + Z(3)*V2(2) )* FD2
    DNORMD(5)=( Z(1)*V2(3) - Z(3)*V2(1) )* FD2
    DNORMD(6)=(-Z(1)*V2(2) + Z(2)*V2(1) )* FD2
    !
    DNORMD(7)=( Y(2)*V2(3) - Y(3)*V2(2) )* FD2 + ( X(3)*V1(2)-X(2)*V1(3))* FD1
    DNORMD(8)=(-Y(1)*V2(3) + Y(3)*V2(1) )* FD2 + ( X(1)*V1(3)-X(3)*V1(1))* FD1
    DNORMD(9)=( Y(1)*V2(2) - Y(2)*V2(1) )* FD2 + (-X(1)*V1(2)+X(2)*V1(1))* FD1
!   COS PART DERIVATIVES
    XNORMQ=1.d0/(TNORM*TNORM)
    DO IK=1,9
       DCOSD(IK)=( DDOTD(IK)*TNORM - PRDT*DNORMD(IK) ) * XNORMQ
    ENDDO
!   CHECK FOR NUMERICAL INCONSISTENCIES
    IF ((ABS(TOR).LT.0.00001D0).OR.(ABS(TOR-PI).LT.0.00001D0).OR. &
        (ABS(TOR-2.D0*PI).LT.0.00001D0)) THEN
       FMOD=   0.D0
    ELSE
       FMOD= - (1.D0 / ABS(SIN(TOR))) * SIGN(1.D0,PRDT_SIGN)
    ENDIF
!   COMPUTE ATOM DERIVATIVES
    DER(1)= - DCOSD(1) * FMOD
    DER(2)= - DCOSD(2) * FMOD
    DER(3)= - DCOSD(3) * FMOD
    DER(4)= ( DCOSD(1) - DCOSD(7) ) * FMOD
    DER(5)= ( DCOSD(2) - DCOSD(8) ) * FMOD
    DER(6)= ( DCOSD(3) - DCOSD(9) ) * FMOD
    DER(7)= ( DCOSD(7) - DCOSD(4) ) * FMOD
    DER(8)= ( DCOSD(8) - DCOSD(5) ) * FMOD
    DER(9)= ( DCOSD(9) - DCOSD(6) ) * FMOD
    DER(10)=  DCOSD(4) * FMOD
    DER(11)=  DCOSD(5) * FMOD
    DER(12)=  DCOSD(6) * FMOD

    F(KK+1) = F(KK+1) + DER(1)    * FACT
    F(KK+2) = F(KK+2) + DER(2)    * FACT
    F(KK+3) = F(KK+3) + DER(3)    * FACT
                                  
    F(JJ+1) = F(JJ+1) + DER(4)    * FACT
    F(JJ+2) = F(JJ+2) + DER(5)    * FACT
    F(JJ+3) = F(JJ+3) + DER(6)    * FACT
                                  
    F(MM+1) = F(MM+1) + DER(7)    * FACT
    F(MM+2) = F(MM+2) + DER(8)    * FACT
    F(MM+3) = F(MM+3) + DER(9)    * FACT
                                  
    F(LL+1) = F(LL+1) + DER(10)   * FACT
    F(LL+2) = F(LL+2) + DER(11)   * FACT
    F(LL+3) = F(LL+3) + DER(12)   * FACT 

    RETURN
  END SUBROUTINE DERI_TORSION

END MODULE DERIVATE_GEO_OBJECTS

MODULE DERIVATE_2_GEO_OBJECTS
  IMPLICIT NONE
CONTAINS
! General use of vector F(:,:)
!
!
!   F(:) CONTAINS THE SECOND DERIVATIVES :
!                                                                     J      K    I
!        |   ...      ....         ....         ....     ...   |   |\                |
!        |                                                     |   |  \              |  J 
!        |   ...     dIx dJx      dIx dJy      dIx dJz   ...   |   |    \            |
!        |                                                     |   |      \          | 
!        |   ...     dIy dJx      dIy dJy      dIy dJz   ...   |   |        \        |  K
!        |                                                     |   |          \      |
!        |   ...     dIz dJx      dIz dJy      dIz dJz   ...   |   |            \    |
!        |                                                     |   |              \  |  I
!        |   ...      ....         ....         ....     ...   |   |                \|
!
!   F can be stored in an UPPER TRIANGULAR matrix applying the Schwarz rule (for 2nd derivatives)
!   The dimension of F is : (3N+1)*3N /2 where N is the number of atomic species 
!       number of variable in input (bond=2, angle=3, torsion=4)

  SUBROUTINE DERI_2_BOND ( I, J, C, F, FACT )
!   CALCOLA LE DERIVATE SECONDE PER UN BOND
!   F IS A DIAGONAL MATRIX OF DIMENSION : 21
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: I, J
    DOUBLE PRECISION, INTENT(IN) :: C(*), FACT
    DOUBLE PRECISION, INTENT(INOUT) :: F(*)
    ! Local Variables
    INTEGER :: II, JJ, K, KX1, KX2
    INTEGER :: ISDERR, INDF
    DOUBLE PRECISION :: XDIST, YDIST, ZDIST
    DOUBLE PRECISION :: XXDIST, XYDIST, XZDIST, YYDIST, YZDIST, ZZDIST
    DOUBLE PRECISION :: R, RINV, RINV3, MODR
    DOUBLE PRECISION, DIMENSION(27) :: SDERR
    DOUBLE PRECISION, DIMENSION(6)  :: FDERR


    II=(I-1)*3
    JJ=(J-1)*3
    XDIST=C(II+1)-C(JJ+1)
    YDIST=C(II+2)-C(JJ+2)
    ZDIST=C(II+3)-C(JJ+3)

    XXDIST = XDIST * XDIST
    XYDIST = XDIST * YDIST
    XZDIST = XDIST * ZDIST
    YYDIST = YDIST * YDIST
    YZDIST = YDIST * ZDIST
    ZZDIST = ZDIST * ZDIST

    R = SQRT( XDIST*XDIST + YDIST*YDIST + ZDIST*ZDIST )
    RINV  = 1.0D0 / R
    RINV3 = 1.0D0 / R**3

    DO K= 1, 27
       SDERR(K) = 0.D0
    END DO
    DO K= 1, 6
       FDERR(K) = 0.D0
    END DO

! FIRST DERIVATIVE OF R
    ! DERIVATIVE WITH RESPECT TO ATOM I
    FDERR( 1) = XDIST * RINV    ! D/DIx
             ! WITH RESPECT TO ATOM I
             SDERR( 1) =   RINV  -  XXDIST * RINV3   ! D/DIx D/DIx 
             SDERR( 2) =         -  XYDIST * RINV3   ! D/DIy D/DIx
             SDERR( 3) =         -  XZDIST * RINV3   ! D/DIz D/DIx
             ! WITH RESPECT TO ATOM J
             SDERR(10) = - SDERR( 1)                 ! D/DJx D/DIx 
             SDERR(16) = - SDERR( 2)                 ! D/DJy D/DIx
             SDERR(22) = - SDERR( 3)                 ! D/DJz D/DIx
    FDERR( 2) = YDIST * RINV    ! D/DIy
             ! WITH RESPECT TO ATOM I
             SDERR( 4) =   SDERR( 2)                 ! D/DIx D/DIy 
             SDERR( 5) =   RINV  -  YYDIST * RINV3   ! D/DIy D/DIy
             SDERR( 6) =         -  YZDIST * RINV3   ! D/DIz D/DIy
             ! WITH RESPECT TO ATOM J
             SDERR(11) = - SDERR( 4)                 ! D/DJx D/DIy 
             SDERR(17) = - SDERR( 5)                 ! D/DJy D/DIy
             SDERR(23) = - SDERR( 6)                 ! D/DJz D/DIy
    FDERR( 3) = ZDIST * RINV    ! D/DIz
             ! WITH RESPECT TO ATOM I
             SDERR( 7) =   SDERR( 3)                 ! D/DIx D/DIz 
             SDERR( 8) =   SDERR( 6)                 ! D/DIy D/DIz
             SDERR( 9) =   RINV  -  ZZDIST * RINV3   ! D/DIz D/DIz
             ! WITH RESPECT TO ATOM J
             SDERR(12) = - SDERR( 7)                 ! D/DJx D/DIz 
             SDERR(18) = - SDERR( 8)                 ! D/DJy D/DIz
             SDERR(24) = - SDERR( 9)                 ! D/DJz D/DIz
    ! DERIVATIVE WITH RESPECT TO ATOM J
    FDERR( 4) = - FDERR( 1)     ! D/DJx
             ! WITH RESPECT TO ATOM J
             SDERR(13) =   SDERR( 1)                 ! D/DJx D/DJx 
             SDERR(14) =   SDERR( 2)                 ! D/DJy D/DJx
             SDERR(15) =   SDERR( 3)                 ! D/DJz D/DJx
    FDERR( 5) = - FDERR( 2)     ! D/DJy
             ! WITH RESPECT TO ATOM J
             SDERR(19) =   SDERR( 4)                 ! D/DJx D/DJy 
             SDERR(20) =   SDERR( 5)                 ! D/DJy D/DJy
             SDERR(21) =   SDERR( 6)                 ! D/DJz D/DJy
    FDERR( 6) = - FDERR( 3)     ! D/DJz
             ! WITH RESPECT TO ATOM J
             SDERR(25) =   SDERR( 7)                 ! D/DJx D/DJz 
             SDERR(26) =   SDERR( 8)                 ! D/DJy D/DJz
             SDERR(27) =   SDERR( 9)                 ! D/DJz D/DJz

!   STORE THE HESSIAN IN F
   ISDERR = 0
   INDF   = 0
   DO KX1 = 1, 6
      DO KX2 = 1, KX1
         INDF   = INDF   + 1
         ISDERR = ISDERR + 1
         F(INDF) = SDERR(ISDERR) * FACT
      END DO
      MODR   = MOD( KX1 , 3 )  
      IF ( MODR .EQ. 0 ) MODR = 3
      ISDERR = ISDERR + ( 3 - MODR )
   END DO

! DEBUG SECTION
#ifdef DEBUG_DERI2
    CALL CHECK_FIRST_DERIVATIVES(LABEL='BOND',I1=I,I2=J,COORD=C,DER1=FDERR,FAC=FACT)
    CALL CHECK_SECOND_DERIVATIVES(LABEL='BOND',I1=I,I2=J,COORD=C,DER2=F,FAC=1.d0)
#endif

   RETURN
 END SUBROUTINE DERI_2_BOND


  SUBROUTINE DERI_2_ANGLE( J, K, I, C, F, FACT )
! CALCOLA LE DERIVATE PER UN ANGLE IN ( J, K, I ) ORDER!!
! F IS A DIAGONAL MATRIX OF DIMENSION : 45
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: I, J, K
    DOUBLE PRECISION, INTENT(IN) :: C(*), FACT
    DOUBLE PRECISION, INTENT(INOUT) :: F(*)
    ! Local Variables
    INTEGER :: II, JJ, KK, KX1, KX2, L
    INTEGER :: IFDERCOS, ISDERCOS, INDF
    DOUBLE PRECISION :: PI, TERM1,    TERM2,    MODR
    DOUBLE PRECISION :: X1DIST,   Y1DIST,   Z1DIST
    DOUBLE PRECISION :: X2DIST,   Y2DIST,   Z2DIST
    DOUBLE PRECISION :: X1X1DIST, X1Y1DIST, X1Z1DIST, Y1Y1DIST, Y1Z1DIST, Z1Z1DIST
    DOUBLE PRECISION :: X2X2DIST, X2Y2DIST, X2Z2DIST, Y2Y2DIST, Y2Z2DIST, Z2Z2DIST
    DOUBLE PRECISION :: X1X2DIST, X1Y2DIST, X1Z2DIST
    DOUBLE PRECISION :: Y1X2DIST, Y1Y2DIST, Y1Z2DIST
    DOUBLE PRECISION :: Z1X2DIST, Z1Y2DIST, Z1Z2DIST
    DOUBLE PRECISION :: R1, R2, R1INV, R2INV, R1R2, R13R23
    DOUBLE PRECISION :: R1SUR2, R2SUR1, R1SUR23, R2SUR13
    DOUBLE PRECISION :: ANG, FMOD1, FMOD2, PRODSC, SQ, SQS
    DOUBLE PRECISION, DIMENSION(54) :: SDERCOS, SDERP, SDERR1R2
    DOUBLE PRECISION, DIMENSION(9)  :: FDERCOS, FDERP, FDERR1R2
    ! Statements execution

    ! Pi value...
    PI=4.D0*ATAN(1.D0)
    II         = (I-1) * 3
    JJ         = (J-1) * 3
    KK         = (K-1) * 3
    X1DIST     = C(KK+1) - C(JJ+1)
    Y1DIST     = C(KK+2) - C(JJ+2)
    Z1DIST     = C(KK+3) - C(JJ+3)
    X1X1DIST   = X1DIST * X1DIST
    X1Y1DIST   = X1DIST * Y1DIST
    X1Z1DIST   = X1DIST * Z1DIST
    Y1Y1DIST   = Y1DIST * Y1DIST
    Y1Z1DIST   = Y1DIST * Z1DIST
    Z1Z1DIST   = Z1DIST * Z1DIST
    R1         = SQRT( X1DIST * X1DIST + Y1DIST * Y1DIST + Z1DIST * Z1DIST )
    R1INV      = 1.D0 / R1
               
    X2DIST     = C(KK+1) - C(II+1)
    Y2DIST     = C(KK+2) - C(II+2)
    Z2DIST     = C(KK+3) - C(II+3)
    X2X2DIST   = X2DIST * X2DIST
    X2Y2DIST   = X2DIST * Y2DIST
    X2Z2DIST   = X2DIST * Z2DIST
    Y2Y2DIST   = Y2DIST * Y2DIST
    Y2Z2DIST   = Y2DIST * Z2DIST
    Z2Z2DIST   = Z2DIST * Z2DIST
    R2         = SQRT( X2DIST * X2DIST + Y2DIST * Y2DIST + Z2DIST * Z2DIST )
    R2INV      = 1.D0 / R2
               
    X1X2DIST   = X1DIST * X2DIST
    X1Y2DIST   = X1DIST * Y2DIST
    X1Z2DIST   = X1DIST * Z2DIST
    Y1X2DIST   = Y1DIST * X2DIST
    Y1Y2DIST   = Y1DIST * Y2DIST
    Y1Z2DIST   = Y1DIST * Z2DIST
    Z1X2DIST   = Z1DIST * X2DIST
    Z1Y2DIST   = Z1DIST * Y2DIST
    Z1Z2DIST   = Z1DIST * Z2DIST
    R1R2       = R1 * R2
    R13R23     = R1**3 * R2**3
    R2SUR1     = R2 * R1INV
    R1SUR2     = R1 * R2INV
    R1SUR23    = R1 * R2INV**3
    R2SUR13    = R2 * R1INV**3               

    PRODSC     = X1DIST * X2DIST + Y1DIST * Y2DIST + Z1DIST * Z2DIST
    SQ         = 1.D0  / ( R1R2 )
    SQS        = SQ * SQ
    ANG        = ACOS( PRODSC * SQ )

!   CHECK FOR NUMERICAL INCONSISTENCIES
    IF ((ABS(ANG).LT.0.00001D0).OR.(ABS(ANG-PI).LT.0.00001D0)) THEN
       FMOD1      = 0.D0
       FMOD2      = 0.D0
    ELSE
       FMOD1      = ( PRODSC * SQ ) / ( SIN(ANG) * SIN(ANG) ) 
       FMOD2      = - 1.D0 / SIN(ANG)
    ENDIF


    DO L = 1, 54
       SDERCOS(L) = 0.D0
       SDERP(L) = 0.D0
       SDERR1R2(L) = 0.D0
    END DO
    DO L= 1, 9
       FDERP(L) = 0.D0
       FDERR1R2(L) = 0.D0
    END DO

!   FIRST DERIVATIVES OF PRODSC
    ! DERIVATIVES RESPECT TO ATOM J
                    ! SECOND DERIVATIVES OF FDERP ( ONLY NON ZERO ONES )
    FDERP(1) = - X2DIST  ! D/DJx
                    SDERP(10) = -1.D0   ! D/DKx D/DJx
                    SDERP(28) =  1.D0   ! D/DIx D/DJx
    FDERP(2) = - Y2DIST  ! D/DJy                     
                    SDERP(17) = -1.D0   ! D/DKy D/DJy
                    SDERP(38) =  1.D0   ! D/DIy D/DJy
    FDERP(3) = - Z2DIST  ! D/DJz                     
                    SDERP(24) = -1.D0   ! D/DKz D/DJz
                    SDERP(48) =  1.D0   ! D/DIz D/DJz
    ! DERIVATIVES RESPECT TO ATOM K                  
    FDERP(4) = X1DIST + X2DIST ! D/DKx               
                    SDERP(13) =  2.D0   ! D/DKx D/DKx
                    SDERP(31) = -1.D0   ! D/DIx D/DKx
    FDERP(5) = Y1DIST + Y2DIST ! D/DKy               
                    SDERP(20) =  2.D0   ! D/DKy D/DKy
                    SDERP(41) = -1.D0   ! D/DIy D/DKy
    FDERP(6) = Z1DIST + Z2DIST ! D/DKz               
                    SDERP(27) =  2.D0   ! D/DKz D/DKz
                    SDERP(51) = -1.D0   ! D/DIz D/DKz
    ! DERIVATIVES RESPECT TO ATOM I
    FDERP(7) = - X1DIST        ! D/DIx
    FDERP(8) = - Y1DIST        ! D/DIy
    FDERP(9) = - Z1DIST        ! D/DIz

!   FIRST DERIVATIVES OF R1*R2

    ! DERIVATIVES RESPECT TO ATOM J  
                    ! SECOND DERIVATIVES OF FDERR1R2 ( ONLY NON ZERO ONES )  
    FDERR1R2(1) = - X1DIST * R2SUR1 ! D/DJx
                    ! WITH RESPECT TO ATOM J (J-J)
                    SDERR1R2( 1) =   R2SUR1 - X1X1DIST * R2SUR13                    ! D/DJx D/DJx
                    SDERR1R2( 2) =          - X1Y1DIST * R2SUR13                    ! D/DJy D/DJx
                    SDERR1R2( 3) =          - X1Z1DIST * R2SUR13                    ! D/DJz D/DJx
                    ! WITH RESPECT TO ATOM K (J-K)
                    SDERR1R2(10) = - R2SUR1 - X1X2DIST * SQ  + X1X1DIST * R2SUR13   ! D/DKx D/DJx
                    SDERR1R2(16) =          - X1Y2DIST * SQ  + X1Y1DIST * R2SUR13   ! D/DKy D/DJx 
                    SDERR1R2(22) =          - X1Z2DIST * SQ  + X1Z1DIST * R2SUR13   ! D/DKz D/DJx
                    ! WITH RESPECT TO ATOM I (J-I)
                    SDERR1R2(28) =            X1X2DIST * SQ                         ! D/DIx D/DJx
                    SDERR1R2(37) =            X1Y2DIST * SQ                         ! D/DIy D/DJx
                    SDERR1R2(46) =            X1Z2DIST * SQ                         ! D/DIz D/DJx
    FDERR1R2(2) = - Y1DIST * R2SUR1 ! D/DJy
                    ! WITH RESPECT TO ATOM J (J-J)
                    SDERR1R2( 4) =   SDERR1R2( 2)                                   ! D/DJx D/DJy
                    SDERR1R2( 5) =   R2SUR1 - Y1Y1DIST * R2SUR13                    ! D/DJy D/DJy
                    SDERR1R2( 6) =          - Y1Z1DIST * R2SUR13                    ! D/DJz D/DJy
                    ! WITH RESPECT TO ATOM K (J-K)
                    SDERR1R2(11) =          - Y1X2DIST * SQ  + X1Y1DIST * R2SUR13   ! D/DKx D/DJy
                    SDERR1R2(17) = - R2SUR1 - Y1Y2DIST * SQ  + Y1Y1DIST * R2SUR13   ! D/DKy D/DJy
                    SDERR1R2(23) =          - Y1Z2DIST * SQ  + Y1Z1DIST * R2SUR13   ! D/DKz D/DJy
                    ! WITH RESPECT TO ATOM I (J-I)
                    SDERR1R2(29) =            Y1X2DIST * SQ                         ! D/DIx D/DJy
                    SDERR1R2(38) =            Y1Y2DIST * SQ                         ! D/DIy D/DJy
                    SDERR1R2(47) =            Y1Z2DIST * SQ                         ! D/DIz D/DJy
    FDERR1R2(3) = - Z1DIST * R2SUR1 ! D/DJz
                    ! WITH RESPECT TO ATOM J (J-J)
                    SDERR1R2( 7) =   SDERR1R2( 3)                                   ! D/DJx D/DJz
                    SDERR1R2( 8) =   SDERR1R2( 6)                                   ! D/DJy D/DJz
                    SDERR1R2( 9) =   R2SUR1 - Z1Z1DIST * R2SUR13                    ! D/DJz D/DJz
                    ! WITH RESPECT TO ATOM K (J-K)
                    SDERR1R2(12) =          - Z1X2DIST * SQ  + X1Z1DIST * R2SUR13   ! D/DKx D/DJz
                    SDERR1R2(18) =          - Z1Y2DIST * SQ  + Y1Z1DIST * R2SUR13   ! D/DKy D/DJz
                    SDERR1R2(24) = - R2SUR1 - Z1Z2DIST * SQ  + Z1Z1DIST * R2SUR13   ! D/DKz D/DJz
                    ! WITH RESPECT TO ATOM I (J-I)
                    SDERR1R2(30) =            Z1X2DIST * SQ                         ! D/DIx D/DJz
                    SDERR1R2(39) =            Z1Y2DIST * SQ                         ! D/DIy D/DJz
                    SDERR1R2(48) =            Z1Z2DIST * SQ                         ! D/DIz D/DJz
    ! DERIVATIVES RESPECT TO ATOM K
    FDERR1R2(4) = X1DIST * R2SUR1 + X2DIST * R1SUR2 ! D/DKx
                    ! WITH RESPECT TO ATOM K (K-K)
                    SDERR1R2(13) =   R2SUR1 +X1X2DIST*SQ*2.D0- X1X1DIST * R2SUR13   & 
                                   + R1SUR2                  - X2X2DIST * R1SUR23   ! D/DKx D/DKx
                    SDERR1R2(14) =            X1Y2DIST * SQ  - X1Y1DIST * R2SUR13   &
                                            + Y1X2DIST * SQ  - X2Y2DIST * R1SUR23   ! D/DKy D/DKx
                    SDERR1R2(15) =            X1Z2DIST * SQ  - X1Z1DIST * R2SUR13   &
                                            + Z1X2DIST * SQ  - X2Z2DIST * R1SUR23   ! D/DKz D/DKx
                    ! WITH RESPECT TO ATOM I (K-I)
                    SDERR1R2(31) =          - X1X2DIST * SQ                         &
                                   - R1SUR2                  + X2X2DIST * R1SUR23   ! D/DIx D/DKx
                    SDERR1R2(40) =          - X1Y2DIST * SQ  + X2Y2DIST * R1SUR23   ! D/DIy D/DKx
                    SDERR1R2(49) =          - X1Z2DIST * SQ  + X2Z2DIST * R1SUR23   ! D/DIz D/DKx
    FDERR1R2(5) = Y1DIST * R2SUR1 + Y2DIST * R1SUR2 ! D/DKy
                    ! WITH RESPECT TO ATOM K (K-K)
                    SDERR1R2(19) =   SDERR1R2(14)                                   ! D/DKx D/DKy
                    SDERR1R2(20) =   R2SUR1 +Y1Y2DIST*SQ*2.D0- Y1Y1DIST * R2SUR13   &
                                   + R1SUR2                  - Y2Y2DIST * R1SUR23   ! D/DKy D/DKy
                    SDERR1R2(21) =            Y1Z2DIST * SQ  - Y1Z1DIST * R2SUR13   &
                                            + Z1Y2DIST * SQ  - Y2Z2DIST * R1SUR23   ! D/DKz D/DKy
                    ! WITH RESPECT TO ATOM I (K-I)
                    SDERR1R2(32) =          - Y1X2DIST * SQ  + X2Y2DIST * R1SUR23   ! D/DIx D/DKy
                    SDERR1R2(41) =          - Y1Y2DIST * SQ                         &
                                   - R1SUR2                  + Y2Y2DIST * R1SUR23   ! D/DIy D/DKy
                    SDERR1R2(50) =          - Y1Z2DIST * SQ  + Y2Z2DIST * R1SUR23   ! D/DIz D/DKy
    FDERR1R2(6) = Z1DIST * R2SUR1 + Z2DIST * R1SUR2 ! D/DKz
                    ! WITH RESPECT TO ATOM K (K-K)
                    SDERR1R2(25) =   SDERR1R2(15)                                   ! D/DKx D/DKz
                    SDERR1R2(26) =   SDERR1R2(21)                                   ! D/DKy D/DKz
                    SDERR1R2(27) =   R2SUR1 +Z1Z2DIST*SQ*2.D0- Z1Z1DIST * R2SUR13   &
                                   + R1SUR2                  - Z2Z2DIST * R1SUR23   ! D/DKz D/DKz
                    ! WITH RESPECT TO ATOM I (K-I)
                    SDERR1R2(33) =          - Z1X2DIST * SQ  + X2Z2DIST * R1SUR23   ! D/DIx D/DKz
                    SDERR1R2(42) =          - Z1Y2DIST * SQ  + Y2Z2DIST * R1SUR23   ! D/DIy D/DKz
                    SDERR1R2(51) =          - Z1Z2DIST * SQ                         &
                                   - R1SUR2                  + Z2Z2DIST * R1SUR23   ! D/DIz D/DKz
    ! DERIVATIVES RESPECT TO ATOM I    
    FDERR1R2(7) = - X2DIST * R1SUR2 ! D/DIx
                    ! WITH RESPECT TO ATOM I (I-I)
                    SDERR1R2(34) =   R1SUR2 - X2X2DIST * R1SUR23                    ! D/DIx D/DIx
                    SDERR1R2(35) =          - X2Y2DIST * R1SUR23                    ! D/DIy D/DIx
                    SDERR1R2(36) =          - X2Z2DIST * R1SUR23                    ! D/DIz D/DIx
    FDERR1R2(8) = - Y2DIST * R1SUR2 ! D/DIy
                    ! WITH RESPECT TO ATOM I (I-I)
                    SDERR1R2(43) =   SDERR1R2(35)                                   ! D/DIx D/DIy
                    SDERR1R2(44) =   R1SUR2 - Y2Y2DIST * R1SUR23                    ! D/DIy D/DIy
                    SDERR1R2(45) =          - Y2Z2DIST * R1SUR23                    ! D/DIz D/DIy
    FDERR1R2(9) = - Z2DIST * R1SUR2 ! D/DIz
                    ! WITH RESPECT TO ATOM I (I-I)
                    SDERR1R2(52) =   SDERR1R2(36)                                   ! D/DIx D/DIz
                    SDERR1R2(53) =   SDERR1R2(45)                                   ! D/DIy D/DIz
                    SDERR1R2(54) =   R1SUR2 - Z2Z2DIST * R1SUR23                    ! D/DIz D/DIz

!   SECOND DERIVATIVES OF COS(ANG) RESPECT TO ATOM  

    IFDERCOS = 0
    ISDERCOS = 0
    DO KX1 = 1, 3          ! Jx, Jy, Jz  
       IFDERCOS = IFDERCOS + 1
       FDERCOS(IFDERCOS) = ( FDERP(KX1) * R1R2 - FDERR1R2(KX1) * PRODSC ) * SQS
       DO KX2 = 1, 3       ! Jx, Jy, Jz  
          ISDERCOS = ISDERCOS + 1
          SDERCOS(ISDERCOS) = LOCFUN(ISDERCOS,KX1,KX2)
       END DO
    END DO

    DO KX1 = 4, 6          ! Kx, Ky, Kz            
       IFDERCOS = IFDERCOS + 1
       FDERCOS(IFDERCOS) = ( FDERP(KX1) * R1R2 - FDERR1R2(KX1) * PRODSC ) * SQS       
       DO KX2 = 1, 6       ! Jx, Jy, Jz, Kx, Ky, Kz 
          ISDERCOS = ISDERCOS + 1
          SDERCOS(ISDERCOS) = LOCFUN(ISDERCOS,KX1,KX2)
       END DO
    END DO

    DO KX1 = 7, 9          ! Ix, Iy, Iz  
       IFDERCOS = IFDERCOS + 1
       FDERCOS(IFDERCOS) = ( FDERP(KX1) * R1R2 - FDERR1R2(KX1) * PRODSC ) * SQS
       DO KX2 = 1, 9       ! Jx, Jy, Jz, Kx, Ky, Kz, Ix, Iy, Iz  
          ISDERCOS = ISDERCOS + 1
          SDERCOS(ISDERCOS) = LOCFUN(ISDERCOS,KX1,KX2)
       END DO
    END DO

       
!   STORE THE FORCES IN F

    ISDERCOS = 0
    INDF     = 0
    DO KX1 = 1, 9
       DO KX2 = 1, KX1
          INDF     = INDF   + 1
          ISDERCOS = ISDERCOS + 1
          TERM1    = FMOD1  * FDERCOS(KX1) * FDERCOS(KX2)  
          TERM2    = SDERCOS(ISDERCOS)
          F(INDF)  = FMOD2 * ( TERM1 + TERM2 )
       END DO
       MODR     = MOD( KX1 , 3 )
       IF ( MODR .EQ. 0 ) MODR = 3
       ISDERCOS = ISDERCOS + ( 3 - MODR )   
    END DO

! DEBUG SECTION
#ifdef DEBUG_DERI2
    CALL CHECK_FIRST_DERIVATIVES(LABEL='ANGLE',I1=J,I2=K,I3=I,COORD=C,DER1=FDERCOS,FAC=FMOD2*FACT)
    CALL CHECK_SECOND_DERIVATIVES(LABEL='ANGLE',I1=J,I2=K,I3=I,COORD=C,DER2=F,FAC=1.d0)
#endif
    RETURN

  CONTAINS 
      ! Local Function
    DOUBLE PRECISION FUNCTION LOCFUN ( I3, I1, I2 )  
      IMPLICIT NONE
      INTEGER :: I1, I2, I3
       LOCFUN = ( ( SDERP(I3) * R1R2                +       & 
                    FDERP(I1) * FDERR1R2(I2)        -       &
                    FDERP(I2) * FDERR1R2(I1)        -       &
                    PRODSC * SDERR1R2(I3)                   &         
                  ) * R1R2 - 2.D0 * FDERR1R2(I2)    *       &
                  ( FDERP(I1) * R1R2                -       &
                    PRODSC * FDERR1R2(I1)                   &
                  )                                         &
                ) / R13R23
      RETURN
    END FUNCTION LOCFUN
  END SUBROUTINE DERI_2_ANGLE



  SUBROUTINE DERI_2_TORSION( K, J, M, L, C, F, FACT )
    USE GEO_OBJECTS, ONLY :  CROSSPRODUCT, &
                             DDOT
! CALCOLA LE DERIVATE PER UN TORSION IN ( K, J, M, L) ORDER!!
! F IS A DIAGONAL MATRIX OF DIMENSION : 78
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: K, J, M, L
    DOUBLE PRECISION, INTENT(IN) :: C(*), FACT
    DOUBLE PRECISION, INTENT(INOUT) :: F(*)
    ! Local Variables
    INTEGER :: MM, JJ, KK, LL, IK, I
    INTEGER :: IFDCOS, ISDCOS, KX1, KX2, INDF
    DOUBLE PRECISION :: V1NORM, V2NORM,   PRDT,      COSPHI,    TOR
    DOUBLE PRECISION :: TNORM,  TNORMINV, TNORMINV2, TNORMINV3, MODR
    DOUBLE PRECISION :: FMOD1,  FMOD2,    TERM1,     TERM2
    DOUBLE PRECISION :: PI, FD1, FD2, FD13, FD23
    DOUBLE PRECISION :: PRDT_SIGN
    DOUBLE PRECISION, DIMENSION(3)  :: X, Y, Z, V1, V2
    DOUBLE PRECISION, DIMENSION(9)  :: DDOTD,  DNORM,   DFD1,   DFD2,   DV1N,   DV2N
    DOUBLE PRECISION, DIMENSION(54) :: DDOTDS, DDNORM,  DDV1N,  DDV2N
    DOUBLE PRECISION, DIMENSION(12) :: FDDOT,  FDNORM,  FDCOS
    DOUBLE PRECISION, DIMENSION(90) :: SDDOT,  SDNORM,  SDCOS

    ! Pi value...
    PI=4.D0*ATAN(1.D0)
    
    KK=(K-1)*3
    JJ=(J-1)*3
    MM=(M-1)*3
    LL=(L-1)*3

    X(1)=C(JJ+1)-C(KK+1)
    X(2)=C(JJ+2)-C(KK+2)    
    X(3)=C(JJ+3)-C(KK+3)    
    
    Y(1)=C(LL+1)-C(MM+1)
    Y(2)=C(LL+2)-C(MM+2)    
    Y(3)=C(LL+3)-C(MM+3)

    Z(1)=C(MM+1)-C(JJ+1)
    Z(2)=C(MM+2)-C(JJ+2)    
    Z(3)=C(MM+3)-C(JJ+3)    
    
    CALL CROSSPRODUCT(X,Z,V1) ! V1(1) = X(2) * Z(3) - Z(2) * X(3)
                              ! V1(2) = X(3) * Z(1) - Z(3) * X(1)
                              ! V1(3) = X(1) * Z(2) - Z(1) * X(2)
                              !
    CALL CROSSPRODUCT(Z,Y,V2) ! V2(1) = Z(2) * Y(3) - Y(2) * Z(3)
                              ! V2(2) = Z(3) * Y(1) - Y(3) * Z(1)
                              ! V2(3) = Z(1) * Y(2) - Y(1) * Z(2)
                              !  
    V1NORM = SQRT( DDOT(V1,V1,3) )
    V2NORM = SQRT( DDOT(V2,V2,3) )

    TNORM  = V1NORM * V2NORM
    PRDT   = DDOT(V1, V2, 3)  ! = V1(1) * V2(1) + V1(2) * V2(2) + V1(3) * V2(3) 
                              ! = (X(2) * Z(3) - Z(2) * X(3))*(Z(2) * Y(3) - Y(2) * Z(3)) +
                              !   (X(3) * Z(1) - Z(3) * X(1))*(Z(3) * Y(1) - Y(3) * Z(1)) +
                              !   (X(1) * Z(2) - Z(1) * X(2))*(Z(1) * Y(2) - Y(1) * Z(2))
    PRDT_SIGN = DDOT(X,V2,3)
    COSPHI = PRDT / TNORM
    IF (ABS(COSPHI).GT.1.D0) COSPHI = COSPHI / ABS(COSPHI)
    TOR    = ACOS ( COSPHI )    

    DO I = 1, 9
       DDOTD ( I)  = 0.D0
       DNORM ( I)  = 0.D0
       DFD1  ( I)  = 0.D0
       DFD2  ( I)  = 0.D0
       DV1N  ( I)  = 0.D0
       DV2N  ( I)  = 0.D0
    END DO
    DO I = 1, 12
       FDDOT ( I)  = 0.D0
       FDNORM( I)  = 0.D0
       FDCOS ( I)  = 0.D0
    END DO
    DO I = 1, 54
       DDOTDS( I)  = 0.D0
       DDNORM( I)  = 0.D0
       DDV1N ( I)  = 0.D0
       DDV2N ( I)  = 0.D0
    END DO
    DO I = 1, 90
       SDDOT ( I)  = 0.D0
       SDNORM( I)  = 0.D0
       SDCOS ( I)  = 0.D0
    END DO

!----------------------------------------------------!
!                                                    !
!   CALCOLO DERIVATE DELL'ANGOLO DI TORSIONE <TOR>   !
!                                                    !
!----------------------------------------------------!

!---------------------------!
!   DDOT PART DERIVATIVES   !
!---------------------------!

    DDOTD( 1) =   Z (2) * V2 (3) - Z (3) * V2 (2)                                     ! D/DX(1)
              DDOTDS( 4) =  - Z (2) * Z (2) - Z (3) * Z (3)      ! D/DY(1) ! D/DX(1)
              DDOTDS( 5) =    Z (2) * Z (1)                      ! D/DY(2) ! D/DX(1)
              DDOTDS( 6) =    Z (3) * Z (1)                      ! D/DY(3) ! D/DX(1)
              DDOTDS( 7) =    Z (2) * Y (2) + Z (3) * Y (3)      ! D/DZ(1) ! D/DX(1)
              DDOTDS( 8) =    V2(3)         - Z (2) * Y (1)      ! D/DZ(2) ! D/DX(1)
              DDOTDS( 9) =  - V2(2)         - Z (3) * Y (1)      ! D/DZ(3) ! D/DX(1)
    DDOTD( 2) = - Z (1) * V2 (3) + Z (3) * V2 (1)                                     ! D/DX(2)
              DDOTDS(13) =    Z (1) * Z (2)                      ! D/DY(1) ! D/DX(2)
              DDOTDS(14) =  - Z (1) * Z (1) - Z (3) * Z (3)      ! D/DY(2) ! D/DX(2)
              DDOTDS(15) =    Z (3) * Z (2)                      ! D/DY(3) ! D/DX(2)
              DDOTDS(16) =  - V2(3)         - Z (1) * Y (2)      ! D/DZ(1) ! D/DX(2)
              DDOTDS(17) =    Z (1) * Y (1) + Z (3) * Y (3)      ! D/DZ(2) ! D/DX(2)
              DDOTDS(18) =    V2(1)         - Z (3) * Y (2)      ! D/DZ(3) ! D/DX(2)
    DDOTD( 3) =   Z (1) * V2 (2) - Z (2) * V2 (1)                                     ! D/DX(3)
              DDOTDS(22) =    Z (1) * Z (3)                      ! D/DY(1) ! D/DX(3)
              DDOTDS(23) =    Z (2) * Z (3)                      ! D/DY(2) ! D/DX(3)
              DDOTDS(24) =  - Z (1) * Z (1) - Z (2) * Z (2)      ! D/DY(3) ! D/DX(3)
              DDOTDS(25) =    V2(2)         - Z (1) * Y (3)      ! D/DZ(1) ! D/DX(3)
              DDOTDS(26) =  - V2(1)         - Z (2) * Y (3)      ! D/DZ(2) ! D/DX(3)
              DDOTDS(27) =    Z (1) * Y (1) + Z (2) * Y (2)      ! D/DZ(3) ! D/DX(3)
    DDOTD( 4) = - Z (2) * V1 (3) + Z (3) * V1 (2)                                     ! D/DY(1)
              DDOTDS(31) =    Z (2) * X (2) + X (3) * Z (3)      ! D/DZ(1) ! D/DY(1)
              DDOTDS(32) =  - V1(3)         - X (1) * Z (2)      ! D/DZ(2) ! D/DY(1)
              DDOTDS(33) =    V1(2)         - X (1) * Z (3)      ! D/DZ(3) ! D/DY(1)
    DDOTD( 5) =   Z (1) * V1 (3) - Z (3) * V1 (1)                                     ! D/DY(2)
              DDOTDS(37) =    V1(3)         - X (2) * Z (1)      ! D/DZ(1) ! D/DY(2)
              DDOTDS(38) =    X (1) * Z (1) + X (3) * Z (3)      ! D/DZ(2) ! D/DY(2)
              DDOTDS(39) =  - V1(1)         - X (2) * Z (3)      ! D/DZ(3) ! D/DY(2)
    DDOTD( 6) = - Z (1) * V1 (2) + Z (2) * V1 (1)                                     ! D/DY(3)
              DDOTDS(43) =  - V1(2)         - X (3) * Z (1)      ! D/DZ(1) ! D/DY(3)
              DDOTDS(44) =    V1(1)         - X (3) * Z (2)      ! D/DZ(2) ! D/DY(3)
              DDOTDS(45) =    X (1) * Z (1) + X (2) * Z (2)      ! D/DZ(3) ! D/DY(3)
    DDOTD( 7) =   Y (2) * V1 (3) - X (2) * V2 (3) - Y (3) * V1 (2) + X (3) * V2 (2)   ! D/DZ(1)
              DDOTDS(46) =  - 2.D0  * X (2) * Y (2) - 2.D0 * Y (3) * X (3)  ! D/DZ(1) ! D/DZ(1)
              DDOTDS(47) =    Y (2) * X (1) + X (2) * Y (1)                 ! D/DZ(2) ! D/DZ(1)
              DDOTDS(48) =    X (1) * Y (3) + Y (1) * X (3)                 ! D/DZ(3) ! D/DZ(1)
    DDOTD( 8) = - Y (1) * V1 (3) + X (1) * V2 (3) + Y (3) * V1 (1) - X (3) * V2 (1)   ! D/DZ(2)
              DDOTDS(49) =    DDOTDS(47)                                    ! D/DZ(1) ! D/DZ(2)
              DDOTDS(50) =  - 2.D0  * X (1) * Y (1) - 2.D0 * X (3) * Y (3)  ! D/DZ(2) ! D/DZ(2)
              DDOTDS(51) =    X (2) * Y (3) + X (3) * Y (2)                 ! D/DZ(3) ! D/DZ(2)
    DDOTD( 9) =   Y (1) * V1 (2) - X (1) * V2 (2) - Y (2) * V1 (1) + X (2) * V2 (1)   ! D/DZ(3)  
              DDOTDS(52) =    DDOTDS(48)                                    ! D/DZ(1) ! D/DZ(3)
              DDOTDS(53) =    DDOTDS(51)                                    ! D/DZ(2) ! D/DZ(3)
              DDOTDS(54) =  - 2.D0  * X (1) * Y (1) - 2.D0 * X (2) * Y (2)  ! D/DZ(3) ! D/DZ(3)
    ! FDDOT DERIVATIVES : First Derivatives of DDOT
              ! SDDOT : Second Derivatives of DDOT
    FDDOT( 1) =  - DDOTD( 1)                      ! D/DKx
              !     D/DX(1)
              SDDOT ( 1) =     DDOTDS( 1)                                         ! D/DKx  D/DKx
              SDDOT ( 2) =     DDOTDS( 2)                                         ! D/DKy  D/DKx
              SDDOT ( 3) =     DDOTDS( 3)                                         ! D/DKz  D/DKx
              SDDOT (10) =     DDOTDS( 7)!- DDOTDS( 1)                            ! D/DJx  D/DKx
              SDDOT (16) =     DDOTDS( 8)!- DDOTDS( 2)                            ! D/DJy  D/DKx
              SDDOT (22) =     DDOTDS( 9)!- DDOTDS( 3)                            ! D/DJz  D/DKx
              SDDOT (28) =     DDOTDS( 4) - DDOTDS( 7)                            ! D/DMx  D/DKx
              SDDOT (37) =     DDOTDS( 5) - DDOTDS( 8)                            ! D/DMy  D/DKx
              SDDOT (46) =     DDOTDS( 6) - DDOTDS( 9)                            ! D/DMz  D/DKx
              SDDOT (55) =   - DDOTDS( 4)                                         ! D/DLx  D/DKx
              SDDOT (67) =   - DDOTDS( 5)                                         ! D/DLy  D/DKx
              SDDOT (79) =   - DDOTDS( 6)                                         ! D/DLz  D/DKx
    FDDOT( 2) =  - DDOTD( 2)                      ! D/DKy                        
              !     D/DX(2)
              SDDOT ( 4) =     SDDOT( 2)                                          ! D/DKx  D/DKy
              SDDOT ( 5) =     DDOTDS(11)                                         ! D/DKy  D/DKy
              SDDOT ( 6) =     DDOTDS(12)                                         ! D/DKz  D/DKy
              SDDOT (11) =     DDOTDS(16)!- DDOTDS( 2)                            ! D/DJx  D/DKy
              SDDOT (17) =     DDOTDS(17)!- DDOTDS(11)                            ! D/DJy  D/DKy
              SDDOT (23) =     DDOTDS(18)!- DDOTDS(12)                            ! D/DJz  D/DKy
              SDDOT (29) =     DDOTDS(13) - DDOTDS(16)                            ! D/DMx  D/DKy
              SDDOT (38) =     DDOTDS(14) - DDOTDS(17)                            ! D/DMy  D/DKy
              SDDOT (47) =     DDOTDS(15) - DDOTDS(18)                            ! D/DMz  D/DKy
              SDDOT (56) =   - DDOTDS(13)                                         ! D/DLx  D/DKy
              SDDOT (68) =   - DDOTDS(14)                                         ! D/DLy  D/DKy
              SDDOT (80) =   - DDOTDS(15)                                         ! D/DLz  D/DKy
    FDDOT( 3) =  - DDOTD( 3)                      ! D/DKz                                 
              !     D/DX(3)
              SDDOT ( 7) =     SDDOT( 3)                                          ! D/DKx  D/DKz
              SDDOT ( 8) =     SDDOT( 6)                                          ! D/DKy  D/DKz
              SDDOT ( 9) =     DDOTDS(21)                                         ! D/DKz  D/DKz
              SDDOT (12) =     DDOTDS(25)!- DDOTDS( 3)                            ! D/DJx  D/DKz
              SDDOT (18) =     DDOTDS(26)!- DDOTDS(12)                            ! D/DJy  D/DKz
              SDDOT (24) =     DDOTDS(27)!- DDOTDS(21)                            ! D/DJz  D/DKz
              SDDOT (30) =     DDOTDS(22) - DDOTDS(25)                            ! D/DMx  D/DKz
              SDDOT (39) =     DDOTDS(23) - DDOTDS(26)                            ! D/DMy  D/DKz
              SDDOT (48) =     DDOTDS(24) - DDOTDS(27)                            ! D/DMz  D/DKz
              SDDOT (57) =   - DDOTDS(22)                                         ! D/DLx  D/DKz
              SDDOT (69) =   - DDOTDS(23)                                         ! D/DLy  D/DKz
              SDDOT (81) =   - DDOTDS(24)                                         ! D/DLz  D/DKz
    FDDOT( 4) =    DDOTD( 1) - DDOTD( 7)          ! D/DJx
              !     D/DX(1)     D/DZ(1)                                
              SDDOT (13) =   - DDOTDS( 7) - DDOTDS( 7) + DDOTDS(46)!+ DDOTDS( 1)  ! D/DJx  D/DJx
              SDDOT (14) =   - DDOTDS( 8) - DDOTDS(16) + DDOTDS(47)!+ DDOTDS( 2)  ! D/DJy  D/DJx
              SDDOT (15) =   - DDOTDS( 9) - DDOTDS(25) + DDOTDS(48)!+ DDOTDS( 3)  ! D/DJz  D/DJx
              SDDOT (31) =     DDOTDS( 7) - DDOTDS(46) - DDOTDS( 4) + DDOTDS(31)  ! D/DMx  D/DJx
              SDDOT (40) =     DDOTDS( 8) - DDOTDS(47) - DDOTDS( 5) + DDOTDS(37)  ! D/DMy  D/DJx
              SDDOT (49) =     DDOTDS( 9) - DDOTDS(48) - DDOTDS( 6) + DDOTDS(43)  ! D/DMz  D/DJx
              SDDOT (58) =     DDOTDS( 4) - DDOTDS(31)                            ! D/DLx  D/DJx
              SDDOT (70) =     DDOTDS( 5) - DDOTDS(37)                            ! D/DLy  D/DJx
              SDDOT (82) =     DDOTDS( 6) - DDOTDS(43)                            ! D/DLz  D/DJx
    FDDOT( 5) =    DDOTD( 2) - DDOTD( 8)          ! D/DJy
              !     D/DX(2)     D/DZ(2)
            ! SDDOT (19) =   - DDOTDS(16) - DDOTDS( 8) + DDOTDS(49)!+ DDOTDS(10)  ! D/DJx  D/DJy
              SDDOT (19) =     SDDOT (14)
              SDDOT (20) =   - DDOTDS(17) - DDOTDS(17) + DDOTDS(50)!+ DDOTDS(11)  ! D/DJy  D/DJy
              SDDOT (21) =   - DDOTDS(18) - DDOTDS(26) + DDOTDS(51)!+ DDOTDS(12)  ! D/DJz  D/DJy
              SDDOT (32) =   - DDOTDS(13) + DDOTDS(32) + DDOTDS(16) - DDOTDS(49)  ! D/DMx  D/DJy
              SDDOT (41) =   - DDOTDS(14) + DDOTDS(38) + DDOTDS(17) - DDOTDS(50)  ! D/DMy  D/DJy
              SDDOT (50) =   - DDOTDS(15) + DDOTDS(44) + DDOTDS(18) - DDOTDS(51)  ! D/DMz  D/DJy
              SDDOT (59) =     DDOTDS(13) - DDOTDS(32)                            ! D/DLx  D/DJy
              SDDOT (71) =     DDOTDS(14) - DDOTDS(38)                            ! D/DLy  D/DJy
              SDDOT (83) =     DDOTDS(15) - DDOTDS(44)                            ! D/DLz  D/DJy
    FDDOT( 6) =    DDOTD( 3) - DDOTD( 9)          ! D/DJz
              !     D/DX(3)     D/DZ(3)
            ! SDDOT (25) =   - DDOTDS(25) + DDOTDS(52) - DDOTDS( 9)!+ DDOTDS(19)  ! D/DJx  D/DJz
            ! SDDOT (26) =   - DDOTDS(26) + DDOTDS(53) - DDOTDS(18)!+ DDOTDS(20)  ! D/DJy  D/DJz
              SDDOT (25) =     SDDOT (15)
              SDDOT (26) =     SDDOT (21)
              SDDOT (27) =   - DDOTDS(27) + DDOTDS(54) - DDOTDS(27)!+ DDOTDS(21)  ! D/DJz  D/DJz
              SDDOT (33) =   - DDOTDS(22) + DDOTDS(33) + DDOTDS(25) - DDOTDS(52)  ! D/DMx  D/DJz
              SDDOT (42) =   - DDOTDS(23) + DDOTDS(39) + DDOTDS(26) - DDOTDS(53)  ! D/DMy  D/DJz
              SDDOT (51) =   - DDOTDS(24) + DDOTDS(45) + DDOTDS(27) - DDOTDS(54)  ! D/DMz  D/DJz
              SDDOT (60) =     DDOTDS(22) - DDOTDS(33)                            ! D/DLx  D/DJz
              SDDOT (72) =     DDOTDS(23) - DDOTDS(39)                            ! D/DLy  D/DJz
              SDDOT (84) =     DDOTDS(24) - DDOTDS(45)                            ! D/DLz  D/DJz
    FDDOT( 7) =    DDOTD( 7) - DDOTD( 4)          ! D/DMx
              !     D/DZ(1)     D/DY(1)
              SDDOT (34) =   - DDOTDS(31) + DDOTDS(46) - DDOTDS(31)!+ DDOTDS(28)  ! D/DMx  D/DMx
              SDDOT (35) =   - DDOTDS(37) + DDOTDS(47) - DDOTDS(32)!+ DDOTDS(29)  ! D/DMy  D/DMx
              SDDOT (36) =   - DDOTDS(43) + DDOTDS(48) - DDOTDS(33)!+ DDOTDS(30)  ! D/DMz  D/DMx
              SDDOT (61) =     DDOTDS(31)!- DDOTDS(28)                            ! D/DLx  D/DMx
              SDDOT (73) =     DDOTDS(37)!- DDOTDS(29)                            ! D/DLy  D/DMx
              SDDOT (85) =     DDOTDS(43)!- DDOTDS(30)                            ! D/DLz  D/DMx
    FDDOT( 8) =    DDOTD( 8) - DDOTD( 5)          ! D/DMy
              !     D/DZ(2)     D/DY(2)
            ! SDDOT (43) =   - DDOTDS(32) + DDOTDS(49) - DDOTDS(37)!+ DDOTDS(34)  ! D/DMx  D/DMy
              SDDOT (43) =     SDDOT (35)
              SDDOT (44) =   - DDOTDS(38) + DDOTDS(50) - DDOTDS(38)!+ DDOTDS(35)  ! D/DMy  D/DMy
              SDDOT (45) =   - DDOTDS(44) + DDOTDS(51) - DDOTDS(39)!+ DDOTDS(36)  ! D/DMz  D/DMy
              SDDOT (62) =     DDOTDS(32)!- DDOTDS(34)                            ! D/DLx  D/DMy
              SDDOT (74) =     DDOTDS(38)!- DDOTDS(35)                            ! D/DLy  D/DMy
              SDDOT (86) =     DDOTDS(44)!- DDOTDS(36)                            ! D/DLz  D/DMy
    FDDOT( 9) =    DDOTD( 9) - DDOTD( 6)          ! D/DMz
              !     D/DZ(3)     D/DY(3)
            ! SDDOT (52) =   - DDOTDS(33) + DDOTDS(52) - DDOTDS(43)!+ DDOTDS(40)  ! D/DMx  D/DMz
            ! SDDOT (53) =   - DDOTDS(39) + DDOTDS(53) - DDOTDS(44)!+ DDOTDS(41)  ! D/DMy  D/DMz
              SDDOT (52) =     SDDOT (36)
              SDDOT (53) =     SDDOT (45)
              SDDOT (54) =   - DDOTDS(45) + DDOTDS(54) - DDOTDS(45)!+ DDOTDS(42)  ! D/DMz  D/DMz
              SDDOT (63) =     DDOTDS(33)!- DDOTDS(40)                            ! D/DLx  D/DMz
              SDDOT (75) =     DDOTDS(39)!- DDOTDS(41)                            ! D/DLy  D/DMz
              SDDOT (87) =     DDOTDS(45)!- DDOTDS(42)                            ! D/DLz  D/DMz
    FDDOT(10) =    DDOTD( 4)                      ! D/DLx
              !     D/DY(1)
              SDDOT (64) =     DDOTDS(28)                                         ! D/DLx  D/DLx
              SDDOT (65) =     DDOTDS(29)                                         ! D/DLy  D/DLx
              SDDOT (66) =     DDOTDS(30)                                         ! D/DLz  D/DLx
    FDDOT(11) =    DDOTD( 5)                      ! D/DLy
              !     D/DY(2)
              SDDOT (76) =     SDDOT (65)                                         ! D/DLx  D/DLy
              SDDOT (77) =     DDOTDS(35)                                         ! D/DLy  D/DLy
              SDDOT (78) =     DDOTDS(36)                                         ! D/DLz  D/DLy
    FDDOT(12) =    DDOTD( 6)                      ! D/DLz
              !     D/DY(3)
              SDDOT (88) =     SDDOT (66)                                         ! D/DLx  D/DLz
              SDDOT (89) =     SDDOT (78)                                         ! D/DLy  D/DLz
              SDDOT (90) =     DDOTDS(42)                                         ! D/DLz  D/DLz

!------------------------------!
!   TNORM PART DERIVATIVES     !
!------------------------------!

     TNORMINV  =   1.D0 / TNORM
     TNORMINV2 =   TNORMINV  * TNORMINV
     TNORMINV3 =   TNORMINV2 * TNORMINV
     FD13      = - V2NORM / (V1NORM)**3
     FD23      = - V1NORM / (V2NORM)**3

    ! V1NORM = SQRT( DDOT(V1,V1,3) )
     DV1N (1) =   Z (2) * V1 (3) - Z (3) * V1 (2)       ! D/DX(1)
              DDV1N( 1) =   Z (2) * Z (2) + Z (3) * Z (3)        ! D/DX(1)   D/DX(1) 
              DDV1N( 2) = - Z (1) * Z (2)                        ! D/DX(2)   D/DX(1) 
              DDV1N( 3) = - Z (1) * Z (3)                        ! D/DX(3)   D/DX(1) 
            ! DDV1N( 4) =   0.D0                                 ! D/DY(1)   D/DX(1) 
            ! DDV1N( 5) =   0.D0                                 ! D/DY(2)   D/DX(1) 
            ! DDV1N( 6) =   0.D0                                 ! D/DY(3)   D/DX(1) 
              DDV1N( 7) = - X (2) * Z (2) - X (3) * Z (3)        ! D/DZ(1)   D/DX(1) 
              DDV1N( 8) =   V1(3)         + X (1) * Z (2)        ! D/DZ(2)   D/DX(1) 
              DDV1N( 9) = - V1(2)         + X (1) * Z (3)        ! D/DZ(3)   D/DX(1) 
     DV1N (2) = - Z (1) * V1 (3) + Z (3) * V1 (1)       ! D/DX(2)
            ! DDV1N(10) = - Z (1) * Z (2)                        ! D/DX(1)   D/DX(2) 
              DDV1N(10) =   DDV1N( 2)
              DDV1N(11) =   Z (1) * Z (1) + Z (3) * Z (3)        ! D/DX(2)   D/DX(2) 
              DDV1N(12) = - Z (2) * Z (3)                        ! D/DX(3)   D/DX(2) 
            ! DDV1N(13) =   0.D0                                 ! D/DY(1)   D/DX(2) 
            ! DDV1N(14) =   0.D0                                 ! D/DY(2)   D/DX(2) 
            ! DDV1N(15) =   0.D0                                 ! D/DY(3)   D/DX(2) 
              DDV1N(16) = - V1(3)         + Z (1) * X (2)        ! D/DZ(1)   D/DX(2) 
              DDV1N(17) = - Z (1) * X (1) - X (3) * Z (3)        ! D/DZ(2)   D/DX(2) 
              DDV1N(18) =   V1(1)         + X (2) * Z (3)        ! D/DZ(3)   D/DX(2) 
     DV1N (3) =   Z (1) * V1 (2) - Z (2) * V1 (1)       ! D/DX(3)
            ! DDV1N(19) = - Z (1) * Z (3)                        ! D/DX(1)   D/DX(3) 
            ! DDV1N(20) = - Z (2) * Z (3)                        ! D/DX(2)   D/DX(3) 
              DDV1N(19) =   DDV1N( 3)
              DDV1N(20) =   DDV1N(12) 
              DDV1N(21) =   Z (1) * Z (1) + Z (2) * Z (2)        ! D/DX(3)   D/DX(3) 
            ! DDV1N(22) =   0.D0                                 ! D/DY(1)   D/DX(3) 
            ! DDV1N(23) =   0.D0                                 ! D/DY(2)   D/DX(3) 
            ! DDV1N(24) =   0.D0                                 ! D/DY(3)   D/DX(3) 
              DDV1N(25) =   V1(2)         + Z (1) * X (3)        ! D/DZ(1)   D/DX(3) 
              DDV1N(26) = - V1(1)         + X (3) * Z (2)        ! D/DZ(2)   D/DX(3) 
              DDV1N(27) = - X (1) * Z (1) - X (2) * Z (2)        ! D/DZ(3)   D/DX(3) 
     DV1N (7) =   X (3) * V1 (2) - X (2) * V1 (3)       ! D/DZ(1)
              DDV1N(46) =   X (3) * X (3) + X (2) * X (2)        ! D/DZ(1)   D/DZ(1) 
              DDV1N(47) = - X (1) * X (2)                        ! D/DZ(2)   D/DZ(1) 
              DDV1N(48) = - X (1) * X (3)                        ! D/DZ(3)   D/DZ(1)
     DV1N (8) =   X (1) * V1 (3) - X (3) * V1 (1)       ! D/DZ(2)
            ! DDV1N(49) = - X (1) * X (2)                        ! D/DZ(1)   D/DZ(2) 
              DDV1N(49) =   DDV1N(47)
              DDV1N(50) =   X (1) * X (1) + X (3) * X (3)        ! D/DZ(2)   D/DZ(2) 
              DDV1N(51) = - X (3) * X (2)                        ! D/DZ(3)   D/DZ(2)
     DV1N (9) = - X (1) * V1 (2) + X (2) * V1 (1)       ! D/DZ(3)
            ! DDV1N(52) = - X (1) * X (3)                        ! D/DZ(1)   D/DZ(3) 
            ! DDV1N(53) = - X (2) * X (3)                        ! D/DZ(2)   D/DZ(3) 
              DDV1N(52) =   DDV1N(48)
              DDV1N(53) =   DDV1N(51)
              DDV1N(54) =   X (1) * X (1) + X (2) * X (2)        ! D/DZ(3)   D/DZ(3)

    ! V2NORM = SQRT( DDOT(V2,V2,3) )
     DV2N (4) = - Z (2) * V2 (3) + Z (3) * V2 (2)       ! D/DY(1)
              DDV2N(28) =   Z (2) * Z (2) + Z (3) * Z (3)        ! D/DY(1)   D/DY(1) 
              DDV2N(29) = - Z (1) * Z (2)                        ! D/DY(2)   D/DY(1) 
              DDV2N(30) = - Z (1) * Z (3)                        ! D/DY(3)   D/DY(1) 
              DDV2N(31) = - Y (2) * Z (2) - Y (3) * Z (3)        ! D/DZ(1)   D/DY(1) 
              DDV2N(32) = - V2(3)         + Y (1) * Z (2)        ! D/DZ(2)   D/DY(1) 
              DDV2N(33) =   V2(2)         + Y (1) * Z (3)        ! D/DZ(3)   D/DY(1) 
     DV2N (5) =   Z (1) * V2 (3) - Z (3) * V2 (1)       ! D/DY(2)
            ! DDV2N(34) = - Z (1) * Z (2)                        ! D/DY(1)   D/DY(2) 
              DDV2N(34) =   DDV2N(29)
              DDV2N(35) =   Z (1) * Z (1) + Z (3) * Z (3)        ! D/DY(2)   D/DY(2) 
              DDV2N(36) = - Z (2) * Z (3)                        ! D/DY(3)   D/DY(2) 
              DDV2N(37) =   V2(3)         + Y (2) * Z (1)        ! D/DZ(1)   D/DY(2) 
              DDV2N(38) = - Y (1) * Z (1) - Y (3) * Z (3)        ! D/DZ(2)   D/DY(2) 
              DDV2N(39) = - V2(1)         + Y (2) * Z (3)        ! D/DZ(3)   D/DY(2) 
     DV2N (6) = - Z (1) * V2 (2) + Z (2) * V2 (1)       ! D/DY(3)
            ! DDV2N(40) = - Z (1) * Z (3)                        ! D/DY(1)   D/DY(3) 
            ! DDV2N(41) = - Z (2) * Z (3)                        ! D/DY(2)   D/DY(3) 
              DDV2N(40) =   DDV2N(30)
              DDV2N(41) =   DDV2N(36)
              DDV2N(42) =   Z (1) * Z (1) + Z (2) * Z (2)        ! D/DY(3)   D/DY(3) 
              DDV2N(43) = - V2(2)         + Y (3) * Z (1)        ! D/DZ(1)   D/DY(3) 
              DDV2N(44) =   V2(1)         + Y (3) * Z (2)        ! D/DZ(2)   D/DY(3) 
              DDV2N(45) = - Y (1) * Z (1) - Y (2) * Z (2)        ! D/DZ(3)   D/DY(3) 
     DV2N (7) =   Y (2) * V2 (3) - Y (3) * V2 (2)       ! D/DZ(1)
              DDV2N(46) =   Y (2) * Y (2) + Y (3) * Y (3)        ! D/DZ(1)   D/DZ(1) 
              DDV2N(47) = - Y (1) * Y (2)                        ! D/DZ(2)   D/DZ(1) 
              DDV2N(48) = - Y (1) * Y (3)                        ! D/DZ(3)   D/DZ(1) 
     DV2N (8) = - Y (1) * V2 (3) + Y (3) * V2 (1)       ! D/DZ(2)
            ! DDV2N(49) = - Y (1) * Y (2)                        ! D/DZ(1)   D/DZ(2) 
              DDV2N(49) =   DDV2N(47)
              DDV2N(50) =   Y (1) * Y (1) + Y (3) * Y (3)        ! D/DZ(2)   D/DZ(2) 
              DDV2N(51) = - Y (2) * Y (3)                        ! D/DZ(3)   D/DZ(2) 
     DV2N (9) =   Y (1) * V2 (2) - Y (2) * V2 (1)       ! D/DZ(3)
            ! DDV2N(52) = - Y (1) * Y (3)                        ! D/DZ(1)   D/DZ(3) 
            ! DDV2N(53) = - Y (2) * Y (3)                        ! D/DZ(2)   D/DZ(3) 
              DDV2N(52) =   DDV2N(48)
              DDV2N(53) =   DDV2N(51)
              DDV2N(54) =   Y (1) * Y (1) + Y (2) * Y (2)        ! D/DZ(3)   D/DZ(3) 



    FD1       =  V2NORM  /  V1NORM 
     DFD1 (1) =  FD13      * DV1N (1)                                 ! D/DX(1)
     DFD1 (2) =  FD13      * DV1N (2)                                 ! D/DX(2)
     DFD1 (3) =  FD13      * DV1N (3)                                 ! D/DX(3)
     DFD1 (4) =  TNORMINV  * DV2N (4)                                 ! D/DY(1)
     DFD1 (5) =  TNORMINV  * DV2N (5)                                 ! D/DY(2)
     DFD1 (6) =  TNORMINV  * DV2N (6)                                 ! D/DY(3)
     DFD1 (7) =  FD13      * DV1N (7)   +   TNORMINV  * DV2N (7)      ! D/DZ(1)
     DFD1 (8) =  FD13      * DV1N (8)   +   TNORMINV  * DV2N (8)      ! D/DZ(2)
     DFD1 (9) =  FD13      * DV1N (9)   +   TNORMINV  * DV2N (9)      ! D/DZ(3)
    FD2       =  V1NORM  /  V2NORM 
     DFD2 (1) =  TNORMINV  * DV1N (1)                                 ! D/DX(1)
     DFD2 (2) =  TNORMINV  * DV1N (2)                                 ! D/DX(2)
     DFD2 (3) =  TNORMINV  * DV1N (3)                                 ! D/DX(3)
     DFD2 (4) =  FD23      * DV2N (4)                                 ! D/DY(1)
     DFD2 (5) =  FD23      * DV2N (5)                                 ! D/DY(2)
     DFD2 (6) =  FD23      * DV2N (6)                                 ! D/DY(3)
     DFD2 (7) =  FD23      * DV2N (7)   +   TNORMINV  * DV1N (7)      ! D/DZ(1)
     DFD2 (8) =  FD23      * DV2N (8)   +   TNORMINV  * DV1N (8)      ! D/DZ(2)
     DFD2 (9) =  FD23      * DV2N (9)   +   TNORMINV  * DV1N (9)      ! D/DZ(3)

    DNORM(1)  =  DV1N (1) * FD1                                              ! D/DX(1)
              DDNORM ( 1) = DV1N( 1) * DFD1( 1) + DDV1N( 1) * FD1  ! D/DX(1) ! D/DX(1)
              DDNORM ( 2) = DV1N( 1) * DFD1( 2) + DDV1N( 2) * FD1  ! D/DX(2) ! D/DX(1)
              DDNORM ( 3) = DV1N( 1) * DFD1( 3) + DDV1N( 3) * FD1  ! D/DX(3) ! D/DX(1)
              DDNORM ( 4) = DV1N( 1) * DFD1( 4)!+ DDV1N( 4) * FD1  ! D/DY(1) ! D/DX(1)
              DDNORM ( 5) = DV1N( 1) * DFD1( 5)!+ DDV1N( 5) * FD1  ! D/DY(2) ! D/DX(1)
              DDNORM ( 6) = DV1N( 1) * DFD1( 6)!+ DDV1N( 6) * FD1  ! D/DY(3) ! D/DX(1)
              DDNORM ( 7) = DV1N( 1) * DFD1( 7) + DDV1N( 7) * FD1  ! D/DZ(1) ! D/DX(1)
              DDNORM ( 8) = DV1N( 1) * DFD1( 8) + DDV1N( 8) * FD1  ! D/DZ(2) ! D/DX(1)
              DDNORM ( 9) = DV1N( 1) * DFD1( 9) + DDV1N( 9) * FD1  ! D/DZ(3) ! D/DX(1)
    DNORM(2)  =  DV1N (2) * FD1                                              ! D/DX(2)
            ! DDNORM (10) = DV1N( 2) * DFD1( 1) + DDV1N(10) * FD1  ! D/DX(1) ! D/DX(2)
              DDNORM (10) = DDNORM ( 2)
              DDNORM (11) = DV1N( 2) * DFD1( 2) + DDV1N(11) * FD1  ! D/DX(2) ! D/DX(2) 
              DDNORM (12) = DV1N( 2) * DFD1( 3) + DDV1N(12) * FD1  ! D/DX(3) ! D/DX(2) 
              DDNORM (13) = DV1N( 2) * DFD1( 4)!+ DDV1N(13) * FD1  ! D/DY(1) ! D/DX(2) 
              DDNORM (14) = DV1N( 2) * DFD1( 5)!+ DDV1N(14) * FD1  ! D/DY(2) ! D/DX(2) 
              DDNORM (15) = DV1N( 2) * DFD1( 6)!+ DDV1N(15) * FD1  ! D/DY(3) ! D/DX(2) 
              DDNORM (16) = DV1N( 2) * DFD1( 7) + DDV1N(16) * FD1  ! D/DZ(1) ! D/DX(2) 
              DDNORM (17) = DV1N( 2) * DFD1( 8) + DDV1N(17) * FD1  ! D/DZ(2) ! D/DX(2) 
              DDNORM (18) = DV1N( 2) * DFD1( 9) + DDV1N(18) * FD1  ! D/DZ(3) ! D/DX(2)  
    DNORM(3)  =  DV1N (3) * FD1                                              ! D/DX(3)
            ! DDNORM (19) = DV1N( 3) * DFD1( 1) + DDV1N(19) * FD1  ! D/DX(1) ! D/DX(3)
            ! DDNORM (20) = DV1N( 3) * DFD1( 2) + DDV1N(20) * FD1  ! D/DX(2) ! D/DX(3)
              DDNORM (19) = DDNORM ( 3)
              DDNORM (20) = DDNORM (12)
              DDNORM (21) = DV1N( 3) * DFD1( 3) + DDV1N(21) * FD1  ! D/DX(3) ! D/DX(3)
              DDNORM (22) = DV1N( 3) * DFD1( 4)!+ DDV1N(22) * FD1  ! D/DY(1) ! D/DX(3)
              DDNORM (23) = DV1N( 3) * DFD1( 5)!+ DDV1N(23) * FD1  ! D/DY(2) ! D/DX(3)
              DDNORM (24) = DV1N( 3) * DFD1( 6)!+ DDV1N(24) * FD1  ! D/DY(3) ! D/DX(3)
              DDNORM (25) = DV1N( 3) * DFD1( 7) + DDV1N(25) * FD1  ! D/DZ(1) ! D/DX(3)
              DDNORM (26) = DV1N( 3) * DFD1( 8) + DDV1N(26) * FD1  ! D/DZ(2) ! D/DX(3)
              DDNORM (27) = DV1N( 3) * DFD1( 9) + DDV1N(27) * FD1  ! D/DZ(3) ! D/DX(3)
    DNORM(4)  =  DV2N (4) * FD2                                              ! D/DY(1)
              DDNORM (28) = DV2N( 4) * DFD2( 4) + DDV2N(28) * FD2  ! D/DY(1) ! D/DY(1)
              DDNORM (29) = DV2N( 4) * DFD2( 5) + DDV2N(29) * FD2  ! D/DY(2) ! D/DY(1)
              DDNORM (30) = DV2N( 4) * DFD2( 6) + DDV2N(30) * FD2  ! D/DY(3) ! D/DY(1)
              DDNORM (31) = DV2N( 4) * DFD2( 7) + DDV2N(31) * FD2  ! D/DZ(1) ! D/DY(1)
              DDNORM (32) = DV2N( 4) * DFD2( 8) + DDV2N(32) * FD2  ! D/DZ(2) ! D/DY(1)
              DDNORM (33) = DV2N( 4) * DFD2( 9) + DDV2N(33) * FD2  ! D/DZ(3) ! D/DY(1)
    DNORM(5)  =  DV2N (5) * FD2                                              ! D/DY(2)
            ! DDNORM (34) = DV2N( 5) * DFD2( 4) + DDV2N(34) * FD2  ! D/DY(1) ! D/DY(2)
              DDNORM (34) = DDNORM (29)
              DDNORM (35) = DV2N( 5) * DFD2( 5) + DDV2N(35) * FD2  ! D/DY(2) ! D/DY(2)
              DDNORM (36) = DV2N( 5) * DFD2( 6) + DDV2N(36) * FD2  ! D/DY(3) ! D/DY(2)
              DDNORM (37) = DV2N( 5) * DFD2( 7) + DDV2N(37) * FD2  ! D/DZ(1) ! D/DY(2)
              DDNORM (38) = DV2N( 5) * DFD2( 8) + DDV2N(38) * FD2  ! D/DZ(2) ! D/DY(2)
              DDNORM (39) = DV2N( 5) * DFD2( 9) + DDV2N(39) * FD2  ! D/DZ(3) ! D/DY(2)
    DNORM(6)  =  DV2N (6) * FD2                                              ! D/DY(3)
            ! DDNORM (40) = DV2N( 6) * DFD2( 4) + DDV2N(40) * FD2  ! D/DY(1) ! D/DY(3)
            ! DDNORM (41) = DV2N( 6) * DFD2( 5) + DDV2N(41) * FD2  ! D/DY(2) ! D/DY(3)
              DDNORM (40) = DDNORM (30)
              DDNORM (41) = DDNORM (36)
              DDNORM (42) = DV2N( 6) * DFD2( 6) + DDV2N(42) * FD2  ! D/DY(3) ! D/DY(3)
              DDNORM (43) = DV2N( 6) * DFD2( 7) + DDV2N(43) * FD2  ! D/DZ(1) ! D/DY(3)
              DDNORM (44) = DV2N( 6) * DFD2( 8) + DDV2N(44) * FD2  ! D/DZ(2) ! D/DY(3)
              DDNORM (45) = DV2N( 6) * DFD2( 9) + DDV2N(45) * FD2  ! D/DZ(3) ! D/DY(3)
    DNORM(7)  =  DV2N (7) * FD2 +  DV1N (7) * FD1                            ! D/DZ(1)
              DDNORM (46) = DV2N( 7) * DFD2( 7) + DDV2N(46) * FD2 +          &
                            DV1N( 7) * DFD1( 7) + DDV1N(46) * FD1  ! D/DZ(1) ! D/DZ(1)
              DDNORM (47) = DV2N( 7) * DFD2( 8) + DDV2N(47) * FD2 +          &
                            DV1N( 7) * DFD1( 8) + DDV1N(47) * FD1  ! D/DZ(2) ! D/DZ(1)
              DDNORM (48) = DV2N( 7) * DFD2( 9) + DDV2N(48) * FD2 +          &
                            DV1N( 7) * DFD1( 9) + DDV1N(48) * FD1  ! D/DZ(3) ! D/DZ(1)
    DNORM(8)  =  DV2N (8) * FD2 +  DV1N (8) * FD1                            ! D/DZ(2)
            ! DDNORM (49) = DV2N( 8) * DFD2( 7) + DDV2N(49) * FD2 +          &
            !               DV1N( 8) * DFD1( 7) + DDV1N(49) * FD1  ! D/DZ(1) ! D/DZ(2)
              DDNORM (49) = DDNORM (47)
              DDNORM (50) = DV2N( 8) * DFD2( 8) + DDV2N(50) * FD2 +          &
                            DV1N( 8) * DFD1( 8) + DDV1N(50) * FD1  ! D/DZ(2) ! D/DZ(2)
              DDNORM (51) = DV2N( 8) * DFD2( 9) + DDV2N(51) * FD2 +          &
                            DV1N( 8) * DFD1( 9) + DDV1N(51) * FD1  ! D/DZ(3) ! D/DZ(2)
    DNORM(9)  =  DV2N (9) * FD2 +  DV1N (9) * FD1                            ! D/DZ(3)
              DDNORM (52) = DV2N( 9) * DFD2( 7) + DDV2N(52) * FD2 +          &
                            DV1N( 9) * DFD1( 7) + DDV1N(52) * FD1  ! D/DZ(1) ! D/DZ(3)
              DDNORM (53) = DV2N( 9) * DFD2( 8) + DDV2N(53) * FD2 +          &
                            DV1N( 9) * DFD1( 8) + DDV1N(53) * FD1  ! D/DZ(2) ! D/DZ(3)
              DDNORM (54) = DV2N( 9) * DFD2( 9) + DDV2N(54) * FD2 +          &
                            DV1N( 9) * DFD1( 9) + DDV1N(54) * FD1  ! D/DZ(3) ! D/DZ(3)

    !First derivatives of TNORM
               !Second derivatives of TNORM
    FDNORM( 1) = - DNORM(1)              ! D/DKx
                   !D/DX(1) 
               SDNORM( 1) =   DDNORM( 1)                                        ! D/DKx  D/DKx
               SDNORM( 2) =   DDNORM( 2)                                        ! D/DKy  D/DKx
               SDNORM( 3) =   DDNORM( 3)                                        ! D/DKz  D/DKx
               SDNORM(10) = - DDNORM( 1) + DDNORM( 7)                           ! D/DJx  D/DKx
               SDNORM(16) = - DDNORM( 2) + DDNORM( 8)                           ! D/DJy  D/DKx
               SDNORM(22) = - DDNORM( 3) + DDNORM( 9)                           ! D/DJz  D/DKx
               SDNORM(28) = - DDNORM( 7) + DDNORM( 4)                           ! D/DMx  D/DKx
               SDNORM(37) = - DDNORM( 8) + DDNORM( 5)                           ! D/DMy  D/DKx
               SDNORM(46) = - DDNORM( 9) + DDNORM( 6)                           ! D/DMz  D/DKx
               SDNORM(55) = - DDNORM( 4)                                        ! D/DLx  D/DKx
               SDNORM(67) = - DDNORM( 5)                                        ! D/DLy  D/DKx
               SDNORM(79) = - DDNORM( 6)                                        ! D/DLz  D/DKx
    FDNORM( 2) = - DNORM(2)              ! D/DKy                               
                   !D/DX(2)                                                    
             ! SDNORM( 4) =   DDNORM(10)                                        ! D/DKx  D/DKy
               SDNORM( 4) =   SDNORM( 2)                                       
               SDNORM( 5) =   DDNORM(11)                                        ! D/DKy  D/DKy
               SDNORM( 6) =   DDNORM(12)                                        ! D/DKz  D/DKy
               SDNORM(11) = - DDNORM(10) + DDNORM(16)                           ! D/DJx  D/DKy
               SDNORM(17) = - DDNORM(11) + DDNORM(17)                           ! D/DJy  D/DKy
               SDNORM(23) = - DDNORM(12) + DDNORM(18)                           ! D/DJz  D/DKy
               SDNORM(29) = - DDNORM(16) + DDNORM(13)                           ! D/DMx  D/DKy
               SDNORM(38) = - DDNORM(17) + DDNORM(14)                           ! D/DMy  D/DKy
               SDNORM(47) = - DDNORM(18) + DDNORM(15)                           ! D/DMz  D/DKy
               SDNORM(56) = - DDNORM(13)                                        ! D/DLx  D/DKy
               SDNORM(68) = - DDNORM(14)                                        ! D/DLy  D/DKy
               SDNORM(80) = - DDNORM(15)                                        ! D/DLz  D/DKy
    FDNORM( 3) = - DNORM(3)              ! D/DKz                               
                   !D/DX(3)                                                    
             ! SDNORM( 7) =   DDNORM(19)                                        ! D/DKx  D/DKz
             ! SDNORM( 8) =   DDNORM(20)                                        ! D/DKy  D/DKz
               SDNORM( 7) =   SDNORM( 3)                                       
               SDNORM( 8) =   SDNORM( 6)                                       
               SDNORM( 9) =   DDNORM(21)                                        ! D/DKz  D/DKz
               SDNORM(12) = - DDNORM(19) + DDNORM(25)                           ! D/DJx  D/DKz
               SDNORM(18) = - DDNORM(20) + DDNORM(26)                           ! D/DJy  D/DKz
               SDNORM(24) = - DDNORM(21) + DDNORM(27)                           ! D/DJz  D/DKz
               SDNORM(30) = - DDNORM(25) + DDNORM(22)                           ! D/DMx  D/DKz
               SDNORM(39) = - DDNORM(26) + DDNORM(23)                           ! D/DMy  D/DKz
               SDNORM(48) = - DDNORM(27) + DDNORM(24)                           ! D/DMz  D/DKz
               SDNORM(57) = - DDNORM(22)                                        ! D/DLx  D/DKz
               SDNORM(69) = - DDNORM(23)                                        ! D/DLy  D/DKz
               SDNORM(81) = - DDNORM(24)                                        ! D/DLz  D/DKz 
    FDNORM( 4) =   DNORM(1) - DNORM(7)   ! D/DJx 
                   !D/DX(1)   !D/DZ(1)
               SDNORM(13) =   DDNORM( 1) - DDNORM( 7) - DDNORM( 7) + DDNORM(46) ! D/DJx  D/DJx
               SDNORM(14) =   DDNORM( 2) - DDNORM(16) - DDNORM( 8) + DDNORM(47) ! D/DJy  D/DJx
               SDNORM(15) =   DDNORM( 3) - DDNORM(25) - DDNORM( 9) + DDNORM(48) ! D/DJz  D/DJx
               SDNORM(31) = - DDNORM( 4) + DDNORM(31) + DDNORM( 7) - DDNORM(46) ! D/DMx  D/DJx
               SDNORM(40) = - DDNORM( 5) + DDNORM(37) + DDNORM( 8) - DDNORM(47) ! D/DMy  D/DJx
               SDNORM(49) = - DDNORM( 6) + DDNORM(43) + DDNORM( 9) - DDNORM(48) ! D/DMz  D/DJx
               SDNORM(58) =   DDNORM( 4) - DDNORM(31)                           ! D/DLx  D/DJx
               SDNORM(70) =   DDNORM( 5) - DDNORM(37)                           ! D/DLy  D/DJx
               SDNORM(82) =   DDNORM( 6) - DDNORM(43)                           ! D/DLz  D/DJx 
    FDNORM( 5) =   DNORM(2) - DNORM(8)   ! D/DJy
                   !D/DX(2)   !D/DZ(2)  
             ! SDNORM(19) =   DDNORM(10) - DDNORM( 8) - DDNORM(16) + DDNORM(49) ! D/DJx  D/DJy
               SDNORM(19) =   SDNORM(14)
               SDNORM(20) =   DDNORM(11) - DDNORM(17) - DDNORM(17) + DDNORM(50) ! D/DJy  D/DJy
               SDNORM(21) =   DDNORM(12) - DDNORM(26) - DDNORM(18) + DDNORM(51) ! D/DJz  D/DJy
               SDNORM(32) = - DDNORM(13) + DDNORM(32) + DDNORM(16) - DDNORM(49) ! D/DMx  D/DJy
               SDNORM(41) = - DDNORM(14) + DDNORM(38) + DDNORM(17) - DDNORM(50) ! D/DMy  D/DJy
               SDNORM(50) = - DDNORM(15) + DDNORM(44) + DDNORM(18) - DDNORM(51) ! D/DMz  D/DJy
               SDNORM(59) =   DDNORM(13) - DDNORM(32)                           ! D/DLx  D/DJy
               SDNORM(71) =   DDNORM(14) - DDNORM(38)                           ! D/DLy  D/DJy
               SDNORM(83) =   DDNORM(15) - DDNORM(44)                           ! D/DLz  D/DJy 
    FDNORM( 6) =   DNORM(3) - DNORM(9)   ! D/DJz
                   !D/DX(3)   !D/DZ(3)
             ! SDNORM(25) =   DDNORM(19) - DDNORM( 9) - DDNORM(25) + DDNORM(52) ! D/DJx  D/DJz
             ! SDNORM(26) =   DDNORM(20) - DDNORM(18) - DDNORM(26) + DDNORM(53) ! D/DJy  D/DJz
               SDNORM(25) =   SDNORM(15)
               SDNORM(26) =   SDNORM(21)
               SDNORM(27) =   DDNORM(21) - DDNORM(27) - DDNORM(27) + DDNORM(54) ! D/DJz  D/DJz
               SDNORM(33) = - DDNORM(22) + DDNORM(33) + DDNORM(25) - DDNORM(52) ! D/DMx  D/DJz
               SDNORM(42) = - DDNORM(23) + DDNORM(39) + DDNORM(26) - DDNORM(53) ! D/DMy  D/DJz
               SDNORM(51) = - DDNORM(24) + DDNORM(45) + DDNORM(27) - DDNORM(54) ! D/DMz  D/DJz
               SDNORM(60) =   DDNORM(22) - DDNORM(33)                           ! D/DLx  D/DJz
               SDNORM(72) =   DDNORM(23) - DDNORM(39)                           ! D/DLy  D/DJz
               SDNORM(84) =   DDNORM(24) - DDNORM(45)                           ! D/DLz  D/DJz 
    FDNORM( 7) =   DNORM(7) - DNORM(4)   ! D/DMx
                   !D/DZ(1)   !D/DY(1)
               SDNORM(34) = - DDNORM(31) + DDNORM(28) + DDNORM(46) - DDNORM(31) ! D/DMx  D/DMx
               SDNORM(35) = - DDNORM(37) + DDNORM(29) + DDNORM(47) - DDNORM(32) ! D/DMy  D/DMx
               SDNORM(36) = - DDNORM(43) + DDNORM(30) + DDNORM(48) - DDNORM(33) ! D/DMz  D/DMx
               SDNORM(61) =   DDNORM(31) - DDNORM(28)                           ! D/DLx  D/DMx
               SDNORM(73) =   DDNORM(37) - DDNORM(29)                           ! D/DLy  D/DMx
               SDNORM(85) =   DDNORM(43) - DDNORM(30)                           ! D/DLz  D/DMx  
    FDNORM( 8) =   DNORM(8) - DNORM(5)   ! D/DMy
                   !D/DZ(2)   !D/DY(2)
             ! SDNORM(43) = - DDNORM(32) + DDNORM(34) + DDNORM(49) - DDNORM(37) ! D/DMx  D/DMy
               SDNORM(43) =   SDNORM(35)
               SDNORM(44) = - DDNORM(38) + DDNORM(35) + DDNORM(50) - DDNORM(38) ! D/DMy  D/DMy
               SDNORM(45) = - DDNORM(44) + DDNORM(36) + DDNORM(51) - DDNORM(39) ! D/DMz  D/DMy
               SDNORM(62) =   DDNORM(32) - DDNORM(34)                           ! D/DLx  D/DMy
               SDNORM(74) =   DDNORM(38) - DDNORM(35)                           ! D/DLy  D/DMy
               SDNORM(86) =   DDNORM(44) - DDNORM(36)                           ! D/DLz  D/DMy
    FDNORM( 9) =   DNORM(9) - DNORM(6)   ! D/DMz
                   !D/DZ(3)   !D/DY(3)
             ! SDNORM(52) = - DDNORM(33) + DDNORM(40) + DDNORM(52) - DDNORM(43) ! D/DMx  D/DMz
             ! SDNORM(53) = - DDNORM(39) + DDNORM(41) + DDNORM(53) - DDNORM(44) ! D/DMy  D/DMz
               SDNORM(52) =   SDNORM(36)
               SDNORM(53) =   SDNORM(45)
               SDNORM(54) = - DDNORM(45) + DDNORM(42) + DDNORM(54) - DDNORM(45) ! D/DMz  D/DMz
               SDNORM(63) =   DDNORM(33) - DDNORM(40)                           ! D/DLx  D/DMz
               SDNORM(75) =   DDNORM(39) - DDNORM(41)                           ! D/DLy  D/DMz
               SDNORM(87) =   DDNORM(45) - DDNORM(42)                           ! D/DLz  D/DMz 
    FDNORM(10) =   DNORM(4)              ! D/DLx                                
                   !D/DY(1)                                                     
               SDNORM(64) =   DDNORM(28)                                        ! D/DLx  D/DLx
               SDNORM(65) =   DDNORM(29)                                        ! D/DLy  D/DLx
               SDNORM(66) =   DDNORM(30)                                        ! D/DLz  D/DLx  
    FDNORM(11) =   DNORM(5)              ! D/DLy                                
                   !D/DY(2)                                                     
             ! SDNORM(76) =   DDNORM(34)                                        ! D/DLx  D/DLy
               SDNORM(76) =   SDNORM(65)                                        
               SDNORM(77) =   DDNORM(35)                                        ! D/DLy  D/DLy
               SDNORM(78) =   DDNORM(36)                                        ! D/DLz  D/DLy
    FDNORM(12) =   DNORM(6)              ! D/DLz                                
                   !D/DY(3)                                                     
             ! SDNORM(88) =   DDNORM(40)                                        ! D/DLx  D/DLz
             ! SDNORM(89) =   DDNORM(41)                                        ! D/DLy  D/DLz
               SDNORM(88) =   SDNORM(66)                                        
               SDNORM(89) =   SDNORM(78)                                        
               SDNORM(90) =   DDNORM(42)                                        ! D/DLz  D/DLz

! FDCOS : First  Derivatives of Cosine function
! SDCOS : Second Derivatives of Cosine function
!
! COSPHI = PRDT / TNORM   
!
    IFDCOS = 0
    ISDCOS = 0
    DO KX1 = 1, 3          ! Kx, Ky, Kz  
       IFDCOS = IFDCOS + 1
       FDCOS(IFDCOS) = ( FDDOT(KX1) * TNORM - PRDT * FDNORM(KX1) ) * TNORMINV2
       DO KX2 = 1, 3       ! Kx, Ky, Kz  
          ISDCOS = ISDCOS + 1
          SDCOS(ISDCOS) = LOCFUN(ISDCOS,KX1,KX2)
       END DO
    END DO

    DO KX1 = 4, 6          ! Jx, Jy, Jz 
       IFDCOS = IFDCOS + 1
       FDCOS(IFDCOS) = ( FDDOT(KX1) * TNORM - PRDT * FDNORM(KX1) ) * TNORMINV2       
       DO KX2 = 1, 6       ! Kx, Ky, Kz, Jx, Jy, Jz  
          ISDCOS = ISDCOS + 1
          SDCOS(ISDCOS) = LOCFUN(ISDCOS,KX1,KX2)
       END DO
    END DO

    DO KX1 = 7, 9          ! Mx, My, Mz 
       IFDCOS = IFDCOS + 1
       FDCOS(IFDCOS) = ( FDDOT(KX1) * TNORM - PRDT * FDNORM(KX1) ) * TNORMINV2       
       DO KX2 = 1, 9       ! Kx, Ky, Kz, Jx, Jy, Jz, Mx, My, Mz
          ISDCOS = ISDCOS + 1
          SDCOS(ISDCOS) = LOCFUN(ISDCOS,KX1,KX2)
       END DO
    END DO

    DO KX1 = 10, 12        ! Lx, Ly, Lz 
       IFDCOS = IFDCOS + 1
       FDCOS(IFDCOS) = ( FDDOT(KX1) * TNORM - PRDT * FDNORM(KX1) ) * TNORMINV2       
       DO KX2 = 1, 12      ! Kx, Ky, Kz, Jx, Jy, Jz, Mx, My, Mz, Lx, Ly, Lz 
          ISDCOS = ISDCOS + 1
          SDCOS(ISDCOS) = LOCFUN(ISDCOS,KX1,KX2)
       END DO
    END DO


!   CHECK FOR NUMERICAL INCONSISTENCIES
    IF ((ABS(TOR).LT.0.00001D0).OR.(ABS(TOR-PI).LT.0.00001D0).OR. &
        (ABS(TOR-2.D0*PI).LT.0.00001D0)) THEN
       FMOD1      = 0.D0
       FMOD2      = 0.D0
    ELSE
       FMOD1      = ( PRDT * TNORMINV ) / ( SIN(TOR) * SIN(TOR) ) 
       FMOD2      = - (1.D0 / ABS(SIN(TOR))) * SIGN(1.D0, PRDT_SIGN)
    ENDIF


!   STORE THE HESSIAN IN F
    ISDCOS     = 0
    INDF       = 0
    DO KX1 = 1, 12
       DO KX2 = 1, KX1
          INDF     = INDF   + 1
          ISDCOS   = ISDCOS + 1
          TERM1    = FMOD1  * FDCOS(KX1) * FDCOS(KX2)  
          TERM2    = SDCOS (ISDCOS)
          F(INDF)  = FMOD2 * ( TERM1 + TERM2 )
       END DO
       MODR     = MOD( KX1 , 3 )
       IF ( MODR .EQ. 0 ) MODR = 3
       ISDCOS = ISDCOS + ( 3 - MODR ) 
    END DO

! DEBUG SECTION
#ifdef DEBUG_DERI2
    CALL CHECK_FIRST_DERIVATIVES(LABEL='TORSION',I1=K,I2=J,I3=M,I4=L,COORD=C,DER1=FDCOS,FAC=FMOD2*FACT)
    CALL CHECK_SECOND_DERIVATIVES(LABEL='TORSION',I1=K,I2=J,I3=M,I4=L,COORD=C,DER2=F,FAC=1.d0)
#endif
    RETURN

  CONTAINS 
      ! Local Function
    DOUBLE PRECISION FUNCTION LOCFUN ( I3, I1, I2 )  
      IMPLICIT NONE
      INTEGER :: I1, I2, I3
       LOCFUN = ( ( SDDOT(I3) * TNORM               +       & 
                    FDDOT(I1) * FDNORM(I2)          -       &
                    FDDOT(I2) * FDNORM(I1)          -       &
                    PRDT      * SDNORM(I3)                  &         
                  ) * TNORM - 2.D0 * FDNORM(I2)     *       &
                  ( FDDOT(I1) * TNORM               -       &
                    PRDT      * FDNORM(I1)                  &
                  )                                         &
                ) * TNORMINV3
      RETURN
    END FUNCTION LOCFUN

  END SUBROUTINE DERI_2_TORSION

#ifdef DEBUG_DERI2
#include "debug_deri_2.DEBUG"
#endif
END MODULE DERIVATE_2_GEO_OBJECTS
