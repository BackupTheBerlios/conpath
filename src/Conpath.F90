!     ===================================================================
!     Authors: T. Laino    -   CSCS (Manno - November 2002)             !
!                          -   NEST - INFM                              !
!                          -   Scuola Normale Superiore di Pisa         !
!              D.Passerone -   CSCS (Manno)                             !
!                                                                       !
!     ===================================================================
!     Preliminary code to perform Conical Intersection search on excited!
!     surfaces.                                                         !
!                                                                       !
!     Copyright (C) 2002  Teodoro Laino                                 !
!                         Daniele Passerone                             !
!                                                                       !
!     This program is free software; you can redistribute it and/or     !
!     modify it under the terms of the GNU General Public License       !
!     as published by the Free Software Foundation; either version 2    !
!     of the License, or (at your option) any later version.            !
!     This program is distributed in the hope that it will be useful,   !
!     but WITHOUT ANY WARRANTY; without even the implied warranty of    !
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !
!     GNU General Public License for more details.                      !
!                                                                       !
!     For information please contact authors at: t.laino@sns.it         !
!                                                                       !
!     ===================================================================
!     Version alpha 0.00                                                !
                                                                        !
!     ===================================================================
      PROGRAM CONPATH                                                   
!     ===================================================================
      USE SYSTEM_UTIL, ONLY: PSYST,           &                         !
                             COMMAND,         &                         !
                             CPSTOP                                     !
      USE START_JOB,   ONLY: CONPATH_START,   &                         !
                             CONPATH_LOGO,    &                         !
                             INPUTFILE,       &                         !
                             NUMAT,           &                         !
                             GEOFILE,         &                         !
                             ABINIT_INP,      &                         !
                             TSTEP,           &                         !
                             ENECONV,         &                         !
                             METHOD,          &                         !
                             COORDINATES,     &                         !
                             DYNTYPE,         &                         !
                             AREWERESTARTING, &                         !
                             TYPEZMAT                                   !
      USE DYNPREPARE,  ONLY: DYN_ALLOCATE,    &                         !
                             DYN_RDGEO,       &                         !
                             C,               &                         !
                             F,               &                         !
                             M,               &                         !
                             V,               &                         !
                             A,               &                         !
                             LAB
      USE DYNACART,          ONLY: DYNAMICS_CART                        !
      USE DYNA_ZMAT,         ONLY: DYNAMICS_ZMAT                        !
      USE DYNACONSTRAINT,    ONLY: DYNAMICS_CONSTR                      !
      USE UNITS_FILE,        ONLY: FLUSH_FILES                          !
      USE PERIODICTABLE,     ONLY: SETPERIODICTABLE                     !
      USE RESTART,           ONLY: READ_RESTART                         !
      USE MARSAGLIAS,        ONLY: AMRSET                               !
      USE READ_ZMAT,         ONLY: READ_ZMAT_INPUT                      !
      USE READ_ZMAT,         ONLY: NATOMLOCZMAT => NUMAT                !
      USE OPTIMIZATION_PROCEDURES, ONLY: INTERNAL_OPTIMIZE_PROC         !
      IMPLICIT NONE                                                     !
!     ------------------------------------------------------------------!
!     Reading Input file and getting all information useful for dynamics!
!     ------------------------------------------------------------------!
      CALL CONPATH_LOGO                                                 !
      CALL AMRSET ( 1000 )                                              !
      CALL SETPERIODICTABLE                                             !      
      CALL CONPATH_START                                                !
      CALL DYN_ALLOCATE                                                 !
      IF (AREWERESTARTING) THEN                                         !
         CALL READ_RESTART                                              !
      ELSE                                                              !
         SELECT CASE (TYPEZMAT)                                         !
         CASE ('MOLPRO')                                                !
            CALL DYN_RDGEO                                              !
            NATOMLOCZMAT = NUMAT                                        ! 
         CASE DEFAULT                                                   !
            CALL READ_ZMAT_INPUT                                        !
         END SELECT                                                     !
      END IF                                                            !
      CALL FLUSH_FILES                                                  !
!     ------------------------------------------------------------------!
!     Here starts the real computational module                         !
!     ------------------------------------------------------------------!
      SELECT CASE (METHOD)                                              !
      CASE ('DYNAMIC   ')                                               !
                                                                        !
         SELECT CASE (COORDINATES)                                      !
         CASE ('CARTESIAN')                                             !
!     ---------------------------- CARTESIAN ---------------------------!
            SELECT CASE(DYNTYPE)                                        !
            CASE('CONICAL   ')                                          !
               CALL DYNAMICS_CART                                       !
            CASE('CONSTRAIN ')                                          !
               CALL DYNAMICS_CONSTR                                     !
            CASE DEFAULT                                                !
               WRITE(6,'(A)')'DYNTYPE: '//DYNTYPE//' UNKNOWN!'          !
               CALL CPSTOP('CONPATH')                                   !
            END SELECT                                                  !
                                                                        !
         CASE ('ZMATCOORD')                                             !
!     ---------------------------- INTERNAL ----------------------------!
            CALL DYNAMICS_ZMAT                                          !
!     ------------------------------------------------------------------!
         CASE DEFAULT                                                   !
            WRITE(6,'(2A)')'CONICALPATH DYNAMICS IS ONLY POSSIBLE IN ' &!
            //'CARTESIAN OR INTERNAL COORDINATES!'                      !
            CALL CPSTOP('CONPATH')                                      !
         END SELECT                                                     !
                                                                        !
      CASE ('OPTIMIZE  ')                                               !
         SELECT CASE (COORDINATES)                                      !
         CASE ('CARTESIAN')                                             !
            WRITE(6,'(A)')'NOT IMPLEMENTED YET!'                        !
            CALL CPSTOP('CONPATH')                                      !
         CASE ('ZMATCOORD')                                             !
            CALL INTERNAL_OPTIMIZE_PROC                                 !
         END SELECT                                                     !
      CASE ('ANALTRAJ  ')                                               !
         WRITE(6,'(A)')'NOT IMPLEMENTED YET!'                           !
         CALL CPSTOP('CONPATH')                                         !         
      CASE DEFAULT                                                      !
                                                                        !
         WRITE(6,'(A)')'METHOD: '//METHOD//' UNKNOWN!'                  !
         CALL CPSTOP('CONPATH')                                         !
                                                                        !
      END SELECT                                                        !
!     ===================================================================
      END PROGRAM CONPATH                                               !
!     ===================================================================
