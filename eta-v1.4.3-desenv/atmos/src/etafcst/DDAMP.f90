    SUBROUTINE DDAMP
!>-------------------------------------------------------------------------------------------------- 
!> SUBROUTINE DDAMP
!> 
!> SUBPROGRAM: DDAMP - DIVERGENCE DAMPING
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 94-03-08       
!>     
!> ABSTRACT:
!> DDAMP MODIFIES THE WIND COMPONENTS SO AS TO REDUCE THE HORIZONTAL DIVERGENCE.
!> A SWITCH PROVIDES THE OPTION OF ALSO MODIFYING THE TEMPERATURE FROM AN ENERGY VIEWPOINT.
!>     
!> PROGRAM HISTORY LOG:
!> 87-08-??  JANJIC     - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 95-03-28  BLACK      - ADDED EXTERNAL EDGE
!> 98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY 
!>
!> INPUT ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>    
!> USE MODULES: CLDWTR
!>              CONTIN
!>              CTLBLK
!>              DYNAM
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              PARMETA
!>              PVRBLS
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!> 
!> DRIVER     : EBU
!>              NEWFLT
!>
!> CALLS      : ZERO2
!>--------------------------------------------------------------------------------------------------                 
    USE CLDWTR
    USE CONTIN
    USE CTLBLK
    USE DYNAM
    USE GLB_TABLE
    USE INDX
    USE LOOPS    , ONLY : JAM
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE PVRBLS
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
#include "sp.h"
!
    LOGICAL(KIND=L4KIND)                                                                        ::&
    & HEAT
! 
    REAL   (KIND=R4KIND), PARAMETER :: RFCP = .25 / 1004.6
!
    REAL   (KIND=R4KIND), PARAMETER :: IMJM = IM * JM - JM / 2
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & RDPDX   , RDPDY   ,                                                                         &
    & UT      , VT      ,                                                                         &
    & CKE     , DPDE
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K
!
    HEAT = .FALSE. 
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
    DO 100 J=MYJS1_P1,MYJE1_P1
        DO 100 I=MYIS_P1,MYIE_P1
            CKE(I,J) = 0.
            CONTINUE
100 END DO
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (CKE    , DPDE    , RDPDX   , RDPDY   , UT      , VT)
!
    DO 150 K=1,LM
!
        CALL ZERO2(CKE)
        CALL ZERO2(DPDE)
!
        DO 110 J=MYJS_P2,MYJE_P2
            DO 110 I=MYIS_P1,MYIE_P1
                DPDE(I,J)   = DETA(K)     * PDSL(I,J)
                 DIV(I,J,K) =  DIV(I,J,K) * HBM2(I,J)
                CONTINUE
    110 END DO
!
        DO 120 J=MYJS2,MYJE2
            DO 120 I=MYIS_P1,MYIE_P1
                RDPDX(I,J) = VTM(I,J,K) / (DPDE(I+IVW(J),J) + DPDE(I+IVE(J),J))
                RDPDY(I,J) = VTM(I,J,K) / (DPDE(I,J-1)      + DPDE(I,J+1))
                CONTINUE
    120 END DO
!
        DO 130 J=MYJS2,MYJE2
            DO 130 I=MYIS1_P1,MYIE1_P1
                 UT(I,J)   = U(I,J,K)
                 VT(I,J)   = V(I,J,K)
                  U(I,J,K) = U(I,J,K) + (DIV(I+IVE(J),J,K) - DIV(I+IVW(J),J,K))                   &
    &                      * RDPDX(I,J) * DDMPU(I,J)
!
                  V(I,J,K) = V(I,J,K) + (DIV(I,J+1,K)      - DIV(I,J-1,K))                        &
    &                      * RDPDY(I,J) * DDMPV(I,J)
!
                CKE(I,J)   = 0.5 * (U(I,J,K) *  U(I,J,K)                                          &
    &                      -       UT(I,J)   * UT(I,J)                                            &
    &                      +        V(I,J,K) *  V(I,J,K)                                          &
    &                      -       VT(I,J)   * VT(I,J))
                CONTINUE
    130 END DO
!
        IF (HEAT) THEN
            DO 140 J=MYJS2,MYJE2
                DO 140 I=MYIS_P1,MYIE_P1
                    T(I,J,K) = T(I,J,K) - RFCP * (CKE(I+IHE(J),J) + CKE(I,J+1)                    &
    &                        +                    CKE(I+IHW(J),J) + CKE(I,J-1))                   &
    &                        * HBM2(I,J)
                    CONTINUE
        140 END DO
        END IF
!
        CONTINUE
!
150 END DO
!
    RETURN
!
    END SUBROUTINE DDAMP
