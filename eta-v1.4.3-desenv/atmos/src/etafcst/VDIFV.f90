    SUBROUTINE VDIFV(LMVK, KTM, DTQ2, UZ0, VZ0, AKMS, U, V, AKM, Z)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE VDIFV
!>
!> SUBPROGRAM: VDIFV - VERTICAL DIFFUSION
!> PROGRAMMER: ?????   
!> ORG: ?????
!> DATE: ??-??-??
!> 
!> ABSTRACT: 
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ?????       - ORIGINATOR
!> 18-01-15  LUCCI       - MODERNIZATION OF THE CODE, INCLUDING:
!>                         * F77 TO F90/F95
!>                         * INDENTATION & UNIFORMIZATION CODE
!>                         * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                         * DOCUMENTATION WITH DOXYGEN
!>                         * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> LMVK  - MASS POINT MODEL SURFACE (LOWEST MODEL HEIGHT)
!> KTM   - 
!> DTQ2  - PHYSICS TINE STEP
!> UZ0   -
!> VZ0   -
!> AKMS  - SURFACE COEF. FOR U AND V DIVIDED BY DELTA Z
!> AKM   - SURFACE EXCHANGE COEF. FOR U AND V
!> Z     - INTERFACE HEIGHT
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> U     - U AT H POINT (TRANSPOSED)
!> V     - V AT H POINT (TRANSPOSED)
!>
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : TURBL
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
#include "sp.h"
! 
    INTEGER(KIND=I4KIND), PARAMETER :: LP1 = LM + 1
    INTEGER(KIND=I4KIND), PARAMETER :: LM1 = LM - 1
!
    REAL   (KIND=R4KIND), DIMENSION(LM)                                   , INTENT(INOUT)       ::&
    & U       ,  V
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                  , INTENT(IN)          ::&
    & AKM
!
    REAL   (KIND=R4KIND), DIMENSION(LM1)                                                        ::&
    & CM      , CR      , RSU     , RSV     , DTOZ
!
    REAL   (KIND=R4KIND), DIMENSION(LP1)                                  , INTENT(IN)          ::&
    & Z
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LMVK    , KTM     
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & LMVM    , LMVP    , KT      , L       , IVI
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & DTQ2    , UZ0     , VZ0     , AKMS
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DTDIF   , DTOZS   , AKMH    , RCMVB   , DTOZAK  , RCML    , DTOZL   , CF
!
    DTDIF = DTQ2 / FLOAT(KTM)
    LMVM = LMVK - 1
    LMVP = LMVK + 1
!
    DO 300 KT=1,KTM
!
        DO 100 L=1,LMVM
            DTOZ(L) = DTDIF / (Z(L) - Z(L+1))
              CR(L) = -DTOZ(L) * AKM(L)
    100 END DO
!    
         CM(1) = DTOZ(1) * AKM(1) + 1.
        RSU(1) =    U(1)
        RSV(1) =    V(1)
!
        DO 110 L=2,LMVM
            DTOZL  = DTOZ(L)
            CF     = -DTOZL * AKM(L-1) / CM(L-1)
             CM(L) = - CR(L-1) * CF + (AKM(L-1) + AKM(L)) * DTOZL + 1.
            RSU(L) = -RSU(L-1) * CF + U(L)
            RSV(L) = -RSV(L-1) * CF + V(L)
    110 END DO
!
        DTOZS = DTDIF / (Z(LMVK) - Z(LMVP))
        AKMH  = AKM(LMVM)
!    
        CF     = -DTOZS * AKMH / CM(LMVM)
        RCMVB  = 1. / ((AKMH + AKMS) * DTOZS - CR(LMVM) * CF + 1.)
        DTOZAK = DTOZS * AKMS
!
        U(LMVK) = (DTOZAK * UZ0 - RSU(LMVM) * CF + U(LMVK)) * RCMVB
        V(LMVK) = (DTOZAK * VZ0 - RSV(LMVM) * CF + V(LMVK)) * RCMVB
!
        DO 120 IVI=1,LMVM
            L    = LMVK - IVI
            RCML = 1. / CM(L)
            U(L) = (-CR(L) * U(L+1) + RSU(L)) * RCML
            V(L) = (-CR(L) * V(L+1) + RSV(L)) * RCML
    120 END DO
!
300 END DO
!
    RETURN
!
    END SUBROUTINE VDIFV


