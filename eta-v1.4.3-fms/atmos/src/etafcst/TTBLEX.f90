    SUBROUTINE TTBLEX(TREF, TTBL, ITB , JTB , KNUM , IARR , JARR, PDSL , AETAL, HTML, PT, PL, QQ, &
    &                   PP, RDP , THE0, STHE, RDTHE, THESP, IPTB, ITHTB)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE TTBLEX
!>
!> SUBPROGRAM: TTBLEX - DESIGNED TO DUPLICATE TIMEF
!> PROGRAMMER: ?????   
!> ORG: ?????
!> DATE: ??-??-??
!> 
!> ABSTRACT: 
!> EXTRACT TEMPERATURE OF THE MOIST ADIABAT FROM THE APPROPRIATE TTBL
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
!> TTBL  - 
!> ITB   - 
!> JTB   - 
!> KNUM  -
!> IARR  -
!> JARR  - 
!> PDSL  -
!> AETAL -
!> HTML  - 
!> PT    - 
!> PL    -
!> RDP   - 
!> THE0  -
!> STHE  - 
!> RDTHE -
!> THESP - 
!> IPTB  - 
!>
!> OUTPUT ARGUMENT LIST:
!> TREF  - 
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> QQ    -
!> PP    -
!> IPTB  - 
!> ITHTB - 

!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : CUCNVC
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
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM_LOC = IDIM2 * JDIM2
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM     = IM    *    JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: IMXJM    = IM    *    JM
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & TREF
!
    REAL   (KIND=R4KIND), DIMENSION(JTB,ITB)                              , INTENT(IN)          ::&
    & TTBL
!
    REAL   (KIND=R4KIND), DIMENSION(ITB)                                  , INTENT(IN)          ::&
    & THE0    , STHE
!
    INTEGER(KIND=I4KIND), DIMENSION(IMJM_LOC)                             , INTENT(IN)          ::&
    & IARR    , JARR
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(IN)          ::&
    & PDSL    , THESP   , HTML
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & PP      , QQ     
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & IPTB    , ITHTB    
!
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & PT      , PL      , RDP     , RDTHE         
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITB     , JTB     , KNUM    , I       , J       , IPTBK   , ITH     , IP      , KK     
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & AETAL   , PK      , TPK     , BTHE00K , STHE00K , BTHE10K , STHE10K , BTHK    , STHK    ,   &
    & TTHK    , T00K    , T10K    , T01K    , T11K 
!
    DO 500 KK=1,KNUM
!---------------------------------- 
! SCALING PRESSURE & TT TABLE INDEX
!---------------------------------- 
        I         = IARR(KK)
        J         = JARR(KK)
         PK       = PDSL(I,J) * AETAL + PT
        TPK       = (PK - PL) * RDP
          QQ(I,J) = TPK - AINT(TPK)
        IPTB(I,J) = INT(TPK) + 1
!--------------------------------- 
! KEEPING INDICES WITHIN THE TABLE 
!--------------------------------- 
        IF (IPTB(I,J) < 1) THEN
            IPTB(I,J) = 1
              QQ(I,J) = 0.
        END IF
!
        IF (IPTB(I,J) >= ITB) THEN
            IPTB(I,J) = ITB - 1
              QQ(I,J) = 0.
        END IF
!-------------------------------- 
! BASE AND SCALING FACTOR FOR THE 
!--------------------------------- 
        IPTBK = IPTB(I,J)
        BTHE00K = THE0(IPTBK)
        STHE00K = STHE(IPTBK)
        BTHE10K = THE0(IPTBK + 1)
        STHE10K = STHE(IPTBK + 1)
!------------------------------- 
! SCALING THE AND TT TABLE INDEX 
!------------------------------- 
        BTHK = (BTHE10K    - BTHE00K) * QQ(I,J) + BTHE00K
        STHK = (STHE10K    - STHE00K) * QQ(I,J) + STHE00K
        TTHK = (THESP(I,J) - BTHK)    / STHK    * RDTHE
!
           PP(I,J) = TTHK - AINT(TTHK)
        ITHTB(I,J) = INT(TTHK) + 1
!--------------------------------- 
! KEEPING INDICES WITHIN THE TABLE
!--------------------------------- 
        IF (ITHTB(I,J) < 1) THEN
            ITHTB(I,J) = 1
               PP(I,J) = 0.
        END IF
!
        IF (ITHTB(I,J) >= JTB) THEN
            ITHTB(I,J) = JTB - 1
               PP(I,J) = 0.
        END IF
!---------------------------------------------- 
! TEMPERATURE AT FOUR SURROUNDING TT TABLE PTS. 
!----------------------------------------------
        ITH  = ITHTB(I,J)
        IP   =  IPTB(I,J)
        T00K = TTBL(ITH  ,IP  )
        T10K = TTBL(ITH+1,IP  )
        T01K = TTBL(ITH  ,IP+1)
        T11K = TTBL(ITH+1,IP+1)
!------------------- 
! PARCEL TEMPERATURE
!-------------------
        TREF(I,J) = (T00K + (T10K - T00K)        * PP(I,J) + (T01K - T00K)    * QQ(I,J)                     &
    &             + (T00K - T10K  - T01K + T11K) * PP(I,J)                    * QQ(I,J)) * HTML(I,J)
!
500 END DO
!
    RETURN
!
    END SUBROUTINE TTBLEX
