    SUBROUTINE OZON2D(LK, POZN, XLAT, RSIN1, RCOS1, RCOS2, QO3)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE OZON2D
!>
!> SUBROUTINE: OZON2D - COMPUTE OZONE IN MODEL COLUMNS
!> PROGRAMMER: KATZ/CAMPANA
!> ORG: W/NP22
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> OZON2D COMPUTES THE OZONE MIXING RATIO IN EACH GRID BOX OF COLUMNS WITHIN THE MODEL DOMAIN
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  KATZ/KC    - ORIGINATOR
!> 96-07-26  BLACK      - MODIFIED FOR ETA MODEL
!> 98-10-28  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> LK    - NUMBER OF LAYERS IN COLUMNS
!> POZM  - ?????
!> XLAT  - GEODETIC LATITUDE OF GRID COLUMNS IN DEGREES
!> RSIN1 - INFORMATION RELATING TO POSITION OF EARTH IN ITS ORBIT
!> RCOS1 - INFORMATION RELATING TO POSITION OF EARTH IN ITS ORBIT
!> RCOS2 - INFORMATION RELATING TO POSITION OF EARTH IN ITS ORBIT
!>
!> OUTPUT ARGUMENT LIST:
!> QO3   - OZONE MIXING RATIO AT MIDLAYERS OF MODEL COLUMNS
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              SEASO3
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : RADTN
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE SEASO3
    USE TEMPCOM
    USE TOPO
! 
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND), PARAMETER :: NLP1  = NL + 1
!
    REAL   (KIND=R4KIND), PARAMETER :: RTD   = 57.2957795
!
#include "sp.h"
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LK)                      , INTENT(INOUT)       ::&
    & QO3     
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, LK)                      , INTENT(IN)          ::&
    & POZN                 
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                          , INTENT(IN)          ::&
    & XLAT   
! 
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & TTHAN
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, NL)                                            ::&
    & QO3O3 
!
    INTEGER(KIND=I4KIND), DIMENSION(IDIM1:IDIM2)                                                ::&
    & JJROW 
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    REAL   (KIND=R4KIND)                                                  , INTENT(IN)          ::&
    & RSIN1   , RCOS1   , RCOS2
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & TH2     , DO3V    , DO3VP   , APLO    , APHI
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & LK
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , K       , NUMITR  , ILOG    , NHALF   , IT  
!
    DO I=MYIS,MYIE
        TH2      = 0.2 * XLAT(I)
        JJROW(I) = 19.001 - TH2
        TTHAN(I) = (19 - JJROW(I)) - TH2
    END DO
!----------------------------------------------- 
! SEASONAL AND SPATIAL INTERPOLATION DONE BELOW.
!----------------------------------------------- 
    DO K=1,NL
        DO I=MYIS,MYIE
            DO3V  = XDUO3N(JJROW(I)  ,K) + RSIN1 * XDO3N2(JJROW(I)  ,K)                           &
    &             +                        RCOS1 * XDO3N3(JJROW(I)  ,K)                           &
    &             +                        RCOS2 * XDO3N4(JJROW(I)  ,K)
!
            DO3VP = XDUO3N(JJROW(I)+1,K) + RSIN1 * XDO3N2(JJROW(I)+1,K)                           &
    &             +                        RCOS1 * XDO3N3(JJROW(I)+1,K)                           &
    &             +                        RCOS2 * XDO3N4(JJROW(I)+1,K)
!--------------------------------------------------------------------------------------------        
! NOW LATITUDINAL INTERPOLATION AND CONVERT O3 INTO MASS MIXING RATIO (ORIG DATA MPY BY 1.E4)
!--------------------------------------------------------------------------------------------         
            QO3O3(I,K) = 1.E-4 * (DO3V + TTHAN(I) * (DO3VP - DO3V))
        END DO
    END DO
!----------------------------------------------------------- 
! VERTICAL INTERPOLATION FOR EACH GRIDPOINT (LINEAR IN LN P)
!----------------------------------------------------------- 
    NUMITR = 0
    ILOG   = NL
 20 CONTINUE
    ILOG = (ILOG + 1) / 2
!
    IF (ILOG == 1) GOTO 25
    NUMITR = NUMITR + 1
    GOTO 20
 25 CONTINUE
!
    DO 60 K=1,LK
!    
        NHALF = (NL + 1) / 2
!
        DO I=MYIS,MYIE
            JJROW(I) = NHALF
        END DO
!    
        DO 40 IT=1,NUMITR
            NHALF = (NHALF + 1) / 2
            DO I=MYIS,MYIE
                IF (POZN(I,K) < PRGFDL(JJROW(I)-1)) THEN
                    JJROW(I) = JJROW(I) - NHALF
                ELSE IF(POZN(I,K) >= PRGFDL(JJROW(I))) THEN
                    JJROW(I) = JJROW(I) + NHALF
                END IF
!
                JJROW(I) = MIN(JJROW(I), NL)
                JJROW(I) = MAX(JJROW(I), 2 )
            END DO
     40 END DO
!    
        DO 50 I=MYIS,MYIE
            IF (POZN(I,K) < PRGFDL(1)) THEN
                QO3(I,K) = QO3O3(I,1 )
            ELSE IF(POZN(I,K) > PRGFDL(NL)) THEN
                QO3(I,K) = QO3O3(I,NL)
            ELSE
                APLO = ALOG(PRGFDL(JJROW(I)-1))
                APHI = ALOG(PRGFDL(JJROW(I)  ))
!
                QO3(I,K) =  QO3O3(I,JJROW(I)       )                                              &
    &                    + (ALOG(POZN(I,K))  - APHI)                                              &
    &                    /       (APLO       - APHI)                                              &
    &                    * (QO3O3(I,JJROW(I) -   1 ) - QO3O3(I,JJROW(I)       ))
            END IF
     50 END DO
!    
 60 END DO
!
    RETURN
!
    END SUBROUTINE OZON2D
