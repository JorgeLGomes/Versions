    SUBROUTINE O3CLIM(MYPE)
!>-------------------------------------------------------------------------------------------------- 
!> SUBROUTINE O3CLIM
!>
!> SUBROUTINE: O3CLIM - GENERATE SEASONAL OZONE DISTRIBUTION
!> PROGRAMMER: GFDL/CAMPANA
!> ORG: W/NP22
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> O3CLIM COMPUTES THE SEASONAL CLIMATOLOGY OF OZONE USING 81-LAYER DATA FROM GFDL.
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  GFDL/KC    - ORIGINATOR
!> 96-07-26  BLACK      - MODIFIED FOR ETA MODEL
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY
!>
!> INPUT  ARGUMENT LIST:
!> MYPE - 
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>
!> USE MODULES: F77KINDS
!>              O3DATA
!>              SEASO3
!>
!> DRIVER     : INIT
!>              INITS       
!>
!> CALLS      : -----     
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
    USE O3DATA
    USE SEASO3
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND), PARAMETER :: NLP1 = NL + 1
    INTEGER(KIND=I4KIND), PARAMETER :: NKK  = 41
    INTEGER(KIND=I4KIND), PARAMETER :: NK   = 81
    INTEGER(KIND=I4KIND), PARAMETER :: NKP  = NK + 1
!
#include "sp.h"
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    REAL   (KIND=R4KIND), DIMENSION(NL)                                                         ::&
    & PSTD    , RBAR    , PHALF
!
    REAL   (KIND=R4KIND), DIMENSION(19)                                                         ::&
    & TEMPN
!
    REAL   (KIND=R4KIND), DIMENSION(37, NL)                                                     ::&
    & O35DEG
!
    REAL   (KIND=R4KIND), DIMENSION(19, NL)                                                     ::& 
    & DDUO3N
!
    REAL   (KIND=R4KIND), DIMENSION(19, 41)                                                     ::&
    & DUO3N 
!
    REAL   (KIND=R4KIND), DIMENSION(10, 41)                                                     ::&
    & RO3     , RO31    , RO32
!
    REAL   (KIND=R4KIND), DIMENSION(10, 40)                                                     ::&
    & RO3M
!
    REAL   (KIND=R4KIND), DIMENSION(81)                                                         ::&
    & RSTD    , RDATA 
!
    INTEGER(KIND=I4KIND)                                                                        ::& 
    & I       , J       , K       , N       , NCASE   , KI      , KK      , KQ      , NKM     ,   &
    & IPLACE  , NKMM    , JJ      , KEN 
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & MYPE
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & AVG     , A1      , B1      , B2 
!--------------------------------------------------------------------------------------------------
! COMPUTE DETAILED O3 PROFILE FROM THE ORIGINAL GFDL PRESSURES WHERE OUTPUT FROM O3INT (PSTD) IS 
! TOP DOWN IN MB*1.E3 AND PSFC=1013.25 MB    
! K.A.C. DEC94
!--------------------------------------------------------------------------------------------------
    DO K=1,NK
        PH(K) =PH(K) * 1013250.
         P(K) = P(K) * 1013250.
    END DO
!
    PH(NKP) = PH(NKP) * 1013250.
!
    DO K=1,NL
        PSTD(K) = P(K)
    END DO
!
    DO K=1,25
        DO N=1,10
            RO31(N,K) = O3HI(N,K)
            RO32(N,K) = O3HI(N,K)
        END DO
    END DO
!
    DO 100 NCASE=1,4
!--------------------------     
! NCASE=1: SPRING (IN N.H.)
! NCASE=2: FALL   (IN N.H.)
! NCASE=3: WINTER (IN N.H.)
! NCASE=4: SUMMER (IN N.H.)
!--------------------------     
        IPLACE = 2
        IF (NCASE == 2) IPLACE = 4
        IF (NCASE == 3) IPLACE = 1
        IF (NCASE == 4) IPLACE = 3
!    
        IF (NCASE == 1 .OR. NCASE == 2) THEN
            DO K=26,41
                DO N=1,10
                    RO31(N,K) = O3LO1(N,K-25)
                    RO32(N,K) = O3LO2(N,K-25)
                END DO
            END DO
        END IF
!    
        IF (NCASE == 3 .OR. NCASE == 4) THEN
            DO K=26,41
                DO N=1,10
                    RO31(N,K) = O3LO3(N,K-25)
                    RO32(N,K) = O3LO4(N,K-25)
                END DO
            END DO
        END IF
!    
        DO 25 KK=1,NKK
            DO N=1,10
                DUO3N(N  ,KK) = RO31(11-N,KK)
                DUO3N(N+9,KK) = RO32(   N,KK)
            END DO
            DUO3N(10,KK) = 0.5 * (RO31(1,KK) + RO32(1,KK))
        25 END DO
!--------------------------------------------------------------------      
! FOR NCASE=2 OR NCASE=4,REVERSE LATITUDE ARRANGEMENT OF CORR. SEASON
!--------------------------------------------------------------------      
        IF (NCASE == 2 .OR. NCASE == 4) THEN
            DO 50 KK=1,NKK
                DO N=1,19
                    TEMPN(N) = DUO3N(20-N,KK)
                END DO
                DO N=1,19
                    DUO3N(N,KK) = TEMPN(N)
                END DO
            50 END DO
        END IF
!----------------------------------------------------------------------      
! DUO3N NOW IS O3 PROFILE FOR APPROPRIATE SEASON AT STD PRESSURE LEVELS
!---------------------------------------------------------------------- 
!
!-----------------------------  
! BEGIN LATITUDE (10 DEG) LOOP
!-----------------------------     
        DO 75 N=1,19
!        
            DO KK=1,NKK
                RSTD(KK) = DUO3N(N,KK)
            END DO
!        
            NKM  = NK-1
            NKMM = NK-3
!-----------------------------------------  
! BESSELS HALF-POINT INTERPOLATION FORMULA
!-----------------------------------------  
            DO K=4,NKMM,2
                KI = K / 2
                RDATA(K) = 0.5 * (RSTD(KI  ) + RSTD(KI+1)) - (RSTD(KI+2)                            &
    &                    -        RSTD(KI+1) - RSTD(KI  )  +  RSTD(KI-1)) / 16.
            END DO
!        
            RDATA(2)   = 0.5 * (RSTD(2)   + RSTD(1))
            RDATA(NKM) = 0.5 * (RSTD(NKK) + RSTD(NKK-1))
!----------------------------------         
! PUT UNCHANGED DATA INTO NEW ARRAY
!----------------------------------         
            DO K=1,NK,2
                KQ = (K+1) / 2
                RDATA(K) = RSTD(KQ)
            END DO
!        
            DO KK=1,NL
                DDUO3N(N,KK) = RDATA(KK) * .01
            END DO
!        
        75 END DO
!---------------------      
! END OF LATITUDE LOOP
!---------------------    
!
!-----------------------------------------------------------------------
! CREATE 5 DEG OZONE QUANTITIES BY LINEAR INTERPOLATION OF 10 DEG VALUES
!-----------------------------------------------------------------------
        DO 90 KK=1,NL
!        
            DO N=1,19
                O35DEG(2*N-1,KK) = DDUO3N(N,KK)
            END DO
!        
            DO N=1,18
                O35DEG(2*N,KK) = 0.5 * (DDUO3N(N,KK) + DDUO3N(N+1,KK))
            END DO
!        
        90 END DO
!    
        DO JJ=1,37
            DO KEN=1,NL
                O3O3(JJ,KEN,IPLACE) = O35DEG(JJ,KEN)
            END DO
        END DO
!    
    100 END DO
!----------------------- 
! END OF LOOP OVER CASES
!----------------------- 
!
!--------------------------------------------------------------------------------------------------
! AVERAGE CLIMATOLOGICAL VALUS OF O3 FROM 5 DEG LAT MEANS, SO THAT TIME AND SPACE INTERPOLATION 
! WILL WORK (SEE SUBR OZON2D)
!--------------------------------------------------------------------------------------------------
    DO I=1,LNGTH
        AVG = 0.25 * ( XRAD1(I) + XRAD2(I)  +  XRAD3(I) + XRAD4(I))
        A1  = 0.5  * ( XRAD2(I) - XRAD4(I))
        B1  = 0.5  * ( XRAD1(I) - XRAD3(I))
        B2  = 0.25 * ((XRAD1(I) + XRAD3(I)) - (XRAD2(I) + XRAD4(I)))
!
        XRAD1(I) = AVG
        XRAD2(I) = A1
        XRAD3(I) = B1
        XRAD4(I) = B2
    END DO
!----------------------------------------
! CONVERT GFDL PRESSURE (MICROBARS) TO PA
!----------------------------------------
    DO N=1,NL
        PRGFDL(N) = PSTD(N) * 1.E-1
    END DO
!
    RETURN
    END SUBROUTINE O3CLIM 	
