    SUBROUTINE FILT25(ZI, TM, IPASS)
!>--------------------------------------------------------------------------------------------------  
!> SUBROUTINE FILT25
!> 
!> SUBPROGRAM: FILT25 - FILTERS THE ARRAY ZI
!> PROGRAMMER: DIMEGO
!> ORG: W/NP22
!> DATE: 86-07-18
!>
!> ABSTRACT:
!> FILTERS AN ARRAY USING A 25PT BLECK FILTER IN THE INTERIOR OF THE DOMAIN.
!>
!> PROGRAM HISTORY LOG:
!> 86-07-18  G DIMEGO  - ORIGINATOR
!> 88-09-23  B SCHMIDT - ADDED THE DOCBLOCK
!> 90-11-27  G DIMEGO  - LEFT Z AS INTERNAL WORK ARRAY ON CRAY
!> 93-06-21  R TREADON - STREAMLINED CODE
!> 95-06-21  T BLACK   - MODIFIED FOR THE E-GRID
!> 99-08-25  T BLACK   - MODIFIED FOR DISTRIBUTED MEMORY
!> 18-01-15  LUCCI     - MODERNIZATION OF THE CODE, INCLUDING:
!>                       * F77 TO F90/F95
!>                       * INDENTATION & UNIFORMIZATION CODE
!>                       * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                       * DOCUMENTATION WITH DOXYGEN
!>                       * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> TM    - ARRAY CONTAINING THE TOPOGRAPHY MASK
!> IPASS - NUMBER OF PASSES THROUGH THE FILTER
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> ZI    - ARRAY CONTAINING THE ARRAY TO BE FILTERED / FILTERED FIELD
!>
!> USE MODULES: F77KINDS
!>              GLB_TABLE
!>              INDX
!>              MAPPINGS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>
!> DRIVER     : NEWFLT
!>
!> CALLS      : EXCH
!>--------------------------------------------------------------------------------------------------
    USE EXCHM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE MAPPINGS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
!
    IMPLICIT NONE
!
    INCLUDE "EXCHM.h"
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(INOUT)       ::&
    & ZI
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & Z
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)             , INTENT(IN)          ::&
    & TM
!
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & IPASS
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    REAL   (KIND=R4KIND)                                                                        ::&
    & CF1     , CF2     , CF3     , CF4     , CF5     , CF6     , HTOT    , SUMH    , RSUMH
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , IP      , IPC     , IIHE    , IIHW   
!
    DATA CF1/0.279372/  , CF2/ 0.171943/    , CF3/-0.006918/    ,                                 &
    &    CF4/0.077458/  , CF5/-0.024693/    , CF6/-0.012940/
!
    IPC = IPASS
!
    DO J=JDIM1,JDIM2
        DO I=IDIM1,IDIM2
            Z(I,J) = 0.
        END DO
    END DO
!---------------------------------------------------
! FILTER THE INTERIOR POINTS WITH 25-PT BLECK FILTER
!---------------------------------------------------
    DO 105 IP=1,IPC
!------- 
! OPENMP
!------- 
!
!$omp parallel do private (HTOT   , I       , IIHE    , IIHW    , RSUMH   , SUMH)
!
        DO 102 J=MYJS4,MYJE4
            DO 101 I=MYIS2,MYIE2
!            
                IIHE = I + IHE(J)
                IIHW = I + IHW(J)
                12345 FORMAT(' MASK=',F2.0,' IHE=',I2,' IHW=',I2,' L=',I2)
                12346 FORMAT(' 2ND ROW=',4(1x,F2.0))
                12347 FORMAT(' 3RD ROW=',4(1x,F2.0))
                12348 FORMAT(' 4TH ROW=',4(1x,F2.0))
                12349 FORMAT(' 5TH ROW=',4(1x,F2.0))
                12350 FORMAT(' 6TH ROW=',4(1x,F2.0))
                12351 FORMAT(' 7TH ROW=',4(1x,F2.0))
!
                HTOT = TM(I     ,J  )                                                             &
    &                + TM(IIHW  ,J-1) + TM(IIHE  ,J-1) + TM(IIHE  ,J+1) + TM(IIHW  ,J+1)          &
    &                + TM(I-1   ,J-2) + TM(I+1   ,J-2) + TM(I+1   ,J+2) + TM(I-1   ,J+2)          &
    &                + TM(I     ,J-2) + TM(I     ,J+2) + TM(I+1   ,J  ) + TM(I-1   ,J  )          &
    &                + TM(IIHE  ,J-3) + TM(IIHE+1,J-1) + TM(IIHW  ,J-3) + TM(IIHE+1,J+1)          &
    &                + TM(IIHW-1,J-1) + TM(IIHE  ,J+3) + TM(IIHW-1,J+1) + TM(IIHW  ,J+3)          &
    &                + TM(I     ,J-4) + TM(I+2   ,J  ) + TM(I-2   ,J  ) + TM(I     ,J+4)
!            
                IF (HTOT > 0.) THEN
                SUMH = CF1 *  TM(I     ,J  )                                                      &
    &                + CF2 * (TM(IIHW  ,J-1) + TM(IIHE  ,J-1) + TM(IIHE  ,J+1) + TM(IIHW  ,J+1))  &
    &                + CF3 * (TM(I-1   ,J-2) + TM(I+1   ,J-2) + TM(I+1   ,J+2) + TM(I-1   ,J+2))  &
    &                + CF4 * (TM(I     ,J-2) + TM(I     ,J+2) + TM(I+1   ,J  ) + TM(I-1   ,J  ))  &
    &                + CF5 * (TM(IIHE  ,J-3) + TM(IIHE+1,J-1) + TM(IIHW  ,J-3) + TM(IIHE+1,J+1)   &
    &                +        TM(IIHW-1,J-1) + TM(IIHE  ,J+3) + TM(IIHW-1,J+1) + TM(IIHW  ,J+3))  &
    &                + CF6 * (TM(I     ,J-4) + TM(I+2   ,J  ) + TM(I-2   ,J  ) + TM(I     ,J+4))
!                
                RSUMH = 1. / SUMH
!                
                Z(I,J) = (CF1  *  ZI(I     ,J  ) * TM(I     ,J  )                                 &
    &                  +  CF2  * (ZI(IIHW  ,J-1) * TM(IIHW  ,J-1)                                 &
    &                  +          ZI(IIHE  ,J-1) * TM(IIHE  ,J-1)                                 &
    &                  +          ZI(IIHE  ,J+1) * TM(IIHE  ,J+1)                                 &
    &                  +          ZI(IIHW  ,J+1) * TM(IIHW  ,J+1))                                &
    &                  +  CF3  * (ZI(I-1   ,J-2) * TM(I-1   ,J-2)                                 &
    &                  +          ZI(I+1   ,J-2) * TM(I+1   ,J-2)                                 &
    &                  +          ZI(I+1   ,J+2) * TM(I+1   ,J+2)                                 &
    &                  +          ZI(I-1   ,J+2) * TM(I-1   ,J+2))                                &
    &                  +  CF4  * (ZI(I     ,J-2) * TM(I     ,J-2)                                 &
    &                  +          ZI(I     ,J+2) * TM(I     ,J+2)                                 &
    &                  +          ZI(I+1   ,J  ) * TM(I+1   ,J  )                                 &
    &                  +          ZI(I-1   ,J  ) * TM(I-1   ,J  ))                                &
    &                  +  CF5  * (ZI(IIHE  ,J-3) * TM(IIHE  ,J-3)                                 &
    &                  +          ZI(IIHE+1,J-1) * TM(IIHE+1,J-1)                                 &
    &                  +          ZI(IIHW  ,J-3) * TM(IIHW  ,J-3)                                 &
    &                  +          ZI(IIHE+1,J+1) * TM(IIHE+1,J+1)                                 &
    &                  +          ZI(IIHW-1,J-1) * TM(IIHW-1,J-1)                                 &
    &                  +          ZI(IIHE  ,J+3) * TM(IIHE  ,J+3)                                 &
    &                  +          ZI(IIHW-1,J+1) * TM(IIHW-1,J+1)                                 &
    &                  +          ZI(IIHW  ,J+3) * TM(IIHW  ,J+3))                                &
    &                  +  CF6  * (ZI(I     ,J-4) * TM(I     ,J-4)                                 &
    &                  +          ZI(I+2   ,J  ) * TM(I+2   ,J  )                                 &
    &                  +          ZI(I-2   ,J  ) * TM(I-2   ,J  )                                 &
    &                  +          ZI(I     ,J+4) * TM(I     ,J+4)))                               &
    &                  * RSUMH
                END IF
!            
        101 END DO
    102 END DO 
!------- 
! OPENMP
!------- 
!
!$omp parallel do
!
        DO J=MYJS4,MYJE4
            DO I=MYIS2,MYIE2
                IF (TM(I,J) > 0.5) THEN
                    ZI(I,J) = Z(I,J)
                END IF
            END DO
        END DO
!    
        CALL EXCH(ZI, 1, 4, 4)
!    
105 END DO
!
    RETURN
!
    END SUBROUTINE FILT25
