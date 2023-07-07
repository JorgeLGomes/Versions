    SUBROUTINE CLTEND(ICLTEND)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE CLTEND
!>
!> SUBPROGRAM: CLTEND - TEMPERATURE CHANGE BY CLOUD PROCESSES
!> PROGRAMMER: FERRIER 
!> ORG: W/NP22
!> DATE: 01-09-26
!>     
!> ABSTRACT: 
!> CLTEND GRADUALLY UPDATES TEMPERATURE TENDENCIES FROM CONVECTION GRID-SCALE MICROPHYSICS AND 
!> PRECIPITATION ASSIMILATION.    
!>
!> PROGRAM HISTORY LOG:
!> 01-09-26  ?????    - ORIGINATOR            
!> 18-01-15  LUCCI    - MODERNIZATION OF THE CODE, INCLUDING:
!>                      * F77 TO F90/F95
!>                      * INDENTATION & UNIFORMIZATION CODE
!>                      * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                      * DOCUMENTATION WITH DOXYGEN
!>                      * OPENMP FUNCTIONALITY
!>
!> INPUT  ARGUMENT LIST:
!> ICLTEND - FLAG SET TO -1 PRIOR TO PHYSICS CALLS, 0 AFTER PHYSICS CALLS, AND 1 FOR UPDATING 
!>           TEMPERATURES EVERY TIME STEP
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!>     
!> USE MODULES: CTLBLK
!>              C_TADJ
!>              F77KINDS
!>              GLB_TABLE
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : EBU
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------                
    USE CTLBLK
    USE C_TADJ
    USE F77KINDS
    USE GLB_TABLE
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE PARMETA
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
!    REAL   (KIND=R4KIND), SAVE, DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2, LM)                         ::&
!    & T_OLD
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                  , INTENT(IN)          ::&
    & ICLTEND
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & DELTPH
!
    IF (ICLTEND < 0) THEN
!--------------------------------------------------------------
! SAVE OLD TEMPERATURE ARRAY BEFORE CUCNVC, GSMDRIVE AND ADJPPT
!--------------------------------------------------------------
        DO K=1,LM
            DO J=MYJS,MYJE
                DO I=MYIS,MYIE
                    T_OLD(I,J,K) = T(I,J,K)
                END DO
            END DO
        END DO
    ELSE IF (ICLTEND == 0) THEN
!------------------------------------------------------------------
! CALCULATE TEMPERATURE TENDENCIES FROM CUCNVC, GSMDRIVE AND ADJPPT 
!------------------------------------------------------------------
        DELTPH = 1. / FLOAT(NPHS)
        DO K=1,LM
            DO J=MYJS,MYJE
                DO I=MYIS,MYIE
                    T_ADJ(I,J,K) =   HTM(I,J,K) * HBM2(I,J) * DELTPH * (T(I,J,K) - T_OLD(I,J,K))
                        T(I,J,K) = T_OLD(I,J,K)
                END DO
            END DO
        END DO
    ELSE
!--------------------------------------------------------------------------------------------------
! GRADUALLY UPDATE TEMPERATURE FROM CUCNVC, GSMDRIVE AND ADJPPT IN SMALL INCREMENTS EVERY DYNAMICS 
! TIME STEP
!--------------------------------------------------------------------------------------------------
        DO K=1,LM
            DO J=MYJS2,MYJE2
                DO I=MYIS,MYIE
                    T(I,J,K) = T(I,J,K) + T_ADJ(I,J,K)
                END DO
            END DO
        END DO
    END IF 
!
    RETURN
!
    END SUBROUTINE CLTEND
