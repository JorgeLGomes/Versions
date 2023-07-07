    SUBROUTINE SGETMO(A, NHRZ1, NHRZ2, LM1, ACOL_T, LM2)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE SGETMO
!>
!> SUBROUTINE: ????? - ?????
!> PROGRAMMER: ?????
!> ORG: ?????
!> DATE: ??-??-??
!>
!> ABSTRACT:
!> ?????
!>
!> PROGRAM HISTORY LOG:
!> ??-??-??  ????? - ORIGINATOR
!> 18-01-15  LUCCI - MODERNIZATION OF THE CODE, INCLUDING:
!>                   * F77 TO F90/F95
!>                   * INDENTATION & UNIFORMIZATION CODE
!>                   * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                   * DOCUMENTATION WITH DOXYGEN
!>                   * OPENMP FUNCTIONALITY
!>
!> INPUT ARGUMENT LIST:
!> A      - 
!> LM1    -
!> LM2    - 
!> NHRZ1  -
!> NHRZ2  -  
!>
!> OUTPUT ARGUMENT LIST:
!> ACOL_T - 
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE 
!>
!> USE MODULES: F77KINDS
!>
!> DRIVER     : GSCOND
!>              PRECPD
!>              TURBL
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    USE F77KINDS
!
    IMPLICIT NONE
!
    REAL   (KIND=R4KIND), DIMENSION(NHRZ1, LM1)                 , INTENT(IN)                    ::&
    & A
!
    REAL   (KIND=R4KIND), DIMENSION(LM1, NHRZ1)                 , INTENT(INOUT)                 ::&
    & ACOL_T
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                        , INTENT(IN)                    ::&
    & LM1     , LM2     , NHRZ1   , NHRZ2
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & L       , N
!
    DO N=1,NHRZ1
        DO L=1,LM1
            ACOL_T(L,N) = A(N,L)
        END DO
    END DO
!
    RETURN
!
    END SUBROUTINE SGETMO 
