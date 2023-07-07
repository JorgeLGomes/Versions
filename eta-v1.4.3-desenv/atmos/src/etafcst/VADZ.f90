    SUBROUTINE VADZ
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE VADZ
!>
!> SUBROUTINE: VADZ - VERTICAL ADVECTION OF HEIGHT
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 93-11-17
!>
!> ABSTRACT:
!> VADV CALCULATES THE CONTRIBUTION OF THE VERTICAL ADVECTION OF HEIGHT IN ORDER TO COMPUTE 
!> W = DZ / DT DIAGNOSTICALLY
!>
!> PROGRAM HISTORY LOG:
!> 93-11-17  JANJIC     - ORIGINATOR
!> 00-01-04  BLACK      - DISTRIBUTED MEMORY AND THREADS
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
!> OUTPUT FILES:
!> NONE
!>
!> USE MODULES: CLDWTR
!>              CONTIN
!>              CTLBLK
!>              DYNAM
!>              F77KINDS
!>              GLB_TABLE
!>              INDX
!>              LOOPS
!>              MAPPINGS
!>              MASKS
!>              MPPCOM
!>              NHYDRO
!>              PARMETA
!>              TEMPCOM
!>              TOPO
!>              VRBLS
!>
!> DRIVER     : EBU
!>
!> CALLS      : -----            
!>--------------------------------------------------------------------------------------------------
    USE CLDWTR
    USE CONTIN
    USE CTLBLK
    USE DYNAM
    USE F77KINDS
    USE GLB_TABLE
    USE INDX
    USE LOOPS
    USE MAPOT
    USE MAPPINGS
    USE MASKS
    USE MPPCOM
    USE NHYDRO
    USE PARMETA
    USE TEMPCOM
    USE TOPO
    USE VRBLS
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM  = IM * JM - JM / 2
    INTEGER(KIND=I4KIND), PARAMETER :: KSMUD = 0
!
    REAL   (KIND=R4KIND), DIMENSION(IDIM1:IDIM2, JDIM1:JDIM2)                                   ::&
    & PRET    , TTB     , ALP1    , FNE     , FSE     , ETADTL
!
    REAL   (KIND=R4KIND)                                                                        ::&
    & G       , RG      , RDT
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & I       , J       , K    

    REAL   (KIND=R4KIND)                                                                        ::&
    & ZETA    , DTLRG   , ALP1P   , DZ      , TTAL
!
    G   = 9.8
    RG  = 1. / 9.8
    RDT = 1. / DT
!
!$omp parallel do 
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            W(I,J,LM+1) = 0.
            Z(I,J,LM+1) = 0.
            IF (SIGMA) Z(I,J,LM+1) = FIS(I,J) * RG
            ALP1(I,J) = ALOG(PINT(I,J,LM+1))
        END DO
    END DO
!
    DO 50 K=LM,1,-1
        ZETA  =  DFL(K) * RG
        DTLRG = DETA(K) * RG
!
!$omp parallel do 
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                ALP1P = ALOG(PINT(I,J,K))
                DZ = (Q(I,J,K) * 0.608 - CWM(I,J,K) + 1.) * T(I,J,K) * (ALP1(I,J) - ALP1P)        &
    &              *  R / (DWDT(I,J,K) * G)
!
                    Z(I,J,K) =   (Z(I,J,K+1) + DZ-ZETA) * HTM(I,J,K) + ZETA
                PDWDT(I,J,K) = DWDT(I,J,K)
                 DWDT(I,J,K) =    W(I,J,K)
                    W(I,J,K) = (DZ - RTOP(I,J,K) * PDSLO(I,J) * DTLRG) * HTM(I,J,K)               &
    &                        * HBM2(I,J) + W(I,J,K+1)
!
                 ALP1(I,J)   = ALP1P
            END DO
        END DO
!
 50 END DO
!
    DO K=1,LM
!
!$omp parallel do 
!
        DO J=MYJS,MYJE
            DO I=MYIS,MYIE
                Z(I,J,K) = (Z(I,J,K) + Z(I,J,K+1)) * 0.5
                W(I,J,K) = (W(I,J,K) + W(I,J,K+1)) * HTM(I,J,K) * HBM2(I,J) * 0.5 * RDT
            END DO
        END DO
    END DO
!
!$omp parallel do 
!
    DO J=MYJS,MYJE
        DO I=MYIS,MYIE
            TTB(I,J) = 0.
        END DO
    END DO
!
    DO K=1,LM-1
!
!$omp parallel do 
!
        DO J=MYJS2,MYJE2
            DO I=MYIS1,MYIE1
                TTAL     = (Z(I,J,K+1) - Z(I,J,K)) * ETADT(I,J,K) * 0.5
                W(I,J,K) = (TTAL+TTB(I,J)) * RDETA(K) + W(I,J,K)
                TTB(I,J) = TTAL
            END DO
        END DO
    END DO
!
!$omp parallel do 
!
    DO J=MYJS2,MYJE2
        DO I=MYIS1,MYIE1
            W(I,J,LM) = TTB(I,J) * RDETA(LM) + W(I,J,LM)
        END DO
    END DO
!
    RETURN
!
    END SUBROUTINE VADZ
