    SUBROUTINE BOCOV
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE BOCOV
!>   
!> SUBPROGRAM: BOCOV - UPDATE WIND POINTS ON BOUNDARY
!> PROGRAMMER: JANJIC
!> ORG: W/NP22
!> DATE: 94-03-08
!>     
!> ABSTRACT:
!> U AND V COMPONENTS OF THE WIND ARE UPDATED ON THE DOMAIN BOUNDARY BY APPLYING THE PRE-COMPUTED 
!> TENDENCIES AT EACH TIME STEP. 
!> AN EXTRAPOLATION FROM INSIDE THE DOMAIN IS USED FOR THE COMPONENT TANGENTIAL TO THE BOUNDARY IF 
!> THE NORMAL COMPONENT IS OUTWARD.
!>     
!> PROGRAM HISTORY LOG:
!> 87-??-??  MESINGER   - ORIGINATOR
!> 95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!> 98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!> 18-01-15  LUCCI      - MODERNIZATION OF THE CODE, INCLUDING:
!>                        * F77 TO F90/F95
!>                        * INDENTATION & UNIFORMIZATION CODE
!>                        * REPLACEMENT OF COMMONS BLOCK FOR MODULES
!>                        * DOCUMENTATION WITH DOXYGEN
!>                        * OPENMP FUNCTIONALITY 
!>
!> INPUT  ARGUMENT LIST:
!> NONE
!>
!> OUTPUT ARGUMENT LIST:
!> NONE
!>
!> INPUT/OUTPUT ARGUMENT LIST:
!> NONE
!> 
!> USE MODULES: BOCO
!>              CTLBLK
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
!> DRIVER     : DIGFLT
!>              EBU
!>              NEWFLT
!>
!> CALLS      : ------
!>--------------------------------------------------------------------------------------------------
    USE BOCO
    USE CTLBLK
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
    REAL   (KIND=R4KIND), PARAMETER :: D06666 = .06666666
!
    INTEGER(KIND=I4KIND), PARAMETER :: IMJM   = IM * JM - JM / 2
!------------------------
! IMPLICIT NONE VARIABLES
!------------------------
    INTEGER(KIND=I4KIND)                                                                        ::&
    & IIM     , JJM     , JL      , N       , I       , II      , J       , JJ
!
#include "sp.h"
!----------------------------------------------------
! TIME INTERPOLATION OF U AND V AT THE OUTER BOUNDARY
!----------------------------------------------------
    IIM = IM - MY_IS_GLB + 1
    JJM = JM - MY_JS_GLB + 1
!
    DO 115 JL=1,LM
!
        N = 1
        DO 111 I=1,IM-1
            UB(N,JL,1) = UB(N,JL,1) + UB(N,JL,2) * DT
            VB(N,JL,1) = VB(N,JL,1) + VB(N,JL,2) * DT
            IF (MY_JS_GLB == 1 .AND. I >= MY_IS_GLB-ILPAD1 .AND. I <= MY_IE_GLB+IRPAD1) THEN
                II = I - MY_IS_GLB + 1
                U(II,1,JL) = UB(N,JL,1)
                V(II,1,JL) = VB(N,JL,1)
            END IF
            N = N + 1
            CONTINUE
    111 END DO
!
        DO 112 I=1,IM-1
            UB(N,JL,1)=UB(N,JL,1)+UB(N,JL,2)*DT
            VB(N,JL,1)=VB(N,JL,1)+VB(N,JL,2)*DT
            IF (MY_JE_GLB == JM .AND. I >= MY_IS_GLB-ILPAD1 .AND. I <= MY_IE_GLB+IRPAD1) THEN
                II = I - MY_IS_GLB + 1
                U(II,JJM,JL) = UB(N,JL,1)
                V(II,JJM,JL) = VB(N,JL,1)
            END IF
            N = N + 1
            CONTINUE
    112 END DO
!
        DO 113 J=2,JM-1,2
            UB(N,JL,1) = UB(N,JL,1) + UB(N,JL,2) * DT
            VB(N,JL,1) = VB(N,JL,1) + VB(N,JL,2) * DT
            IF (MY_IS_GLB == 1 .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                JJ = J - MY_JS_GLB + 1
                U(1,JJ,JL) = UB(N,JL,1)
                V(1,JJ,JL) = VB(N,JL,1)
            END IF
            N = N + 1
            CONTINUE
    113 END DO
!
        DO 114 J=2,JM-1,2
            UB(N,JL,1) = UB(N,JL,1) + UB(N,JL,2) * DT
            VB(N,JL,1) = VB(N,JL,1) + VB(N,JL,2) * DT
            IF (MY_IE_GLB == IM .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                JJ = J - MY_JS_GLB + 1
                U(IIM,JJ,JL) = UB(N,JL,1)
                V(IIM,JJ,JL) = VB(N,JL,1)
            END IF
            N = N + 1
            CONTINUE
    114 END DO
        CONTINUE
!
115 END DO
!-------------------------------------------------------
! EXTRAPOLATION OF TANGENTIAL VELOCITY AT OUTFLOW POINTS
!-------------------------------------------------------
    DO 125 JL=1,LM
!
        IF (IBROW == 1) THEN
            DO 121 I=MYIS1_P1,MYIE2_P1
                IF (V(I,1,JL) < 0.) U(I,1,JL) = (VTM(I,5,JL) + 1.) * U(I,3,JL)                    &
    &                                         -  VTM(I,5,JL)       * U(I,5,JL)
                CONTINUE
        121 END DO
        END IF
!
        IF (ITROW == 1) THEN
            DO 122 I=MYIS1_P1,MYIE2_P1
                IF (V(I,JJM,JL) > 0.) U(I,JJM,JL) = (VTM(I,JJM-4,JL) + 1.) * U(I,JJM-2,JL)        &
    &                                             -  VTM(I,JJM-4,JL)       * U(I,JJM-4,JL)                         
                CONTINUE
        122 END DO
        END IF
!
        DO 123 J=4,JM-3,2
            IF (ILCOL == 1) THEN
                IF (MY_IS_GLB == 1 .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                    JJ = J - MY_JS_GLB + 1
                    IF (U(1,JJ,JL) < 0.) V(1,JJ,JL) = (VTM(3,JJ,JL) + 1.) * V(2,JJ,JL)            & 
    &                                               -  VTM(3,JJ,JL)       * V(3,JJ,JL)
                END IF
            END IF
            CONTINUE
    123 END DO
!
        DO 124 J=4,JM-3,2
            IF (IRCOL == 1) THEN
                IF (MY_IE_GLB == IM .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                    JJ = J - MY_JS_GLB + 1
                    IF (U(IIM,JJ,JL) > 0.)  V(IIM,JJ,JL) = (VTM(IIM-2,JJ,JL)+1.) * V(IIM-1,JJ,JL) &
    &                                                    -  VTM(IIM-2,JJ,JL)     * V(IIM-2,JJ,JL)
                END IF
            END IF
            CONTINUE
    124 END DO
        CONTINUE
!
125 END DO
!-----------------------------------------------------
! SPACE INTERPOLATION OF U AND V AT THE INNER BOUNDARY
!-----------------------------------------------------
    DO 140 JL=1,LM
!
        IF (IBROW == 1 .AND. ILCOL == 1) THEN
            U(2,2,JL) = D06666 * (4. * (U(1,1,JL) + U(2,1,JL) + U(2,3,JL))                        &
    &                 +                 U(1,2,JL) + U(1,4,JL) + U(2,4,JL))
!
            V(2,2,JL) = D06666 * (4. * (V(1,1,JL) + V(2,1,JL) + V(2,3,JL))                        &
    &                 +                 V(1,2,JL) + V(1,4,JL) + V(2,4,JL))
        END IF
!
        IF (IBROW == 1 .AND. IRCOL == 1) THEN
            U(IIM-1,2,JL) = D06666 * (4. * (U(IIM-2,1,JL) + U(IIM-1,1,JL) + U(IIM-2,3,JL))        &
    &                     +                 U(IIM  ,2,JL) + U(IIM  ,4,JL) + U(IIM-1,4,JL))
!
            V(IIM-1,2,JL) = D06666 * (4. * (V(IIM-2,1,JL) + V(IIM-1,1,JL) + V(IIM-2,3,JL))        &
    &                     +                 V(IIM  ,2,JL) + V(IIM  ,4,JL) + V(IIM-1,4,JL))
        END IF
!
        IF (ITROW == 1 .AND. ILCOL == 1) THEN
            U(2,JJM-1,JL) = D06666 * (4. * (U(1,JJM  ,JL) + U(2,JJM  ,JL) + U(2,JJM-2,JL))        &
    &                     +                 U(1,JJM-1,JL) + U(1,JJM-3,JL) + U(2,JJM-3,JL))
!
            V(2,JJM-1,JL) = D06666 * (4. * (V(1,JJM  ,JL) + V(2,JJM  ,JL) + V(2,JJM-2,JL))        &
    &                     +                 V(1,JJM-1,JL) + V(1,JJM-3,JL) + V(2,JJM-3,JL))
        END IF
!
        IF (ITROW == 1 .AND. IRCOL == 1) THEN
            U(IIM-1,JJM-1,JL) = D06666 * (4. * (U(IIM-2,JJM  ,JL)  + U(IIM-1,JJM  ,JL)            &
    &                         +                 U(IIM-2,JJM-2,JL)) + U(IIM,JJM-1  ,JL)            &
    &                         +                 U(IIM,JJM-3  ,JL)  + U(IIM-1,JJM-3,JL))
!
            V(IIM-1,JJM-1,JL) = D06666 * (4. * (V(IIM-2,JJM  ,JL)  + V(IIM-1,JJM  ,JL)            &
    &                         +                 V(IIM-2,JJM-2,JL)) + V(IIM,JJM-1  ,JL)            &
    &                         +                 V(IIM  ,JJM-3,JL)  + V(IIM-1,JJM-3,JL))
        END IF
!-----------------------------------------------------
! SPACE INTERPOLATION OF U AND V AT THE INNER BOUNDARY
!-----------------------------------------------------
        IF (IBROW == 1) THEN
            DO 131 I=MYIS2,MYIE2
                U(I,2,JL) = (U(I-1,1,JL) + U(I,1,JL) + U(I-1,3,JL) + U(I,3,JL)) * 0.25
                V(I,2,JL) = (V(I-1,1,JL) + V(I,1,JL) + V(I-1,3,JL) + V(I,3,JL)) * 0.25
                CONTINUE
        131 END DO
        END IF
!
        IF (ITROW == 1) THEN
            DO 132 I=MYIS2,MYIE2
                U(I,JJM-1,JL) = (U(I-1,JJM-2,JL) + U(I,JJM-2,JL)                                  &
    &                         +  U(I-1,JJM  ,JL) + U(I,JJM  ,JL)) * 0.25
!
                V(I,JJM-1,JL) = (V(I-1,JJM-2,JL) + V(I,JJM-2,JL)                                  &
    &                         +  V(I-1,JJM  ,JL) + V(I,JJM  ,JL)) * 0.25
                CONTINUE
        132 END DO
        END IF
!
        DO 133 J=3,JM-2,2
            IF (ILCOL == 1) THEN
                IF (MY_IS_GLB == 1 .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                    JJ = J - MY_JS_GLB + 1
                    U(1,JJ,JL) = (U(1,JJ-1,JL) + U(2,JJ-1,JL) + U(1,JJ+1,JL) + U(2,JJ+1,JL)) * 0.25
                    V(1,JJ,JL) = (V(1,JJ-1,JL) + V(2,JJ-1,JL) + V(1,JJ+1,JL) + V(2,JJ+1,JL)) * 0.25
                END IF
            END IF
            CONTINUE
    133 END DO
!
        IF (IRCOL == 1) THEN
            DO 134 J=3,JM-2,2
                IF (MY_IE_GLB == IM .AND. J >= MY_JS_GLB-JBPAD1 .AND. J <= MY_JE_GLB+JTPAD1) THEN
                    JJ = J - MY_JS_GLB + 1
                    U(IIM-1,JJ,JL) = 0.25 * (U(IIM-1,JJ-1,JL) + U(IIM,JJ-1,JL)                    &
    &                              +         U(IIM-1,JJ+1,JL) + U(IIM,JJ+1,JL))
!
                    V(IIM-1,JJ,JL) = 0.25 * (V(IIM-1,JJ-1,JL) + V(IIM,JJ-1,JL)                    &
    &                              +         V(IIM-1,JJ+1,JL) + V(IIM,JJ+1,JL))
                END IF
                CONTINUE
        134 END DO
        END IF
!
        CONTINUE
140 END DO
!
    RETURN
!
    END SUBROUTINE BOCOV
