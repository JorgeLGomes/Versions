!>-------------------------------------------------------------------------------------------------- 
!> INTERFACE EXCHM.h
!>
!> USE MODULES: -----
!>  
!> DRIVER     : DIGFLT
!>              EBU
!>              EPS
!>              FILT25
!>              GOSSIP
!>              INITS
!>              NEWFLT
!>              PDTEDT
!>              TURBL
!>
!> CALLS      : -----
!>--------------------------------------------------------------------------------------------------
    INTERFACE EXCH
!
    MODULE PROCEDURE EXCH0   , EXCH00  , EXCH01  , EXCH011 , EXCH0001111  ,                       &
    &                EXCH1   , EXCH11  , EXCH111 , EXCH1111, EXCH11111    , EXCH111111   ,        &
    &                IEXCH
!
    END INTERFACE

