C---------------------------------------------------------------------- 
C  ENTRY POINT BORT IS REQUIRED FOR NON-CRAY SYSTEMS                   
C---------------------------------------------------------------------- 
      SUBROUTINE BORT(STR)                                             
      CHARACTER*(*) STR                                                 
      PRINT*
      PRINT*,'************************ABORT**************************'
      PRINT*,STR                                                        
      PRINT*,'************************ABORT**************************'
      PRINT*
      CALLEXIT(49)                                                          
      END                                                               
