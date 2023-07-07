C---------------------------------------------------------------------- 
C  UPDATE THE MESSAGE BUFFER WITH NEW SUBSET                            
C---------------------------------------------------------------------- 
      SUBROUTINE MSGUPD(LUNIT,LUN)                                      
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /MSGPTR/ NBY0,NBY1,NBY2,NBY3,NBY4,NBY5                     
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(32),MBAY(5000,32)     
                                                                        
      CHARACTER*8 CBAY                                                  
      EQUIVALENCE (CBAY,JBAY)                                           
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  PAD THE SUBSET BUFFER                                                
C  ---------------------                                                
                                                                        
      CALL PAD(IBAY,IBIT,IBYT,8)                                        
      GOTO 1                                                            
                                                                        
C  SPECIAL ENTRY POINT FOR COPYSB                                       
C  ------------------------------                                       
                                                                        
      ENTRY SUBUPD(LUNIT,LUN,JBYT)                                      
      IBYT = JBYT                                                       
                                                                        
C  SEE IF THE NEW SUBSET FITS                                           
C  --------------------------                                           
                                                                        
1     IF(MBYT(LUN)+IBYT.GT.MAXBYT) THEN                                 
         CALL MSGWRT(LUNIT,MBAY(1,LUN),MBYT(LUN))                       
         CALL MSGINI(LUN)                                               
      ENDIF                                                             
                                                                        
      IF(IBYT.GT.MAXBYT-MBYT(LUN)) GOTO 900                             
                                                                        
C  SET A BYTE COUNT AND TRANSFER THE SUBSET BUFFER INTO THE MESSAGE     
C  ----------------------------------------------------------------     
                                                                        
      LBIT = 0                                                          
      CALL PKB(IBYT,16,IBAY,LBIT)                                       
      CALL MVB(IBAY,1,MBAY(1,LUN),MBYT(LUN)-3,IBYT)                     
                                                                        
C  UPDATE THE SUBSET AND BYTE COUNTERS                                  
C  --------------------------------------                               
                                                                        
      MBYT(LUN)   = MBYT(LUN)   + IBYT                                  
      NSUB(LUN)   = NSUB(LUN)   + 1                                     
                                                                        
      LBIT = (NBY0+NBY1+NBY2+4)*8                                       
      CALL PKB(NSUB(LUN),16,MBAY(1,LUN),LBIT)                           
                                                                        
      LBYT = NBY0+NBY1+NBY2+NBY3                                        
      NBYT = IUPB(MBAY(1,LUN),LBYT+1,24)                                
      LBIT = LBYT*8                                                     
      CALL PKB(NBYT+IBYT,24,MBAY(1,LUN),LBIT)                           
                                                                        
C  RESET THE USER ARRAYS AND EXIT NORMALLY                              
C  ---------------------------------------
                                                                        
      CALL USRTPL(LUN,1,1)                                              
      RETURN                                                            

C  ON ENCOUTERING OVERLARGE REPORTS RESET THE USER ARRAYS AND EXIT 'GRACEFULLY'
C  ----------------------------------------------------------------------------

900   PRINT*,'MSGUPD - SUBSET LONGER THAN ANY POSSIBLE MESSAGE   '
      PRINT*,'>>>>>>>OVERLARGE SUBSET DISCARDED FROM FILE<<<<<<<<'
      PRINT*,'MSGUPD - SUBSET LONGER THAN ANY POSSIBLE MESSAGE   '
      RETURN                                                            

      END                                                               
