C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE JSTIFY                                                 
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*1  SIGN                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ENTRY JSTCHR(STR)                                                 
                                                                        
      LSTR = LEN(STR)                                                   
                                                                        
      IF(STR.EQ.' ') GOTO 900                                           
1     IF(STR(1:1).EQ.' ') THEN                                          
         STR  = STR(2:LSTR)                                             
         GOTO 1                                                         
      ENDIF                                                             
      RETURN                                                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      ENTRY JSTNUM(STR,SIGN,IRET)                                       
                                                                        
      IRET = 0                                                          
      LSTR = LEN(STR)                                                   
                                                                        
      IF(STR.EQ.' ') GOTO 900                                           
2     IF(STR(1:1).EQ.' ') THEN                                          
         STR  = STR(2:LSTR)                                             
         GOTO 2                                                         
      ENDIF                                                             
      IF(STR(1:1).EQ.'+') THEN                                          
         STR  = STR(2:LSTR)                                             
         SIGN = '+'                                                     
      ELSEIF(STR(1:1).EQ.'-') THEN                                      
         STR  = STR(2:LSTR)                                             
         SIGN = '-'                                                     
      ELSE                                                              
         SIGN = '+'                                                     
      ENDIF                                                             
                                                                        
      CALL STRNUM(STR,NUM)                                              
      IF(NUM.LT.0) IRET = -1                                            
      RETURN                                                            
                                                                        
900   CALL BORT('JSTIFY - BLANK STRING NOT ALLOWED')                   
      END                                                               
