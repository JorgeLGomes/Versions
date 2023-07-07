C-----------------------------------------------------------------------
C  INTEGER FROM A STRING                                                
C-----------------------------------------------------------------------
      SUBROUTINE STRNUM(STR,NUM)                                        
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*20  STR2                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NUM = 0                                                           
      K = 0                                                             
                                                                        
      CALL STRSUC(STR,STR2,NUM)                                         
                                                                        
      DO I=1,NUM                                                        
      READ(STR(I:I),'(I1)',ERR=99) J                                    
      IF(J.EQ.0 .AND. STR(I:I).NE.'0') GOTO 99                          
      K = K*10+J                                                        
      ENDDO                                                             
                                                                        
      NUM = K                                                           
      RETURN                                                            
                                                                        
99    NUM = -1                                                          
      RETURN                                                            
      END                                                               
