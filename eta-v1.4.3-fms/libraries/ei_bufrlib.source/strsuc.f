C-----------------------------------------------------------------------
C  DEAD SPACE FORM A STRING                                             
C-----------------------------------------------------------------------
      SUBROUTINE STRSUC(STR1,STR2,LENS)                                 
                                                                        
      CHARACTER*(*) STR1,STR2                                           
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LENS = 0                                                          
      LSTR = LEN(STR1)                                                  
                                                                        
      DO I=1,LSTR                                                       
      IF(STR1(I:I).NE.' ') GOTO 2                                       
      ENDDO                                                             
      RETURN                                                            
                                                                        
2     DO J=I,LSTR                                                       
      IF(STR1(J:J).EQ.' ') GOTO 3                                       
      LENS = LENS+1                                                     
      STR2(LENS:LENS) = STR1(J:J)                                       
      ENDDO                                                             
      RETURN                                                            
                                                                        
3     DO I=J,LSTR                                                       
      IF(STR1(I:I).NE.' ') LENS = -1                                    
      ENDDO                                                             
      RETURN                                                            
                                                                        
      END                                                               
