C---------------------------------------------------------------------- 
C  CONVERT A SIX CHARACTER (FXY) ASCII DESCRIPTOR TO AN INTEGER         
C---------------------------------------------------------------------- 
      FUNCTION IFXY(ADSC)                                               
                                                                        
      CHARACTER*6 ADSC                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      READ(ADSC,'(I1,I2,I3)') IF,IX,IY                                  
      IFXY = IF*2**14 + IX*2**8 + IY                                    
      RETURN                                                            
      END                                                               
