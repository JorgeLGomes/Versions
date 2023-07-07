C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION NEMOCK(NEMO)                                             
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*38  CHRSET                                              
                                                                        
      DATA CHRSET /'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.'/            
      DATA NCHR   /38/                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  GET THE LENGTH OF NEMO                                               
C  ----------------------                                               
                                                                        
      LNEMO = 0                                                         
                                                                        
      DO I=LEN(NEMO),1,-1                                               
      IF(NEMO(I:I).NE.' ') THEN                                         
         LNEMO = I                                                      
         GOTO 1                                                         
      ENDIF                                                             
      ENDDO                                                             
                                                                        
1     IF(LNEMO.LT.1 .OR. LNEMO.GT.8) THEN                               
         NEMOCK = -1                                                    
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  SCAN NEMO FOR ALLOWABLE CHARACTERS                                   
C  ----------------------------------                                   
                                                                        
      DO 10 I=1,LNEMO                                                   
      DO J=1,NCHR                                                       
      IF(NEMO(I:I).EQ.CHRSET(J:J)) GOTO 10                              
      ENDDO                                                             
      NEMOCK = -1                                                       
      RETURN                                                            
10    ENDDO                                                             
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      NEMOCK = 0                                                        
      RETURN                                                            
      END                                                               
