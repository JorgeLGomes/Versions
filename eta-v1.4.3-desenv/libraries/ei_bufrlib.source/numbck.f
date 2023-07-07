C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION NUMBCK(NUMB)                                             
                                                                        
      CHARACTER*10 CHRSET                                               
      CHARACTER*6  NUMB                                                 
      CHARACTER*1  FC                                                   
      LOGICAL      DIGIT
                                                                        
      DATA CHRSET /'0123456789'/                                        
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NUMBCK = 0                                                        
      LNUMB  = 0                                                        
      FC     = NUMB(1:1)                                                
                                                                        
C  CHECK THE FIRST CHARACTER OF NUMB                                    
C  ---------------------------------                                    
                                                                        
      IF(.NOT.(FC.EQ.'A' .OR. FC.EQ.'0' .OR. FC.EQ.'3')) GOTO 900       
                                                                        
C  CHECK THE REST OF NUMB                                               
C  ----------------------                                               
                                                                        
      DO 10 I=2,6                                                       
      DO J=1,10                                                         
      IF(NUMB(I:I).EQ.CHRSET(J:J)) GOTO 10                              
      ENDDO                                                             
      GOTO 900                                                          
10    ENDDO                                                             
                                                                        
C  CHECK FOR A VALID DESCRIPTOR                                         
C  ----------------------------                                         
                                                                        
      IF(DIGIT(NUMB(2:6))) THEN
         READ(NUMB,'(1X,I2,I3)') IX,IY                             
      ELSE
         GOTO 900
      ENDIF
      IF(IX.LT.0 .OR. IX.GT. 63) GOTO 900                               
      IF(IY.LT.0 .OR. IY.GT.255) GOTO 900                               
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      NUMBCK = 0                                                        
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  ----------                                                           
                                                                        
900   NUMBCK = -1                                                       
      RETURN                                                            
      END                                                               
