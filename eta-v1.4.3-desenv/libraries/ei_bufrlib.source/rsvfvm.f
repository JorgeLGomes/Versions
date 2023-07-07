C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RSVFVM(NEM1,NEM2)                                      
                                                                        
      CHARACTER*8 NEM1,NEM2                                             
                                                                        
      DO I=1,LEN(NEM1)                                                  
      IF(I.EQ.1) THEN                                                   
         J = 1                                                          
      ELSE                                                              
         IF(NEM1(I:I).EQ.'.') THEN                                      
            NEM1(I:I) = NEM2(J:J)                                       
            J = J+1                                                     
         ENDIF                                                          
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
