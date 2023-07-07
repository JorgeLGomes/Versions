C---------------------------------------------------------------------- 
C  CHARACTER TRANSFER TO A STRING WITH EBCDIC TRANSLATION               
C---------------------------------------------------------------------- 
      SUBROUTINE CHRTRNA(STR,CHR,N)                                     
                                                                        
      COMMON /CHARAC/ IASCII,IATOE(0:255),IETOA(0:255)                  
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*1   CHR(N)                                              
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      DO I=1,N                                                          
      STR(I:I) = CHR(I)                                                 
      IF(IASCII.EQ.0) CALL IPKM(STR(I:I),1,IATOE(IUPM(STR(I:I),8)))     
      ENDDO                                                             
      RETURN                                                            
      END                                                               
