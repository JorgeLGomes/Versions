C---------------------------------------------------------------------- 
C  REAL NUMBER FROM A STRING                                            
C---------------------------------------------------------------------- 
      FUNCTION VALX(STR)                                                
                                                                        
      CHARACTER*(*) STR
      CHARACTER*99  BSTR
      CHARACTER*8   FMT
      REAL*8        BMISS                                           
                                                                        
      DATA BMISS /10E10/
      data noinline /0/                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      LENS = LEN(STR)
      IF(LENS.GT.99) CALL BORT('VALX - ARG TOO LONG')
      BSTR(1:LENS) = STR            
      RJ = RJUST(BSTR(1:LENS))
      WRITE(FMT,'(''(F'',I2,''.0)'')') LENS                             
      READ(BSTR,FMT,ERR=900) VAL                                         
      VALX = VAL                                                        
      RETURN                                                            
900   VALX = BMISS                                                      
      RETURN                                                            
      END                                                               
