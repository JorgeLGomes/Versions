C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE STATUS(LUNIT,LUN,IL,IM)                                

      PARAMETER (NFILES=32)
                                                                        
      COMMON /STBFR/ IOLUN(NFILES),IOMSG(NFILES)                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IF(LUNIT.LE.0 .OR. LUNIT.GT.99) GOTO 900                          
                                                                        
C  CLEAR THE STATUS INDICATORS                                          
C  ---------------------------                                          
                                                                        
      LUN = 0                                                           
      IL  = 0                                                           
      IM  = 0                                                           
                                                                        
C  SEE IF THE UNIT IS DEFINED                                           
C  --------------------------                                           
                                                                        
      DO I=1,NFILES                                                     
      IF(ABS(IOLUN(I)).EQ.LUNIT) LUN = I                                
      ENDDO                                                             
                                                                        
C  IF NOT, CHECK FOR FILE SPACE - RETURN LUN=0 IF NO FILE SPACE         
C  ------------------------------------------------------------         
                                                                        
      IF(LUN.EQ.0) THEN                                                 
         DO I=1,NFILES                                                  
         IF(IOLUN(I).EQ.0) LUN = I                                      
         IF(IOLUN(I).EQ.0) RETURN                                       
         ENDDO                                                          
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  IF FILE DEFINED RETURN STATUSES                                      
C  -------------------------------                                      
                                                                        
      IL = SIGN(1,IOLUN(LUN))                                           
      IM = IOMSG(LUN)                                                   
                                                                        
      RETURN                                                            
900   CALL BORT('STATUS - ILLEGAL UNIT GIVEN')                         
      END                                                               
