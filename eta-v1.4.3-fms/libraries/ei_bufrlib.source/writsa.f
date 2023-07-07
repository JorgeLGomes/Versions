C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
      SUBROUTINE WRITSA(LUNXX,MSGT,MSGL)                                
                                                                        
      COMMON /BUFRMG/ MSGLEN,MSGTXT(5000)                               
                                                                        
      DIMENSION MSGT(*)                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      LUNIT = ABS(LUNXX)                                                
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.EQ.0) GOTO 902                                              
                                                                        
C  SEE IF A MEMORY MESSAGE IS WAITING OR FORCED                         
C  --------------------------------------------                         
                                                                        
      IF(LUNXX.LT.0) CALL CLOSMG(LUNIT)                                 
                                                                        
      IF(MSGLEN.GT.0) THEN                                              
         MSGL = MSGLEN                                                  
         DO N=1,MSGL                                                    
         MSGT(N) = MSGTXT(N)                                            
         ENDDO                                                          
         MSGLEN = 0                                                     
      ELSE                                                              
         MSGL = 0                                                       
      ENDIF                                                             
                                                                        
      IF(LUNXX.LT.0) RETURN                                             
                                                                        
C  PACK UP THE SUBSET AND PUT IT INTO THE MESSAGE                       
C  ----------------------------------------------                       
                                                                        
      CALL WRTREE(LUN)                                                  
      CALL MSGUPD(LUNIT,LUN)                                            
                                                                        
      RETURN                                                            
900   CALL BORT('WRITSA - FILE IS CLOSED                     ')        
901   CALL BORT('WRITSA - FILE IS OPEN FOR INPUT             ')        
902   CALL BORT('WRITSA - NO MESSAGE OPEN                    ')        
      END                                                               
