C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE OPENMG(LUNIT,SUBSET,JDATE)                             
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
      CHARACTER*(*) SUBSET                                              
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.NE.0) CALL CLOSMG(LUNIT)                                    
      CALL WTSTAT(LUNIT,LUN,IL, 1)                                      
                                                                        
C  GET SOME SUBSET PARTICULARS                                          
C  ---------------------------                                          
                                                                        
      CALL NEMTBA(LUN,SUBSET,MTYP,MSTB,INOD)                            
      INODE(LUN) = INOD                                                 
      IDATE(LUN) = I4DY(JDATE)
                                                                        
C  INITIALIZE THE OPEN MESSAGE                                          
C  ---------------------------                                          
                                                                        
      CALL MSGINI(LUN)                                                  
      CALL USRTPL(LUN,1,1)                                              
                                                                        
      RETURN                                                            
900   CALL BORT('OPENMG - FILE IS CLOSED            ')                 
901   CALL BORT('OPENMG - FILE IS OPEN FOR INPUT    ')                 
      END                                                               
