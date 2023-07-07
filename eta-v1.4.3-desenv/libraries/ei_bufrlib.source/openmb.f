C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE OPENMB(LUNIT,SUBSET,JDATE)                             
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
                                                                        
      CHARACTER*(*) SUBSET                                              
      LOGICAL       OPEN                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
                                                                        
C  GET SOME SUBSET PARTICULARS                                          
C  ---------------------------                                          
                                                                         
      CALL NEMTBA(LUN,SUBSET,MTYP,MSTB,INOD)                            
      OPEN = IM.EQ.0.OR.INOD.NE.INODE(LUN).OR.I4DY(JDATE).NE.IDATE(LUN)   
                                                                        
C  MAYBE OPEN A NEW OR DIFFERENT TYPE OF MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      IF(OPEN) THEN                                                     
         CALL CLOSMG(LUNIT)                                             
         CALL WTSTAT(LUNIT,LUN,IL, 1)                                   
         INODE(LUN) = INOD                                              
         IDATE(LUN) = I4DY(JDATE)
         CALL MSGINI(LUN)                                               
         CALL USRTPL(LUN,1,1)                                           
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CALL BORT('OPENMB - FILE IS CLOSED            ')                 
901   CALL BORT('OPENMB - FILE IS OPEN FOR INPUT    ')                 
      END                                                               
