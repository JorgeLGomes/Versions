C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CLOSBF(LUNIT)                                          
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.GT.0 .AND. IM.NE.0) CALL CLOSMG(LUNIT)                      
      CALL WTSTAT(LUNIT,LUN,0,0)                                        
      CLOSE(LUNIT)                                                      
      RETURN                                                            
      END                                                               
