C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE NENUCK(NEMO,NUMB,LUN)                                  
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*8   NEMO                                                
      CHARACTER*6   NUMB                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK TABLE A                                                        
C  -------------                                                        
                                                                        
      ENTRY NENUAA(NEMO,NUMB,LUN)                                       
                                                                        
      DO N=1,NTBA(LUN)                                                  
      IF(NUMB(4:6).EQ.TABA(N,LUN)(1: 3)) GOTO 900                       
      IF(NEMO     .EQ.TABA(N,LUN)(4:11)) GOTO 900                       
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  CHECK TABLE B AND D                                                  
C  -------------------                                                  
                                                                        
      ENTRY NENUBD(NEMO,NUMB,LUN)                                       
                                                                        
      DO N=1,NTBB(LUN)                                                  
      IF(NUMB.EQ.TABB(N,LUN)(1: 6)) GOTO 900                            
      IF(NEMO.EQ.TABB(N,LUN)(7:14)) GOTO 900                            
      ENDDO                                                             
                                                                        
      DO N=1,NTBD(LUN)                                                  
      IF(NUMB.EQ.TABD(N,LUN)(1: 6)) GOTO 900                            
      IF(NEMO.EQ.TABD(N,LUN)(7:14)) GOTO 900                            
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  ----------                                                           
                                                                        
900   CALL BORT('NENUCK - DUPLICATE NEM/NUM '//NEMO//' '//NUMB)        
      END                                                               
