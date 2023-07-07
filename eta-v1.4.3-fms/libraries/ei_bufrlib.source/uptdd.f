C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UPTDD(ID,LUN,IENT,IRET)                                
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*56  DXSTR                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LDD = LDXD(IDXV+1)+1                                              
                                                                        
C  CHECK IF IENT IS IN BOUNDS                                           
C  --------------------------                                           
                                                                        
      NDSC = IUPM(TABD(ID,LUN)(LDD:LDD),8)                              
                                                                        
      IF(IENT.EQ.0) THEN                                                
         IRET = NDSC                                                    
         RETURN                                                         
      ELSEIF(IENT.LT.0 .OR. IENT.GT.NDSC) THEN                          
         CALL BORT('UPTDD - IENT OUT OF RANGE')                        
      ENDIF                                                             
                                                                        
C  RETURN THE DESCRIPTOR INDICATED BY IENT                              
C  ---------------------------------------                              
                                                                        
      IDSC = LDD+1 + (IENT-1)*2                                         
      IRET = IUPM(TABD(ID,LUN)(IDSC:IDSC),16)                           
                                                                        
      RETURN                                                            
      END                                                               
