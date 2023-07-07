C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PKTDD(ID,LUN,IDN,IRET)                                 
                                                                        
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
                                                                        
C  ZERO THE COUNTER IF IDN IS ZERO                                      
C  -------------------------------                                      
                                                                        
      IF(IDN.EQ.0) THEN                                                 
         CALL IPKM(TABD(ID,LUN)(LDD:LDD),1,0)                           
         IRET = 0                                                       
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  UPDATE THE STORED DESCRIPTOR COUNT FOR THIS TABLE D ENTRY            
C  ---------------------------------------------------------            
                                                                        
      ND = IUPM(TABD(ID,LUN)(LDD:LDD),8)                                
                                                                        
      IF(ND.LT.0 .OR. ND.EQ.250) THEN                                   
         IRET = -1                                                      
         RETURN                                                         
      ELSE                                                              
         ND = ND+1                                                      
         CALL IPKM(TABD(ID,LUN)(LDD:LDD),1,ND)                          
         IRET = ND                                                      
      ENDIF                                                             
                                                                        
C  PACK AND STORE THE DESCRIPTOR                                        
C  -----------------------------                                        
                                                                        
      IDM = LDD+1 + (ND-1)*2                                            
      CALL IPKM(TABD(ID,LUN)(IDM:IDM),2,IDN)                            
                                                                        
      RETURN                                                            
      END                                                               
