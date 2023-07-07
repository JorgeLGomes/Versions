C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE INCTAB(ATAG,ATYP,NODE)                                 
                                                                        
      COMMON /TABLES/ MAXTAB,NTAB,TAG(15000),TYP(15000),KNT(15000),     
     .                JUMP(15000),LINK(15000),JMPB(15000),              
     .                IBT(15000),IRF(15000),ISC(15000),                 
     .                ITP(15000),VALI(15000),KNTI(15000),               
     .                ISEQ(15000,2),JSEQ(15000)                         
                                                                        
      CHARACTER*(*) ATAG,ATYP                                           
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NTAB = NTAB+1                                                     
      IF(NTAB.GT.MAXTAB) CALL BORT('INCTAB - TOO MANY ENTRIES')        
                                                                        
      TAG(NTAB) = ATAG                                                  
      TYP(NTAB) = ATYP                                                  
      NODE = NTAB                                                       
                                                                        
      RETURN                                                            
      END                                                               
