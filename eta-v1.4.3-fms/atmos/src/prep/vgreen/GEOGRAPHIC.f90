    SUBROUTINE Geographic(GLATR, GLONR,IM,JM,DLMD,DPHD,TLM0D,TPH0D) 
    USE CONSTANTS
    USE DIAGNOSTIC
!
    IMPLICIT NONE
!
    INTEGER(KIND=R4KIND)                          , INTENT(IN)              :: IM
    INTEGER(KIND=R4KIND)                          , INTENT(IN)              :: JM  
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: DLMD
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: DPHD
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: TLM0D
    REAL   (KIND=R4KIND)                          , INTENT(IN)              :: TPH0D
    REAL   (KIND=R4KIND), DIMENSION(IM,JM)        , INTENT(INOUT)           :: GLATR
    REAL   (KIND=R4KIND), DIMENSION(IM,JM)        , INTENT(INOUT)           :: GLONR    
    INTEGER(KIND=R4KIND)                                                    :: ICTPH0
    REAL   (KIND=R4KIND)                                                    :: SINPHI
    REAL   (KIND=R4KIND)                                                    :: CTPH0
    REAL   (KIND=R4KIND)                                                    :: DPH
    REAL   (KIND=R4KIND)                                                    :: TPH
    REAL   (KIND=R4KIND)                                                    :: SB
    REAL   (KIND=R4KIND)                                                    :: SBD
    REAL   (KIND=R4KIND)                                                    :: TLM
    REAL   (KIND=R4KIND)                                                    :: WB
    REAL   (KIND=R4KIND)                                                    :: WBD
    REAL   (KIND=R4KIND)                                                    :: TDLM
    REAL   (KIND=R4KIND)                                                    :: TDPH
    REAL   (KIND=R4KIND)                                                    :: DLM
    REAL   (KIND=R4KIND)                                                    :: STPH
    REAL   (KIND=R4KIND)                                                    :: CTPH
    REAL   (KIND=R4KIND)                                                    :: SINPH
    REAL   (KIND=R4KIND)                                                    :: STPH0
    REAL   (KIND=R4KIND)                                                    :: COSLAM
    REAL   (KIND=R4KIND)                                                    :: TPH0
    REAL   (KIND=R4KIND)                                                    :: FACT
    REAL   (KIND=R4KIND)                                                    :: ELAT
    REAL   (KIND=R4KIND)                                                    :: ELON
    REAL   (KIND=R4KIND)                                                    :: DIF
    INTEGER(KIND=R4KIND)                                                    :: I
    INTEGER(KIND=R4KIND)                                                    :: J
!
!
    WBD=-(IM-1)*DLMD
    SBD=-((JM-1)*0.5)*DPHD
    tph0=tph0d*dtr                                                    
    wb=wbd*dtr                                                        
    sb=sbd*dtr                                                        
    dlm=dlmd*dtr                                                      
    dph=dphd*dtr                                                      
    tdlm=dlm+dlm                                                      
    tdph=dph+dph                                                      

    tph=sb-dph

    stph0=sin(tph0)                                                   
    ctph0=cos(tph0)  
    IF (DIAG) write(6,*) WBD,SBD, WB, SB, DLM, DPH,  dph ,  tph                                      
    DO J=1,JM                                            
      tlm=wb-tdlm+mod(J+1,2)*dlm                                
      tph=tph+dph                                               
      stph=sin(tph)                                             
      ctph=cos(tph)                                           
      DO I=1,IM
        tlm=tlm+tdlm                                                      
        sinphi=ctph0*stph+stph0*ctph*cos(tlm)                             
        glatr(I,J)=asin(sinphi)
        coslam=ctph*cos(tlm)/(cos(glatr(I,J))*ctph0)-tan(glatr(I,J))*tan(tph0) 
        coslam=min(coslam,1.)
        fact=1.                                                           
        if (tlm .gt. D00) fact=-1.
        glonr(I,J)=-tlm0d*dtr+fact*acos(coslam) 
        IF (DIAG) write(6,*)I,J, H90+GLATR(I,J)/DTR,H360-GLONR(I,J)/DTR
       enddo
     enddo

     END SUBROUTINE Geographic





