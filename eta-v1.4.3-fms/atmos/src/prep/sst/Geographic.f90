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
    REAL   (KIND=R4KIND)                                                    :: STLM
    REAL   (KIND=R4KIND)                                                    :: CTLM
    REAL   (KIND=R4KIND)                                                    :: STPH0
    REAL   (KIND=R4KIND)                                                    :: TPH0
    REAL   (KIND=R4KIND)                                                    :: anum
    REAL   (KIND=R4KIND)                                                    :: denom
    REAL   (KIND=R4KIND)                                                    :: relm
    INTEGER(KIND=R4KIND)                                                    :: I
    INTEGER(KIND=R4KIND)                                                    :: J
!
!
    WBD=-(IM-1)*DLMD
    SBD=-((JM-1)*0.5)*DPHD
    tph0=tph0d*dtr                                                    
    stph0=sin(tph0)                                                   
    ctph0=cos(tph0)  
    wb=wbd*dtr                                                        
    sb=sbd*dtr                                                        
    dlm=dlmd*dtr                                                      
    dph=dphd*dtr                                                      
    tdlm=dlm+dlm                                                      
    tdph=dph+dph                                                      

    tph=sb-dph

    IF (DIAG) write(6,*) WBD,SBD, WB, SB, DLM, DPH,  dph ,  tph                                      
    DO J=1,JM                                            
      tlm=wb-tdlm+mod(J+1,2)*dlm                                
      tph=tph+dph                                               
      stph=sin(tph)                                             
      ctph=cos(tph)                                           
      DO I=1,IM
        tlm=tlm+tdlm
        stlm=sin(tlm)                                             
        ctlm=cos(tlm)                                                     
        sinphi=ctph0*stph+stph0*ctph*ctlm                             
        glatr(I,J)=asin(sinphi)/DTR
        anum=ctph*stlm
        denom=(ctlm*ctph-stph0*sinphi)/ctph0
        relm=atan2(anum,denom)
        glonr(I,J)=relm/dtr+tlm0d
!
        if(glonr(I,J).gt. 180.)    glonr(I,J)=glonr(I,J)-360.
        if(glonr(I,J).lt.-180.)    glonr(I,J)=glonr(I,J)+360.
!                      
        IF (DIAG) write(6,*)I,J, GLATR(I,J),GLONR(I,J)
       enddo
     enddo

     END SUBROUTINE Geographic





