C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MSGWRT(LUNIT,MBAY,MBYT)                                
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /BUFRMG/ MSGLEN,MSGTXT(5000)                               
                                                                        
      CHARACTER*4 BUFR,SEVN                                             
      DIMENSION   MBAY(*)                                               
                                                                        
      DATA BUFR/'BUFR'/                                                 
      DATA SEVN/'7777'/                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  MAKE SURE ALL SECTIONS HAVE EVEN NUMBER OF BYTES                     
C  ------------------------------------------------                     
                                                                        
      IAD1 = 8                                                          
      LEN1 = IUPB(MBAY,IAD1+1,24)                                       
      LEN2 = IUPB(MBAY,IAD1+8, 1)                                       
      MTYP = IUPB(MBAY,IAD1+9, 8)                                       
      IAD2 = IAD1+LEN1                                                  
      LEN2 = IUPB(MBAY,IAD2+1,24)*LEN2                                  
      IAD3 = IAD2+LEN2                                                  
      LEN3 = IUPB(MBAY,IAD3+1,24)                                       
      IAD4 = IAD3+LEN3                                                  
      LEN4 = IUPB(MBAY,IAD4+1,24)                                       
                                                                        
      IF(MOD(LEN1,2).NE.0) GOTO 901                                     
      IF(MOD(LEN2,2).NE.0) GOTO 902                                     
      IF(MOD(LEN3,2).NE.0) GOTO 903                                     
      IF(MOD(LEN4,2).NE.0) THEN                                         
         IAD5 = IAD4+LEN4                                               
         IBIT = IAD4*8                                                  
         LEN4 = LEN4+1                                                  
         CALL PKB(LEN4,24,MBAY,IBIT)                                    
         IBIT = IAD5*8                                                  
         CALL PKB(0,8,MBAY,IBIT)                                        
         MBYX = MBYT+1                                                  
      ELSE                                                              
         MBYX = MBYT                                                    
      ENDIF                                                             
                                                                        
C  WRITE SECTION 0 BYTE COUNT AND SECTION 5                             
C  ----------------------------------------                             
                                                                        
      IBIT = 0                                                          
      KBIT = (MBYX-4)*8                                                 
                                                                        
      CALL PKC(BUFR, 4,MBAY,IBIT)                                       
      CALL PKB(MBYX,24,MBAY,IBIT)                                       
      CALL PKC(SEVN, 4,MBAY,KBIT)                                       
                                                                        
C  ZERO OUT THE EXTRA BYTES WHICH WILL BE WRITTEN                       
C  ----------------------------------------------                       
                                                                        
      IMSG = 8/NBYTW                                                    
      MWRD = (MBYX/8+1)*IMSG                                            
      MBZZ = MWRD*NBYTW-MBYX
      DO I=1,MBZZ
      CALL PKB(0,8,MBAY,KBIT)                                       
      ENDDO
                                                                        
C  WRITE THE MESSAGE PLUS PADDING TO A WORD BOUNDARY                    
C  -------------------------------------------------                    
                                                                        
      IMSG = 8/NBYTW                                                    
      MWRD = (MBYX/8+1)*IMSG                                            
      WRITE(LUNIT) (MBAY(I),I=1,MWRD)                                   
C     PRINT*,'MSGWRT - LUNIT=',LUNIT,' BYTES=',MBYX                     
                                                                        
C  save a memory copy of this message - no bufr tables though           
C  ----------------------------------------------------------           
                                                                        
      IF(MTYP.NE.11) THEN                                               
         MSGLEN = MWRD                                                  
         DO I=1,MSGLEN                                                  
         MSGTXT(I) = MBAY(I)                                            
         ENDDO                                                          
      ENDIF                                                             
                                                                        
      RETURN                                                            
901   CALL BORT('MSGWRT - UNEVEN SECTION 1')                           
902   CALL BORT('MSGWRT - UNEVEN SECTION 2')                           
903   CALL BORT('MSGWRT - UNEVEN SECTION 3')                           
      END                                                               
