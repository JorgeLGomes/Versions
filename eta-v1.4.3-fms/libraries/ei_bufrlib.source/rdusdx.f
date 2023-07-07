C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RDUSDX(LUNDX,LUN)                                      
                                                                        
      COMMON /TABABD/ NTBA(0:32),NTBB(0:32),NTBD(0:32),MTAB(50,32),     
     .                IDNA(50,32,2),IDNB(250,32),IDND(250,32),          
     .                TABA(50,32),TABB(250,32),TABD(250,32)             
                                                                        
      CHARACTER*600 TABD                                                
      CHARACTER*128 TABB                                                
      CHARACTER*128 TABA                                                
      CHARACTER*80  CARD                                                
      CHARACTER*8   NEMO                                  
      CHARACTER*6   NUMB                                                
      LOGICAL       DIGIT
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE THE DX-TABLE PARTITION AND SOURCE FILE                    
C  -------------------------------------------------                    
                                                                        
      CALL DXINIT(LUN,1)                                                
      REWIND LUNDX                                                      
                                                                        
C  READ USER CARDS UNTIL THERE ARE NO MORE                              
C  ---------------------------------------                              
                                                                        
1     READ(LUNDX,'(A80)',END=100) CARD                                  
                                                                        
C  REREAD IF NOT A DEFINITION CARD                                      
C  -------------------------------                                      
                                                                        
      IF(CARD(1: 1).EQ.       '*') GOTO 1                               
      IF(CARD(3:10).EQ.'--------') GOTO 1                               
      IF(CARD(3:10).EQ.'        ') GOTO 1                               
      IF(CARD(3:10).EQ.'MNEMONIC') GOTO 1                               
      IF(CARD(3:10).EQ.'TABLE  D') GOTO 1                               
      IF(CARD(3:10).EQ.'TABLE  B') GOTO 1                               
                                                                        
C  PARSE A DESCRIPTOR DEFINITION CARD                                   
C  ----------------------------------                                   
                                                                        
      IF(CARD(12:12).EQ.'|' .AND. CARD(21:21).EQ.'|') THEN              
                                                                        
         NEMO = CARD(3:10)                                              
         NUMB = CARD(14:19)                                             
         IF(NEMOCK(NEMO).NE.0) GOTO 900                                 
         IF(NUMBCK(NUMB).NE.0) GOTO 900                                 
                                                                        
         IF(NUMB(1:1).EQ.'A') THEN                                      
            N = NTBA(LUN)+1                                             
            IF(N.GT.NTBA(0)) GOTO 901                                   
            CALL NENUAA(NEMO,NUMB,LUN)                             
            TABA(N,LUN)( 1: 3) = NUMB(4:6)                              
            TABA(N,LUN)( 4:11) = NEMO                                   
            TABA(N,LUN)(13:67) = CARD(23:77)                            
            NTBA(LUN) = N                                               
                                                                        
            IF(DIGIT(NEMO(3:8))) THEN
               READ(NEMO,'(2X,2I3)') MTYP,MSBT                            
               IDNA(N,LUN,1) = MTYP                                     
               IDNA(N,LUN,2) = MSBT              
            ELSE
               READ(NUMB(4:6),'(I3)') IDNA(N,LUN,1)                           
               IDNA(N,LUN,2) = 0                                              
            ENDIF                                                             
                                                                        
            NUMB(1:1) = '3'                                             
         ENDIF                                                          
                                                                        
         IF(NUMB(1:1).EQ.'0') THEN                                      
            N = NTBB(LUN)+1                                             
            IF(N.GT.NTBB(0)) GOTO 902                                   
            CALL NENUBD(NEMO,NUMB,LUN)                                  
            IDNB(N,LUN) = IFXY(NUMB)                                    
            TABB(N,LUN)( 1: 6) = NUMB                                   
            TABB(N,LUN)( 7:14) = NEMO                                   
            TABB(N,LUN)(16:70) = CARD(23:77)                            
            NTBB(LUN) = N                                               
            GOTO 1                                                      
         ENDIF                                                          
                                                                        
         IF(NUMB(1:1).EQ.'3') THEN                                      
            N = NTBD(LUN)+1                                             
            IF(N.GT.NTBD(0)) GOTO 903                                   
            CALL NENUBD(NEMO,NUMB,LUN)                                  
            IDND(N,LUN) = IFXY(NUMB)                                    
            TABD(N,LUN)( 1: 6) = NUMB                                   
            TABD(N,LUN)( 7:14) = NEMO                                   
            TABD(N,LUN)(16:70) = CARD(23:77)                            
            NTBD(LUN) = N                                               
            GOTO 1                                                      
         ENDIF                                                          
                                                                        
         GOTO 904                                                       
      ENDIF                                                             
                                                                        
C  PARSE A SEQUENCE DEFINITION CARD                                     
C  --------------------------------                                     
                                                                        
      IF(CARD(12:12).EQ.'|' .AND. CARD(19:19).NE.'|') THEN              
         CALL SEQSDX(CARD,LUN)                                          
         GOTO 1                                                         
      ENDIF                                                             
                                                                        
C  PARSE AN ELEMENT DEFINITION CARD                                     
C  --------------------------------                                     
                                                                        
      IF(CARD(12:12).EQ.'|' .AND. CARD(19:19).EQ.'|') THEN              
         CALL ELEMDX(CARD,LUN)                                          
         GOTO 1                                                         
      ENDIF                                                             
                                                                        
C  CAN'T FIGURE OUT WHAT KIND OF CARD IT IS                             
C  ----------------------------------------                             
                                                                        
      GOTO 905                                                          
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
100   CALL MAKESTAB                                                     
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  ----------                                                           
                                                                        
900   PRINT*,CARD                                                       
      CALL BORT('RDUSDX - NEMO OR NUMB ERROR             '//CARD)      
901   CALL BORT('RDUSDX - TOO MANY TABLE A ENTRIES       '//CARD)      
902   CALL BORT('RDUSDX - TOO MANY TABLE B ENTRIES       '//CARD)      
903   CALL BORT('RDUSDX - TOO MANY TABLE D ENTRIES       '//CARD)      
904   CALL BORT('RDUSDX - BAD DESCRIPTOR NUMBER          '//CARD)      
905   CALL BORT('RDUSDX - BAD CARD FORMAT                '//CARD)      
      END                                                               
