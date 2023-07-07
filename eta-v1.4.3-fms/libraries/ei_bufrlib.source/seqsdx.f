C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SEQSDX(CARD,LUN)                                       
                                                                        
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
                                                                        
      CHARACTER*80  CARD,SEQS                                           
      CHARACTER*12  ATAG,TAGS(250)                                      
      CHARACTER*8   NEMO,NEMA,NEMB                                      
      CHARACTER*3   TYPS                                                
      CHARACTER*1   REPS,TAB                                            
                                                                        
      DATA MAXTGS /250/                                                 
      DATA MAXTAG /12/                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  FIND THE SEQUENCE TAG IN TABLE D AND PARSE THE SEQUENCE STRING       
C  --------------------------------------------------------------       
                                                                        
      NEMO = CARD( 3:10)                                                
      SEQS = CARD(14:78)                                                
                                                                        
      CALL NEMTAB(LUN,NEMO,IDN,TAB,ISEQ)                                
      CALL PARSEQ(SEQS,TAGS,MAXTGS,NTAG)                                
      IF(TAB.NE.'D') GOTO 900                                           
      IF(NTAG.EQ.0 ) GOTO 900                                           
                                                                        
      DO N=1,NTAG                                                       
      ATAG = TAGS(N)                                                    
      IREP = 0                                                          
                                                                        
C  CHECK FOR REPLICATOR                                                 
C  --------------------                                                 
                                                                        
      DO I=1,5                                                          
      IF(ATAG(1:1).EQ.REPS(I,1)) THEN                                   
         DO J=2,MAXTAG                                                  
         IF(ATAG(J:J).EQ.REPS(I,2)) THEN                                
            IF(J.EQ.MAXTAG) GOTO 901                                    
            CALL STRNUM(ATAG(J+1:MAXTAG),NUMR)                          
            IF(I.EQ.1 .AND. NUMR.LE.0  ) GOTO 901                       
            IF(I.EQ.1 .AND. NUMR.GT.255) GOTO 901                       
            IF(I.NE.1 .AND. NUMR.NE.0  ) GOTO 901                       
            ATAG = ATAG(2:J-1)                                          
            IREP = I                                                    
            GOTO 1                                                      
         ENDIF                                                          
         ENDDO                                                          
         GOTO 901                                                       
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  CHECK FOR VALID TAG                                                  
C  -------------------                                                  
                                                                        
1     IF(NEMOCK(ATAG).NE.0) GOTO 901                                    
      CALL NEMTAB(LUN,ATAG,IDN,TAB,IRET)                                
      IF(IRET.GT.0) THEN                                                
         IF(TAB.EQ.'B' .AND. IREP.NE.0) GOTO 902                        
         IF(ATAG(1:1).EQ.'.') THEN                                      
            NEMB = TAGS(N+1)                                            
            CALL NUMTAB(LUN,IDN,NEMA,TAB,ITAB)                          
            CALL NEMTAB(LUN,NEMB,JDN,TAB,IRET)                          
            CALL RSVFVM(NEMA,NEMB)                                      
            IF(NEMA.NE.ATAG) GOTO 903                                   
            IF(N.GT.NTAG ) GOTO 905                                     
            IF(TAB.NE.'B') GOTO 906                                     
         ENDIF                                                          
      ELSE                                                              
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
C  WRITE THE DESCRIPTOR STRING INTO TABD ARRAY                          
C  -------------------------------------------                          
                                                                        
10    IF(IREP.GT.0) CALL PKTDD(ISEQ,LUN,IDNR(IREP,1)+NUMR,IRET)         
      IF(IRET.LT.0) GOTO 904                                            
      CALL PKTDD(ISEQ,LUN,IDN,IRET)                                     
      IF(IRET.LT.0) GOTO 904                                            
                                                                        
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('SEQSDX - UNDEFINED SEQUENCE: '             //   NEMO) 
901   CALL BORT('SEQSDX - BAD TAG IN SEQUENCE: '            //TAGS(N)) 
902   CALL BORT('SEQSDX - REPLICATED ELEMENTS NOT ALLOWED:' //TAGS(N)) 
903   CALL BORT('SEQSDX - UNDEFINED TAG: '                  //TAGS(N)) 
904   CALL BORT('SEQSDX - TOO MANY DESCRIPTORS IN STRING:'  //   NEMO) 
905   CALL BORT('SEQSDX - FOLLOWING-VALUE LAST IN STRING:'  //   NEMA) 
906   CALL BORT('SEQSDX - FOLLOWING VALUE NOT FROM TABLEB:' //   NEMB) 
      END                                                               
