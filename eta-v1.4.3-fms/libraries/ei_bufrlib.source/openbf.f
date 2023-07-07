C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE OPENBF(LUNIT,IO,LUNDX)                                 
                                                                        
      COMMON /MSGCWD/ NMSG(32),NSUB(32),MSUB(32),INODE(32),IDATE(32)    
      COMMON /STBFR / IOLUN(32),IOMSG(32)                               
      COMMON /QUIET / IPRT                                              
                                                                        
      CHARACTER*(*) IO                                                  
      CHARACTER*4   BUFR,MSTR                                           
      LOGICAL       SKIPDX,APPEND                                       
                                                                        
      DATA IFIRST/0/                                                    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

                                                                        
      IF(IFIRST.EQ.0) THEN                                              
         CALL WRDLEN                                                    
         CALL BFRINI                                                    
         IFIRST = 1                                                     
      ENDIF                                                             
                                                                        
      IF(IO.EQ.'QUIET') THEN                                            
         IPRT = LUNDX                                                   
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  SEE IF A FILE CAN BE OPENED                                          
C  ---------------------------                                          
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(LUN.EQ.0) GOTO 900                                             
      IF(IL .NE.0) GOTO 901                                             
                                                                        
C  CHECK FOR NO BUFR DATA OR NO DATA AT ALL IN AN "IN" FILE             
C  --------------------------------------------------------             
                                                                        
      IF(IO.EQ.'IN' .AND. LUNIT.EQ.LUNDX) THEN                          
         REWIND LUNIT                                                   
         READ(LUNIT,END=100,ERR=902) MSTR                               
         IBIT = 0                                                       
         CALL UPC(BUFR,4,MSTR,IBIT)                                     
         IF(BUFR.NE.'BUFR') GOTO 902                                    
      ENDIF                                                             
                                                                        
C  SET INITIAL OPEN DEFAULTS                                            
C  -------------------------                                            
                                                                        
      REWIND LUNIT                                                      
      NMSG (LUN) = 0                                                    
      NSUB (LUN) = 0                                                    
      MSUB (LUN) = 0                                                    
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
      SKIPDX = .FALSE.                                                  
      APPEND = .FALSE.                                                  
                                                                        
C  DECIDE HOW TO SETUP THE DICTIONARY                                   
C  ----------------------------------                                   
                                                                        
      IF(IO.EQ.'IN') THEN                                               
         CALL WTSTAT(LUNIT,LUN,-1,0)                                    
         CALL READDX(LUNIT,LUN,LUNDX)                                   
      ELSE IF(IO.EQ.'OUT') THEN                                         
         CALL WTSTAT(LUNIT,LUN, 1,0)                                    
         CALL WRITDX(LUNIT,LUN,LUNDX)                                   
      ELSE IF(IO.EQ.'APN' .OR. IO.EQ.'APX') THEN                        
         CALL WTSTAT(LUNIT,LUN, 1,0)                                    
         CALL READDX(LUNIT,LUN,LUNDX)                                   
         IF(IO.EQ.'APN') CALL POSAPN(LUNIT)                             
         IF(IO.EQ.'APX') CALL POSAPX(LUNIT)                             
      ELSE                                                              
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
                                                                        
C  FILE OPENED FOR INPUT IS EMPTY - LET READMG GIVE THE BAD NEWS        
C  -------------------------------------------------------------        
                                                                        
100   REWIND LUNIT                                                      
      CALL WTSTAT(LUNIT,LUN,-1,0)                                       
      CALL DXINIT(LUN,0)                                                
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   CALL BORT('OPENBF - TOO MANY FILES OPENED ALREADY       ')       
901   CALL BORT('OPENBF - FILE ALREADY OPEN                   ')       
902   CALL BORT('OPENBF - INPUT FILE HAS NON-BUFR DATA        ')       
903   CALL BORT('OPENBF - IO MUST BE ONE OF "IN" "OUT" "APN"  ')       
      END                                                               
