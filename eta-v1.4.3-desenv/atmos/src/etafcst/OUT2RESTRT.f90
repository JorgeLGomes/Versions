    SUBROUTINE OUT2RESTRT(ITAG)
!>--------------------------------------------------------------------------------------------------
!> SUBROUTINE OUT2RESTRT
!>
!> SUBPROGRAM: OUT2RESTRT - 
!>
!> DRIVER     : CHKOUT
!>              
!>
!>--------------------------------------------------------------------------------------------------
!    USE ABCI        !andre check
    USE ACMCLD      !andre check
    USE ACMCLH      !andre check
    USE ACMPRE      !andre check
    USE ACMRDL      !andre check
    USE ACMRDS      !andre check
    USE ACMSFC      !andre check
    USE ASTSAV      !andre check
    USE BANDTA      !andre check
    USE BDCOMB      !andre check
    USE BDWIDE      !andre check
    USE BOCO        !andre check
    USE BUFFER      !andre check
    USE C_FRACN     !andre check
    USE CLDWTR      !andre check
    USE CMICRO_CONS !andre check
    USE CMICRO_START!andre check
    USE CMICRO_STATS!andre check
    USE CMY600      !andre check  
    USE CNVCLD      !andre check
    USE CO2BD2      !andre check
    USE CO2BD3      !andre check
    USE CO2BD4      !andre check
    USE CO2BD5      !andre check
    USE COMPVS0     !andre check    
    USE COMPVS      !andre check
    USE CONTIN      !andre check   
    USE C_TADJ      !andre check    
    USE CTLBLK      !andre check
    USE CUINIT      !andre check
!    USE DIUCON      !andre check
    USE DYNAM       !andre check
!    USE MODULE EXCH_BUF_REAL!andre check
    USE F77KINDS
    USE GLB_TABLE   !andre check
    USE HCON
    USE IACCR_TABLES!andre check
    USE IMASS_TABLES!andre check
    USE INDX        !andre check
    USE INPUT       !andre check
    USE IRATE_TABLES!andre check
    USE IRIME_TABLES!andre check
    USE IVENT_TABLES!andre check
    USE LOOPS       !andre check
    USE MAPOT       !andre check
    USE MAPPINGS    !andre check
    USE MASKS       !andre check
    USE MOMENTO     !andre check
    USE MPPCOM
    USE NHYDRO      !andre check
    USE NSOILTYPE   !andre check
    USE OPTIONS     !andre check
!    USE OUTFIL      !andre check  
!    USE PHYCON      !andre check: all variables now are parameter 
    USE PHYS        !andre check 
    USE PPTASM      !andre check
    USE PRFHLD      !andre check
    USE PVRBLS      !andre check
    USE RACCR_TABLES!andre check
    USE RD1TIM      !andre check
    USE RDFSAV      !andre check
!    USE RITE        !andre check
    USE RMASS_TABLES!andre check
    USE RRATE_TABLES!andre check   
    USE RVELR_TABLES!andre check   
    USE RVENT_TABLES!andre check
!    USE SAVMEM      !andre check
    USE SCRTCH      !andre check
    USE SDENS_TABLES!andre check
    USE SEASO3      !andre check
    USE SLOPES      !andre check
    USE SOIL        !andre check
    USE SSALB       !andre check
    USE SWRSAV      !andre check
    USE TABCOM
    USE TEMPV       !andre check
    USE TOPO        !andre check
    USE VRBLS       !andre check
    USE UPDT
    USE Z0EFFT      !andre check 
!
    IMPLICIT NONE
!
    INTEGER(KIND=I4KIND)                                                                        ::&
    & ITAG,  LISTOUT 
!
!
!
!
    CHARACTER(LEN=6)                                                                            ::&
    & C_NHRS     
!   
    CHARACTER(LEN=4)                                                                            ::&
    & C_MYPE     
!   
    CHARACTER(LEN=150)                                                                          ::&
    & OUTFILE 
!
!
    WRITE(C_NHRS,'(I6.6)') ITAG  
!
    WRITE(C_MYPE,'(I4.4)') MYPE
!            
    OUTFILE = "OUT2RESTRT/OUT2RESTRT_"//C_MYPE//"."//C_NHRS
!
    LISTOUT=1001
!
    OPEN (UNIT=LISTOUT, FILE= TRIM(OUTFILE), FORM = 'UNFORMATTED')
!
!MODULE_ABCI
!         WRITE(LISTOUT) NSOLD
!         WRITE(LISTOUT) AI
!	 WRITE(LISTOUT) BI
!	 WRITE(LISTOUT) CI
!MODULE ACMCLD
	 WRITE(LISTOUT) NCLOD    
	 WRITE(LISTOUT) TCLOD      
	 WRITE(LISTOUT) NCFRCV
	 WRITE(LISTOUT) NCFRST
	 WRITE(LISTOUT) ACFRCV
	 WRITE(LISTOUT) ACFRST  
!MODULE_ACMCLH 
         WRITE(LISTOUT) NHEAT
         WRITE(LISTOUT) THEAT
         WRITE(LISTOUT) AVRAIN
         WRITE(LISTOUT) AVCNVC
         WRITE(LISTOUT) ARATIM
         WRITE(LISTOUT) ACUTIM
         WRITE(LISTOUT) TRAIN
         WRITE(LISTOUT) TCUCN 
!MODULE_ACMPRE
	 WRITE(LISTOUT) TPREC
         WRITE(LISTOUT) ACSNOW
         WRITE(LISTOUT) ACSNOM 
         WRITE(LISTOUT) SSROFF
         WRITE(LISTOUT) BGROFF 
!MODULE_ACMRDL
	 WRITE(LISTOUT) NRDLW
	 WRITE(LISTOUT) ARDLW
	 WRITE(LISTOUT) TRDLW
         WRITE(LISTOUT) RLWIN
         WRITE(LISTOUT) RLWOUT 
         WRITE(LISTOUT) RLWTOA 
         WRITE(LISTOUT) ALWIN
         WRITE(LISTOUT) ALWOUT
         WRITE(LISTOUT) ALWTOA 
         WRITE(LISTOUT) RLWTT 
!MODULE_ACMRDS
	 WRITE(LISTOUT) NRDSW 
	 WRITE(LISTOUT) ARDSW
	 WRITE(LISTOUT) TRDSW
         WRITE(LISTOUT) RSWIN 
         WRITE(LISTOUT) RSWOUT
         WRITE(LISTOUT) RSWTOA
         WRITE(LISTOUT) ASWIN
         WRITE(LISTOUT) ASWOUT
         WRITE(LISTOUT) ASWTOA
         WRITE(LISTOUT) RSWTT
!MODULE_ACMSFC
	 WRITE(LISTOUT) NSRFC
	 WRITE(LISTOUT) ASRFC
	 WRITE(LISTOUT) TSRFC
	 WRITE(LISTOUT) APHTIM
	 WRITE(LISTOUT) SFCSHX
	 WRITE(LISTOUT) SFCLHX  								 
	 WRITE(LISTOUT) SUBSHX 
	 WRITE(LISTOUT) SNOPCX  								 
	 WRITE(LISTOUT) SFCUVX 
	 WRITE(LISTOUT) SFCEVP  								 
	 WRITE(LISTOUT) POTEVP  
	 WRITE(LISTOUT) POTFLX
!MODULE_ASTSAV
	 WRITE(LISTOUT) SOLC    
	 WRITE(LISTOUT) RSIN1
	 WRITE(LISTOUT) RCOS1
	 WRITE(LISTOUT) RCOS2
!MODULE_BANDTA
         WRITE(LISTOUT) ARNDM
	 WRITE(LISTOUT) BRNDM  
	 WRITE(LISTOUT) BETAD   
	 WRITE(LISTOUT) AP     
	 WRITE(LISTOUT) BP      
	 WRITE(LISTOUT) ATP    
	 WRITE(LISTOUT) BTP    
	 WRITE(LISTOUT) BANDLO
	 WRITE(LISTOUT) BANDHI
         WRITE(LISTOUT) AO3RND  
	 WRITE(LISTOUT) BO3RND
	 WRITE(LISTOUT) AB15
         WRITE(LISTOUT) ARNDM1
	 WRITE(LISTOUT) ARNDM2  							
         WRITE(LISTOUT) BRNDM1  
         WRITE(LISTOUT) BRNDM2  							
         WRITE(LISTOUT) AP1    
         WRITE(LISTOUT) AP2     							
         WRITE(LISTOUT) BP1    
         WRITE(LISTOUT) BP2     							
         WRITE(LISTOUT) ATP1   
         WRITE(LISTOUT) ATP2    							
         WRITE(LISTOUT) BTP1   
         WRITE(LISTOUT) BTP2    							
         WRITE(LISTOUT) BETAD1  
         WRITE(LISTOUT) BETAD2  							
         WRITE(LISTOUT) BANDL1  
         WRITE(LISTOUT) BANDL2  							
         WRITE(LISTOUT) BANDH1 
         WRITE(LISTOUT) BANDH2   
         WRITE(LISTOUT) ARNDM3  							
         WRITE(LISTOUT) BRNDM3  							
         WRITE(LISTOUT) AP3     							
         WRITE(LISTOUT) BP3     							
         WRITE(LISTOUT) ATP3    							
         WRITE(LISTOUT) BTP3    							
         WRITE(LISTOUT) BETAD3  							
         WRITE(LISTOUT) BANDL3  							
         WRITE(LISTOUT) BANDH3
!MODULE_BDCOMB
         WRITE(LISTOUT) IBAND
         WRITE(LISTOUT) ACOMB
         WRITE(LISTOUT) BCOMB
         WRITE(LISTOUT) BETACM
         WRITE(LISTOUT) APCM
         WRITE(LISTOUT) BPCM
         WRITE(LISTOUT) ATPCM
         WRITE(LISTOUT) BTPCM
         WRITE(LISTOUT) BDLOCM
         WRITE(LISTOUT) BDHICM
         WRITE(LISTOUT) AB15CM
         WRITE(LISTOUT) AO3CM 
         WRITE(LISTOUT) BO3CM
         WRITE(LISTOUT) BETINC
!MODULE_BDWIDE
         WRITE(LISTOUT) AWIDE  
         WRITE(LISTOUT) BWIDE  
         WRITE(LISTOUT) BETAWD 
         WRITE(LISTOUT) APWD   
         WRITE(LISTOUT) BPWD   
         WRITE(LISTOUT) ATPWD  
         WRITE(LISTOUT) BTPWD  
         WRITE(LISTOUT) BDLOWD 
         WRITE(LISTOUT) BDHIWD 
         WRITE(LISTOUT) BETINW 
         WRITE(LISTOUT) AB15WD 
         WRITE(LISTOUT) SKO2D  
         WRITE(LISTOUT) SKC1R	
         WRITE(LISTOUT) SKO3R
!MODULE_BOCO	 
	 WRITE(LISTOUT) PDB
	 WRITE(LISTOUT) TB
	 WRITE(LISTOUT) QB
	 WRITE(LISTOUT) UB
	 WRITE(LISTOUT) VB
	 WRITE(LISTOUT) Q2B
	 WRITE(LISTOUT) CWMB     
!MODULE_BUFFER
	 WRITE(LISTOUT) IP
	 WRITE(LISTOUT) BUF
!MODULE_C_FRACN
	 WRITE(LISTOUT) F_ICE	   									  
	 WRITE(LISTOUT) F_RAIN     									  
	 WRITE(LISTOUT) F_RIMEF
!MODULE_CLDWTR
	 WRITE(LISTOUT) CWM
	 WRITE(LISTOUT) LC
	 WRITE(LISTOUT) U00
	 WRITE(LISTOUT) SR
	 WRITE(LISTOUT) UL
!MODULE_CMICRO_CONS
	 WRITE(LISTOUT)	ABFR		  
	 WRITE(LISTOUT)	CBFR	
	 WRITE(LISTOUT) CIACW 
	 WRITE(LISTOUT)	CIACR	      
	 WRITE(LISTOUT) C_N0R0
	 WRITE(LISTOUT)	CN0R0	      
	 WRITE(LISTOUT) CN0R_DMRMIN
	 WRITE(LISTOUT)	CN0R_DMRMAX   
	 WRITE(LISTOUT) CRACW 
	 WRITE(LISTOUT)	CRAUT	      
	 WRITE(LISTOUT) ESW0  
	 WRITE(LISTOUT)	QAUT0	
	 WRITE(LISTOUT)	RFMAX	
	 WRITE(LISTOUT)	RHGRD	      
	 WRITE(LISTOUT) RQR_DR1
	 WRITE(LISTOUT)	RQR_DR2 
	 WRITE(LISTOUT)	RQR_DR3 
	 WRITE(LISTOUT)	RQR_DRMIN
	 WRITE(LISTOUT)	RQR_DRMAX
	 WRITE(LISTOUT) RR_DR1
	 WRITE(LISTOUT)	RR_DR2  
	 WRITE(LISTOUT)	RR_DR3  
	 WRITE(LISTOUT)	RR_DRMIN	  
	 WRITE(LISTOUT)	RR_DRMAX
!MODULE_CMICRO_START
         WRITE(LISTOUT)	MICRO_START
!MODULE_CMICRO_STATS
         WRITE(LISTOUT)	NSTATS
         WRITE(LISTOUT)	QMAX
         WRITE(LISTOUT)	QTOT
!MODULE_CMY600
         WRITE(LISTOUT)	MY_GROWTH
!MODULE_CNVCLD
         WRITE(LISTOUT)	CUPPT
	 WRITE(LISTOUT) CFRACL  
	 WRITE(LISTOUT) CFRACM 
	 WRITE(LISTOUT) CFRACH
!MODULE_CO2BD2
	 WRITE(LISTOUT)CO231
	 WRITE(LISTOUT)CO238
	 WRITE(LISTOUT)CDT31
	 WRITE(LISTOUT)CDT38
	 WRITE(LISTOUT)C2D31
	 WRITE(LISTOUT)C2D38
!MODULE_CO2BD3
	 WRITE(LISTOUT)CO251
	 WRITE(LISTOUT)CO258
	 WRITE(LISTOUT)CDT51
	 WRITE(LISTOUT)CDT58
	 WRITE(LISTOUT)C2D51
	 WRITE(LISTOUT)C2D58
	 WRITE(LISTOUT)CO2M51
	 WRITE(LISTOUT)CO2M58
	 WRITE(LISTOUT)CDTM51
	 WRITE(LISTOUT)CDTM58
	 WRITE(LISTOUT)C2DM51
	 WRITE(LISTOUT)C2DM58 
	 WRITE(LISTOUT)STEMP
	 WRITE(LISTOUT)GTEMP
!MODULE_CO2BD4
	 WRITE(LISTOUT) CO271
	 WRITE(LISTOUT) CO278
	 WRITE(LISTOUT) CDT71
	 WRITE(LISTOUT) CDT78
 	 WRITE(LISTOUT) C2D71
 	 WRITE(LISTOUT) C2D78
!MODULE_CO2BD5
         WRITE(LISTOUT) CO211
	 WRITE(LISTOUT) CO218
!MODULE_COMPVS0
         WRITE(LISTOUT) C1XPVS0
         WRITE(LISTOUT) C2XPVS0
         WRITE(LISTOUT) TBPVS0
!MODULE_COMPVS
         WRITE(LISTOUT) C1XPVS
         WRITE(LISTOUT) C2XPVS
         WRITE(LISTOUT) TBPVS
!MODULE_CONTIN
	 WRITE(LISTOUT) PDSL
         WRITE(LISTOUT) PDSLO
	 WRITE(LISTOUT) PSDT
	 WRITE(LISTOUT) RTOP
	 WRITE(LISTOUT) OMGALF
	 WRITE(LISTOUT) DIV
	 WRITE(LISTOUT) ETADT
!MODULE_C_TADJ
 	 WRITE(LISTOUT) T_ADJ
	 WRITE(LISTOUT) T_OLD
!MODULE_CTLBLK
 	 WRITE(LISTOUT) RUN 
 	 WRITE(LISTOUT) FIRST   
 	 WRITE(LISTOUT) RESTRT  
 	 WRITE(LISTOUT) SIGMA   
 	 WRITE(LISTOUT) NEST    
 	 WRITE(LISTOUT) SINGLRST
 	 WRITE(LISTOUT) SUBPOST
 	 WRITE(LISTOUT) IHRST   
 	 WRITE(LISTOUT) NFCST   
 	 WRITE(LISTOUT) NBC     
 	 WRITE(LISTOUT) LIST    
 	 WRITE(LISTOUT) IOUT    
 	 WRITE(LISTOUT) NTSD    
 	 WRITE(LISTOUT) NTSTM   
 	 WRITE(LISTOUT) NSTART  
 	 WRITE(LISTOUT) NTDDMP  
 	 WRITE(LISTOUT) NPREC	
 	 WRITE(LISTOUT) IDTAD	
 	 WRITE(LISTOUT) NBOCO	
 	 WRITE(LISTOUT) NSHDE	
 	 WRITE(LISTOUT) NCP	
 	 WRITE(LISTOUT) NPHS	
 	 WRITE(LISTOUT) NCNVC	
 	 WRITE(LISTOUT) NRADS	
 	 WRITE(LISTOUT) NRADL 
 	 WRITE(LISTOUT) IDAT
 	 WRITE(LISTOUT) DT
!MODULE_CUINIT
         WRITE(LISTOUT) CURAD
!MODULE_CUPARM, 
!         ONLY PARAMETER
!MODULE_DIUCON
!         WRITE(LISTOUT) SEASON
!         WRITE(LISTOUT) XXXX
!         WRITE(LISTOUT) JDNMC
!         WRITE(LISTOUT) FCSTDA
!         WRITE(LISTOUT) FJDNMC
!         WRITE(LISTOUT) TSLAG
!         WRITE(LISTOUT) RLAG
!         WRITE(LISTOUT) TIMIN
!         WRITE(LISTOUT) TPI 
!         WRITE(LISTOUT) HPI   
!         WRITE(LISTOUT) YEAR
!         WRITE(LISTOUT) DAY    
!         WRITE(LISTOUT) DHR
!         WRITE(LISTOUT) JTIME
!         WRITE(LISTOUT) DAZ
!MODULE_DYNAM
         WRITE(LISTOUT) PT
         WRITE(LISTOUT) R      
         WRITE(LISTOUT) DY
         WRITE(LISTOUT) CPGFV
         WRITE(LISTOUT) EN
         WRITE(LISTOUT) ENT   
         WRITE(LISTOUT) F4D
         WRITE(LISTOUT) F4Q   
         WRITE(LISTOUT) EF4T  
         WRITE(LISTOUT) AETA
         WRITE(LISTOUT) DETA  
         WRITE(LISTOUT) RDETA
         WRITE(LISTOUT) F4Q2  
         WRITE(LISTOUT) ETA   
         WRITE(LISTOUT) DFL
         WRITE(LISTOUT) EM 
         WRITE(LISTOUT) EMT
         WRITE(LISTOUT) DX  
         WRITE(LISTOUT) WPDAR 
         WRITE(LISTOUT) CPGFU
         WRITE(LISTOUT) CURV  
         WRITE(LISTOUT) FCP   
         WRITE(LISTOUT) FDIV  
         WRITE(LISTOUT) F     
         WRITE(LISTOUT) DDMPU 
         WRITE(LISTOUT) DDMPV 
         WRITE(LISTOUT) FAD
!MODULE_EXCH_BUF_REAL
!         WRITE(LISTOUT) BUF0
!         WRITE(LISTOUT) BUF1
!         WRITE(LISTOUT) BUF2
!         WRITE(LISTOUT) BUF3
!MODULE_GLB_TABLE
         WRITE(LISTOUT) IS_GLB_TABLE
         WRITE(LISTOUT) IE_GLB_TABLE    
         WRITE(LISTOUT) JS_GLB_TABLE
         WRITE(LISTOUT) JE_GLB_TABLE
!MODULE_HCON
!MODULE_IACCR_TABLES
         WRITE(LISTOUT) ACCRI
!MODULE_IMASS_TABLES
         WRITE(LISTOUT) MASSI
!MODULE_INDX
!         WRITE(LISTOUT)IHE
!         WRITE(LISTOUT)IVE
!         WRITE(LISTOUT)IHW
!         WRITE(LISTOUT)IHW
!         WRITE(LISTOUT)IRAD
!         WRITE(LISTOUT)IHEG
!         WRITE(LISTOUT)IVEG
!         WRITE(LISTOUT)IHWG
!         WRITE(LISTOUT)IVWG
!         WRITE(LISTOUT)IRADG
!MODULE_INPUT
         WRITE(LISTOUT) ICH     
         WRITE(LISTOUT) ICM     
         WRITE(LISTOUT) ICT     
         WRITE(LISTOUT) ICB
         WRITE(LISTOUT) CH      
         WRITE(LISTOUT) CM      
         WRITE(LISTOUT) CL      
         WRITE(LISTOUT) EMCH    
         WRITE(LISTOUT) EMCM    
         WRITE(LISTOUT) EMCL
         WRITE(LISTOUT) RR      
         WRITE(LISTOUT) QQO3
         WRITE(LISTOUT) DTEMP   
         WRITE(LISTOUT) PPRESS 
!MODULE_IRATE_TABLES
         WRITE(LISTOUT) VSNOWI
!MODULE_IRIME_TABLES
         WRITE(LISTOUT) VEL_RF
!MODULE_IVENT_TABLES	 
         WRITE(LISTOUT) VENTI1
         WRITE(LISTOUT) VENTI2	 
!MODULE_LOOPS	 
	 WRITE(LISTOUT) IHLA
	 WRITE(LISTOUT) IHHA
	 WRITE(LISTOUT) IVLA
	 WRITE(LISTOUT) IVHA
	 WRITE(LISTOUT) JRA	 
	 WRITE(LISTOUT) LMH
	 WRITE(LISTOUT) LMV
!MODULE MAPOT
         WRITE(LISTOUT) LSL
         WRITE(LISTOUT) IXM
	 WRITE(LISTOUT) IYM
         WRITE(LISTOUT) ISHDE
         WRITE(LISTOUT) TSPH
         WRITE(LISTOUT) WBD
         WRITE(LISTOUT) SBD
         WRITE(LISTOUT) TLM0D
         WRITE(LISTOUT) TPH0D
         WRITE(LISTOUT) DLMD
	 WRITE(LISTOUT) DPHD
         WRITE(LISTOUT) CMLD
         WRITE(LISTOUT) DP30
         WRITE(LISTOUT) X1P
         WRITE(LISTOUT) Y1P
         WRITE(LISTOUT) DISLP
         WRITE(LISTOUT) Z0SLP
 	 WRITE(LISTOUT) SPL		 
	 WRITE(LISTOUT) ALSL 
	 WRITE(LISTOUT) TSHDE
	 WRITE(LISTOUT) ERLAM0 
	 WRITE(LISTOUT) CPHI0	 	 
	 WRITE(LISTOUT) SPHI0
!MODULE_MAPPINGS
	 WRITE(LISTOUT) G2LI    
	 WRITE(LISTOUT) L2GI
	 WRITE(LISTOUT) G2LJ
	 WRITE(LISTOUT) L2GJ
!MODULE_MASKS
	 WRITE(LISTOUT) VBM2    
	 WRITE(LISTOUT) VBM3    
	 WRITE(LISTOUT) SM      
	 WRITE(LISTOUT) SICE    
	 WRITE(LISTOUT) HBM2    
	 WRITE(LISTOUT) HBM3
	 WRITE(LISTOUT) HTM     
	 WRITE(LISTOUT) VTM
!MODULE_MOMENTO
	 WRITE(LISTOUT) UMFLX
	 WRITE(LISTOUT) VMFLX  
	 WRITE(LISTOUT) AKMS10
!MODULE_MPPCOM
         !colocar depois
!MODULE_NHYDRO
 	 WRITE(LISTOUT) HYDRO
 	 WRITE(LISTOUT) SPLINE	 
 	 WRITE(LISTOUT) DWDT
 	 WRITE(LISTOUT) PDWDT
 	 WRITE(LISTOUT) PINT	 
 	 WRITE(LISTOUT) W
 	 WRITE(LISTOUT) Z
!MODULE_NSOILTYPE
 	 WRITE(LISTOUT) NSOTYP
!MODULE_OPTIONS
 	 WRITE(LISTOUT) SPVAL
 	 WRITE(LISTOUT) IBESSL
 	 WRITE(LISTOUT) KSB     
 	 WRITE(LISTOUT) IOFFS   
 	 WRITE(LISTOUT) IFLAG   
 	 WRITE(LISTOUT) SATDEL
!MODULE_OUTFIL
! 	 WRITE(LISTOUT) RSTFIL
! 	 WRITE(LISTOUT) ITAG   
! 	 WRITE(LISTOUT) LRSTRT
!MODULE_PHYCON
! 	 WRITE(LISTOUT) AMOLWT  
! 	 WRITE(LISTOUT) CSUBP   
! 	 WRITE(LISTOUT) DIFFCTR 
! 	 WRITE(LISTOUT) G       
! 	 WRITE(LISTOUT) GRAVDR  
! 	 WRITE(LISTOUT) O3DIFCTR
! 	 WRITE(LISTOUT) P0      
! 	 WRITE(LISTOUT) P0XZP2  
! 	 WRITE(LISTOUT) P0XZP8  
! 	 WRITE(LISTOUT) P0X2    
! 	 WRITE(LISTOUT) RADCON  
! 	 WRITE(LISTOUT) RGAS    
! 	 WRITE(LISTOUT) RGASSP  
! 	 WRITE(LISTOUT) SECPDA  
! 	 WRITE(LISTOUT) RATCO2MW
! 	 WRITE(LISTOUT) RATH2OMW                                                                       
! 	 WRITE(LISTOUT) RADCON1                                                                        
! 	 WRITE(LISTOUT) GINV    
! 	 WRITE(LISTOUT) P0INV   
! 	 WRITE(LISTOUT) GP0INV
!MODULE_PHYS
 	 WRITE(LISTOUT) KTM
 	 WRITE(LISTOUT) DTQ2    
 	 WRITE(LISTOUT) TDTQ2   
 	 WRITE(LISTOUT) DTD     
 	 WRITE(LISTOUT) TDTD    		   
 	 WRITE(LISTOUT) ROS     
 	 WRITE(LISTOUT) CS      
 	 WRITE(LISTOUT) DS      
 	 WRITE(LISTOUT) ROI     
 	 WRITE(LISTOUT) CI      
 	 WRITE(LISTOUT) DI      		   
 	 WRITE(LISTOUT) PL      
 	 WRITE(LISTOUT) THL     
 	 WRITE(LISTOUT) RDQ     
 	 WRITE(LISTOUT) RDTH    
 	 WRITE(LISTOUT) RDP     
 	 WRITE(LISTOUT) RDTHE   		   
 	 WRITE(LISTOUT) PLQ     
 	 WRITE(LISTOUT) RDPQ    
 	 WRITE(LISTOUT) RDTHEQ
 	 WRITE(LISTOUT) DFRLG
 	 WRITE(LISTOUT) QS0     
 	 WRITE(LISTOUT) SQS
 	 WRITE(LISTOUT) THE0    
 	 WRITE(LISTOUT) STHE
 	 WRITE(LISTOUT) THE0Q   
 	 WRITE(LISTOUT) STHEQ
 	 WRITE(LISTOUT) MXSNAL  
 	 WRITE(LISTOUT) EPSR	    				 
 	 WRITE(LISTOUT) RADIN
 	 WRITE(LISTOUT) RADOT	    				 
 	 WRITE(LISTOUT) GLAT
 	 WRITE(LISTOUT) GLON	    				 
 	 WRITE(LISTOUT) CZEN	    				 
 	 WRITE(LISTOUT) HTOP
 	 WRITE(LISTOUT) HBOT	    				 
 	 WRITE(LISTOUT) CNVTOP 
 	 WRITE(LISTOUT) CNVBOT      				 
 	 WRITE(LISTOUT) TG
 	 WRITE(LISTOUT) GFFC	    				 
 	 WRITE(LISTOUT) SST
 	 WRITE(LISTOUT) ALBASE      				 
 	 WRITE(LISTOUT) ALBEDO      				 
 	 WRITE(LISTOUT) HDAC
 	 WRITE(LISTOUT) HDACV	    				 
 	 WRITE(LISTOUT) CZMEAN 
 	 WRITE(LISTOUT) SIGT4
 	 WRITE(LISTOUT) PTBL
 	 WRITE(LISTOUT) TTBL
 	 WRITE(LISTOUT) TTBLQ
!MODULE_PPTASM
 	 WRITE(LISTOUT) MTSTPE  
 	 WRITE(LISTOUT) ITSTLOC 
 	 WRITE(LISTOUT) JTSTLOC 
 	 WRITE(LISTOUT) PHOUR   								     
 	 WRITE(LISTOUT) APREC
 	 WRITE(LISTOUT) TLATCU  								     
 	 WRITE(LISTOUT) TLATGS
 	 WRITE(LISTOUT) PPTDAT
!MODULE_PRFHLD
 	 WRITE(LISTOUT) TLMIN  
 	 WRITE(LISTOUT) TLMAX 
!MODULE_PVRBLS  
 	 WRITE(LISTOUT) Z0
 	 WRITE(LISTOUT) USTAR
 	 WRITE(LISTOUT) UZ0 
 	 WRITE(LISTOUT) VZ0 
 	 WRITE(LISTOUT) THZ0   
 	 WRITE(LISTOUT) QZ0   
 	 WRITE(LISTOUT) THS
 	 WRITE(LISTOUT) QS
 	 WRITE(LISTOUT) AKMS 
 	 WRITE(LISTOUT) AKHS 
 	 WRITE(LISTOUT) RF
 	 WRITE(LISTOUT) TWBS 
 	 WRITE(LISTOUT) QWBS 
 	 WRITE(LISTOUT) SNO 
 	 WRITE(LISTOUT) SI
 	 WRITE(LISTOUT) CLDEFI 
 	 WRITE(LISTOUT) PREC
 	 WRITE(LISTOUT) ACPREC    
 	 WRITE(LISTOUT) ACCLIQ
 	 WRITE(LISTOUT) CUPREC
 	 WRITE(LISTOUT) TH100
 	 WRITE(LISTOUT) Q100
 	 WRITE(LISTOUT) U100
 	 WRITE(LISTOUT) V100
 	 WRITE(LISTOUT) TH10    
 	 WRITE(LISTOUT) Q10     
 	 WRITE(LISTOUT) U10
 	 WRITE(LISTOUT) V10
 	 WRITE(LISTOUT) TSHLTR
 	 WRITE(LISTOUT) TSHLTR1
 	 WRITE(LISTOUT) QSHLTR
 	 WRITE(LISTOUT) PSHLTR    
 	 WRITE(LISTOUT) XMOMFLUX
 	 WRITE(LISTOUT) YMOMFLUX
	 WRITE(LISTOUT) AKM10 
	 WRITE(LISTOUT) AKM10V 
!	 WRITE(LISTOUT) PLM     
         WRITE(LISTOUT) Q2   
!MODULE_RACCR_TABLES    
         WRITE(LISTOUT) ACCRR
!MODULE_RD1TIM
         WRITE(LISTOUT) K400
         WRITE(LISTOUT) RAD1
         WRITE(LISTOUT) CTHK
         WRITE(LISTOUT) TAUCV
         WRITE(LISTOUT) PTOPC
         WRITE(LISTOUT) LTOP
         WRITE(LISTOUT) LVL
!MODULE_RDFSAV
         WRITE(LISTOUT) EMISP   
         WRITE(LISTOUT) EMIST   
         WRITE(LISTOUT) XLATT   
         WRITE(LISTOUT) XLATP   
         WRITE(LISTOUT) Q19001  
         WRITE(LISTOUT) HP98    
         WRITE(LISTOUT) H3M6    
         WRITE(LISTOUT) HP75    
         WRITE(LISTOUT) H6M2    
         WRITE(LISTOUT) HP537   
         WRITE(LISTOUT) H74E1   
         WRITE(LISTOUT) H15E1   
         WRITE(LISTOUT) Q14330  
         WRITE(LISTOUT) HP2     
         WRITE(LISTOUT) TWENTY  
         WRITE(LISTOUT) HNINE   
         WRITE(LISTOUT) DEGRAD  
         WRITE(LISTOUT) HSIGMA 
         WRITE(LISTOUT) DAYSEC  
         WRITE(LISTOUT) RCO2
         WRITE(LISTOUT) CAO3SW  
         WRITE(LISTOUT) CAH2SW  
         WRITE(LISTOUT) CBSW
!MODULE_RITE
!         WRITE(LISTOUT) BETA    
!         WRITE(LISTOUT) DRIP    
!         WRITE(LISTOUT) EC      
!         WRITE(LISTOUT) EDIR    
!         WRITE(LISTOUT) ETT     
!         WRITE(LISTOUT) FLX1    
!         WRITE(LISTOUT) FLX2    
!         WRITE(LISTOUT) FLX3    
!         WRITE(LISTOUT) RUNOF   
!         WRITE(LISTOUT) DEW     
!         WRITE(LISTOUT) RIB     
!         WRITE(LISTOUT) RUNOFF3
!MODULE_RMASS_TABLES
         WRITE(LISTOUT) MASSR
!MODULE_RRATE_TABLES
         WRITE(LISTOUT) RRATE
!MODULE_RVELR_TABLES
         WRITE(LISTOUT) VRAIN
!MODULE_RVENT_TABLES
         WRITE(LISTOUT) VENTR1
         WRITE(LISTOUT) VENTR2
!MODULE_SAVEMEM
!         WRITE(LISTOUT) DDUO3N
!         WRITE(LISTOUT) DDO3N2  
!         WRITE(LISTOUT) DDO3N3  
!         WRITE(LISTOUT) DDO3N4
!         WRITE(LISTOUT) RAD1   ! aparece tambem no modulo RD1TIM 
!         WRITE(LISTOUT) RAD2    
!         WRITE(LISTOUT) RAD3    
!         WRITE(LISTOUT) RAD4
!MODULE_SCRTCH
         WRITE(LISTOUT) SUM    
	 WRITE(LISTOUT) PERTSM 
	 WRITE(LISTOUT) SUM3   
	 WRITE(LISTOUT) SUMWDE
         WRITE(LISTOUT) SRCWD
         WRITE(LISTOUT) SRC1NB 
         WRITE(LISTOUT) DBDTNB
         WRITE(LISTOUT) ZMASS  
	 WRITE(LISTOUT) ZROOT
         WRITE(LISTOUT) SC     
         WRITE(LISTOUT) DSC    
         WRITE(LISTOUT) XTEMV
         WRITE(LISTOUT) TFOUR
         WRITE(LISTOUT) FORTCU
         WRITE(LISTOUT) X     
         WRITE(LISTOUT) X1
         WRITE(LISTOUT) SRCS 
         WRITE(LISTOUT) SUM4 
         WRITE(LISTOUT) SUM6
         WRITE(LISTOUT) SUM7	
	 WRITE(LISTOUT) SUM8	
	 WRITE(LISTOUT) SUM4WD  
	 WRITE(LISTOUT) R1  	 
         WRITE(LISTOUT) R2	
         WRITE(LISTOUT) S2	
         WRITE(LISTOUT) T3   
         WRITE(LISTOUT) R1WD
	 WRITE(LISTOUT) X2	
         WRITE(LISTOUT) EXPO	
         WRITE(LISTOUT) FAC
         WRITE(LISTOUT) CNUSB	
         WRITE(LISTOUT) DNUSB
         WRITE(LISTOUT) ALFANB  
         WRITE(LISTOUT) AROTNB
         WRITE(LISTOUT) ANB	
         WRITE(LISTOUT) BNB	
         WRITE(LISTOUT) CENTNB  
         WRITE(LISTOUT) DELNB	
         WRITE(LISTOUT) BETANB
!MODULE_SDENS_TABLES
         WRITE(LISTOUT) SDENS
!MODULE_SEASO3
	 WRITE(LISTOUT) XDUO3N
	 WRITE(LISTOUT) XDO3N2
	 WRITE(LISTOUT) XDO3N3
	 WRITE(LISTOUT) XDO3N4
	 WRITE(LISTOUT) PRGFDL
	 WRITE(LISTOUT) O3O3
	 WRITE(LISTOUT) XRAD1
	 WRITE(LISTOUT) XRAD2
	 WRITE(LISTOUT) XRAD3
	 WRITE(LISTOUT) XRAD4
!MODULE_SLOPES
         WRITE(LISTOUT) ISLD
         WRITE(LISTOUT) VTMS 
!MODULE_SOIL	  
         WRITE(LISTOUT) IVGTYP 
         WRITE(LISTOUT) ISLTYP 
	 WRITE(LISTOUT) ISLOPE
         WRITE(LISTOUT) SOILTB
         WRITE(LISTOUT) SFCEXC
	 WRITE(LISTOUT) SMSTAV
	 WRITE(LISTOUT) SMSTOT
         WRITE(LISTOUT) GRNFLX 
         WRITE(LISTOUT) PCTSNO
         WRITE(LISTOUT) VEGFRC 
         WRITE(LISTOUT) CMC
	 WRITE(LISTOUT) SMC
	 WRITE(LISTOUT) STC
	 WRITE(LISTOUT) SH2O
	 WRITE(LISTOUT) SLDPTH
	 WRITE(LISTOUT) RTDPTH	 
!MODULE_SSALB	 
	 WRITE(LISTOUT) ZA
	 WRITE(LISTOUT) TRN
	 WRITE(LISTOUT) DZA
	 WRITE(LISTOUT) ALBD
	 WRITE(LISTOUT) ALB1
	 WRITE(LISTOUT) ALB2
	 WRITE(LISTOUT) ALB3	 
!MODULE_SWRSAV 
	 WRITE(LISTOUT) ABCFF
	 WRITE(LISTOUT) PWTS
	 WRITE(LISTOUT) CFCO2  
	 WRITE(LISTOUT) CFO3   
	 WRITE(LISTOUT) REFLO3 
	 WRITE(LISTOUT) RRAYAV
!MODULE_TABCOM
         !colocar depois
!MODULE_TEMPV
         WRITE(LISTOUT) P0
         WRITE(LISTOUT) T0
	 WRITE(LISTOUT) Q0
!MODULE_TOPO
         WRITE(LISTOUT) TEMP2X
         WRITE(LISTOUT) TTVG
	 WRITE(LISTOUT) HTMG
!MODULE_VRBLS
         WRITE(LISTOUT) PD
         WRITE(LISTOUT) FIS
         WRITE(LISTOUT) RES
         WRITE(LISTOUT) T
         WRITE(LISTOUT) U
         WRITE(LISTOUT) V
         WRITE(LISTOUT) Q
!MODULE_UPDT
         !colocar depois
!MODULE_Z0EFFT
	 WRITE(LISTOUT) ZEFFIJ
!
	 CLOSE(LISTOUT)
	 
	 IF (MYPE == 0) THEN
            WRITE(*,*) ' OU2RESTRT: END OF WRITING'
            WRITE(*,*) ' '
         END IF
   
    END SUBROUTINE OUT2RESTRT
