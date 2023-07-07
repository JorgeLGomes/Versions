      subroutine eta_const
c
c *** Original code received from U. of Athens and modified at FSL.
c
c *** Two-dimensional horizontal indexing into one-dimensional,     
c     defines dummy initial and boundary moisture values        
c     and calculates constants needed for the one-dimensional  
c     version of the ub/nmc model.                            

C	revisions adapted from Z. Janjic version of this code.  Now 
C	stays in im,jm indexing, and seems to produce better results

c
      implicit none
c
      include 'ecommons.h'
c
      real*4 pt,aeta(lm),eta(lmp1),dfl(lmp1),detac(lm)

	REAL,ALLOCATABLE:: RES(:,:),SNO(:,:)
	REAL,ALLOCATABLE:: SST(:,:),SI(:,:),CMC(:,:),ALBEDO(:,:)

	REAL,ALLOCATABLE:: HBM2(:,:),VBM2(:,:),VBM3(:,:)
	REAL,ALLOCATABLE:: SM(:,:),SICE(:,:),HTM(:,:,:),VTM(:,:,:)
c
      real*4 pdb(kb,2)
     .      ,tb(kb,lm,2),qb(kb,lm,2)
     .      ,ub(kb,lm,2),vb(kb,lm,2)
     .      ,pd(im,jm)
     .      ,pdt(im,jm),fist(im,jm)  

	REAL,ALLOCATABLE:: TT(:,:,:),UT(:,:,:),VT(:,:,:),QT(:,:,:)
	REAL,ALLOCATABLE:: SFCGRID(:,:,:),SMC(:,:,:),STC(:,:,:)
	REAL,ALLOCATABLE:: SH2O(:,:,:)

       real dum2(im,jm)

c
C	real smc,stc,tg
     	real	tg(im,jm)
      integer*4 idat(3),ihrst,ntsd,i,j,k,l,len,II,JJ
	common /mytime/idat
	real*4 pdmin
c
      logical run
c_______________________________________________________________________________
c
      print*,'im,jm=',im,jm

c
      len=index(init_out//' ',' ')-1
      if (init_out(len:len) .ne. '/') then
         len=len+1
         init_out(len:len)='/'
      endif
      open(1,file=init_out(1:len)//'preproc.init'
     .    ,status='old',form='unformatted')

	ALLOCATE(UT(IM,JM,LM),VT(IM,JM,LM),TT(IM,JM,LM),QT(IM,JM,LM))
	ALLOCATE(SM(IM,JM),RES(IM,JM))
	write(6,*) 'to preproc.init read'
C
C	       L   I3    I    
      read(1) run,idat,ihrst,ntsd,ut,vt
	write(6,*) 'middle of preproc.init read'
      read(1) tt,qt,pdt,fist,sm,res,eta,pt,detac,aeta,dfl            
	write(6,*) 'past preproc.init read'

	write(6,*) 'some res values: ', (res(i,jm/2),i=1,im)
	

C        do J=jm,1,-1
C        write(6,666) (pdt(i,j)/100.,i=im-18,im)
C        enddo
  666   format(20(f5.0,1x))

	ALLOCATE(SFCGRID(IM,JM,12))
	read(1) sfcgrid
	
	ALLOCATE(STC(IM,JM,NSOIL),SMC(IM,JM,NSOIL))

	
	if (GRIBSOIL .and. .NOT. REAN) then

	write(6,*) 'pulling soil from sfcgrid'

	do K=1,4
          do j=1,jm
            do i=1,im
              stc(i,j,K)=sfcgrid(i,j,2*K+3)
              smc(i,j,K)=sfcgrid(i,j,2*K+4)
            enddo
          enddo
	enddo

C Isabel Jul2016 Include soil layers
        do K=5,8
           do j=1,jm
             do i=1,im
               stc(i,j,K)=sfcgrid(i,j,11)
               smc(i,j,K)=sfcgrid(i,j,12)
             enddo
           enddo
        enddo    
C Isabel


	endif

      close(1)
c
	ALLOCATE(HBM2(IM,JM),VBM2(IM,JM),VBM3(IM,JM),HTM(IM,JM,LM))
	ALLOCATE(VTM(IM,JM,LM),SICE(IM,JM),SNO(IM,JM))
	ALLOCATE(SST(IM,JM),SI(IM,JM),ALBEDO(IM,JM),SH2O(IM,JM,NSOIL))
      call cnsts(pt,aeta,eta,dfl,detac
     .          ,res,fist,sno,sst,si,albedo
     .          ,hbm2,vbm2,vbm3,sm,sice,htm,vtm,smc,stc,tg,sh2o)
c
	do j=1,jm
	do i=1,im
         res(i,j)=1./res(i,j)                                                  
        enddo
	enddo
C
      print *,' '
      print *,'Creating eta initial condition file...'
      print *,'Write to file : ',init_out(1:len)//'init.file'
      open(1,file=init_out(1:len)//'init.file'
     .    ,status='unknown',form='unformatted')
      write(1) run,idat,ihrst,ntsd

Cmp----------------------------
	write(6,*) 'in CONST: pd values'
	pdmin=9999999.
	do JJ=jm,1,-1
	do II=1,IM
	if (pdt(II,JJ) .lt. pdmin) pdmin=pd(II,JJ)
	enddo
	enddo
	write(6,*) 'pdmin in CONST is: ', pdmin
      write(1) pdt
      write(6,*) 'pdt: ',pdt(1,1),pdt(21,21),pdt(im,jm)
Cmp---------------------------
      write(1) res
C-----------------------------
	write(6,*) 'fist: ', fist(1,1),fist(21,21),fist(im,jm)
      write(1) fist
C-----------------------------
      do l=1,lm
	do j=1,jm
	do i=1,im
	dum2(i,j)=ut(i,j,l)
	enddo
	enddo
         write(1) dum2
        print *,'u: ',L,dum2(1,1),dum2(21,21),dum2(im,jm)
      enddo

	DEALLOCATE(UT)
C-----------------------------
      do l=1,lm
	do j=1,jm
        do i=1,im
        dum2(i,j)=vt(i,j,l)
        enddo
        enddo
         write(1) dum2
        if (l.eq.20) print *,'v:',dum2(1,1),dum2(21,21),dum2(im,jm)
      enddo
	DEALLOCATE(VT)
C----------------------------
      do l=1,lm
	do J=1,jm
	do I=1,im
	dum2(i,j)=tt(i,j,l)
	enddo
	enddo
         write(1) dum2
        if (l.eq.20) print *,'t:',dum2(1,1),dum2(21,21),dum2(im,jm)
      enddo
	DEALLOCATE(TT)
C------------------------------
      do l=1,lm
        do J=1,jm
        do I=1,im
        dum2(i,j)=qt(i,j,l)
        enddo
        enddo
         write(1) dum2
        if (l.eq.20) print *,'q:',dum2(1,1),dum2(21,21),dum2(im,jm)
      enddo
	DEALLOCATE(QT)
C------------------------------
C	wet --> si
	write(6,*) 'SI VALUES'
	do JJ=JM,1,-JM/30
	write(6,969) (si(II,JJ),II=1,IM,IM/12)
  969	format(25(f4.2,1x))
	enddo
      write(1) si
	write(6,*) 'SNO VALUES'
	do JJ=JM,1,-JM/30
	write(6,969) (sno(II,JJ),II=1,IM,IM/12)
	enddo

	do J=1,jm
	do I=1,im
	if (smc(i,j,1) .eq. 0. .and. sm(i,j) .eq. 0) then
!	write(6,*) 'LAND POINT WITH ZERO SMC VAL: ', I,J,SICE(I,J)
	do II=1,nsoil
	SMC(I,J,II)=0.14
	STC(I,J,II)=271.
	enddo
	endif
	enddo
	enddo
      write(1) sno
      write(1) smc

	DEALLOCATE(SNO,SMC)
C-----------------------------------
C  SET VEGETATION CANOPY WATER CONTENT EQUAL TO ZERO EVERYWHERE FOR NOW
	ALLOCATE(CMC(IM,JM))
	cmc=0.
      write(1) cmc
	DEALLOCATE(CMC)
C-----------------------------------
	do JJ=1,JM
	do II=1,IM
	if (STC(II,JJ,1) .eq. 0 .and. sm(II,JJ) .eq. 0) then
	write(6,*) 'STC zero at land point!!! ', II,JJ
	endif
	enddo
	enddo
      write(1) stc
	write(6,*) ' soil tmps (layer 1) '
	do JJ=JM,1,-JM/30
	write(6,922) (stc(II,JJ,1),II=1,IM,IM/12)
  922	format(25(f4.0,1x))
	enddo

	DEALLOCATE(STC)

	write(6,*) 'sh2o values: '
        do JJ=JM,1,-JM/30
        write(6,923) (sh2o(II,JJ,1),II=1,IM,IM/12)
	enddo

	write(1) sh2o
	deallocate (sh2o)
	write(6,*) 'albedo values: '
        do JJ=JM,1,-JM/30
        write(6,923) (albedo(II,JJ),II=1,IM,IM/12)
	enddo
  923	format(25(f4.2,1x))

	write(1) albedo
	deallocate (albedo)
C-------------------------------------
      close(1)
      end                                        
c
c===============================================================================
c
      subroutine cnsts(pt,aeta,eta,dfl,detac
     .                ,res,fis,sno,sst,si,albedo
     .                ,hbm2,vbm2,vbm3,sm,sice,htm,vtm
     .		      ,smc,stc,tg,sh2o)
c
c *** Routine for initialization of constants and variables         
c
      implicit none
c
      include 'ecommons.h'
      include 'econstants.h'
c
      integer*4 itb,jtb,itbq,jtbq,JJ,NSOTYP,NFLT

      parameter (itb=76,jtb=134
     .          ,itbq=152,jtbq=440)

      real frh2o,rsnow,snofac

c
c     common/pteta/                       
      real*4 pt,aeta(lm),eta(lmp1),dfl(lmp1),detac(lm)
     .      ,res(im,jm),fis(im,jm),sno(im,jm)              
     .      ,sst(im,jm),si(im,jm),cmc(im,jm)
c
c     common/masks/                        
      real*4 hbm2(im,jm),vbm2(im,jm),vbm3(im,jm)
     .      ,sm(im,jm),sice(im,jm)  
     .      ,htm(im,jm,lm),vtm(im,jm,lm)

C
C       COMMON /OLDALB/ ALBASE(IM,JM)
       REAL ALBASE(IM,JM)
C

c
      real*4 dxj(jm),wpdarj(jm),cpgfuj(jm),curvj(jm),fcpj(jm)      
     .      ,fdivj(jm),emj(jm),emtj(jm),fadj(jm)                  
     .      ,ddmpuj(jm),ddmpvj(jm),hdacj(jm)                              
     .      ,qsold(jtb),pold(jtb),qsnew(jtb),pnew(jtb)                  
     .      ,y2p(jtb),app(jtb),aqp(jtb)                              
     .      ,theold(jtb),told(jtb),thenew(jtb),tnew(jtb)                  
     .      ,y2t(jtb),apt(jtb),aqt(jtb)                              

Cmp
C	real smc,stc,tg
	real smc(im,jm,nsoil),stc(im,jm,nsoil),
     +			tg(im,jm),sh2o(im,jm,nsoil),bx,fk
Cmp
c
      integer*4 khla(jam),khha(jam),kvla(jam),kvha(jam)                  
     .         ,khl2(jm),khh2(jm),kvl2(jm),kvh2(jm),kvl                  
Cmp     .	       ,ihw(jm),ihe(jm),ivw(jm),ive(jm),ihl,ihh
     .	       ,ihl,ihh,ihe,ihw
c
      integer lmh(im,jm),lmv(im,jm)                                        
c
      real*4 rdeta(lm),f4q2(lm),diff
     .      ,em(jam),emt(jam)                                          
     .      ,dx(im,jm),wpdar(im,jm),cpgfu(im,jm),curv(im,jm),fcp(im,jm) 
     .      ,fdiv(im,jm),fad(im,jm),f(im,jm),ddmpu(im,jm),ddmpv(im,jm) 
     .      ,glat(im,jm),glon(im,jm),glatr(im,jm),glonr(im,jm)        
     .      ,deta2(lm),aeta2(lm),eta2(lmp1),dfrlg(lmp1)                         
     .      ,qs0(jtb),sqs(jtb),the0(itb),sthe(itb)                  
     .      ,epsr(im,jm),gffc(im,jm)         
     .      ,hdac(im,jm),hdacv(im,jm)                                        
     .      ,ptbl(itb,jtb),ttbl(jtb,itb)                                  
cds new convection
     .      ,ttblq(jtbq,itbq),the0q(itbq),stheq(itbq)
c
      integer*4 idum2(im,jm),idat(3),ierr
c
      real*4 deta1(lm),aeta1(lm),eta1(lmp1)                                
     .      ,dum2(im,jm),dum2b(im,jm),dum2a(im,jm)
     .      ,dum2sm(im,jm)
c
      integer*4 kpm,kpm1,kthm,kthm1,kvh,khl,khh
     .         ,nddamp
     .         ,list,len,i,j,k,l,ja,nfcst,nbc,nx,ny
c
      real*4 vegfrc(im,jm),sldpth(nsoil),rtdpth(nsoil)
     .      ,f4d,f4q,ctlm,ef4t,stlm,aph,stph,ctph,tlm
     .      ,sinphi,coslam,fp,fact,thl,en,ent,dy,cpgfv
     .      ,dtad,tsph,acdt,dxp,wbi
     .      ,rdph,rdlm,sbi,anbi,ebi
     .      ,rdq,rdth,rdthe,rdtheq,rdp,rdpq
     .      ,pl,pt1,pt2,r1,ph,tph0,tph,tdph,tdlm,cddamp
     .      ,thh,ti0,wb,dph,dlm,sb

C Soil parameter arrays used in cold start SH2O initialization
       real,allocatable:: PSIS(:),BETA(:),SMCMAX(:)
       
       REAL  PSISG(15),BETAG(15),SMCMAXG(15)
       REAL  PSISD(18),BETAD(18),SMCMAXD(18)
       REAl  PSISZ(9),BETAZ(9),SMCMAXZ(9)
C--------------------------------------------------------------
c
	INTEGER*4 I1D(360,180),IONETA(im,jm), JULM(13),JULD,
     .  ILON1,ILON2,ILAT1,ILAT2,ILONX,ILATX,GDS(200)
     .  
       REAL*4 ALBC1(361,180),ALBC2 (361,180),ALBC3(361,180),
     .	  ALBC4(361,180),ALBC  (361,180),ALBEDO(IM,JM),
     .		SNALMX(361,180),MXSNAL(IM,JM),
     .   S1,S2,W1,W2,AR1,AR2,AR3,AR4,H1,D5,D01,HM1,ELON,ELAT,
     .    H90,H360,WGT1,WGT2,D00,DIF,alb1d(IMJM),radfac

        real,allocatable:: GGSTC1(:,:),GGSTC2(:,:)
        real,allocatable:: GGSMC1(:,:),GGSMC2(:,:)
        real,allocatable:: SSTRAW(:,:),SNOWRAW(:,:)
        real,allocatable:: GGLAND(:,:),GGICE(:,:)
        
        parameter(nx=192,ny=94)


	parameter (H1=1.0,D5=5.E-1,D01=1.00E-2,HM1=-1.0,
     .	H90=90.0,H360=360.0,D00=0.0,radfac=57.295777951)
    
Cmp
	common /mytime/idat
	DATA JULM/0,31,59,90,120,151,181,212,243,273,304,334,365/
c
c *** Universal constants.
c
      real*4 a,twom,SNUP,SALP,ARX
      data a/6376000./,twom/.00014584/
      DATA SNUP,SALP /0.040, 2.6/

c
c *** Dissipation & turbulence.
c
      real*4 coac,codamp,tddamp,dfc,ddfc
	data codamp/0150.00/,w/0.25/
CGSM     .    ,tddamp/48.00/,dfc/01./,ddfc/8.0/
     .    ,tddamp/72.00/,dfc/01./,ddfc/1.0/
     
c
c
c *** Surface data.
c
Cmp	these should be same as operational now.

      real*4 ros,cs,ds,aks,dzg,tg0,tga,roi,ci,di,aki,dzi,elwv,plq
      data ros/1500./,cs/1339.2/,ds/.100/,aks/.0000005  /,dzg/02.50/       
     .,tg0/258.16/,tga/30./                                             
     .,roi/916.6/,ci/2060.0/,di/.100/,aki/.000001075/,dzi/2.0/          
     .,ti0/271.16/                                                      
     .,elwv/2.50e6/                                                     
     .,plq/70000./

C-----------------------------------------------------------------------
C Constants/data used in cold start SH2O initialization
	REAL HLICE, GRAV, BLIM
C      DATA HLICE/3.335E5/,GRAV/9.81/,T0/273.15/
      DATA HLICE/3.335E5/,GRAV/9.81/
      DATA BLIM/5.5/
      
      logical vegflag
      logical noz,lusaf
      
C Parameter used in cold start SH2O initialization
C      write(6,*) 'lendo a QUANTIDADE DE TIPOS DE SOLO '
      open(1,file='TYPSOLO',form='formatted',status='old')
      read(1,10)NSOTYP
10    format(I2)    
C      write(6,*) 'leu a QUANTIDADE DE TIPOS DE SOLO '
      close(1)
C      NSOTYP=15
      print*,"QUANTIDADE DE TIPOS DE SOLO",NSOTYP
      
      allocate (PSIS(NSOTYP),BETA(NSOTYP),SMCMAX(NSOTYP)) 

CPARAMETROS ZOGBLER
      DATA PSISZ /0.04,0.62,0.47,0.14,0.10,0.26,0.14,0.36,0.04/
      DATA BETAZ /4.26,8.72,11.55,4.74,10.73,8.17,6.77,5.25,4.26/
      DATA SMCMAXZ /0.421,0.464,0.468,0.434,0.406,
     &             0.465,0.404,0.439,0.421/
CPARAMETROS ZOGBLER

CDANIEL PARAMETROS CLASIFICACION DE MOIRA
       DATA PSISD /0.265,0.103,0.252,0.106,0.088,0.052,
     &           0.109,0.042,0.018,0.424,0.400,0.220, 
     &           0.601,0.500,0.132,0.929,1.294,0.8500/
       DATA BETAD /2.02,2.78,3.03,2.25,2.62,3.02,
     &           1.97,2.53,2.98,2.13,2.44,2.90,
     &           2.01,2.42,2.97,2.81,1.83,2.29/
       DATA SMCMAXD /0.240,0.354,0.514,0.402,0.477,0.564,
     &             0.440,0.566,0.631,0.443,0.517,0.566,
     &             0.478,0.562,0.675,0.505,0.534,0.580/
CDANIEL

CAROL PARAMETROS NEW GLOBAL SOIL 
      DATA PSISG /0.069,0.036,0.141,0.759,0.759,0.355,
     &           0.135,0.617,0.263,0.098,0.324,0.468, 
     &           0.355,0.069,0.036/
      DATA BETAG /2.79,4.26,4.74,5.33,5.33,5.25,
     &           6.66,8.72,8.17,10.73,10.39,11.55,
     &           5.25,2.79,4.26/
      DATA SMCMAXG /0.339, 0.421, 0.434, 0.476, 0.476, 0.439,
     &            0.404, 0.464, 0.465, 0.406, 0.468, 0.468,
     &            0.439, 0.200, 0.421/
CAROL      
C-----------------------------------------------------------------------

c_______________________________________________________________________________
c
      
      IF (NSOTYP.eq.15) THEN
	   PSIS=PSISG
	   BETA=BETAG
       	   SMCMAX=SMCMAXG
      ENDIF
      
      IF (NSOTYP.eq.9) THEN	
       	   PSIS=PSISZ
	   BETA=BETAZ
	   SMCMAX=SMCMAXZ
      ENDIF
      
      IF (NSOTYP.eq.18) THEN	
	    PSIS=PSISD
	    BETA=BETAD
	    SMCMAX=SMCMAXD
      ENDIF
      
      namelist/VEG/vegflag

CJLG     
	write(6,*) 'lendo VEGETATION '
       open(1,file='VEGETATION',form='formatted',status='old')
       read(1,VEG)
       close(1)
	write(6,*) 'leu VEGETATION '
CJLG      

c------------------------------------------------------------------------------------------------------
Cmp
	tg=-99.
Cmp
      noz=.true.
      list=3
	nfcst=13
	nbc=16

	write(6,*) 'executing cnsts....GRIBSOIL is: ', GRIBSOIL

	write(6,*) 'const believes there are ', nsoil, ' soil layers'

Cmp
Cmp	compute coac based on resolution!
Cmp	
Cmp	coac=.04*(.53846154/DPHD)**2.


C	removed squaring per t. black recommendation

	   
CGSM 1km	tddamp=nhour

	coac=.04*(.53846154/DPHD)

C	cap at 0.2
CGSM	coac=amin1(coac,0.2)
CGSM        coac=amin1(coac,0.35)

CGSM 1km
      IF (DPHD.LT.0.108) THEN

!GSM        coac=(0.04+0.08*(0.108-DPHD)/0.099)*(.53846154/DPHD)
!GSM          coac=0.8
!GSM          coac=1
!GSM          coac=0.35
!GSM          coac=0.583
!GSM          coac=1.01 
!GSM          coac=0.6
          coac=0.6
      ENDIF

	write(6,*) 'computed coac based on grid resolution...'
	write(6,*) 'will use coac, dfc, ddfc, tddamp= ', coac, dfc, ddfc,tddamp
Cmp

c *** Derived vertical constants.
c
      do l=1,lmp1                             
         dfrlg(l)=dfl(l)                            
         dfl(l)=dfl(l)*g                            
      enddo
c
c *** Dummy constants.
c
      do l=1,lmp1                             
         eta1(l)=eta(l)                                                    
         eta2(l)=eta(l)                                                    
      enddo
      do l=1,lm                              
         deta1(l)=detac(l)                                                  
         aeta1(l)=aeta(l)                                                  
         deta2(l)=detac(l)                                                  
         aeta2(l)=aeta(l)                                                  
      enddo
c
c *** Derived geometrical constants.
c
      tph0=tph0d*dtr                                                    
      wb=wbd*dtr                                                        
      sb=sbd*dtr                                                        
      dlm=dlmd*dtr                                                      
      dph=dphd*dtr                                                      
      tdlm=dlm+dlm                                                      
      tdph=dph+dph                                                      
      rdlm=1./dlm                                                       
      rdph=1./dph                                                       
c                                                                       
      wbi=wb+tdlm                                                       
      sbi=sb+tdph                                                       
      ebi=wb+im2*tdlm                                                   
      anbi=sb+jm3*dph                                                   
c                                                                       
c     stph0=sin(tph0)                                                   
c     ctph0=cos(tph0)                                                   
c
c *** Time stepping related constants.
c
      tsph=3600./dt                                                     
      nddamp=tddamp*tsph+.5                                             
c                                                                       
      dtad=idtad                                                        
c
c *** Derived horizontal grid constants.
c
      dy=a*dph                                                          
      cpgfv=-dt/(48.*dy)                                                
      en= dt/( 4.*dy)*dtad                                              
      ent=dt/(16.*dy)*dtad                                              
c                                                                       
      tph=sb-dph                                                        
      do j=1,jm                                            
C      ihw(j)=-mod(j,2)
C      ihe(j)=ihw(j)+1
C      ivw(j)=-mod(j+1,2)
C      ive(j)=ivw(j)+1
c                                                                       
         tph=tph+dph                                               
         dxp=a*dlm*cos(tph)                                        
         dxj(j)=dxp                                                
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
csd      wpdarj(j)=-dt*w/(32.*dxp*dy)                      
csd      wpdarj(j)=-dt*w*100000./(32.*dxp*dy)                      
         wpdarj(j)=-w*((a*dlm*amin1(cos(anbi),cos(sbi)))**2+dy**2)
     .            /(dt*32.*dxp*dy)*.88
caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
         cpgfuj(j)=-dt/(48.*dxp)                                   
         curvj(j)=.5*dt*tan(tph)/a                                 
         fcpj(j)=dt/(cp*192.*dxp*dy)                               
         fdivj(j)=1./(12.*dxp*dy)                                  
         emj(j)= dt/( 4.*dxp)*dtad                                 
         emtj(j)=dt/(16.*dxp)*dtad                                 
         fadj(j)=-dt/(48.*dxp*dy)*dtad                             
Cmp         acdt=coac*dt                                              
Cmp     .       *sqrt((a*dlm*amin1(cos(anbi),cos(sbi)))**2+dy**2)     
        ACDT=DT
     2            *SQRT((A*DLM*AMIN1(COS(ANBI),COS(SBI)))**2+DY**2)
Cmp         cddamp=codamp*acdt                                        
         cddamp=.04*codamp*acdt                                        
Cmp         hdacj(j)=acdt/(4.*dxp*dy)                                 
         hdacj(j)=coac*acdt/(4.*dxp*dy)                                 
         ddmpuj(j)=cddamp/dxp                                      
         ddmpvj(j)=cddamp/dy                                       
      enddo
c
c *** Spreading of upstream height-point advection factor.
c
      ja=0                                                              
      do j=3,5                                              
         ja=ja+1                                                       
      khla(ja)=2
      khha(ja)=im-1-mod(j+1,2)
         emt(ja)=emtj(j)                                               
      enddo

      do j=jm4,jm2                                          
         ja=ja+1                                                       
         khla(ja)=2                                 
         khha(ja)=im-1-mod(j+1,2)                                            
         emt(ja)=emtj(j)                                               
      enddo

      do j=6,jm5                                            
         ja=ja+1                                                       
	khla(ja)=2
        khha(ja)=2+mod(j,2)
        emt(ja)=emtj(j)                                               
      enddo

      do j=6,jm5                                            
         ja=ja+1                                                       
      khla(ja)=im-2
      khha(ja)=im-1-mod(j+1,2)
      emt(ja)=emtj(j)
      enddo
      print*,'ja=',ja,' jam=',jam
c
c *** Spreading of upstream velocity-point advection factor.
c
      ja=0                                                              
      do j=3,5                                              
         ja=ja+1                                                       
	 kvla(ja)=2
         kvha(ja)=im-1-mod(j,2)
         em(ja)=emj(j)                                                 
      enddo

      do j=jm4,jm2                                          
         ja=ja+1                                                       
         kvla(ja)=2     
         kvha(ja)=im-1-mod(j,2)                                           
         em(ja)=emj(j)                                                 
      enddo

      do j=6,jm5                                            
         ja=ja+1                                                       
         kvla(ja)=2                                          
         kvha(ja)=2+mod(j+1,2)                                              
         em(ja)=emj(j)                                                 
      enddo

      do j=6,jm5                                            
         ja=ja+1                                                       
         kvla(ja)=im-2
         kvha(ja)=im-1-mod(j,2)
         em(ja)=emj(j)                                                 
      enddo
c
c *** Coriolis parameter in tll system & related constants.
c
      tph=sb-dph                                                
      do j=1,jm                                            
c                                                                       
         tlm=wb-tdlm+mod(j,2)*dlm                                  
         tph=tph+dph                                               
         stph=sin(tph)                                             
         ctph=cos(tph)                                             
c                                                                       
         do i=1,im                                             
            tlm=tlm+tdlm                                                      
            fp=twom*(ctph0*stph+stph0*ctph*cos(tlm))                          
            f(i,j)=0.5*dt*fp                                                    
         enddo
      enddo
c
c *** Geographic lat and long of tll grid points.
c
      tph=sb-dph                                                
      do j=1,jm                                            
         tlm=wb-tdlm+mod(j+1,2)*dlm                                
         tph=tph+dph                                               
         stph=sin(tph)                                             
         ctph=cos(tph)                                             
c                                                                       
	do i=1,im
            tlm=tlm+tdlm                                                      
            sinphi=ctph0*stph+stph0*ctph*cos(tlm)                             
            glatr(i,j)=asin(sinphi)
            coslam=ctph*cos(tlm)/(cos(glatr(i,j))*ctph0)
     .            -tan(glatr(i,j))*tan(tph0) 
            coslam=min(coslam,1.)
            fact=1.                                                           
            if (tlm .gt. 0.0) fact=-1.
            glonr(i,j)=-tlm0d*dtr+fact*acos(coslam)                            
         enddo
      enddo
c
c *** Derived vertical grid constants.
c
      ef4t=.5*dt/cp                                                     
      f4q =   -dt*dtad                                                  
      f4d =-.5*dt*dtad                                                  
      do l=1,lm                          
         rdeta(l)=1./detac(l)                        
         f4q2(l)=-.25*dt*dtad/detac(l)               
      enddo
c
c *** Boundary masks.
c
      do j=1,jm                                          
	do i=1,im
         hbm2(i,j)=0.                                                        
         vbm2(i,j)=0.                                                        
         vbm3(i,j)=0.                                                        
	enddo
      enddo
	
      do j=3,jm2                                            
         do i=2,im-1-mod(j+1,2)                                      
            hbm2(i,j)=1.                                                        
         enddo
	 do i=2,im-1-mod(j,2)
            vbm2(i,j)=1.
         enddo
      enddo

      do j=4,jm3                                            
          do i=2+mod(j+1,2),im-2
            vbm3(i,j)=1.                                                        
          enddo
      enddo
c
c *** Topography masks & maximum vertical indices.
c
	do j=1,jm
         do i=1,im                                          
         lmh(i,j)=lm                                                         
         lmv(i,j)=lm                                                         
         enddo
	enddo
      do l=1,lm                              
      do j=1,jm                                          
	do i=1,im
         htm(i,j,l)=1.                                                       
         vtm(i,j,l)=1.                                                       
        enddo
      enddo
      enddo
c
c *** Topography masks & maximum vertical indices.
c
      if (.not. sigma) then                                            
c
      DO L=1,LM
       DO J=1,JM
        DO I=1,IM
         IF(LMH(I,J).EQ.LM.AND.ETA(L+1).GE.RES(I,J)) LMH(I,J)=L
        ENDDO
       ENDDO
      ENDDO

c                                                                       
         do l=1,lm                              
         do j=1,jm                                          
	 do i=1,im
            if (eta(l+1) .gt. res(i,j)) htm(i,j,l)=0.
              IF(I.EQ.33.AND.J.EQ.33) THEN
               DIFF=ETA(L+1)-RES(I,J)
               WRITE(6,28000)I,J,L,ETA(L+1),RES(I,J),DIFF
28000          FORMAT(1X,3I6,3(1X,E18.11))
              ENDIF
         enddo
         enddo
	 enddo
c                                                                       
      DO L=1,LM
       DO J=2,JM1
        IHE=MOD(J+1,2)
        IHW=IHE-1
        KHL=1+MOD(J,2)
        KHH=IM-1
        DO I=KHL,KHH
           IF(ETA(L+1).GT.RES(I,J)) THEN
                  VTM(I+IHE,J,L)=0.
                  VTM(I+IHW,J,L)=0.
                  VTM(I,J-1,L)=0.
                  VTM(I,J+1,L)=0.
           ENDIF
        ENDDO
       ENDDO
      ENDDO

c                                                                       

      DO L=1,LM
       DO J=2,JM1
       KVL=2-MOD(J,2)
         DO I=KVL,IM-1
           IF(LMV(I,J).EQ.LM.AND.VTM(I,J,L).LT.0.1) LMV(I,J)=L-1
         ENDDO
       ENDDO
      ENDDO

c
      endif
c

C next 80 lines or so from operational CNSTS.f
C--------------SPREADING OF LATITUDE DEPENDENT CONSTANTS----------------

      DO J=1,JM
        DO I=1,IM
c       KHH=IM-1+MOD(J,2)
c       DO I=1,KHH
          DX(I,J)=DXJ(J)
Cmp          WPDAR(I,J)=WPDARJ(J)*HBM2(I,J)
	if (.not. SIGMA) then
          WPDAR(I,J)=WPDARJ(J)*HBM2(I,J)
	else
	  if (I .eq. 1 .and. J .eq. 1) then
	 	write(6,*) 'REDUCING WPDAR BECAUSE RUNNING IN SIGMA!!!'
	  endif
	WPDAR(I,J)=WPDARJ(J)*HBM2(I,J)*0.5
	endif
          FCP(I,J)=FCPJ(J)
          FDIV(I,J)=FDIVJ(J)
          FAD(I,J)=FADJ(J)
          HDAC(I,J)=HDACJ(J)*1.25*HBM2(I,J)
c       ENDDO
c       KVH=IM-MOD(J,2)
c       DO I=1,KVH
          HDACV(I,J)=HDACJ(J)
          CPGFU(I,J)=CPGFUJ(J)
          CURV(I,J)=CURVJ(J)
        ENDDO
      ENDDO
C--------------INCREASING DIFFUSION ALONG THE BOUNDARIES----------------
      DO J=3,JM2
        IF (J.LE.5.OR.J.GE.JM4) THEN
          KHH=IM-2+MOD(J,2)
          DO I=2,KHH
           HDAC(I,J)=HDAC(I,J)* DFC
          ENDDO
        ELSE
          KHH=2+MOD(J,2)
          DO I=2,KHH
           HDAC(I,J)=HDAC(I,J)* DFC
          ENDDO
          KHH=IM-2+MOD(J,2)
          DO I=IM-2,KHH
           HDAC(I,J)=HDAC(I,J)* DFC
          ENDDO
        ENDIF
      ENDDO
CMEB          KHL2(J)=IM*(J-1)-(J-1)/2+2 ! 2 if j odd or j even
CMEB          KHL3(J)=IM*(J-1)-J/2+3     ! 3 if j odd, 2 if j even
CMEB          KHH2(J)=IM*J-J/2-1         ! IM-1 if j odd, IM-2 if j even
CMEB          KHH3(J)=IM*J-(J+1)/2-1     ! IM-2 if j odd, or even

C-----------------------------------------------------------------------
      DO J=1,JM
CMEB          KVL0(J)=IM*(J-1)-J/2+1  1 if j odd, 1 if J even
CMEB          KVH0(J)=IM*J-(J+1)/2    IM-1 if j odd, IM if J even
        KVH=IM-MOD(J,2)
        DO I=1,KVH
         DDMPU(I,J)=DDMPUJ(J)*VBM2(I,J)
         DDMPV(I,J)=DDMPVJ(J)*VBM2(I,J)
         HDACV(I,J)=HDACV(I,J)*VBM2(I,J)
	 write(1201,*)I,J,DDMPU(I,J),VBM2(I,J)
        ENDDO
      ENDDO
C--------------INCREASING DIFFUSION ALONG THE BOUNDARIES----------------
      DO J=3,JM2
CMEB          KVL3(J)=IM*(J-1)-(J-1)/2+2 2 if j odd, 3 if j even
CMEB          KVL2(J)=IM*(J-1)-J/2+2     2 if j odd or even
CMEB          KVH2(J)=IM*J-(J+1)/2-1     IM-2 if j odd, IM-1 if j even
CMEB          KVH3(J)=IM*J-J/2-2         IM-2 if j odd or even
        IF (J.LE.5.OR.J.GE.JM4) THEN
          KVH=IM-1-MOD(J,2)
          DO I=2,KVH
            DDMPU(I,J)=DDMPU(I,J)*DDFC
            DDMPV(I,J)=DDMPV(I,J)*DDFC
            HDACV(I,J)=HDACV(I,J)* DFC
          ENDDO
        ELSE
          KVH=3-MOD(J,2)
          DO I=2,KVH
            DDMPU(I,J)=DDMPU(I,J)*DDFC
            DDMPV(I,J)=DDMPV(I,J)*DDFC
            HDACV(I,J)=HDACV(I,J)* DFC
          ENDDO
          KVH=IM-1-MOD(J,2)
          DO I=IM-2,KHH
            DDMPU(I,J)=DDMPU(I,J)*DDFC
            DDMPV(I,J)=DDMPV(I,J)*DDFC
            HDACV(I,J)=HDACV(I,J)* DFC
          ENDDO
        ENDIF
      ENDDO

c
c *** Surface parameters.


	if (REAN) then
	write (6,*) 'PROCESSING THE REANALYSIS FILE!'

C       special section for reanalysis data
C

        allocate(GGSTC1(nx,ny),GGSTC2(nx,ny))
        allocate(GGSMC1(nx,ny),GGSMC2(nx,ny))
        allocate(SSTRAW(nx,ny),SNOWRAW(nx,ny))
        allocate(GGLAND(nx,ny),GGICE(nx,ny))

        call snowsoilsst(im,jm,tlm0d,tph0d,dlmd,dphd,rean_sfc,
     +  GGSTC1,GGSTC2,GGSMC1,GGSMC2,SSTRAW,SNOWRAW,GGLAND,GGICE,GDS)

        call process_gaus(im,jm,tlm0d,tph0d,dlmd,dphd,sm,stc,smc,sst,
     +  sice,si,GGSTC1,GGSTC2,GGSMC1,GGSMC2,SSTRAW,SNOWRAW,GGLAND,GGICE,
     +  gds)

	do J=1,JM
	do I=1,IM
	if (smc(i,j,1) .eq. 0. .and. sm(i,j) .eq. 0.) then
	write(6,*) 'post process_gaus, bad smc ', i,j
	endif
	enddo
	enddo

        deallocate(GGSTC1,GGSTC2,GGSMC1,GGSMC2,SSTRAW,
     +		SNOWRAW,GGLAND,GGICE)


C       end special section for reanalysis data
	endif
C----------------------------------------
c
C READ MAX SNOW ALBEDO FILE
C 90N-20N VIA DAVE ROBINSON, JCAM, 1985, VOL. 24, 402-411
C 20N-90S VIA LAND-SFC TYPE 'CORRELATED' TO ROBINSON DATABASE
C SNALMX = GLOBAL 1-DEG x 1-DEG MAXIMUM SNOW ALBEDO
C Units in percent, later converted to fraction
C values over sea=0.0 (non-land mass)
C values over land between 21.0 and 80.0
	write(6,*) 'snalmx values read in:'
      READ(20)SNALMX
	do J=150,10,-10
	write(6,663) (.01*snalmx(i,j),i=1,361,20)
	enddo
  663	format(20(f4.2,1x))
C-----------------------------------------------------------------------
C READ ALBEDO FILES
C GLOBAL 1-DEG x 1-DEG (4) SEASONAL SNOWFREE ALBEDO
C 90N-90S VIA MATTHEWS ***get reference
C ALBC1 = 3-MONTH AVERAGE CENTERED ON 31 Jan
C ALBC2 = 3-MONTH AVERAGE CENTERED ON 30 Apr
C ALBC3 = 3-MONTH AVERAGE CENTERED ON 31 Jul
C ALBC4 = 3-MONTH AVERAGE CENTERED ON 31 Oct
C Units in percent, later converted to fraction
C values over sea=6.0 (non-land mass)
C values over land between 11.0 and 32.0 (exception, glacier=0.75)
C *NOTE: in future, replace w/high-res 0.144-deg albedo NESDIS database
C   ...and follow similar spatial/temporal averaging for greenness frac
C

	
      READ(21)ALBC1
      READ(22)ALBC2
      READ(23)ALBC3
      READ(24)ALBC4

C fIND jULIAN DAY OF YEAR TO DO THE TEMPORAL INTERPOLATION
      JULD = JULM(IDAT(1)) + IDAT(2)
      IF(JULD.LE.32) THEN
        S1 = 32 - JULD
        S2 = JULD + 30
        WGT1 = S2/(S1+S2)
        WGT2 = S1/(S1+S2)
        DO J = 1,180
        DO I = 1,361
        ALBC(I,J) = WGT1 * ALBC1(I,J) + WGT2 * ALBC4(I,J)
        END DO
        END DO
      ELSE IF(JULD.LE.121.AND.JULD.GT.32) THEN
        S1 = 121 - JULD
        S2 = JULD - 32
        WGT1 = S2/(S1+S2)
        WGT2 = S1/(S1+S2)
        DO J = 1,180
        DO I = 1,361
        ALBC(I,J) = WGT1 * ALBC2(I,J) + WGT2 * ALBC1(I,J)
        END DO
        END DO
      ELSE IF(JULD.LE.213.AND.JULD.GT.121) THEN
        S1 = 213 - JULD
        S2 = JULD - 121
        WGT1 = S2/(S1+S2)
        WGT2 = S1/(S1+S2)
        DO J = 1,180
        DO I = 1,361
        ALBC(I,J) = WGT1 * ALBC3(I,J) + WGT2 * ALBC2(I,J)
        END DO
        END DO
      ELSE IF(JULD.LE.305.AND.JULD.GT.213) then
        S1 = 305 - JULD
        S2 = JULD - 213
        WGT1 = S2/(S1+S2)
        WGT2 = S1/(S1+S2)
        DO J = 1,180
        DO I = 1,361
        ALBC(I,J) = WGT1 * ALBC4(I,J) + WGT2 * ALBC3(I,J)
        END DO
        END DO
      ELSE
        S1 = 365 - JULD + 32
        S2 = JULD - 305
        WGT1 = S2/(S1+S2)
        WGT2 = S1/(S1+S2)
        DO J = 1,180
        DO I = 1,361
        ALBC(I,J) = WGT1 * ALBC1(I,J) + WGT2 * ALBC4(I,J)
        END DO
        END DO
      END IF
C
C BEGIN SPATIAL INTERPOLATION FOR BASELINE SNOWFREE ALBEDO AND MAX SNOW
C ALBEDOS
C
      DO J=1,JM
        DO I=1,IM
Cnew
                ALBASE(I,J)=20.
                MXSNAL(I,J)=55.
Cnew
          IF (SM(I,J) .LT. 0.9) THEN
C
C-----------------------------------------------------------------------
C IF LANDMASS, DETERMINE LAT,LON AND WEIGHTS
C
      ELAT=90.+glatr(I,J)/DTR
      ELON=360.-glonr(I,J)/DTR
      IF(ELON.GT.360.)ELON=ELON-360.
      ILON1=INT(ELON)
      DIF=ELON-ILON1
      IF(DIF.GT.D5)ILON1=ILON1+1
      IF(ILON1.EQ.D00)ILON1=360
      ILON2=ILON1+1
      ILAT1=INT(ELAT)
      DIF=ELAT-ILAT1
      IF(DIF.GT.D5)ILAT1=MIN(ILAT1+1,179)
      ILAT2=ILAT1+1
      W1=ELON-ILON1+D5
      IF(W1.LT.D00)W1=W1+360.
      W2=ELAT-ILAT1+D5
C
      AR1=W1*W2
      AR2=W1*(H1-W2)
      AR3=(H1-W1)*(H1-W2)
      AR4=(H1-W1)*W2
C-----------------------------------------------------------------------
C INTERPOLATE BASELINE SNOWFREE ALBEDO TO E GRID
C INPUT:  ALBC (GLOBAL 1-DEG x 1-DEG)
C OUTPUT: ALBASE (SPECIFIED ETA GRID)
C
C INTERPOLATE MAX SNOW ALBEDO TO E GRID
C INPUT:  SNALMX (GLOBAL 1-DEG x 1-DEG)
C OUTPUT: MXSNAL (SPECIFIED ETA GRID)
C
C DOING BOTH BASELINE SNOWFREE ALBEDO AND MAX SNOW ALBEDO INTERPOLATIONS
C IN THE SAME BLOCK IS POSSIBLE SINCE THEY HAVE IDENTICAL LAND-SEA MASKS
c            IF ( (ALBC(ILON2,ILAT2) .NE. 0.) .AND.
c     .           (ALBC(ILON2,ILAT1) .NE. 0.) .AND.
c     .           (ALBC(ILON1,ILAT1) .NE. 0.) .AND.
c     .           (ALBC(ILON1,ILAT2) .NE. 0.) ) THEN
c the quarterly mathews albedo data base sea values = 6.0
C beginning lat/long=-90,0 (southpole, prime meredian)
c lowest land value = 11.0
c max non-glacial value about 32.0 +/-
c glacial value = 75.0
            IF ( (ALBC(ILON2,ILAT2) .GT. 7.) .AND.
     .           (ALBC(ILON2,ILAT1) .GT. 7.) .AND.
     .           (ALBC(ILON1,ILAT1) .GT. 7.) .AND.
     .           (ALBC(ILON1,ILAT2) .GT. 7.) ) THEN
C-----------------------------------------------------------------------
C ALL 4 SURROUNDING POINTS ARE LAND POINTS
              ALBASE(I,J)=AR1*ALBC(ILON2,ILAT2)+
     .                    AR2*ALBC(ILON2,ILAT1)+
     .                    AR3*ALBC(ILON1,ILAT1)+
     .                    AR4*ALBC(ILON1,ILAT2)
              MXSNAL(I,J)=AR1*SNALMX(ILON2,ILAT2)+
     .                    AR2*SNALMX(ILON2,ILAT1)+
     .                    AR3*SNALMX(ILON1,ILAT1)+
     .                    AR4*SNALMX(ILON1,ILAT2)
C
c            if(i.eq.303 .and .j.eq.143 .or.
c     &         i.eq.309 .and. j.eq.143) then
c           if(mod(i,20).eq.0 .and. mod(j,20).eq.0) then
c              print *,'albase ',i,j,ar1,ar2,ar3,ar4,ALBC(ILON2,ILAT2),
c     &            ALBC(ILON2,ILAT1),ALBC(ILON1,ILAT1),
c     &            ALBC(ILON1,ILAT2),ALBASE(I,J)
c              print *,'mxsnal ',i,j,ar1,ar2,ar3,ar4,SNALMX(ILON2,ILAT2),
c     &            SNALMX(ILON2,ILAT1),SNALMX(ILON1,ILAT1),
c     &            SNALMX(ILON1,ILAT2),MXSNAL(I,J)
c            endif
C
            ELSE
C-----------------------------------------------------------------------
C ONE OR MORE OF THE 4 SURROUNDING POINT ARE NOT LAND POINTS
C TAKE NEAREST NEIGHBOR LAND POINT IN THE FOLLOWING ORDER:
C (ILON2,ILAT2),(ILON2,ILAT1),(ILON1,ILAT1),(ILON1,ILAT2)
C
              ARX=-999.
              IF (ALBC(ILON2,ILAT2) .GT. 7.) THEN
                IF (AR1 .GT. ARX) THEN
                  ARX=AR1
                  ILONX=ILON2
                  ILATX=ILAT2
                ENDIF
              ENDIF
              IF (ALBC(ILON2,ILAT1) .GT. 7.) THEN
                IF (AR2 .GT. ARX) THEN
                  ARX=AR2
                  ILONX=ILON2
                  ILATX=ILAT1
                ENDIF
              ENDIF
              IF (ALBC(ILON1,ILAT1) .GT. 7.) THEN
                IF (AR3 .GT. ARX) THEN
                  ARX=AR3
                  ILONX=ILON1
                  ILATX=ILAT1
                ENDIF
              ENDIF
              IF (ALBC(ILON1,ILAT2) .GT. 7.) THEN
                IF (AR4 .GT. ARX) THEN
                  ARX=AR4
                  ILONX=ILON1
                  ILATX=ILAT2
                ENDIF
              ENDIF
C-----------------------------------------------------------------------
C Use nearest land neighbor:
              IF (ARX .GT. 0.) THEN
                ALBASE(I,J)=ALBC(ILONX,ILATX)
                MXSNAL(I,J)=SNALMX(ILONX,ILATX)
              ELSE
C-----------------------------------------------------------------------
C NO SURROUNDING POINTS ARE LAND (ARX=-999):
C SET DEFAULT SNOWFREE ALBEDO=20 PERCENT, DEFAULT SNOWALB=55 PERCENT
                ALBASE(I,J)=20.
                MXSNAL(I,J)=55.
                PRINT *,'AT: LAT,LON',ELAT,ELON
                PRINT *,'SNOWFREE ALBEDO SET TO DEFAULT VALUE OF 20%'
                PRINT *,'MAX SNOW ALBEDO SET TO DEFAULT VALUE OF 55%'
C-----------------------------------------------------------------------
C end of ALBASE,MXSNAL land points<4 block
              ENDIF
C-----------------------------------------------------------------------
C end of ALBASE,MXSNAL interpolation block
            ENDIF
C-----------------------------------------------------------------------
C end of land (SM=0) block
          ENDIF
C-----------------------------------------------------------------------
C CONVERT ALBEDO UNITS: PERCENT TO FRACTION
          ALBASE(I,J)=ALBASE(I,J)*D01
          MXSNAL(I,J)=MXSNAL(I,J)*D01
        ENDDO
      ENDDO

Cmp    END ALBEDO SECTION **************************************************

      tph=sb-dph                                                
      print *,' '
CJLG     moved to external module
CJLG        if (.NOT. REAN) then
CJLG        CALL SSTHIRES(sst,sm,glatr,glonr,IDAT,6,DTR)
CJLG        write(6,*) 'back from SSTHIRES call'
CJLG        endif

	lusaf=.true.
	if (.NOT. REAN) then
	CALL SNOHIRES(si,sm,glatr,glonr,6,LUSAF,DTR)
	write(6,*) 'back from SNOHIRES call'
	endif

	write(6,*) 'STPH0= ', stph0
      do j=1,jm                                            
	TLM=WB-TDLM+MOD(J+1,2)*DLM
        TPH=TPH+DPH
        STPH=SIN(TPH)
        CTPH=COS(TPH)
        KHH=IM-MOD(J,2)
	do i=1,im
          TLM=TLM+TDLM
          CTLM=COS(TLM)
          STLM=SIN(TLM)
          APH=ASIN(STPH0*CTPH*CTLM+CTPH0*STPH)

Cmp	operational line only used if dont have snowdepth data
        if (.NOT.lusaf.and.sm(i,j).lt.0.9) then
        write(6,*) 'NOT LUSAF!'
	SI(I,J)=0.333*SI(I,J)*SIN(APH)
	endif

        enddo
      enddo
	do J=1,JM
	do I=1,IM
	if (smc(i,j,1) .eq. 0. .and. sm(i,j) .eq. 0.) then
	write(6,*) 'pre important doloop, bad smc ', i,j
	endif
	enddo
	enddo


C------------------------------------------------
c
C  DO-LOOP BELOW DEFINES SOME MASKS AND CHANGES THE MEANING OF OTHERS.
C  MASKS DEFINITIONS FOLLOWING THIS DO-LOOP ARE THOSE CONVENTIONS
C    EXPECTED BY ETA MODEL.
C  --
C  PRIOR TO THIS DO-LOOP:
C
C    SM      LANDMASS (=0), OPEN SEA OR SEAICE (=1)
C    SI      SM=0
C              SI=SNOW DEPTH (0 TO 10.9 M)
C            SM=1
C              SI=1, SEAICE
C              SI=0, OPEN SEA
C    SNO     UNDEFINED
C    SICE    UNDEFINED
C    ALBEDO  UNDEFINED
C  --
C  AFTER THIS DO-LOOP:
C
C    SM      LANDMASS OR SEAICE (=0), OPEN SEA (=1)
C    SI      SNOW DEPTH (M) (ALWAYS ZERO OVER OPEN SEA OR SEAICE)
C    SNO     SNOW WATER EQUIVALENT (M), = 0.2*SI
C    SICE    SICE=0
C              LANDMASS OR OPEN SEA
C            SICE=1
C              SEAICE
C    ALBEDO  SM=1
C              ALBEDO=0.06 OVER OPEN SEA
C            SM=0
C              SICE=0, ALBEDO=DYNAMIC LAND ALBEDO, A FUNCTION OF THE
C                BASELINE SNOWFREE ALBEDO (ALBASE, FROM MATTHEWS) AND
C                INCLUDING A SNOW EFFECT (FROM ROBINSON)
C                SO, DYNAMIC LAND ALBEDO=FCT(ALBASE,MXSNAL,SNO)
C              SICE=1, ALBEDO=0.60 OVER SEAICE
C-----------------------------------------------------------------------
C
      DO J=1,JM
        DO I=1,IM

   	  IF (REAN .and. SICE(I,J) .eq. 1) SI(I,J)=1

	  SICE(I,J)=0.
          SNO (I,J)=0.
          IF(ALBEDO(I,J).GT.0.99999999)    ALBEDO(I,J)=0.99999999

        IF(SM(I,J).GT.0.9) THEN
C  SEA
          EPSR(I,J)=.97
          GFFC(I,J)=0.
          ALBEDO(I,J)=.06
	  ALBASE(I,J)=.06

          IF(SI (I,J).GT.0.    ) THEN
C  SEA-ICE
            SM(I,J)=0.
	    SI(I,J)=0.
            SICE(I,J)=1.
             TG(I,J)=TI0
            GFFC(I,J)=ROI*CI*AKI/DZI
            ALBEDO(I,J)=.60
            ALBASE(I,J)=.60
          END IF

        ELSE
C  LAND
          EPSR(I,J)=1.0
          GFFC(I,J)=ROS*CS*AKS/DZG
          SICE(I,J)=0.
          SNO(I,J)=SI(I,J)*.20
        END IF
        ENDDO
      ENDDO
Cmp
C	convert lat/lon from rad to deg
	do J=1,JM
	do I=1,IM
	glon(I,J)=360.-glonr(I,J)*radfac
	glat(I,J)=glatr(I,J)*radfac
	if (glon(i,j) .ge. 360.) glon(I,J)=glon(I,J)-360.
	enddo
	enddo

	if (.NOT. GRIBSOIL) then
	CALL READSFC(glat,glon,sm,sice,smc,stc,tg)
	endif

C-----------------------------------------------------------------------
      DO J=1,JM
       DO I=1,IM
        GFFC(I,J)=GFFC(I,J)*HBM2(I,J)
        EPSR(I,J)=EPSR(I,J)*HBM2(I,J)
       ENDDO
      ENDDO

c
c *** Coarse look-up table for saturation point.
c
      kthm=jtb                                                          
      kpm=itb                                                           
      kthm1=kthm-1                                                      
      kpm1=kpm-1                                                        
c                                                                       
      thl=210.                                                          
      thh=350.                                                          
      pl=pt                                                             
      ph=105000.                                                        
      r1=r
      pt1=pt
c-----------------------------------------------------------------------
csd new convection
      call table(ptbl,ttbl,pt
     .          ,rdq,rdth,rdp,rdthe,pl,thl,qs0,sqs,sthe,the0)
      call tableq(ttblq,rdpq,rdtheq,plq,thl,stheq,the0q)
c-----------------------------------------------------------------------
      len=index(init_out//' ',' ')-1

      open(1,file=init_out(1:len)//'cnst.file'
     .    ,status='unknown',form='unformatted')
c
      write(1)                                                      
     . nfcst,nbc,list
     .,dt,idtad,sigma                                                   
     .,khla,khha,kvla,kvha,khl2,khh2,kvl2,kvh2
c
      write(1) lmh
      write(1) lmv
      write(1) hbm2 
      write(1) vbm2 
      write(1) vbm3 
      write(1) sm
Cmp
        write(6,*) 'sea ice being written'
        do J=jm,1,-JM/5
        write(6,372) (sice(I,J),I=1,im,2)
  372   format(50(f2.0))
        enddo

      write(1) sice
      do l=1,lm
         write(1) ((htm(i,j,l),i=1,im),j=1,jm)
      enddo
      do l=1,lm
         write(1) ((vtm(i,j,l),i=1,im),j=1,jm)
      enddo
      write(1) dy,cpgfv,en,ent,r,pt,tddamp
     .        ,f4d,f4q,ef4t,detac,rdeta,aeta,f4q2,eta,dfl
     .        ,em,emt
      write(1) dx
      write(1) wpdar
      write(1) cpgfu
      write(1) curv
      write(1) fcp
      write(1) fdiv
      write(1) fad
      write(1) f
      write(1) ddmpu
      write(1) ddmpv
      pt2=pt
      write(1) pt2,glatr
      write(1) glonr
      write(1) plq,rdpq,rdtheq,stheq,the0q
      write(1) ros,cs,ds,roi,ci,di
     .        ,pl,thl,rdq,rdth,rdp,rdthe
     .        ,deta2,aeta2,dfrlg
     .        ,qs0,sqs,sthe,the0
      write(1) mxsnal
      write(1) epsr
Cmp
	do J=1,jm
	 do I=1,im
	  tg(I,J)=amax1(stc(I,J,4),273.15)
	 enddo
	enddo
Cmp
      write(1) tg
      write(1) gffc 

C  20020531
C
C       9-point filter SST data
CJLG        do NFLT=1,4
CJLG        do J=2,JM-1
CJLG        do I=2,IM-1
CJLG        SST(I,J)=0.25*SST(I,J)+0.125*(SST(I+1,J)+SST(I-1,J)+
CJLG     &                               SST(I,J+1)+SST(I,J-1))
CJLG     &               +0.0625*(SST(I+1,J+1)+SST(I+1,J-1)+
CJLG     &                        SST(I-1,J+1)+SST(I-1,J-1))
CJLG        enddo
CJLG        enddo
CJLG        enddo
      write(6,*) 'SST WITH ZERO VAL IN CNST.file'
      sst=0.
      write(1) sst
C
      write(6,*) 'sample SST values '
      do J=jm,1,-JM/30
        write(6,369) (sst(I,J),I=1,im,IM/12)
  369	format(31(f4.0,x))
      enddo
C
      write(1) albase
      write(1) hdac
      write(1) hdacv
      write(1) ttblq
      write(1) ptbl,ttbl
     .        ,r1,pt1,tsph
     .        ,wbd,sbd,tlm0d,tph0d,dlmd,dphd,tlm0d,dp30
     .        ,x1p,y1p,ixm,iym
     .        ,deta1,aeta1,eta1

C NEW CONTENT *********************************************

C    TIME TO WRITE THE SOIL MODEL STUFF TO NHIBU (1)
C
	
C **   FIRST THE VEGETATION TYPES ************************

       REWIND 30
       READ (30) I1D
C        SET DEFAULT FOR ETA LAND FAR FROM LAND IN GLOBAL
       IONETA=7
              CALL PUTEM(glat,glon,I1D,sm,si,IONETA)
	write(6,*) 'writing vegetation type'
	if (vegflag) then
          read(45)IONETA
          write(6,*) 'read new vegetation map'
        endif
	write(6,*)'Done vegetation type'
       WRITE(1) IONETA

	write(6,*) 'sample veg type values '
	do J=jm,1,-jm/30
	write(6,431) (IONETA(I,J),I=1,im,im/10)
  431	format(60(I2,x))
	enddo

C **   SECOND THE SOIL TYPES TYPES **********************

      IF ((NSOTYP.eq.15).or.(NSOTYP.eq.18)) THEN
          read(49)IONETA
      ELSEIF (NSOTYP.eq.9) THEN	
	 REWIND 31
         READ (31) I1D
C        SET DEFAULT FOR ETA LAND FAR FROM LAND IN GLOBAL
       IONETA=2
          CALL PUTEM(glat,glon,I1D,sm,si,IONETA)
       DO J = 1,JM
         DO I = 1,IM
         IF(IONETA(I,J).EQ.13) IONETA(I,J) = 9
         END DO
       END DO
      ENDIF
      
       write(6,*) 'writing soil type'
       WRITE(1) IONETA
	write(6,*) 'sample soil type values'
	do J=jm,1,-jm/30
	write(6,397) (IONETA(I,J),I=1,im,im/10)
  397	format(60(i2,x))
	enddo
       
C Add sh2o stuff right here
	DO 200 J=1,JM
	DO 100 I=1,IM
	DO 70 K=1,NSOIL
C ----------------------------------------------------------------------
C cold start:  determine liquid soil water content (SH2O)
C SH2O <= SMC for T < 273.149K (-0.001C)
C
Cnew
C
	if (SM(I,J) .eq. 1) then
	STC(I,J,K)=amax1(273.15,STC(I,J,K))
	endif
C
Cnew
C
            IF (STC(I,J,K) .LT. 273.149) THEN
C ----------------------------------------------------------------------
C first guess following explicit solution for Flerchinger Eqn from Koren
C et al, JGR, 1999, Eqn 17 (KCOUNT=0 in FUNCTION FRH2O).
              BX = BETA(IONETA(I,J))
              IF ( BETA(IONETA(I,J)) .GT. BLIM ) BX = BLIM
              FK = (((HLICE/(GRAV*(-PSIS(IONETA(I,J)))))*
     1             ((STC(I,J,K)-T0)/STC(I,J,K)))**
     1             (-1/BX))*SMCMAX(IONETA(I,J))
              IF (FK .LT. 0.02) FK = 0.02
              SH2O(I,J,K) = MIN ( FK, SMC(I,J,K) )
C ----------------------------------------------------------------------
C now use iterative solution for liquid soil water content using
C FUNCTION FRH2O (from the Eta "NOAH" land-surface model) with the
C initial guess for SH2O from above explicit first guess.
              SH2O(I,J,K)=FRH2O(STC(I,J,K),SMC(I,J,K),SH2O(I,J,K),
     .                    SMCMAX(IONETA(I,J)),BETA(IONETA(I,J)),
     .                    PSIS(IONETA(I,J)))
            ELSE
C ----------------------------------------------------------------------
C SH2O = SMC for T => 273.149K (-0.001C)
              SH2O(I,J,K)=SMC(I,J,K)
c
            ENDIF

   70	CONTINUE
  100	CONTINUE
  200	CONTINUE

C **   THIRD THE SURFACE SLOPE TYPES ***********************

       REWIND 32
       READ (32) I1D
C        SET DEFAULT FOR ETA LAND FAR FROM LAND IN GLOBAL
       IONETA=1
              CALL PUTEM(glat,glon,I1D,sm,si,IONETA)
	write(6,*) 'writing slope type'
       WRITE(1) IONETA

Cmp END NEW CONTENT ***************************************

C-----------------------------------------------------------------------

        write(6,*) 'calling VFRAC'
        CALL VFRAC(im,jm,glatr,glonr,TLM0D,TPH0D,DLMD,DPHD,sm,si,
     +  vegfrc)

        write(6,*) 'sample VEGFRC values (* 1000) '
        do J=jm,1,-jm/30
        write(6,396) (vegfrc(I,J)*1000,I=1,im,im/12)
  396   format(31(f4.0,x))
        enddo

C-----------------------------------------------------------------------
C DETERMINE ALBEDO OVER LAND
      DO J=1,JM
        DO I=1,IM
          IF(SM(I,J).LT.0.9.AND.SICE(I,J).LT.0.9) THEN
C SNOWFREE ALBEDO
            IF ( (SNO(I,J) .EQ. 0.0) .OR.
     .           (ALBASE(I,J) .GE. MXSNAL(I,J) ) ) THEN
CC        write(6,*) 'set albedo to albase: ', ALBASE(I,J)
              ALBEDO(I,J) = ALBASE(I,J)
            if(i.eq.303 .and .j.eq.143 .or.
     &         i.eq.309 .and. j.eq.143) then
                print *,i,j,sno(i,j),albase(i,j),mxsnal(i,j)
            endif
            ELSE
C MODIFY ALBEDO IF SNOWCOVER:
C BELOW SNOWDEPTH THRESHOLD...
              IF (SNO(I,J) .LT. SNUP) THEN
                RSNOW = SNO(I,J)/SNUP
                SNOFAC = 1. - ( EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))
C ABOVE SNOWDEPTH THRESHOLD...
              ELSE
                SNOFAC = 1.0
              ENDIF
C CALCULATE ALBEDO ACCOUNTING FOR SNOWDEPTH AND VEGFRC
              ALBEDO(I,J) = ALBASE(I,J)
     .          + (1.0-VEGFRC(I,J))*SNOFAC*(MXSNAL(I,J)-ALBASE(I,J))
            ENDIF
          END IF
        ENDDO
      ENDDO


      write(1) vegfrc

CGSM	sldpth(1)=0.1
CGSM	sldpth(2)=0.3
CGSM	sldpth(3)=0.6
CGSM	sldpth(4)=1.0
CGSM	rtdpth(1)=0.1
CGSM	rtdpth(2)=0.3
CGSM	rtdpth(3)=0.6
CGSM	rtdpth(4)=0.0

CIsabel
       sldpth(1)=0.1
       sldpth(2)=0.3
       sldpth(3)=0.6
       sldpth(4)=1.0
       sldpth(5)=1.5
       sldpth(6)=2.1
       sldpth(7)=2.8
       sldpth(8)=3.6
C    
       rtdpth(1)=0.1
       rtdpth(2)=0.3
       rtdpth(3)=0.6
       rtdpth(4)=0.0
       rtdpth(5)=0.0
       rtdpth(6)=0.0
       rtdpth(7)=0.0
       rtdpth(8)=0.0
CIsabel

C	write(6,*) 'writing sldpth ... ', sldpth
      write(1) sldpth  
C	write(6,*) 'writing rtdpth... ', rtdpth
      write(1) rtdpth   
      close(1)
c
      return                                     
      end
c
c===============================================================================
c
      subroutine spline(jtb,nold,xold,yold,y2,nnew,xnew,ynew,p,q)           
c
c *** This is a one-dimensional cubic spline fitting routine 
c        programed for a small scalar machine.                       
c
c *** Programer: Z. Janjic, Yugoslav Fed. Hydromet. Inst., Beograd 
c
c *** nold - number of given values of the function.  must be ge 3. 
c     xold - locations of the points at which the values of the     
c            function are given.  must be in ascending order.       
c     yold - the given values of the function at the points xold.   
c     y2   - the second derivatives at the points xold.  if natural 
c            spline is fitted y2(1)=0. and y2(nold)=0. must be      
c            specified.                                             
c     nnew - number of values of the function to be calculated.     
c     xnew - locations of the points at which the values of the     
c            function are calculated.  xnew(k) must be ge xold(1)   
c            and le xold(nold).                                     
c     ynew - the values of the function to be calculated.           
c     p, q - auxiliary vectors of the length nold-2.                
c                                                                   
      real*4 xold(jtb),yold(jtb),y2(jtb),p(jtb),q(jtb)                        
     .      ,xnew(jtb),ynew(jtb)                                              
c_______________________________________________________________________________
c
      noldm1=nold-1                                                     
c                                                                       
      dxl=xold(2)-xold(1)                                               
      dxr=xold(3)-xold(2)                                               
      dydxl=(yold(2)-yold(1))/dxl                                       
      dydxr=(yold(3)-yold(2))/dxr                                       
      rtdxc=.5/(dxl+dxr)                                                
c                                                                       
      p(1)= rtdxc*(6.*(dydxr-dydxl)-dxl*y2(1))                          
      q(1)=-rtdxc*dxr                                                   
c                                                                       
      if(nold.eq.3) go to 700                                           
c-----------------------------------------------------------------------
      k=3                                                               
c                                                                       
 100  dxl=dxr                                                           
      dydxl=dydxr                                                       
      dxr=xold(k+1)-xold(k)                                             
      dydxr=(yold(k+1)-yold(k))/dxr                                     
      dxc=dxl+dxr                                                       
      den=1./(dxl*q(k-2)+dxc+dxc)                                       
c                                                                       
      p(k-1)= den*(6.*(dydxr-dydxl)-dxl*p(k-2))                         
      q(k-1)=-den*dxr                                                   
c                                                                       
      k=k+1                                                             
      if(k.lt.nold) go to 100                                           
c-----------------------------------------------------------------------
 700  k=noldm1                                                          
c                                                                       
 200  y2(k)=p(k-1)+q(k-1)*y2(k+1)                                       
c                                                                       
      k=k-1                                                             
      if(k.gt.1) go to 200                                              
c-----------------------------------------------------------------------
      k1=1                                                              
c                                                                       
 300  xk=xnew(k1)                                                       
c                                                                       
      do 400 k2=2,nold                                                  
      if(xold(k2).le.xk) go to 400                                      
      kold=k2-1                                                         
      go to 450                                                         
 400  continue                                                          
      ynew(k1)=yold(nold)                                               
      go to 600                                                         
c                                                                       
 450  if(k1.eq.1)   go to 500                                           
      if(k.eq.kold) go to 550                                           
c                                                                       
 500  k=kold                                                            
c                                                                       
      y2k=y2(k)                                                         
      y2kp1=y2(k+1)                                                     
      dx=xold(k+1)-xold(k)                                              
      rdx=1./dx                                                         
c                                                                       
      ak=.1666667*rdx*(y2kp1-y2k)                                       
      bk=.5*y2k                                                         
      ck=rdx*(yold(k+1)-yold(k))-.1666667*dx*(y2kp1+y2k+y2k)            
c                                                                       
 550  x=xk-xold(k)                                                      
      xsq=x*x                                                           
c                                                                       
      ynew(k1)=ak*xsq*x+bk*xsq+ck*x+yold(k)                             
c                                                                       
 600  k1=k1+1                                                           
      if(k1.le.nnew) go to 300                                          
c
      return                                    
      end                                       
c
c===============================================================================
c
      subroutine albedo_vege(ivegetype,albedo)
c
c *** ssib vegetation types (dorman and sellers, 1989; jam)
c
c      1:   broadleaf-evergreen trees  (tropical forest)
c      2:   broadleaf-deciduous tress
c      3:   broadleaf and needleleaf tress (mixed forest)
c      4:   needleleaf-evergreen trees
c      5:   needleleaf-deciduous tress (larch)
c      6:   broadleaf tress with groundcover (savanna)
c      7:   groundcover only (perennial)
c      8:   broadleaf shrubs with perennial groundcover
c      9:   broadleaf shrubs with bare soil
c     10:   dwarf trees and shrubs with groundcover (tundra)
c     11:   bare soil
c     12:   cultivations (the same parameters for the type 7)
c     13:   glacial
c     14:   water, according to eta sm
c
      dimension albedo_veg(14)
c
      data albedo_veg/0.11, 0.19, 0.16, 0.13, 0.19, 0.19, 0.19
     .               ,0.29, 0.29, 0.14, 0.15, 0.19, 0.15, 0.10/
c_______________________________________________________________________________
c
      albedo=albedo_veg(ivegetype)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C	BLOCK DATA GSOIL

C	logical GRIBSOIL
C	DATA GRIBSOIL /.FALSE./

C	END


C-------------------------------------------------------

      FUNCTION FRH2O(TKELV,SMC,SH2O,SMCMAX,B,PSIS)

      IMPLICIT NONE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  PURPOSE:  CALCULATE AMOUNT OF SUPERCOOLED LIQUID SOIL WATER CONTENT
CC  IF TEMPERATURE IS BELOW 273.15K (T0).  REQUIRES NEWTON-TYPE ITERATION
CC  TO SOLVE THE NONLINEAR IMPLICIT EQUATION GIVEN IN EQN 17 OF
CC  KOREN ET AL. (1999, JGR, VOL 104(D16), 19569-19585).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C New version (FEB 2001): much faster and more accurate newton iteration
c achieved by first taking log of eqn cited above -- less than 4
c (typically 1 or 2) iterations achieves convergence.  Also, explicit
c 1-step solution option for special case of parameter Ck=0, which reduces
c the original implicit equation to a simpler explicit form, known as the
c ""Flerchinger Eqn". Improved handling of solution in the limit of
c freezing point temperature T0.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C INPUT:
C
C   TKELV.........Temperature (Kelvin)
C   SMC...........Total soil moisture content (volumetric)
C   SH2O..........Liquid soil moisture content (volumetric)
C   SMCMAX........Saturation soil moisture content (from REDPRM)
C   B.............Soil type "B" parameter (from REDPRM)
C   PSIS..........Saturated soil matric potential (from REDPRM)
C
C OUTPUT:
C   FRH2O.........supercooled liquid water content.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL B
      REAL BLIM
      REAL BX
      REAL CK
      REAL DENOM
      REAL DF
      REAL DH2O
      REAL DICE
      REAL DSWL
      REAL ERROR
      REAL FK
      REAL FRH2O
      REAL GS
      REAL HLICE
      REAL PSIS
      REAL SH2O
      REAL SMC
      REAL SMCMAX
      REAL SWL
      REAL SWLK
      REAL TKELV
      REAL T0

      INTEGER NLOG
      INTEGER KCOUNT

      PARAMETER (CK=8.0)
C      PARAMETER (CK=0.0)
      PARAMETER (BLIM=5.5)
C      PARAMETER (BLIM=7.0)
      PARAMETER (ERROR=0.005)

      PARAMETER (HLICE=3.335E5)
      PARAMETER (GS = 9.81)
      PARAMETER (DICE=920.0)
      PARAMETER (DH2O=1000.0)
      PARAMETER (T0=273.15)

C  ###   LIMITS ON PARAMETER B: B < 5.5  (use parameter BLIM)  ####
C  ###   SIMULATIONS SHOWED IF B > 5.5 UNFROZEN WATER CONTENT  ####
C  ###   IS NON-REALISTICALLY HIGH AT VERY LOW TEMPERATURES    ####
C##################################################################
C
      BX = B
      IF ( B .GT. BLIM ) BX = BLIM
C------------------------------------------------------------------

C INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
      NLOG=0
      KCOUNT=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  IF TEMPERATURE NOT SIGNIFICANTLY BELOW FREEZING (T0), SH2O = SMC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF (TKELV .GT. (T0 - 1.E-3)) THEN

        FRH2O=SMC

      ELSE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (CK .NE. 0.0) THEN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCC OPTION 1: ITERATED SOLUTION FOR NONZERO CK CCCCCCCCCCC
CCCCCCCCCCCC IN KOREN ET AL, JGR, 1999, EQN 17 CCCCCCCCCCCCCCCCC
C
C INITIAL GUESS FOR SWL (frozen content)
        SWL = SMC-SH2O
C KEEP WITHIN BOUNDS.
        IF(SWL .GT. (SMC-0.02)) THEN
          SWL=SMC-0.02
c         SWL=MAX(0.,SMC-0.02)
        ENDIF
        IF(SWL .LT. 0.)  THEN
          SWL=0.
        ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  START OF ITERATIONS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DO WHILE (NLOG .LT. 10 .AND. KCOUNT .EQ. 0)
         NLOG = NLOG+1
!	write(6,*) 'PSIS,GS,HLICE: ', PSIS,GS,HLICE
!	write(6,*) 'CK,SWL,SMCMAX,SMC,SWL,BX: ', 
!     +		CK,SWL,SMCMAX,SMC,SWL,BX
!	write(6,*) 'TKELV,T0: ', TKELV,T0
         DF = ALOG(( PSIS*GS/HLICE ) * ( ( 1.+CK*SWL )**2. ) *
     &        ( SMCMAX/(SMC-SWL) )**BX) - ALOG(-(TKELV-T0)/TKELV)
         DENOM = 2. * CK / ( 1.+CK*SWL ) + BX / ( SMC - SWL )
         SWLK = SWL - DF/DENOM
!	write(6,*) 'DF,DENOM,SWLK: ', DF,DENOM,SWLK
         DSWL=ABS(SWLK-SWL)
         SWL=SWLK
C KEEP WITHIN BOUNDS.
         IF(SWL .GT. (SMC-0.02)) THEN
           SWL=SMC-0.02
c          SWL=MAX(0.,SMC-0.02)
         ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)
CC WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF ( DSWL .LE. ERROR )  THEN
           KCOUNT=KCOUNT+1
         END IF
        END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  END OF ITERATIONS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C LOWER BOUND FOR SWL (ICE CONTENT)
        IF(SWL .LT. 0.)  THEN
          SWL=0.
        ENDIF
C UPPER BOUND FOR SWL ALREADY APPLIED WITHIN DO-BLOCK
        FRH2O = SMC - SWL
C
CCCCCCCCCCCCCCCCCCCCCCCC END OPTION 1 CCCCCCCCCCCCCCCCCCCCCCCCCCC

       ENDIF

       IF (KCOUNT .EQ. 0) THEN
!         Print*,'Flerchinger used. Iterations=',NLOG

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0 CCCCCCCC
CCCCCCCCCCCCC IN KOREN ET AL., JGR, 1999, EQN 17  CCCCCCCCCCCCCCC
C
        FK=(((HLICE/(GS*(-PSIS)))*((TKELV-T0)/TKELV))**(-1/BX))*SMCMAX
        IF (FK .LT. 0.02) FK = 0.02
        FRH2O = MIN ( FK, SMC )
C
CCCCCCCCCCCCCCCCCCCCCCCCC END OPTION 2 CCCCCCCCCCCCCCCCCCCCCCCCCC

       ENDIF

      ENDIF

      RETURN
      END 
