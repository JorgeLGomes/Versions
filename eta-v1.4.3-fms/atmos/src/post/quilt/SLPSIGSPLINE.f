      SUBROUTINE SLPSIGSPLINE(PD,FIS,TSIG,QSIG
     +,            SPL,LSL,DETA,PT,PSLP)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C   SUBROUTINE:  SLPSIGSPLINE      SLP REDUCTION
C   02-1981 TO 03 1998 JANJIC, D. JOVIC, S. NICKOVIC
C   09-2000  H CHUANG: MODIFIED TO BE INCLUDED IN OPERATIONAL
C                      MPI QUILT                 
C   NOTE: CURRENT POST ONLY PROCESS UP TO 1000 MB PRESSURE LEVEL  
C         2 ADDITIONAL LEVELS (1025 1050) ARE ADDED DURING SPLINE 
C         FITTING COMPUTATION BUT THE FIELDS ON THESE TWO LEVELS  
C         ARE NOT OUTPUT TO GRID FILES 
C
C ABSTRACT:  THIS ROUTINE COMPUTES THE SEA LEVEL PRESSURE
C            REDUCTION USING CUBIC SPLINE FITTING
C            METHOD OR THE STANDARD NCEP REDUCTION FOR
C            SIGMA COORDINATES.  
C
C USAGE:  CALL SLPSIG FROM SUBROUITNE QUILT
C
C   INPUT ARGUMENT LIST:
C     PD   - SFC PRESSURE MINUS PTOP
C     FIS  - SURFACE GEOPOTENTIAL
C     TSIG - TEMPERATURE 
C     QSIG - SPECIFIC HUMIDITY
C     FI   - GEOPOTENTIAL
C     PT   - TOP PRESSURE OF DOMAIN
C
C   OUTPUT ARGUMENT LIST:
C     PSLP - THE FINAL REDUCED SEA LEVEL PRESSURE ARRAY
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:
C             NONE
C
C-----------------------------------------------------------------------
      INCLUDE "parmeta"
      INCLUDE "PARA.comm"
      INCLUDE "mpif.h"
      INCLUDE "mpp.h"
      PARAMETER (LMP1=LM+1)
cmp      parameter(nsmud=100,lp2=lm+2,lsmp2=lsm+2)
      parameter(nsmud=25,lp2=lm+2,lsmp2=lsm+2)
      parameter(iii=22,jjj=293)
C-----------------------------------------------------------------------
                            P A R A M E T E R
     & (RD=287.04,ROG=RD/9.8,PQ0=379.90516,A2=17.2693882
     &, A3=273.16,A4=35.86,GAMMA=6.5E-3,RGAMOG=GAMMA*ROG
     &, H1M12=1.E-12,g=9.8)
C-----------------------------------------------------------------------
                            R E A L
     & PD(IM,MY_JSD:MY_JED),FIS(IM,MY_JSD:MY_JED)
     &,tmask(IM,MY_JSD:MY_JED)
     &,FI(IM,MY_JSD:MY_JED,LSMP2),TSIG(IM,MY_JSD:MY_JED,LM)
     &,QSIG(IM,MY_JSD:MY_JED,LM)
                            R E A L
     & PSLP(IM,MY_JSD:MY_JED),fsll(IM,MY_JSD:MY_JED)
                            R E A L
     & SPL(LSM)
                            R E A L
     & DETA(LM),AETA(LM),ETA(LMP1),zth(lm+2),hcol  (lm+2)
     &,y2(lm+2),ztsl  (lsmp2),ovrlx (lsmp2),hcolsl(lsmp2)
     &,phld  (lm+2),qhld  (lm+2)
     &,hcol3(3),pcol3(3),y23(3)
!-----------------------------------------------------------------------
                              d a t a
     & y2/lp2*0./
!-----------------------------------------------------------------------
                              d a t a
     & ovrlx/lsmp2*0.175/
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                            I N T E G E R
     & IHE(JM),IHW(JM)
C-----------------------------------------------------------------------
                            L O G I C A L
     & SIGMA,STDRD
C-----------------------------------------------------------------------
      STDRD=.FALSE.
C-----------------------------------------------------------------------
C***
C***  CALCULATE THE I-INDEX EAST-WEST INCREMENTS
C***
c      print*,'in slpsig'
      DO J=1,JM
        IHE(J)=MOD(J+1,2)
        IHW(J)=IHE(J)-1
      ENDDO
C-----------------------------------------------------------------------
C***
C***  INITIALIZE ARRAYS.  LOAD SLP ARRAY WITH SURFACE PRESSURE.
C***
CC!$omp parallel do 
	write(6,*) 'believe j bounds are: ', jsta_i,jend_i
      DO J=JSTA_I,JEND_I
      DO I=1,IM
        PSLP(I,J)=PD(I,J)+PT
      ENDDO
      ENDDO
c       	  
      ztsl(1:lsm)=alog(spl(1:lsm))**2
      ztsl(lsm+1)=alog(102500.)**2  
      ztsl(lsmp2)=alog(105000.)**2
      ztbot=ztsl(lsmp2)
C      
C*** COMPUTE ETA AND AETA 
      ETA(1)=0.0
      DO L=2,LM+1
       ETA(L)=ETA(L-1)+DETA(L-1)
      END DO 
      DO L=1,LM
       AETA(L)=0.5*(ETA(L)+ETA(L+1))
      ENDDO
c

       do 200 j=jsta_i,jend_i
        ihl=1
        ihh=im-1+mod(j,2)
        do 201 i=ihl,ihh
	 pdp=pd(i,j)
	 hcol(lm+1)=fis(i,j)/g
         alp1l=alog(pt+pdp)
         zth(lm+1)=(alog(pdp+pt))**2	 
         do 220 ivi=1,lm
          l=lm+1-ivi
	  alp1u=alog(eta(l)*pdp+pt)
          alp2u=(alog(eta(l)*pdp+pt))**2
	  dh=(qsig(i,j,l)*0.608+1.)*(alp1l-alp1u)*tsig(i,j,l)*ROG
	  zth(l)=alp2u
	  hcol(l)=hcol(l+1)+dh	 
	  alp1l=alp1u
 220     continue	  
         if(zth(lm+1).ge.ztbot)then
          lmd=lm+1      
         else
          lmd=lm+2
          zth(lm+2)=ztbot
          x=ztbot-zth(lm+1)
          dum=(hcol(lm+1)-hcol(lm-5))/(zth(lm+1)-zth(lm-5))
          d2=0.
          hcol(lm+2)=d2*x*x+dum*x+hcol(lm+1)
         endif
         y2(lmd)=0.
c	 if(i.eq.iii.and.j.eq.jjj)then
c	  print*,'sample input height'
c	  do l=1,lmd
c	   print*,l,exp(zth(l)**0.5),hcol(l)
c	  end do
c	 end if  
         call splinef(lmd,zth,hcol,y2,lsmp2,ztsl,hcolsl,phld,qhld)
         do 240 l=1,lsmp2
          fi(i,j,l)=hcolsl(l)
c	  if(i.eq.iii.and.j.eq.jjj)print*,'sample H'
c     +    ,exp(ZTSL(L)**0.5),hcolsl(l)
 240     continue
 201    continue  
 200   continue
c--------------filling remaining undefined mass point boundaries-------
       do 300 l=1,lsmp2
        if(mod(jsta_i,2) .lt. 1)then   !even jsta_i
         do 351 j=jsta_i,jend_i,2
          fi(im,j,l)=fi(im-1,j,l)
 351     continue
        else
	 do 352 j=jsta_i+1,jend_i,2
          fi(im,j,l)=fi(im-1,j,l)
 352     continue
        end if  
	rlx=ovrlx(l)
        do 310 j=jsta_i,jend_i
         ihl=1
         ihh=im-1+mod(j,2)
         do 311 i=ihl,ihh
          href=fis(i,j)/g
          if(href.gt.300.)then
           href=href+1500.
          endif
          if(href.ge.fi(i,j,l))then
           tmask(i,j)=rlx
          else
           tmask(i,j)=0.
          endif
 311     continue
 310    continue	  
c        print*,'before relaxation at',l
        do 320 n=1,nsmud
c	 if(l.eq.41)print*,'calling update at n=',n
	 call update(fi(1,my_jsd,l))
c         if(l.eq.41)print*,'after calling update at n=',n
c	 print*,'jsta_im jend_im=',jsta_im,jend_im
         do 330 j=jsta_im,jend_im
          ihl=1+mod(j,2)
          ihh=im-1
          do 331 i=ihl,ihh
           if(tmask(i,j).gt.0.05) then
             fsll(i,j)=fi(i+ihw(j),j-1,l)+fi(i+ihe(j),j-1,l)
     2         +fi(i+ihw(j),j+1,l)+fi(i+ihe(j),j+1,l)-fi(i,j,l)*4.
           endif
 331      continue
 330     continue
c         PRINT*,'AFTER SMOOTHING'
         do 340 j=jsta_im,jend_im
          ihl=1+mod(j,2)
          ihh=im-1
          do 341 i=ihl,ihh
           if(tmask(i,j).gt.0.05) then
            fi(i,j,l)=fsll(i,j)*tmask(i,j)+fi(i,j,l)
           endif
 341      continue
 340     continue
 320    continue
c        if(i.eq.iii.and.j.eq.jjj)print*,'H after relaxation',i,j
c     +,  fi(i,j,l)	
 300   continue  	           
c--------------sea level pressure---------------------------------------
       do 600 j=jsta_i,jend_i
        ihl=1
        ihh=im-1+mod(j,2)
        do 601 i=ihl,ihh
         if(fis(i,j).gt.-1..and.fis(i,j).lt.1.) then
          pslp(i,j)=pt+pd(i,j) 
         else
	  if(fi(i,j,lsmp2).gt.0)then
	   gar=exp(ztsl(lsmp2)**0.5)+5000.     !assume pslp < 1100 mb
	   hcol3(1)=fi(i,j,lsmp2)+(fi(i,j,lsmp2)-fi(i,j,lsmp2-5))
     +     *(alog(gar)**2.-ztsl(lsmp2))/(ztsl(lsmp2)-ztsl(lsmp2-5))
           hcol3(2)=fi(i,j,lsmp2)
	   hcol3(3)=fi(i,j,lsmp2-1) 
	   pcol3(1)=gar      
           pcol3(2)=exp(ztsl(lsmp2)**0.5)       
           pcol3(3)=exp(ztsl(lsmp2-1)**0.5)	 
	  else   
	   do l=lsmp2,2,-1 
	    if(fi(i,j,l).le.0. .and. fi(i,j,l-1).gt.0.)then 
	     lll=l
	     go to 603
	    end if 
	   end do  
 603       hcol3(1)=fi(i,j,lll)      
           hcol3(2)=fi(i,j,lll-1)    
           hcol3(3)=fi(i,j,lll-2)    
           pcol3(1)=exp(ztsl(lll)**0.5)       
           pcol3(2)=exp(ztsl(lll-1)**0.5)       
           pcol3(3)=exp(ztsl(lll-2)**0.5)        
          end if
          do l=1,3
           y23(l)=0.0
          end do
          if(hcol3(1).eq.hcol3(2) .or. hcol3(2).eq.hcol3(3))
     +    print*,'problem slp',i,j,hcol3(1),hcol3(2),hcol3(3)	  
          call splinef(3,hcol3,pcol3,y23,1,0.0,slp1,phld,qhld)
          pslp(i,j)=slp1
         endif
c         if(i.eq.iii.and.j.eq.jjj)then
c	  print*,'on P after relaxation at ',iii,jjj
c	  do l=1,lsmp2
c	   print*,'l,H,P= ',l,fi(i,j,l),exp(ztsl(l)**0.5)	
c          end do
c	  print*,'pslp = ',pslp(i,j)
c 	 end if    
c         if(pslp(i,j).gt.100200. .and. pslp(i,j).lt.100400.)
c     +   print*,'target pslp',i,j,pslp(i,j)
 601    continue
 600   continue
C--------------------------------------------------------------------
C     SKIP THE STANDARD SCHEME.
C--------------------------------------------------------------------
      GO TO 430
C--------------------------------------------------------------------
C***
C***  IF YOU WANT THE "STANDARD" ETA/SIGMA REDUCTION
C***  THIS IS WHERE IT IS DONE.
C***
C
      doout410: DO J=JSTA_I,JEND_I
        doin410: DO I=1,IM
      IF(FIS(I,J).GE.1.)THEN
        LMA=LM
        ALPP1=ALOG(PD(I,J)+PT)
        SLOP=0.0065*ROG*TSIG(I,J,LM)
        IF(SLOP.LT.0.50)THEN
          SLPP=ALPP1+FIS(I,J)/(RD*TSIG(I,J,LMA))
        ELSE
          TTT=-(ALOG(PD(I,J)+PT)+ALPP1)
     1         *SLOP*0.50+TSIG(I,J,LMA)
          SLPP=(-TTT+SQRT(TTT*TTT+2.*SLOP*
     1          (FIS(I,J)/RD+
     2          (TTT+0.50*SLOP*ALPP1)*ALPP1)))/SLOP
        ENDIF
        PSLP(I,J)=EXP(SLPP)
      ENDIF
      END DO doin410
      END DO doout410
C
C****************************************************************
C     AT THIS POINT WE HAVE A SEA LEVEL PRESSURE FIELD BY
C     EITHER METHOD.  5-POINT AVERAGE THE FIELD ON THE E-GRID.
C****************************************************************
C
  430 CONTINUE
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SPLINEF(NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)           
C                                                                       
C     ******************************************************************
C     *                                                                *
C     *  THIS IS A ONE-DIMENSIONAL CUBIC SPLINE FITTING ROUTINE        *
C     *  PROGRAMED FOR A SMALL SCALAR MACHINE.                         *
C     *                                                                *
C     *  PROGRAMER[ Z. JANJIC, YUGOSLAV FED. HYDROMET. INST., BEOGRAD  *
C     *                                                                *
C     *                                                                *
C     *                                                                *
C     *  NOLD - NUMBER OF GIVEN VALUES OF THE FUNCTION.  MUST BE GE 3. *
C     *  XOLD - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
C     *         FUNCTION ARE GIVEN.  MUST BE IN ASCENDING ORDER.       *
C     *  YOLD - THE GIVEN VALUES OF THE FUNCTION AT THE POINTS XOLD.   *
C     *  Y2   - THE SECOND DERIVATIVES AT THE POINTS XOLD.  IF NATURAL *
C     *         SPLINE IS FITTED Y2(1)=0. AND Y2(NOLD)=0. MUST BE      *
C     *         SPECIFIED.                                             *
C     *  NNEW - NUMBER OF VALUES OF THE FUNCTION TO BE CALCULATED.     *
C     *  XNEW - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
C     *         FUNCTION ARE CALCULATED.  XNEW(K) MUST BE GE XOLD(1)   *
C     *         AND LE XOLD(NOLD).                                     *
C     *  YNEW - THE VALUES OF THE FUNCTION TO BE CALCULATED.           *
C     *  P, Q - AUXILIARY VECTORS OF THE LENGTH NOLD-2.                *
C     *                                                                *
C     ******************************************************************
C                                                                       
                             D I M E N S I O N                          
     2 XOLD(NOLD),YOLD(NOLD),Y2(NOLD),P(NOLD),Q(NOLD)
     3,XNEW(NNEW),YNEW(NNEW)
C-----------------------------------------------------------------------
      NOLDM1=NOLD-1                                                     
C                                                                       
      DXL=XOLD(2)-XOLD(1)                                               
      DXR=XOLD(3)-XOLD(2)                                               
      DYDXL=(YOLD(2)-YOLD(1))/DXL                                       
      DYDXR=(YOLD(3)-YOLD(2))/DXR                                       
      RTDXC=.5/(DXL+DXR)                                                
C                                                                       
      P(1)= RTDXC*(6.*(DYDXR-DYDXL)-DXL*Y2(1))                          
      Q(1)=-RTDXC*DXR                                                   
C                                                                       
      IF(NOLD.EQ.3) GO TO 700                                           
C-----------------------------------------------------------------------
      K=3                                                               
C                                                                       
 100  DXL=DXR                                                           
      DYDXL=DYDXR                                                       
      DXR=XOLD(K+1)-XOLD(K)  
c      if(i.eq.120.and.j.eq.279)then
c       print*,'in spline',k,XOLD(K),XOLD(K+1),dxr 
c      end if                                           
      DYDXR=(YOLD(K+1)-YOLD(K))/DXR                                     
      DXC=DXL+DXR                                                       
      DEN=1./(DXL*Q(K-2)+DXC+DXC)                                       
C                                                                       
      P(K-1)= DEN*(6.*(DYDXR-DYDXL)-DXL*P(K-2))                         
      Q(K-1)=-DEN*DXR                                                   
C                                                                       
      K=K+1                                                             
      IF(K.LT.NOLD) GO TO 100                                           
C-----------------------------------------------------------------------
 700  K=NOLDM1                                                          
C                                                                       
 200  Y2(K)=P(K-1)+Q(K-1)*Y2(K+1)                                       
C                                                                       
      K=K-1                                                             
      IF(K.GT.1) GO TO 200                                              
C-----------------------------------------------------------------------
      K1=1                                                              
C                                                                       
 300  XK=XNEW(K1)                                                       
C                                                                       
      DO 400 K2=2,NOLD                                                  
      IF(XOLD(K2).LE.XK) GO TO 400                                      
      KOLD=K2-1                                                         
      GO TO 450                                                         
 400  CONTINUE                                                          
      YNEW(K1)=YOLD(NOLD)                                               
      GO TO 600                                                         
C                                                                       
 450  IF(K1.EQ.1)   GO TO 500                                           
      IF(K.EQ.KOLD) GO TO 550                                           
C                                                                       
 500  K=KOLD                                                            
C                                                                       
      Y2K=Y2(K)                                                         
      Y2KP1=Y2(K+1)                                                     
      DX=XOLD(K+1)-XOLD(K)                                              
      RDX=1./DX                                                         
C                                                                       
      AK=.1666667*RDX*(Y2KP1-Y2K)                                       
      BK=.5*Y2K                                                         
      CK=RDX*(YOLD(K+1)-YOLD(K))-.1666667*DX*(Y2KP1+Y2K+Y2K)            
C                                                                       
 550  X=XK-XOLD(K)                                                      
      XSQ=X*X                                                           
C                                                                       
      YNEW(K1)=AK*XSQ*X+BK*XSQ+CK*X+YOLD(K)                             
C                                                                       
 600  K1=K1+1                                                           
      IF(K1.LE.NNEW) GO TO 300       
C-----------------------------------------------------------------------
      RETURN                                                            
      END
