      subroutine conh12t(h1,imjm,lm,h2,imt,jmt)                       
c
c *** Routine for reordering thibue-height-point-1-dimensional      
c        matrices for 2-dimensional indexing (imt,jmt).
c
      implicit none
c
      integer*4 imt,jmt,imjm,lm,i,j,k,l
c
      real*4 h1(imjm,lm),h2(imt,jmt,lm)                               
c_______________________________________________________________________________
c
      do l=1,lm                              
         k=0                                        
         do j=1,jmt                                            
            do i=1,imt,2                                              
               k=k+1          
               h2(i,j,l)=h1(k,l)    
            enddo
c                                                                       
            if (j .lt. jmt) then
               do i=2,imt,2                                              
                  k=k+1     
                  h2(i,j,l)=h1(k,l)  
               enddo
            endif 
         enddo
      enddo
c
      return                                     
      end                                        
c
c===============================================================================
      subroutine conh12t_sfc(h1,imjm,h2,imt,jmt)                       
c
c *** Routine for reordering thibue-height-point-1-dimensional      
c        matrices for 2-dimensional indexing (imt,jmt).
c
      implicit none
c
      integer*4 imt,jmt,imjm,i,j,k
c
      real*4 h1(imjm),h2(imt,jmt)                               
c_______________________________________________________________________________
c
         k=0                                        
        do j=1,jmt                                            

            do i=1,imt,2                                              
               k=k+1          
               h2(i,j)=h1(k)    
            enddo
c                                                                       
            if (j .lt. jmt) then
               do i=2,imt,2                                              
                  k=k+1     
                  h2(i,j)=h1(k)  
               enddo
            endif 

        enddo
c
      return                                     
      end                                        
c
c===============================================================================
c
      subroutine conv12t(v1,imjm,lm,v2,imt,jmt)                      
c
c *** Routine for reordering thibue-velocity-point-1-dimensional 
c        matrices for 2-dimensional indexing (imt,jmt).
c
      implicit none
c
      integer*4 imt,jmt,imjm,lm,i,j,k,l
c
      real*4 v1(imjm,lm),v2(imt,jmt,lm)
c_______________________________________________________________________________
c
      do l=1,lm                              
         k=0                                        
         do j=1,jmt                                            
            do i=2,imt,2                                              
               k=k+1        
               v2(i,j,l)=v1(k,l) 
            enddo
c                                                                       
            if (j .lt. jmt) then
c                                                                       
               do i=1,imt,2                                              
                  k=k+1 
                  v2(i,j,l)=v1(k,l)  
               enddo
            endif
c                                                                       
         enddo
      enddo
c
      return                                     
      end                                        
c
c===============================================================================
c
      subroutine conh21t(h2,imt,jmt,lm,h1,imjm)
c
c *** Routine for reordering thibue-height-point-two-dimensional  
c     matrices (imt,jmt) for one-dimenional indexing.                       
c
      implicit none
c
      integer*4 imt,jmt,imjm,lm,i,j,k,l
c
      real*4 h2(imt,jmt,lm),h1(imjm,lm)
c_______________________________________________________________________________
c
      do l=1,lm
         k=0
         do j=1,jmt
            do i=1,imt,2
               k=k+1
               h1(k,l)=h2(i,j,l)
            enddo
c
            if (j .lt. jmt) then
               do i=2,imt,2
                  k=k+1
                  h1(k,l)=h2(i,j,l) 
               enddo
            endif
         enddo
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine conv21t(v2,imt,jmt,lm,v1,imjm)
c
c *** Routine for reordering thibue-height-point-two-dimensional  
c     matrices (imt,jmt) for one-dimenional indexing.                       
c
      implicit none
c
      integer*4 imt,jmt,imjm,lm,i,j,k,l
c
      real*4 v2(imt,jmt,lm),v1(imjm,lm)
c_______________________________________________________________________________
c
      do l=1,lm
         k=0
         do j=1,jmt
            do i=2,imt,2
               k=k+1
               v1(k,l)=v2(i,j,l)
            enddo
            if (j .ne. jmt) then
               do i=1,imt,2
                  k=k+1
                  v1(k,l)=v2(i,j,l)
               enddo
            endif
         enddo
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine conh12(h1,imjm,lm,h2,im,jm)                       
c
c *** Routine for reordering thibue-height-point-1-dimensional
c        matrices (im,jm) for 2-dimensional indexing.
c
      implicit none
c
      integer*4 im,jm,imjm,lm,i,j,k,l,imm1
c
      real*4 h1(imjm,lm),h2(im,jm,lm)
c_______________________________________________________________________________
c
      imm1=im-1
      do l=1,lm
         k=0                                        
         do j=1,jm                                            
         do i=1,imm1+mod(j,2)
            k=k+1                                                             
            h2(i,j,l)=h1(k,l)                                                 
         enddo
         enddo
         do j=2,jm-1,2
            h2(im,j,l)=h2(imm1,j,l)
         enddo
      enddo
c
      return                                     
      end                                        
c
c===============================================================================

	subroutine conh12_TEST(h1,imjm,lm,h2,im,jm)
	real h2(im,jm,lm),h1(imjm,lm)

	do l=1,lm

      DO J=1,JM
      DO I=1,IM
        H2(I,J,l)=-999.
      ENDDO
      ENDDO

      IMT=2*IM-1
      DO K=1,IMJM
        JX=(K-1)/IMT+1
        IX=K-(JX-1)*IMT
        IF(IX.LE.IM)THEN
          I=IX
          J=2*JX-1
        ELSE
          I=IX-IM
          J=2*JX
        ENDIF
        H2(I,J,l)=H1(K,l)
        IF(IX.EQ.IMT)H2(I+1,J,l)=H2(I,J,l)
      ENDDO
	enddo
C
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
c
      subroutine conv12(v1,imjm,lm,v2,im,jm)                      
c
c *** Routine for reordering thibue-velocity-point-1-dimensional 
c        matrices (im,jm) for 2-dimensional indexing.
c
      implicit none
c
      integer*4 im,jm,imjm,lm,i,j,k,l,imm1
c
      real*4 v1(imjm,lm),v2(im,jm,lm)
c_______________________________________________________________________________
c
      imm1=im-1
      do l=1,lm
         k=0                                        
         do j=1,jm                                            
         do i=1,im-mod(j,2)
            k=k+1                                                             
            v2(i,j,l)=v1(k,l)                                                
         enddo
         enddo
         do j=1,jm,2
            v2(im,j,l)=v2(imm1,j,l)
         enddo
      enddo
c
      return                                     
      end                                        
c
c===============================================================================
c
      subroutine conh21(h2,im,jm,lm,h1,imjm)
c
c *** Routine for reordering thibue-height-point-two-dimensional  
c     matrices (im,jm) for one-dimenional indexing.                       
c
      implicit none
c
      integer*4 im,jm,imjm,lm,i,j,k,l
c
      real*4 h2(im,jm,lm),h1(imjm,lm)
c_______________________________________________________________________________
c
      
      do l=1,lm
         k=0
         do j=1,jm
         do i=1,im
            if (.not. (mod(j,2) .eq. 0 .and. i .eq. im)) then
               k=k+1
               h1(k,l)=h2(i,j,l)
            endif
         enddo
         enddo
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine conv21(v2,im,jm,lm,v1,imjm)
c
c *** Routine for reordering thibue-height-point-two-dimensional  
c     matrices (im,jm) for one-dimenional indexing.                       
c
      implicit none
c
      integer*4 im,jm,imjm,lm,i,j,k,l
c
      real*4 v2(im,jm,lm),v1(imjm,lm)
c_______________________________________________________________________________
c
      do l=1,lm
         k=0
         do j=1,jm
         do i=1,im
            if (.not. (mod(j,2) .eq. 1 .and. i .eq. im)) then
               k=k+1
               v1(k,l)=v2(i,j,l)
            endif
         enddo
         enddo
      enddo
c
      return
      end
c
c===============================================================================
c
      subroutine iconh12(h1,imjm,lm,h2,im,jm)                       
c
c *** Routine for reordering thibue-height-point-1-dimensional
c        matrices (im,jm) for 2-dimensional indexing.
c
      implicit none
c
      integer*4 im,jm,imjm,lm,i,j,k,l,imm1
c
      integer*4 h1(imjm,lm),h2(im,jm,lm)
c_______________________________________________________________________________
c
      imm1=im-1
      do l=1,lm
         k=0                                        
         do j=1,jm                                            
         do i=1,imm1+mod(j,2)
            k=k+1                                                             
            h2(i,j,l)=h1(k,l)                                                 
         enddo
         enddo
         do j=2,jm-1,2
            h2(im,j,l)=h2(imm1,j,l)
         enddo
      enddo
c
      return                                     
      end                                        
c
c===============================================================================
c
      subroutine iconv12(v1,imjm,lm,v2,im,jm)                      
c
c *** Routine for reordering thibue-velocity-point-1-dimensional 
c        matrices (im,jm) for 2-dimensional indexing.
c
      implicit none
c
      integer*4 im,jm,imjm,lm,i,j,k,l,imm1
c
      integer*4 v1(imjm,lm),v2(im,jm,lm)
c_______________________________________________________________________________
c
      imm1=im-1
      do l=1,lm
         k=0                                        
         do j=1,jm                                            
         do i=1,im-mod(j,2)
            k=k+1                                                             
            v2(i,j,l)=v1(k,l)                                                
         enddo
         enddo
         do j=1,jm,2
            v2(im,j,l)=v2(imm1,j,l)
         enddo
      enddo
c
      return                                     
      end                                        
