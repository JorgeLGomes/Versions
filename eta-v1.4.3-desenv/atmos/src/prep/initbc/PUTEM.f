C
       subroutine putem(hlat,hlon,i1d,sm,sice,ieta)
C
       INCLUDE "parmeta"
C
       real*4 hlat(im,jm),hlon(im,jm)
	integer*4 i1d(360,180), ieta(im,jm)
       real*4 sm(im,jm),sice(im,jm)
C
C   This subroutine moves the value stored in i1d (integer 1 by 1 deg global
C   grid) to the nearest eta land point in the ieta integer array.
C   Special effort is made to ensure that islands are taken care of when
C   not resolved on the 1 x 1 deg. grid.
c
        do j = 1,jm
        do i = 1,im

          if((sm(i,j).gt.0.5).or.(sice(i,j).eq.1.0)) then
            ieta(i,j) = 0
          else
            indx = nint(hlon(i,j) - 179.5)
            if(indx.lt.1) indx = nint(hlon(i,j) + 180.5)
            indy = nint(hlat(i,j) + 90.5)
            if(i1d(indx,indy).ne.0) then
Cmp	     write(6,*) 'assigned a value'
                ieta(i,j) = i1d(indx,indy)
            else
c We have an isolated point and need to look at surrounding
c 1 degree grid points, at this latitude.
                do krad=1,40
		do jj = indy-krad,indy+krad
                do ii=indx-krad,indx+krad
                  itt = min(ii,360)
                  itt = max(itt,1)
                  jtt = min(jj,180)
                  jtt = max(jtt,1)
                  if(i1d(itt,jtt).ne.0) then
                    ieta(i,j) = i1d(itt,jtt)
                    go to 33
                  end if
                end do
                end do
                end do
           print *,"oh oh ",hlat(i,j),hlon(i,j),indx,indy
33              continue
            end if
          end if
        end do
        end do
          return
          end
