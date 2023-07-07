      program ConvInput2ModelData
!
      implicit none
!
!
      character(LEN=255):: gribfile,outdir
!
!_______________________________________________________________________________
!
! *** Initialize tables.
!
      call es_ini
!
! *** Read file names.
!
      read(5,'(a)') gribfile
      read(5,'(a)') outdir
!
      call conv2model(gribfile,outdir)

!
      end program ConvInput2ModelData
