      subroutine mpi_isend(buf,count,datatype,source,
     & tag,comm,request,ierror)
      integer buf(*), count,datatype,source,tag,comm,
     & request,ierror
      call mpi_error()
      return
      end  

      subroutine mpi_irecv(buf,count,datatype,source,
     & tag,comm,request,ierror)
      integer buf(*), count,datatype,source,tag,comm,
     & request,ierror
      call mpi_error()
      return
      end

      subroutine mpi_send(buf,count,datatype,dest,tag,comm,ierror)
      integer buf(*), count,datatype,dest,tag,comm,ierror
      call mpi_error()
      return
      end
      
      subroutine mpi_recv(buf,count,datatype,source,
     & tag,comm,status,ierror)
      integer buf(*), count,datatype,source,tag,comm,
     & status(*),ierror
      call mpi_error()
      return
      end

      subroutine mpi_comm_split(comm,color,key,newcomm,ierror)
      integer comm,color,key,newcomm,ierror
      return
      end

      subroutine mpi_comm_rank(comm, rank,ierr)
      implicit none
      integer comm, rank,ierr
      rank = 0
      return
      end

      subroutine mpi_comm_size(comm, size, ierr)
      implicit none
      integer comm, size, ierr
      size = 1
      return
      end

Cmp	additions for internal quilt Eta model
	
	subroutine mpi_group_excl(intin,n,ranks,ngroup,ierr)
	integer intin,n,ranks(n),ngroup,ierr
	ngroup=intin
	return
	end

	subroutine mpi_group_free (intin,ierr )
	integer intin,ierr
	return
	end

	subroutine mpi_intercomm_create (lcom,llead,ipeer,irem,itag,
     &			newcom,ierr )
	integer lcom,llead,ipeer,irem,itag,newcom,ierr
	newcom=lcom
	write(6,*) 'shouldnt be calling this!'
	return
	end

        subroutine mpi_comm_create (com,group,ncom,ierr ) 
        integer com,group,ncom,ierr
        ncom=com
        return
        end

        subroutine mpi_comm_group (com,group,ierr )
        include 'mpif.h'
        integer com,group,ierr
        write(6,*) 'dont want to see this'
        group=com
        return
        end


Cmp	end additions


      double precision function mpi_wtime()
      implicit none
      double precision t
c This function must measure wall clock time, not CPU time. 
c Since there is no portable timer in Fortran (77)
c we call a routine compiled in C (though the C source may have
c to be tweaked). 
      call wtime(t)
c The following is not ok for "official" results because it reports
c CPU time not wall clock time. It may be useful for developing/testing
c on timeshared Crays, though. 
c     call second(t)

      mpi_wtime = t

      return
      end


c may be valid to call this in single processor case
      subroutine mpi_barrier(comm,ierror)
      return
      end

c may be valid to call this in single processor case
      subroutine mpi_bcast(buf, nitems, type, root, comm, ierr)
      implicit none
      integer buf(*), nitems, type, root, comm, ierr
      return
      end

      subroutine mpi_comm_dup(oldcomm, newcomm,ierror)
      integer oldcomm, newcomm,ierror
      newcomm= oldcomm
      return
      end

      subroutine mpi_error()
      print *, 'mpi_error called'
      stop
      end 

      subroutine mpi_abort(comm, errcode, ierr)
      implicit none
      integer comm, errcode, ierr
      print *, 'mpi_abort called'
      stop
      end

      subroutine mpi_finalize(ierr)
      return
      end

      subroutine mpi_init(ierr)
      return
      end


c assume double precision, which is all SP uses 
      subroutine mpi_reduce(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      implicit none
      include 'mpif.h'
      integer nitems, type, op, root, comm, ierr
      double precision inbuf(*), outbuf(*)

      if (type .eq. mpi_real8) then
         call mpi_reduce_dp(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      else if (type .eq.  mpi_double_complex) then
         call mpi_reduce_dc(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      else if (type .eq.  mpi_complex) then
         call mpi_reduce_complex(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      else if (type .eq.  mpi_real) then
         call mpi_reduce_real(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      else if (type .eq.  mpi_integer) then
         call mpi_reduce_int(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      else 
         print *, 'mpi_reduce: unknown type ', type
      end if
      return
      end


      subroutine mpi_reduce_real(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      implicit none
      integer nitems, type, op, root, comm, ierr, i
      real inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end

      subroutine mpi_reduce_dp(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      implicit none
      integer nitems, type, op, root, comm, ierr, i
      double precision inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end

      subroutine mpi_reduce_dc(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      implicit none
      integer nitems, type, op, root, comm, ierr, i
      double complex inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end


      subroutine mpi_reduce_complex(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      implicit none
      integer nitems, type, op, root, comm, ierr, i
      complex inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end

      subroutine mpi_reduce_int(inbuf, outbuf, nitems, 
     $                      type, op, root, comm, ierr)
      implicit none
      integer nitems, type, op, root, comm, ierr, i
      integer inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end

      subroutine mpi_allreduce(inbuf, outbuf, nitems, 
     $                      type, op, comm, ierr)
      implicit none
      integer nitems, type, op, comm, ierr
      double precision inbuf(*), outbuf(*)

      call mpi_reduce(inbuf, outbuf, nitems, 
     $                      type, op, 0, comm, ierr)
      return
      end

      subroutine mpi_alltoall(inbuf, nitems, type, outbuf, nitems_dum, 
     $                        type_dum, comm, ierr)
      implicit none
      include 'mpif.h'
      integer nitems, type, comm, ierr, nitems_dum, type_dum
      double precision inbuf(*), outbuf(*)
      if (type .eq. mpi_real8) then
         call mpi_alltoall_dp(inbuf, outbuf, nitems, 
     $                      type, comm, ierr)
      else if (type .eq.  mpi_double_complex) then
         call mpi_alltoall_dc(inbuf, outbuf, nitems, 
     $                      type, comm, ierr)
      else if (type .eq.  mpi_complex) then
         call mpi_alltoall_complex(inbuf, outbuf, nitems, 
     $                      type, comm, ierr)
      else if (type .eq.  mpi_real) then
         call mpi_alltoall_real(inbuf, outbuf, nitems, 
     $                      type, comm, ierr)
      else if (type .eq.  mpi_integer) then
         call mpi_alltoall_int(inbuf, outbuf, nitems, 
     $                      type, comm, ierr)
      else 
         print *, 'mpi_alltoall: unknown type ', type
      end if
      return
      end

      subroutine mpi_alltoall_dc(inbuf, outbuf, nitems, 
     $                           type, comm, ierr)
      implicit none
      integer nitems, type, comm, ierr, i
      double complex inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end


      subroutine mpi_alltoall_complex(inbuf, outbuf, nitems, 
     $                           type, comm, ierr)
      implicit none
      integer nitems, type, comm, ierr, i
      double complex inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end

      subroutine mpi_alltoall_dp(inbuf, outbuf, nitems, 
     $                           type, comm, ierr)
      implicit none
      integer nitems, type, comm, ierr, i
      double precision inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end

      subroutine mpi_alltoall_real(inbuf, outbuf, nitems, 
     $                             type, comm, ierr)
      implicit none
      integer nitems, type, comm, ierr, i
      real inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end

      subroutine mpi_alltoall_int(inbuf, outbuf, nitems, 
     $                            type, comm, ierr)
      implicit none
      integer nitems, type, comm, ierr, i
      integer inbuf(*), outbuf(*)
      do i = 1, nitems
         outbuf(i) = inbuf(i)
      end do
      
      return
      end

      subroutine mpi_wait(request,status,ierror)
      integer request,status,ierror
Cmp      call mpi_error()
	write(6,*) 'not quitting with mpi_wait'
      return
      end

      subroutine mpi_waitall(count,requests,status,ierror)
      integer count,requests(*),status(*),ierror
      call mpi_error()
      return
      end

CC      three more added 2000/11/14 to deal with MPI post


        subroutine mpi_gatherv(sendbuf,sendcount,sendtype,
     +  recvbuf,recvcounts,displs,recvtype,root,comm,ierr)

        integer sendcount,sendtype,recvcounts(*),displs(*)
        integer recvtype,root,comm,ierr
        write(6,*) 'should never execute this when running with one CPU'

        return
        end


        subroutine mpi_sendrecv(sendbuf,sendcount,sendtype,
     +  dest,sendtag,recvbuf,recvcount,recvtype,source,recvtag,
     +  comm,status,ierr)

        integer sendcount,sendtype
        integer dest,sendtag,recvcount
        integer recvtype,source,recvtag,comm,ierr,status(*)

        write(6,*) 'should never execute this when running with one CPU'

        return
        end

        subroutine mpi_scatterv(sendbuf,sendcounts,displs,
     +  sendtype,recvbuf,recvcount,recvtype,root,comm,ierr)

        integer sendcounts(*),sendtype,recvcount,displs(*)
        integer recvtype,root,comm,ierr

        write(6,*) 'should never execute this when running with one CPU'

        return
        end

