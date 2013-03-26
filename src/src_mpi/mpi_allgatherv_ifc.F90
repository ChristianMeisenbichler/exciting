 module mpi_allgatherv_ifc_module
 #ifdef MPI
 use mpi
 #endif
 contains
 subroutine mpi_allgatherv_ifc(set,rlen,ibuf,rbuf,zbuf)
        use modmpi , only: MPIglobal
        use setpatrtitioning, only:firstofset,nofset
        implicit none
        integer, intent(in) :: set,rlen
        integer, intent(inout), optional :: ibuf(*)
        real(8), intent(inout), optional :: rbuf(*)
        complex(8), intent(inout), optional :: zbuf(*)
        ! local variables
#ifndef MPI1
#define BUFFER mpi_in_place
#endif
#ifdef MPI1
        complex(8), allocatable :: bufz(:)
        real(8), allocatable :: bufr(:)
        integer, allocatable :: bufi(:)
#endif
        integer, allocatable :: buf_n(:),buf_dspls(:)
        integer :: j
        logical :: ti,tr,tz
        ti=present(ibuf); tr=present(rbuf); tz=present(zbuf)
        if (count((/ti,tr,tz/)).ne.1) then
          write(*,*)
          write(*,'("Error(mpi_allgatherv_ifc): exactly one array must be defined.")')
          write(*,*)
          stop
        end if
#ifdef MPI
        allocate(buf_n(MPIglobal%procs),buf_dspls(MPIglobal%procs))
        ! displacements within receive buffer (flattened array)
        buf_dspls=(/(rlen*(firstofset(j,set)-1),j=0,MPIglobal%procs-1)/)
        ! number of elements in send buffer (flattened array)
        buf_n=(/(rlen*nofset(j,set),j=0,MPIglobal%procs-1)/)
        ! use recieve buffer as sendbuffer by specifying mpi_in_place
        if (ti) then
#ifdef MPI1
#define BUFFER bufi
          allocate(bufi(buf_n(MPIglobal%rank+1)))
          bufi(:)=ibuf(buf_dspls(MPIglobal%rank+1)+1:buf_dspls(MPIglobal%rank+1)+buf_n(MPIglobal%rank+1))
#endif
          call mpi_allgatherv(BUFFER, &
            buf_n(MPIglobal%rank+1), &
            mpi_integer, &
            ibuf, &
            buf_n, &
            buf_dspls, &
            mpi_integer, &
            mpi_comm_world, &
            MPIglobal%ierr)
#ifdef MPI1
          deallocate(bufi)
#undef BUFFER
#endif
        end if
        if (tr) then
#ifdef MPI1
#define BUFFER bufr
          allocate(bufr(buf_n(MPIglobal%rank+1)))
          bufr(:)=rbuf(buf_dspls(MPIglobal%rank+1)+1:buf_dspls(MPIglobal%rank+1)+buf_n(MPIglobal%rank+1))
#endif
          call mpi_allgatherv(BUFFER, &
            buf_n(MPIglobal%rank+1), &
            mpi_double_precision, &
            rbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_precision, &
            mpi_comm_world, &
            MPIglobal%ierr)
#ifdef MPI1
          deallocate(bufr)
#undef BUFFER
#endif
        end if
        if (tz) then
#ifdef MPI1
#define BUFFER bufz
          allocate(bufz(buf_n(MPIglobal%rank+1)))
          bufz(:)=zbuf(buf_dspls(MPIglobal%rank+1)+1:buf_dspls(MPIglobal%rank+1)+buf_n(MPIglobal%rank+1))
#endif
          call mpi_allgatherv(BUFFER, &
            buf_n(MPIglobal%rank+1), &
            mpi_double_complex, &
            zbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_complex, &
            mpi_comm_world, &
            MPIglobal%ierr)
#ifdef MPI1
          deallocate(bufz)
#undef BUFFER
#endif
        end if
        deallocate(buf_n,buf_dspls)
#endif
      end subroutine
end module
