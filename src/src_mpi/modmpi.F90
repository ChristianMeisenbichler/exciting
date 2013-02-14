! Copyright (C) 2006-2008 C. Ambrosch-Draxl. C. Meisenbichler S. Sagmeister
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!
! !MODULE:  modmpi
! !DESCRIPTION:
!   MPI variables and interface functions
!   In case of compiled without MPI support it defines the
!   mpi specific variables such that the code behaves exactly as
!   the unmodified scalar version
!
!
! !REVISION HISTORY:
!   Created October 2006 (CHM)
!   Added wrapper routines, 2007-2008 (S. Sagmeister)
!   Added allgatherv interface, August 2010 (S. Sagmeister)
!
!
!
Module modmpi
#ifdef MPI
use mpi
#endif
implicit none

      Type MPIinfo

        Integer :: rank, myprocrow, myproccol
        Integer :: procs, nprocrows, nproccols
        Integer :: blocksize = 64
        Integer :: ierr
        Integer :: comm, context

      End Type MPIinfo

      Logical :: splittfile
      Type(MPIinfo) :: MPIglobal
!
!!$  character(256)::filename
Contains
      Subroutine initMPI
#ifdef MPI
    !        mpi init
         Call mpi_init (MPIglobal%ierr)
         MPIglobal%comm = mpi_comm_world
         Call mpi_comm_size (mpi_comm_world, MPIglobal%procs, MPIglobal%ierr)
         Call mpi_comm_rank (mpi_comm_world, MPIglobal%rank, MPIglobal%ierr)
         splittfile = .True.
!TODO: decent processor grid initialization
         MPIglobal%nprocrows = 1
         MPIglobal%nproccols = 1
#endif
#ifndef MPI
         MPIglobal%comm  = 0
         MPIglobal%procs = 1
         MPIglobal%rank  = 0
         splittfile = .False.
#endif
      End Subroutine initMPI
!
      Subroutine finitMPI
#ifdef MPI
         Call MPI_Finalize (MPIglobal%ierr)
#endif
      End Subroutine finitMPI
!
! service functions to partition k points
! still kept around but should be replaced by generig functions
! firstofset nofset lastofset ....
      Function nofk (process)
         Use modmain, Only: nkpt
         Integer :: nofk
         Integer, Intent (In) :: process
         nofk = nkpt / MPIglobal%procs
         If ((Mod(nkpt, MPIglobal%procs) .Gt. process)) nofk = nofk + 1
      End Function nofk
!
      Function firstk (process)
         Integer :: firstk
         Integer, Intent (In) :: process
         integer::i
         firstk = 1
         Do i = 0, process - 1
            firstk = firstk + nofk (i)
         End Do
      End Function firstk
!
      Function lastk (process)
         Integer :: lastk,i
         Integer, Intent (In) :: process
         lastk = 0
         Do i = 0, process
            lastk = lastk + nofk (i)
         End Do
      End Function lastk
!
      Function procofk (k)
         Integer :: procofk
         Integer, Intent (In) :: k
         Integer :: iproc
         procofk = 0
         Do iproc = 0, MPIglobal%procs - 1
            If (k .Gt. lastk(iproc)) procofk = procofk + 1
         End Do
      End Function procofk
!
!-------------generalized partition!
      Function nofset (process, set)
         Integer :: nofset
         Integer, Intent (In) :: process, set
         nofset = set / MPIglobal%procs
         If ((Mod(set, MPIglobal%procs) .Gt. process)) nofset = nofset + 1
      End Function nofset
!
      Function firstofset (process, set)
         Integer :: firstofset
         Integer, Intent (In) :: process, set
         integer::i
         firstofset = 1
         Do i = 0, process - 1
            firstofset = firstofset + nofset (i, set)
         End Do
      End Function firstofset
!
      Function lastofset (process, set)
         Integer :: lastofset
         Integer, Intent (In) :: process, set
         integer::i
         lastofset = 0
         Do i = 0, process
            lastofset = lastofset + nofset (i, set)
         End Do
      End Function lastofset
!
      Function procofindex (k, set)
         Integer :: procofindex
         Integer, Intent (In) :: k, set
         Integer :: iproc
         procofindex = 0
         Do iproc = 0, MPIglobal%procs - 1
            If (k .Gt. lastofset(iproc, set)) procofindex = procofindex &
           & + 1
         End Do
      End Function procofindex
!
      Function lastproc (row, set)
         Implicit None
         Integer :: lastproc
         Integer, Intent (In) :: row, set
         If (row .Ne. nofset(0, set)) Then
            lastproc = MPIglobal%procs
         Else
            lastproc = modulo (set, MPIglobal%procs)
            If (lastproc .Eq. 0) lastproc = MPIglobal%procs
         End If
         lastproc = lastproc - 1
      End Function lastproc
!
!
!
!------------------interface to MPI_barrier for xs-part
      Subroutine barrier
         Implicit None
  ! do nothing if only one process
#ifndef MPI
         If (MPIglobal%procs .Eq. 1) Return
#endif
  ! call the MPI barrier
#ifdef MPI
         Call MPI_barrier (mpi_comm_world, MPIglobal%ierr)
#endif
      End Subroutine barrier


      Subroutine endloopbarrier (set, mult)
         Implicit None
         Integer, Intent (In) :: set, mult
         Integer :: i
         Do i = 1, (nofset(0, set)-nofset(MPIglobal%rank, set)) * mult
            Call barrier
         End Do
      End Subroutine endloopbarrier
!
!------------------wrappers for MPI communication

      subroutine mpi_allgatherv_ifc(set,rlen,ibuf,rbuf,zbuf)
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


      ! taken from SIESTA
      subroutine set_processorYdefault(Nodes,procYdefault)

      ! Finds a sensible default value for the processor grid in the Y
      ! direction to try to get an equal split in X and Y. Currently only
      ! looks for factors of 2, 3 and 5.
      !
      ! Input :
      !
      ! integer Nodes        : total number of processors
      !
      ! Output :
      !
      ! integer procYdefault : default value of Y grid parameter
      !
      ! Written by Julian Gale, November 1999
      !
      ! Passed arguments
            integer :: Nodes, procYdefault
      ! Local variables
            integer :: Nx, Ny, Nrem
            logical :: factor

      ! Initialise values
            Nx = 1
            Ny = 1
            Nrem = Nodes
            factor = .true.

      ! Loop looking for factors
            do while (factor.and.Nrem.gt.1)
              factor = .false.
              if (mod(Nrem,2).eq.0) then
                Nrem = Nrem/2
                factor = .true.
                if (Nx.gt.Ny) then
                  Ny = 2*Ny
                else
                  Nx = 2*Nx
                endif
              endif
              if (mod(Nrem,3).eq.0) then
                Nrem = Nrem/3
                factor = .true.
                if (Nx.gt.Ny) then
                  Ny = 3*Ny
                else
                  Nx = 3*Nx
                endif
              endif
              if (mod(Nrem,5).eq.0) then
                Nrem = Nrem/5
                factor = .true.
                if (Nx.gt.Ny) then
                  Ny = 5*Ny
                else
                  Nx = 5*Nx
                endif
              endif
            enddo

      ! Choose default value as lowest of Nx and Ny
            if (Nx.gt.Ny) then
              procYdefault = Ny
            else
              procYdefault = Nx
            endif

      end subroutine set_processorYdefault


End Module modmpi
