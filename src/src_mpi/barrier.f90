
!------------------interface to MPI_barrier for xs-part
      Subroutine barrier
      use modmpi
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
         use setpatrtitioning
         use modmpi

         Implicit None
         Integer, Intent (In) :: set, mult
         Integer :: i
         Do i = 1, (nofset(0, set)-nofset(MPIglobal%rank, set)) * mult
            Call barrier
         End Do
      End Subroutine endloopbarrier
!

