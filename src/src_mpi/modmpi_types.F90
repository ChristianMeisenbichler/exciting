
Module modmpi_types
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

! not needed, so not tested and not implemented
!       Integer, Parameter :: DISTRIBUTE_ROWS = 1
      Integer, Parameter :: DISTRIBUTE_COLS = 2
      Integer, Parameter :: DISTRIBUTE_2D   = 3

      Logical       :: splittfile
      Type(MPIinfo) :: MPIglobal
      Type(MPIinfo) :: MPIkpt_1D, MPIkpt_2D



End Module modmpi_types
