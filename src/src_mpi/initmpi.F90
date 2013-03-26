  Subroutine initMPI
#ifdef MPI
    Use modmpi_types
    !        mpi init
         Call mpi_init(MPIglobal%ierr)

         MPIglobal%comm = mpi_comm_world
         Call mpi_comm_size (mpi_comm_world, MPIglobal%procs, MPIglobal%ierr)
         Call mpi_comm_rank (mpi_comm_world, MPIglobal%rank,  MPIglobal%ierr)

         splittfile = .True.

!TODO: decent processor grid initialization based on splitting procs over k-points
         MPIkpt_2D%comm  = mpi_comm_world
         MPIkpt_2D%procs = MPIglobal%procs
         MPIkpt_2D%rank  = MPIglobal%rank
         Call set_processorYdefault(MPIkpt_2D%procs, MPIkpt_2D%nprocrows)
         MPIkpt_2D%nproccols = MPIkpt_2D%procs/MPIkpt_2D%nprocrows
         Call BLACS_GET(0, 0, MPIkpt_2D%context)
         Call BLACS_GRIDINIT(MPIkpt_2D%context,'r', MPIkpt_2D%nprocrows, MPIkpt_2D%nproccols)
         Call BLACS_GRIDINFO(MPIkpt_2D%context, MPIkpt_2D%nprocrows, MPIkpt_2D%nproccols, MPIkpt_2D%myprocrow, MPIkpt_2D%myproccol)

         MPIkpt_1D%comm  = mpi_comm_world
         MPIkpt_1D%procs = MPIglobal%procs
         MPIkpt_1D%rank  = MPIglobal%rank
         MPIkpt_1D%nprocrows = 1
         MPIkpt_1D%nproccols = MPIglobal%procs
         Call BLACS_GET(0, 0, MPIkpt_1D%context)
         Call BLACS_GRIDINIT(MPIkpt_1D%context,'r', MPIkpt_1D%nprocrows, MPIkpt_1D%nproccols)
         Call BLACS_GRIDINFO(MPIkpt_1D%context, MPIkpt_1D%nprocrows, MPIkpt_1D%nproccols, MPIkpt_1D%myprocrow, MPIkpt_1D%myproccol)
#else
         MPIglobal%comm  = 0
         MPIglobal%procs = 1
         MPIglobal%rank  = 0
         splittfile = .False.
#endif
      End Subroutine initMPI
!
