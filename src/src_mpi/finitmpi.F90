



      Subroutine finitMPI
      use modmpi_types
#ifdef MPI
         Call BLACS_GRIDEXIT(MPIkpt_2D%context)
         Call BLACS_GRIDEXIT(MPIkpt_1D%context)
         Call MPI_Finalize(MPIglobal%ierr)
#endif
      End Subroutine finitMPI

!
