! !MODULE:  test_helpers
! !DESCRIPTION:
! General functions and tools for the unittests
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC)
!
module test_helpers
#ifdef MPI
    use modmpi
#endif

    implicit none

    real(8), parameter :: tol  = 1E-12
    real(8), parameter :: zero = 0.0D0

    contains

#ifdef MPI
    subroutine setupProcGrid(nprocrows_sub, nproccols_sub, comm_sub, context_sub, ierror)
      implicit none

      integer, intent(in)  :: nprocrows_sub, nproccols_sub
      integer, intent(out) :: comm_sub, context_sub, ierror

      integer              :: nprocs, i, world_group, procs_group
      integer, dimension(:), allocatable :: proclist

      nprocs = nprocrows_sub*nproccols_sub
      allocate(proclist(nprocs))
      proclist = (/(i,i=0,nprocs-1)/)
      Call MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierror)
      Call MPI_GROUP_INCL(world_group, nprocs, proclist, procs_group, ierror)
      Call MPI_COMM_CREATE(MPI_COMM_WORLD, procs_group, comm_sub, ierror)

      Call BLACS_GET(0, 0, context_sub)
      Call BLACS_GRIDINIT(context_sub,'r', nprocrows_sub, nproccols_sub)

      Call MPI_GROUP_FREE(world_group, ierror) ! this just frees the handle, which is not needed any further

    end subroutine setupProcGrid
#endif


#ifdef MPI
    subroutine finalizeProcGrid(comm, context, ierror)
      implicit none

      integer, intent(in)   :: comm, context, ierror

      Call BLACS_GRIDEXIT(context)

      Call MPI_COMM_FREE(comm, ierror)

    end subroutine finalizeProcGrid
#endif


#ifdef MPI
    subroutine getBlacsGridInfo(MPIdata)
      implicit none
      Type(MPIinfo), intent(inout) :: MPIdata

      Call BLACS_GRIDINFO(MPIdata%context, MPIdata%nprocrows, MPIdata%nproccols, MPIdata%myprocrow, MPIdata%myproccol)
      MPIdata%procs = MPIdata%nprocrows * MPIdata%nproccols

    end subroutine getBlacsGridInfo
#endif


    subroutine fillIdentity(mat, size)
      implicit none

      complex(8), dimension(:,:), intent(out) :: mat
      integer, intent(in)                     :: size

      integer :: row
      mat = 0
      do row=1,size
        mat(row, row) = cmplx(1.0D0, 0.0D0,8)
      end do

    end subroutine fillIdentity


    subroutine fillHermitian(mat, size)
      implicit none

      complex(8), dimension(:,:), intent(out) :: mat
      integer, intent(in)                     :: size

      integer :: row, col

      do row=1,size
        do col=1,row  
          if (row .ne. col) then 
            mat(row, col) = cmplx(row, col,8)
            mat(col, row) = cmplx(row,-col,8)
          else
            mat(col, row) = cmplx(row, 0.0D0,8)
          end if
        end do
      end do
    end subroutine fillHermitian


    subroutine fillComplex(mat, nrows, ncols)
      implicit none

      complex(8), dimension(:,:), intent(out) :: mat
      integer, intent(in)                     :: nrows, ncols

      integer :: row, col

      do row=1,nrows
        do col=1,ncols  
          mat(row, col) = cmplx(row,col+10,8)
        end do
      end do
    end subroutine fillComplex
    

    ! extract the local data for the specified process in a grid
    ! no communication!! the global matrix has to be present locally!!
#ifdef MPI
    subroutine getBlockDistributedLoc(global, local, MPIdata)
      implicit none

      Integer, External   :: NUMROC

      ! arguments
      complex(8), dimension(:,:), intent(in)  :: global
      complex(8), dimension(:,:), intent(out) :: local
      Type(MPIinfo), intent(in)               :: MPIdata

      ! local variables
      integer :: n_procs_x, n_procs_y, proc_x, proc_y, blocksize
      integer :: n_rows_glob, n_cols_glob, n_rows_loc, n_cols_loc, &
                 n_blocks_x_loc, i_block_x, glob_x_start, glob_x_end, loc_x_start, loc_x_end, &
                 n_blocks_y_loc, i_block_y, glob_y_start, glob_y_end, loc_y_start, loc_y_end

      n_procs_x = MPIdata%nproccols
      n_procs_y = MPIdata%nprocrows
      proc_x    = MPIdata%myproccol
      proc_y    = MPIdata%myprocrow
      blocksize = MPIdata%blocksize

      n_rows_glob = SIZE(global, 1)
      n_cols_glob = SIZE(global, 2)

      n_rows_loc = NUMROC(n_rows_glob, blocksize, proc_y, 0, n_procs_y)
      n_cols_loc = NUMROC(n_cols_glob, blocksize, proc_x, 0, n_procs_x)

      n_blocks_x_loc = CEILING(FLOAT(n_cols_loc)/blocksize)
      n_blocks_y_loc = CEILING(FLOAT(n_rows_loc)/blocksize)

      do i_block_x=0,n_blocks_x_loc-1
        glob_x_start = (proc_x + i_block_x*n_procs_x)*blocksize + 1
        glob_x_end   = MIN(glob_x_start + blocksize - 1, n_cols_glob)
        loc_x_start  = i_block_x*blocksize + 1
        loc_x_end    = MIN(loc_x_start + blocksize - 1, n_cols_loc)

        do i_block_y=0,n_blocks_y_loc-1
          glob_y_start = (proc_y + i_block_y*n_procs_y)*blocksize + 1
          glob_y_end   = MIN(glob_y_start + blocksize - 1, n_rows_glob)
          loc_y_start  = i_block_y*blocksize + 1
          loc_y_end    = MIN(loc_y_start + blocksize - 1, n_rows_loc)
          local(loc_y_start:loc_y_end,loc_x_start:loc_x_end) = global(glob_y_start:glob_y_end,glob_x_start:glob_x_end)
        end do
      end do

    end subroutine getBlockDistributedLoc
#endif


#ifdef MPI
    subroutine getLocalIndices(n_rows_glob, n_cols_glob, idx_rows, idx_cols, MPIdata)
      implicit none

      Integer, External   :: NUMROC

      ! arguments
      integer,               intent(in)  :: n_rows_glob, n_cols_glob
      integer, dimension(:), intent(out) :: idx_rows, idx_cols
      Type(MPIinfo),         intent(in)  :: MPIdata

      ! local variables
      integer :: n_procs_x, n_procs_y, proc_x, proc_y, blocksize, i
      integer :: n_rows_loc, n_cols_loc, &
                 n_blocks_x_loc, i_block_x, glob_x_start, glob_x_end, loc_x_start, loc_x_end, &
                 n_blocks_y_loc, i_block_y, glob_y_start, glob_y_end, loc_y_start, loc_y_end

      n_procs_x = MPIdata%nproccols
      n_procs_y = MPIdata%nprocrows
      proc_x    = MPIdata%myproccol
      proc_y    = MPIdata%myprocrow
      blocksize = MPIdata%blocksize

      n_rows_loc = NUMROC(n_rows_glob, blocksize, proc_y, 0, n_procs_y)
      n_cols_loc = NUMROC(n_cols_glob, blocksize, proc_x, 0, n_procs_x)

      n_blocks_x_loc = CEILING(FLOAT(n_cols_loc)/blocksize)
      n_blocks_y_loc = CEILING(FLOAT(n_rows_loc)/blocksize)

      do i_block_x=0,n_blocks_x_loc-1
        glob_x_start = (proc_x + i_block_x*n_procs_x)*blocksize + 1
        glob_x_end   = MIN(glob_x_start + blocksize - 1, n_cols_glob)
        loc_x_start  = i_block_x*blocksize + 1
        loc_x_end    = MIN(loc_x_start + blocksize - 1, n_cols_loc)

        idx_cols(loc_x_start:loc_x_end) = (/(i,i=glob_x_start,glob_x_end)/)
      end do

      do i_block_y=0,n_blocks_y_loc-1
        glob_y_start = (proc_y + i_block_y*n_procs_y)*blocksize + 1
        glob_y_end   = MIN(glob_y_start + blocksize - 1, n_rows_glob)
        loc_y_start  = i_block_y*blocksize + 1
        loc_y_end    = MIN(loc_y_start + blocksize - 1, n_rows_loc)
          
        idx_rows(loc_y_start:loc_y_end) = (/(i,i=glob_y_start,glob_y_end)/)
      end do

    end subroutine getLocalIndices
#endif

#ifdef MPI
    subroutine printLocalMatrix(matrix, n_rows, n_cols, name, MPIdata)
      implicit none

      ! arguments
      real(8), dimension(:,:), Intent(in)  :: matrix
      integer,                 Intent(in)  :: n_rows, n_cols
      Character(len=*),        Intent(in)  :: name
      Type(MPIinfo),           Intent(in)  :: MPIdata

      ! local variables
      Integer :: i, n_cols_cut
      Complex :: dummyv
      real(8), dimension(100) :: row

      dummyv = 1
      Do i=1,50000*MPIdata%rank
        dummyv=dummyv*3.33D0
      End Do

      write (*,*) ' '
      write (*,*) name
      write (*,*) 'local matrix at proc ', MPIdata%rank
      n_cols_cut = min(100,n_cols)
      Do i=1,n_rows
        row(1:n_cols_cut) = matrix(i,1:n_cols_cut)
        write (*,"(100f9.3)") row(1:n_cols_cut)
      End Do

    end subroutine printLocalMatrix
#endif

#ifdef MPI
    subroutine printGlobalMatrix(matrix, n_rows, n_cols, name, MPIdata)
#else 
    subroutine printGlobalMatrix(matrix, n_rows, n_cols, name)
#endif
      implicit none

      ! arguments
      real(8), dimension(:,:), Intent(in)  :: matrix
      integer,                 Intent(in)  :: n_rows, n_cols
      Character(len=*),        Intent(in)  :: name
#ifdef MPI
      Type(MPIinfo),           Intent(in)  :: MPIdata
#endif

      ! local variables
      Integer :: i, n_cols_cut
      real(8), dimension(100) :: row

#ifdef MPI
      If (MPIdata%rank == 0) then
#endif
        write (*,*) 'global matrix', name
        n_cols_cut = min(100,n_cols)
        Do i=1,n_rows
          row(1:n_cols_cut) = matrix(i,1:n_cols_cut)
          write (*,"(100f9.3)") row(1:n_cols_cut)
        End Do
#ifdef MPI
      End If
#endif

    end subroutine printGlobalMatrix


end module test_helpers
