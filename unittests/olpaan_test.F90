! !MODULE:  modOlpaan_test
! !DESCRIPTION:
! Modules with test and helper routines for testing olpaan
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
module modOlpaan_test
    Use fruit
    Use test_helpers
    Use modinput,        Only: input
    Use mod_muffin_tin,  Only: idxlm, lmmaxapw,lmmaxmat
    Use mod_APW_LO,      Only: apword, apwordmax
    Use mod_Gkvector,    Only: ngkmax
    Use mod_atoms,       Only: natmtot, idxas
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix
#ifdef MPI
    use modmpi
#endif

    Double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains

    Subroutine testcaseOlpaanSerial
      Implicit None

      Call set_test_name ('overlap matrix calculation')
      Call testOlpaan_Serial

    End Subroutine testcaseOlpaanSerial


#ifdef MPI
    Subroutine testcaseOlpaan1Proc
      Implicit None

      Call set_test_name ('overlap matrix calculation')
      Call testOlpaan_1Proc

    End Subroutine testcaseOlpaan1Proc
#endif


#ifdef MPI
    Subroutine testcaseOlpaan4Proc
      Implicit None

      Call set_test_name ('overlap matrix calculation')
      Call testOlpaan_4Proc

    End Subroutine testcaseOlpaan4Proc
#endif

! initialisation of global variables
    Subroutine InitOverlapGlobals(lmaxmat,lmaxapw,gsize)
      Implicit None

      Integer, Intent(in) :: lmaxmat,lmaxapw,gsize

      Integer l,m,lm

      If (.not.associated(input%groundstate)) Allocate(input%groundstate)
      input%groundstate%lmaxmat=lmaxmat
      input%groundstate%lmaxapw=lmaxapw
      lmmaxapw=(lmaxapw+1)**2
      lmmaxmat=(lmaxmat+1)**2
! copied from init0.f90     
      If (allocated(idxlm)) Deallocate (idxlm)
      Allocate (idxlm(0:input%groundstate%lmaxapw,-input%groundstate%lmaxapw:input%groundstate%lmaxapw))
      lm = 0
      Do l = 0, input%groundstate%lmaxapw
         Do m = - l, l
            lm = lm + 1
            idxlm (l, m) = lm
         End Do
      End Do
  
      ngkmax=gsize
      natmtot=1
      apwordmax=1
      apword(:,:)=1

      idxas (:, :) = 1

    End Subroutine InitOverlapGlobals


! deallocation of global variables       
    Subroutine FreeOverlapGlobals
      Implicit None

      Deallocate(idxlm,input%groundstate)

    End Subroutine FreeOverlapGlobals




!------------------------------------------------------------------------------
! !TEST: testOlpaan_Serial
!------------------------------------------------------------------------------
! !DESCRIPTION:
! 1st test, serial
! The purpose is to test whether the matching coefficients (apwalm) are treated properly.
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
    Subroutine testOlpaan_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: overlap,overlap_ref

! initialisation of global variables
      Call InitOverlapGlobals(lmaxmat,lmaxapw,gsize)

      Call newmatrix(overlap,nmatp)
      Call newmatrix(overlap_ref,nmatp)
! initialisation is finished

      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      apwalm(:,:,:,:)=0d0
      Do l1=0,lmaxapw
         Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
               apwalm(g1,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
            Enddo 
         Enddo
      Enddo

      Do g2=1,gsize
         Do g1=1,gsize      
            overlap_ref%za(g1,g2)=cmplx(g1*g2,0,8)
         Enddo
      Enddo

      Call olpaan(overlap,1,1,gsize,apwalm)

      Call assert_equals(nmatp, overlap%size, 'checking result rank')
      Call assert_equals(nmatp, size(overlap%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(overlap%za,2), 'checking result size cols')
      Call assert_equals(overlap_ref%za, overlap%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(overlap)
      Call deletematrix(overlap_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call FreeOverlapGlobals    
    End Subroutine testOlpaan_Serial



!------------------------------------------------------------------------------
! !TEST: testOlpaan_1Proc
!------------------------------------------------------------------------------
! !DESCRIPTION:
! 1st test, MPI with 1 proc
! The purpose is to test whether the matching coefficients (apwalm) are treated properly.
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
#ifdef MPI
    Subroutine testOlpaan_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: overlap,overlap_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      If (MPIglobal_1D%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
         Call InitOverlapGlobals(lmaxmat,lmaxapw,gsize)

         Call newmatrix(overlap,nmatp)
         Call newmatrix(overlap_ref,nmatp)
! initialisation is finished

         allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
         apwalm(:,:,:,:)=0d0
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1=1,gsize
                  apwalm(g1,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
               Enddo 
            Enddo
         Enddo

         Do g2=1,gsize
            Do g1=1,gsize      
               overlap_ref%za(g1,g2)=cmplx(g1*g2,0,8)
            Enddo
         Enddo

         Call olpaan(overlap,1,1,gsize,apwalm,gsize)

         Call assert_equals(nmatp, overlap%size, 'checking result rank')
         Call assert_equals(nmatp, size(overlap%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(overlap%za,2), 'checking result size cols')
         Call assert_equals(overlap_ref%za, overlap%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
         Call deletematrix(overlap)
         Call deletematrix(overlap_ref)
         Deallocate(apwalm)
! deallocation of global variables   
         Call FreeOverlapGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testOlpaan_1Proc
#endif


!------------------------------------------------------------------------------
! !TEST: testOlpaan_4Proc
!------------------------------------------------------------------------------
! !DESCRIPTION:
! 1st test, MPI with 4 procs
! The purpose is to test whether the matching coefficients (apwalm) are treated properly.
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
#ifdef MPI
    Subroutine testOlpaan_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g1_loc,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: overlap,overlap_ref
      Complex (8), Dimension(nmatp,nmatp) :: overlap_ref_global

! MPI related variables
      Integer :: n_procs_test, ierror_t
      Integer :: nrows_loc, ncols_loc
      Integer :: gsize_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc2glob, dummy

! Externals
      Integer, External   :: NUMROC

      n_procs_test = 4
      Call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      If (MPIglobal_1D%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
         Call InitOverlapGlobals(lmaxmat,lmaxapw,gsize)

         Call newmatrix(overlap, nmatp, DISTRIBUTE_COLS)
         Call newmatrix(overlap_ref, nmatp, DISTRIBUTE_COLS)

! init datastructures for splitting apwalm
        gsize_loc = NUMROC(gsize, MPIglobal_1D%blocksize, MPIglobal_1D%myproccol, 0, MPIglobal_1D%nproccols)
        allocate(dummy(gsize), apwalm1_loc2glob(gsize_loc))
        Call getLocalIndices(gsize, gsize, dummy, apwalm1_loc2glob, MPIglobal_1D)

! initialisation is finished

         allocate(apwalm (gsize_loc, apwordmax, lmmaxapw, natmtot))      
         apwalm(:,:,:,:)=0d0
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1_loc=1,gsize_loc !splitting of apwalm along first dimension!
                  g1=apwalm1_loc2glob(g1_loc) 
                  apwalm(g1_loc,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
               Enddo 
            Enddo
         Enddo

         overlap_ref_global = Cmplx(0,0,8)
         Do g2=1,gsize
            Do g1=1,gsize      
               overlap_ref_global(g1,g2)=cmplx(g1*g2,0,8)
            Enddo
         Enddo
         Call getBlockDistributedLoc(overlap_ref_global, overlap_ref%za, MPIglobal_1D)

         nrows_loc = nmatp
         select case (MPIglobal_1D%rank)
            case (0)
               ncols_loc = 4
            case (1)
               ncols_loc = 3
            case (2)
               ncols_loc = 2
            case (3)
               ncols_loc = 2
         End select

         Call olpaan(overlap,1,1,gsize,apwalm,gsize_loc)

         Call assert_equals(nmatp, overlap%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(overlap%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(overlap%za,2), 'checking result size cols')
         Call assert_equals(overlap_ref%za, overlap%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
         Call deletematrix(overlap)
         Call deletematrix(overlap_ref)
         Deallocate(apwalm, apwalm1_loc2glob, dummy)
! deallocation of global variables   
         Call FreeOverlapGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testOlpaan_4Proc
#endif

end module modOlpaan_test
