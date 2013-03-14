! !MODULE:  modOlpalon_test
! !DESCRIPTION:
! Modules with test and helper routines for testing olpalon
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
module modOlpalon_test
    Use fruit
    Use test_helpers
    Use modinput,        Only: input
    Use mod_muffin_tin,  Only: idxlm, lmmaxapw,lmmaxmat
    Use mod_APW_LO,      Only: apword, apwordmax, nlomax, lorbl, nlorb, lolmmax
    Use mod_Gkvector,    Only: ngkmax
    Use mod_atoms,       Only: natmtot
    Use mod_eigensystem, Only: oalo, idxlo
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix,evsystem,newsystem,deletesystem
#ifdef MPI
    use modmpi
#endif

    Double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains

    Subroutine testcaseOlpalonSerial
      Implicit None

      Call set_test_name ('Matching coefficients, apwalm')
      Call testOlpalon_APW_Serial
      Call set_test_name ('Radial integrals, oalo')
      Call testOlpalon_APWLO_Serial

    End Subroutine testcaseOlpalonSerial


#ifdef MPI
    Subroutine testcaseOlpalon1Proc
      Implicit None

      Call set_test_name ('Matching coefficients, apwalm')
      Call testOlpalon_APW_1Proc
      Call set_test_name ('Radial integrals, oalo')
      Call testOlpalon_APWLO_1Proc

    End Subroutine testcaseOlpalon1Proc
#endif


#ifdef MPI
    Subroutine testcaseOlpalon4Proc
      Implicit None

      Call set_test_name ('Matching coefficients, apwalm')
      Call testOlpalon_APW_4Proc
      Call set_test_name ('Radial integrals, oalo')
      Call testOlpalon_APWLO_4Proc

    End Subroutine testcaseOlpalon4Proc
#endif


! initialisation of global variables
    Subroutine initGlobals(lmaxmat,lmaxapw,gsize)
      Implicit None

      Integer, Intent(in) :: lmaxmat,lmaxapw,gsize

      Integer l,m,lm,i,ilo

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

      nlomax=10

      If (allocated(oalo)) Deallocate (oalo)
      Allocate (oalo(apwordmax, nlomax, natmtot))
      
      lolmmax=25
      nlorb(1)=7
      lorbl(1:7,1)=(/1,0,1,0,0,2,0/)
      If (allocated(idxlo)) Deallocate (idxlo)
      Allocate (idxlo(lolmmax, nlomax, natmtot)) 
! copied from genidxlo
      i=0
      Do ilo = 1, nlorb (1)
        l = lorbl (ilo, 1)
          Do m = - l, l
            i = i + 1
            lm = idxlm (l, m) 
            idxlo (lm, ilo, 1) = i
          End Do
      End Do
    End Subroutine initGlobals


! deallocation of global variables       
    Subroutine freeGlobals
      Implicit None

      Deallocate(oalo)
      Deallocate(idxlm,idxlo,input%groundstate)

    End Subroutine freeGlobals



!------------------------------------------------------------------------------
! test testOlpalon_APW_Serial
!------------------------------------------------------------------------------
! 1st test, serial
! The purpose is to test whether the matching coefficients apwalm are handled properly.
! The matching coefficients apwalm depend on g1.
! The radial integrals oalo are constant. 
    Subroutine testOlpalon_APW_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,gsize=9,nmatp=27)

      Integer :: l1,m1,lm1,lm2,g1
      integer :: ilo,last
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type(evsystem)           :: system
      Type(HermitianMatrix)    :: overlap_ref


! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,gsize)

      Call newsystem(system,nmatp,gsize)
      Call newmatrix(overlap_ref,nmatp)
! initialisation is finished

! The condition for the spherically symmetric evaluation
      oalo(:,:,:)=1d0

      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
      Do l1=0,lmaxapw
         Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
               apwalm(g1,:,lm1,1)=cmplx(g1,0,8)
            Enddo
         Enddo
      Enddo

      last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
      do g1=1,gsize
         do ilo=1,nlorb (1)
            l1=lorbl(ilo,1)
            lm1=idxlm(l1,-l1)
            lm2=idxlm(l1,l1)
            overlap_ref%za(g1,gsize+idxlo(lm1,ilo,1):gsize+idxlo(lm2,ilo,1))=cmplx(g1,0,8)
         enddo
      enddo

      Call olpalon(system,1,1,apwalm)

      Call assert_equals(nmatp, system%overlap%size, 'checking result rank')
      Call assert_equals(nmatp, size(system%overlap%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(system%overlap%za,2), 'checking result size cols')
      Call assert_equals(overlap_ref%za, system%overlap%za, nmatp, nmatp, tol, 'checking result numbers')
      Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking hamilton=0')

! finalisation
      Call deletesystem(system)
      Call deletematrix(overlap_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testOlpalon_APW_Serial


!------------------------------------------------------------------------------
! test testOlpalon_APW_Serial
!------------------------------------------------------------------------------
! 2nd test, serial
! The purpose is to test whether the radial integrals oalo are handled properly.
! The matching coefficients apwalm depend on g1.
! The radial integrals oalo depend on the LO index. 

    Subroutine testOlpalon_APWLO_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,gsize=9,nmatp=27)

      Integer :: l1,m1,lm1,lm2,g1
      integer :: ilo,last
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type(evsystem)           :: system
      Type(HermitianMatrix)    :: overlap_ref


! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,gsize)

      Call newsystem(system,nmatp,gsize)
      Call newmatrix(overlap_ref,nmatp)
! initialisation is finished

! The condition for the spherically symmetric evaluation
      oalo(:,:,:)=1d8
      oalo(:,:,1)=0d0

      system%overlap%za(:,:)=cmplx(0,0,8)
      overlap_ref%za(:,:)=cmplx(0,0,8)

      Do ilo= 1,nlorb (1)
         oalo(1,ilo,1)=dble(ilo)
      End Do

      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
      Do l1=0,lmaxapw
         Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
               apwalm(g1,:,lm1,1)=cmplx(g1,0,8)
            Enddo
         Enddo
      Enddo

      last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
      do g1=1,gsize
         do ilo=1,nlorb (1)
            l1=lorbl(ilo,1)
            lm1=idxlm(l1,-l1)
            lm2=idxlm(l1,l1)
            overlap_ref%za(g1,gsize+idxlo(lm1,ilo,1):gsize+idxlo(lm2,ilo,1))=cmplx(g1*ilo,0,8)
         enddo
      enddo

      Call olpalon(system,1,1,apwalm)

      Call assert_equals(nmatp, system%overlap%size, 'checking result rank')
      Call assert_equals(nmatp, size(system%overlap%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(system%overlap%za,2), 'checking result size cols')
      Call assert_equals(overlap_ref%za, system%overlap%za, nmatp, nmatp, tol, 'checking result numbers')
      Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking hamilton=0')

! finalisation
      Call deletesystem(system)
      Call deletematrix(overlap_ref)
      Deallocate(apwalm)
! deallocation of global variables 
      Call freeGlobals
    End Subroutine testOlpalon_APWLO_Serial


!------------------------------------------------------------------------------
! test testOlpalon_APW_1Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 1 proc
! The purpose is to test whether the matching coefficients apwalm are handled properly.
! The matching coefficients apwalm depend on g1.
! The radial integrals oalo are constant. 
#ifdef MPI
    Subroutine testOlpalon_APW_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,gsize=9,nmatp=27)

      Integer :: l1,m1,lm1,lm2,g1
      integer :: ilo
      Complex (8), allocatable  :: apwalm (:, :, :, :)
      Type(evsystem)            :: system
      Type(HermitianMatrix)     :: overlap_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

! initialisation of global variables
         Call initGlobals(lmaxmat,lmaxapw,gsize)

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(overlap_ref,nmatp)
! initialisation is finished

! The condition for the spherically symmetric evaluation
         oalo(:,:,:)=1d0

         allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1=1,gsize
                  apwalm(g1,:,lm1,1)=cmplx(g1,0,8)
               Enddo
            Enddo
         Enddo

         do g1=1,gsize
            do ilo=1,nlorb (1)
               l1=lorbl(ilo,1)
               lm1=idxlm(l1,-l1)
               lm2=idxlm(l1,l1)
               overlap_ref%za(g1,gsize+idxlo(lm1,ilo,1):gsize+idxlo(lm2,ilo,1))=cmplx(g1,0,8)
            enddo
         enddo

         Call olpalon(system,1,1,apwalm)

         Call assert_equals(nmatp, system%overlap%size, 'checking result rank')
         Call assert_equals(nmatp, size(system%overlap%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(system%overlap%za,2), 'checking result size cols')
         Call assert_equals(overlap_ref%za, system%overlap%za, nmatp, nmatp, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking hamilton=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(overlap_ref)
         Deallocate(apwalm)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testOlpalon_APW_1Proc
#endif


!------------------------------------------------------------------------------
! test testOlpalon_APW_1Proc
!------------------------------------------------------------------------------
! 2nd test, MPI with 1 proc
! The purpose is to test whether the radial integrals oalo are handled properly.
! The matching coefficients apwalm depend on g1.
! The radial integrals oalo depend on the LO index. 
#ifdef MPI
    Subroutine testOlpalon_APWLO_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,gsize=9,nmatp=27)

      Integer :: l1,m1,lm1,lm2,g1
      integer :: ilo
      Complex (8), allocatable  :: apwalm (:, :, :, :)
      Type(evsystem)        :: system
      Type(HermitianMatrix) :: overlap_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

! initialisation of global variables
         Call initGlobals(lmaxmat,lmaxapw,gsize)

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(overlap_ref,nmatp)
! initialisation is finished

! The condition for the spherically symmetric evaluation
         oalo(:,:,:)=1d8
         oalo(:,:,1)=0d0

         Do ilo= 1,nlorb (1)
            oalo(1,ilo,1)=dble(ilo)
         End Do

         allocate(apwalm(gsize, apwordmax, lmmaxapw, natmtot))
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1=1,gsize
                  apwalm(g1,:,lm1,1)=cmplx(g1,0,8)
               Enddo
            Enddo
         Enddo

         do g1=1,gsize
            do ilo=1,nlorb (1)
               l1=lorbl(ilo,1)
               lm1=idxlm(l1,-l1)
               lm2=idxlm(l1,l1)
               overlap_ref%za(g1,gsize+idxlo(lm1,ilo,1):gsize+idxlo(lm2,ilo,1))=cmplx(g1*ilo,0,8)
            enddo
         enddo

         Call olpalon(system,1,1,apwalm)

         Call assert_equals(nmatp, system%overlap%size, 'checking result rank')
         Call assert_equals(nmatp, size(system%overlap%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(system%overlap%za,2), 'checking result size cols')
         Call assert_equals(overlap_ref%za, system%overlap%za, nmatp, nmatp, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking hamilton=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(overlap_ref)
         Deallocate(apwalm)
! deallocation of global variables 
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testOlpalon_APWLO_1Proc
#endif


!------------------------------------------------------------------------------
! test testOlpalon_APW_4Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 4 procs
! The purpose is to test whether the matching coefficients apwalm are handled properly.
! The matching coefficients apwalm depend on g1.
! The radial integrals oalo are constant. 
#ifdef MPI
    Subroutine testOlpalon_APW_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,gsize=9,nmatp=27)

      Integer :: l1,m1,lm1,lm2,g1,g1_loc
      integer :: ilo
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type(evsystem)           :: system
      Type(HermitianMatrix)    :: overlap_ref

      Complex (8), Dimension(nmatp,nmatp) :: overlap_ref_global

! Externals
      Integer,    External :: NUMROC

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, ncols_loc, gsize_nrows_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc2glob, dummy

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

! initialisation of global variables
         Call initGlobals(lmaxmat,lmaxapw,gsize)

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(overlap_ref,nmatp, DISTRIBUTE_2D)

! init datastructures for splitting apwalm
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         allocate(apwalm1_loc2glob(gsize_nrows_loc), dummy(MPIglobal%blocksize))
         Call getLocalIndices(gsize, MPIglobal%nproccols*MPIglobal%blocksize, apwalm1_loc2glob, dummy, MPIglobal)

! initialisation is finished

! The condition for the spherically symmetric evaluation
         oalo(:,:,:)=1d0

         allocate(apwalm (gsize_nrows_loc, apwordmax, lmmaxapw, natmtot))
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1_loc=1,gsize_nrows_loc !splitting of apwalm along first dimension
                                          ! according to the processor rows
                  g1=apwalm1_loc2glob(g1_loc) 
                  apwalm(g1_loc,:,lm1,1)=cmplx(g1,0,8)
               Enddo
            Enddo
         Enddo

         overlap_ref_global = Cmplx(0,0,8)
         do g1=1,gsize
            do ilo=1,nlorb (1)
               l1=lorbl(ilo,1)
               lm1=idxlm(l1,-l1)
               lm2=idxlm(l1,l1)
               overlap_ref_global(g1,gsize+idxlo(lm1,ilo,1):gsize+idxlo(lm2,ilo,1))=cmplx(g1,0,8)
            enddo
         enddo
         Call getBlockDistributedLoc(overlap_ref_global, overlap_ref%za, MPIglobal)

         Select Case (MPIglobal%myprocrow)
            Case (0)
               nrows_loc = 14
            Case (1)
               nrows_loc = 13
         End Select
         Select Case (MPIglobal%myproccol)
            Case (0)
               ncols_loc = 14
            Case (1)
               ncols_loc = 13
         End Select

         Call olpalon(system,1,1,apwalm)

         Call assert_equals(nmatp, system%overlap%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(system%overlap%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(system%overlap%za,2), 'checking result size cols')
         Call assert_equals(overlap_ref%za, system%overlap%za, nrows_loc, ncols_loc, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking hamilton=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(overlap_ref)
         Deallocate(apwalm)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testOlpalon_APW_4Proc
#endif


!------------------------------------------------------------------------------
! test testOlpalon_APW_4Proc
!------------------------------------------------------------------------------
! 2nd test, MPI with 4 procs
! The purpose is to test whether the radial integrals oalo are handled properly.
! The matching coefficients apwalm depend on g1.
! The radial integrals oalo depend on the LO index. 
#ifdef MPI
    Subroutine testOlpalon_APWLO_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,gsize=9,nmatp=27)

      Integer :: l1,m1,lm1,lm2,g1,g1_loc
      integer :: ilo
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type(evsystem)        :: system
      Type(HermitianMatrix) :: overlap_ref
      Complex (8), Dimension(nmatp,nmatp) :: overlap_ref_global

! Externals
      Integer,    External :: NUMROC

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, ncols_loc, gsize_nrows_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc2glob, dummy

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

! initialisation of global variables
         Call initGlobals(lmaxmat,lmaxapw,gsize)

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(overlap_ref,nmatp, DISTRIBUTE_2D)

! init datastructures for splitting apwalm
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         allocate(apwalm1_loc2glob(gsize_nrows_loc), dummy(MPIglobal%blocksize))
         Call getLocalIndices(gsize, MPIglobal%nproccols*MPIglobal%blocksize, apwalm1_loc2glob, dummy, MPIglobal)

! initialisation is finished

! The condition for the spherically symmetric evaluation
         oalo(:,:,:)=1d8
         oalo(:,:,1)=0d0

         Do ilo= 1,nlorb (1)
            oalo(1,ilo,1)=dble(ilo)
         End Do

         allocate(apwalm (gsize_nrows_loc, apwordmax, lmmaxapw, natmtot))
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1_loc=1,gsize_nrows_loc !splitting of apwalm along first dimension
                                          ! according to the processor rows
                  g1=apwalm1_loc2glob(g1_loc) 
                  apwalm(g1_loc,:,lm1,1)=cmplx(g1,0,8)
               Enddo
            Enddo
         Enddo

         overlap_ref_global = Cmplx(0,0,8)
         do g1=1,gsize
            do ilo=1,nlorb (1)
               l1=lorbl(ilo,1)
               lm1=idxlm(l1,-l1)
               lm2=idxlm(l1,l1)
               overlap_ref_global(g1,gsize+idxlo(lm1,ilo,1):gsize+idxlo(lm2,ilo,1))=cmplx(g1*ilo,0,8)
            enddo
         enddo
         Call getBlockDistributedLoc(overlap_ref_global, overlap_ref%za, MPIglobal)

         Select Case (MPIglobal%myprocrow)
            Case (0)
               nrows_loc = 14
            Case (1)
               nrows_loc = 13
         End Select
         Select Case (MPIglobal%myproccol)
            Case (0)
               ncols_loc = 14
            Case (1)
               ncols_loc = 13
         End Select

         Call olpalon(system,1,1,apwalm)

         Call assert_equals(nmatp, system%overlap%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(system%overlap%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(system%overlap%za,2), 'checking result size cols')
         Call assert_equals(overlap_ref%za, system%overlap%za, nrows_loc, ncols_loc, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking hamilton=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(overlap_ref)
         Deallocate(apwalm)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testOlpalon_APWLO_4Proc
#endif

end module modOlpalon_test
