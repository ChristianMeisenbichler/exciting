! !MODULE:  modHmllolon_test
! !DESCRIPTION:
! Modules with test and helper routines for testing hmllolon
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
module modHmllolon_test
    Use fruit
    Use test_helpers
    Use modinput,        Only: input
    Use mod_muffin_tin,  Only: idxlm, lmmaxapw,lmmaxvr,lmmaxmat
    Use mod_APW_LO,      Only: nlomax, lorbl, nlorb, lolmmax
    Use mod_Gkvector,    Only: ngkmax
    Use mod_atoms,       Only: natmtot
    Use mod_eigensystem, Only: gntyry, hlolo, idxlo
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix,evsystem,newsystem,deletesystem
#ifdef MPI
    use modmpi
#endif

    Double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains


    Subroutine testcaseHmllolonSerial
      Implicit None

      Call set_test_name ('Gaunt coefficients')
      Call testHmllolon_Gnt_Serial
      Call set_test_name ('Spherically symmetric contributions, hlolo')
      Call testHmllolon_SphSym_Serial
      Call set_test_name ('Spherically symmetric and asymmetric contributions, hlolo')
      Call testHmllolon_SphSymAsym_Serial

    End Subroutine testcaseHmllolonSerial


#ifdef MPI
    Subroutine testcaseHmllolon1Proc
      Implicit None

      Call set_test_name ('Gaunt coefficients')
      Call testHmllolon_Gnt_1Proc
      Call set_test_name ('Spherically symmetric contributions, hlolo')
      Call testHmllolon_SphSym_1Proc
      Call set_test_name ('Spherically symmetric and asymmetric contributions, hlolo')
      Call testHmllolon_SphSymAsym_1Proc

    End Subroutine testcaseHmllolon1Proc
#endif


#ifdef MPI
    Subroutine testcaseHmllolon4Proc
      Implicit None

      Call set_test_name ('Gaunt coefficients')
      Call testHmllolon_Gnt_4Proc
      Call set_test_name ('Spherically symmetric contributions, hlolo')
      Call testHmllolon_SphSym_4Proc
      Call set_test_name ('Spherically symmetric and asymmetric contributions, hlolo')
      Call testHmllolon_SphSymAsym_4Proc

    End Subroutine testcaseHmllolon4Proc
#endif


! initialisation of global variables
    Subroutine initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)
      Implicit None

      Integer, Intent(in) :: lmaxmat,lmaxapw,lmaxvr,gsize

      Integer l,m,lm,i,ilo

      If (.not.associated(input%groundstate)) Allocate(input%groundstate)
      input%groundstate%lmaxmat=lmaxmat
      input%groundstate%lmaxapw=lmaxapw
      input%groundstate%lmaxvr=lmaxvr
      lmmaxapw=(lmaxapw+1)**2
      lmmaxvr=(lmaxvr+1)**2
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

      nlomax=10

      If (allocated(hlolo)) Deallocate (hlolo)
      Allocate (hlolo(nlomax, nlomax, lmmaxvr, natmtot))
      
      lolmmax=25
      nlorb(1)=3
      lorbl(1:3,1)=(/0,1,2/)
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

      Deallocate(hlolo)
      Deallocate(gntyry)
      Deallocate(idxlm,idxlo,input%groundstate)

    End Subroutine freeGlobals



! allocate and generate complex Gaunt coefficient array
! assumes that the global variables are set
! copied from init1.f90
    Subroutine initGntyry
      Implicit None

      Integer l1,m1,lm1,l2,m2,lm2,l3,m3,lm3

      Complex (8) gauntyry
      External gauntyry

      If (allocated(gntyry)) Deallocate (gntyry)
      Allocate (gntyry(lmmaxmat, lmmaxvr, lmmaxapw))
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do l2 = 0, input%groundstate%lmaxvr
               Do m2 = - l2, l2
                  lm2 = idxlm (l2, m2)
                  Do l3 = 0, input%groundstate%lmaxapw
                     Do m3 = - l3, l3
                        lm3 = idxlm (l3, m3)
                        gntyry (lm1, lm2, lm3) = gauntyry (l1, l2, l3, m1, m2, m3)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do


    End Subroutine initGntyry


!------------------------------------------------------------------------------
! test testHmllolon_Gnt_Serial
!------------------------------------------------------------------------------
! 1st test, serial
! The purpose is to test whether gaunt coefficients are handled properly.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmllolon_Gnt_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=43)
      
      Integer :: l1,m1,l2,m2,l3,m3
      integer :: i,j,ilo,ilo2
      complex(8)               :: test
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newsystem(system,nmatp,gsize)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      hlolo(:,:,:,:)=1d0

      i=0
      Do ilo= 1,nlorb (1)
        l1 = lorbl (ilo, 1)
        Do m1 = - l1, l1
          i=i+1
          j=0
          Do ilo2= 1,nlorb (1)
            l3 = lorbl (ilo2, 1)
            Do m3 = - l3, l3
              j=j+1
              if (i.le.j) then
                test=cmplx(0,0,8)
                Do l2 = 0, input%groundstate%lmaxvr
                  Do m2 = - l2, l2
                    test = test + gauntyry (l1, l2, l3, m1, m2, m3)
                  End Do
                End Do
                hamilton_ref%za(i+gsize,j+gsize)=test
              endif
            End Do
          End Do
        End Do
      End Do

      Call hmllolon(system,1,1)

      Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(system%hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(system%hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, system%hamilton%za, nmatp, nmatp, tol, 'checking result numbers')
      Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
      Call deletesystem(system)
      Call deletematrix(hamilton_ref)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmllolon_Gnt_Serial


!------------------------------------------------------------------------------
! test testHmllolon_SphSym_Serial
!------------------------------------------------------------------------------
! 2nd test, serial
! The purpose is to test whether the radial integrals are handled properly.
! The radial integrals hlolo depend on the local orbital index ilo.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmllolon_SphSym_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=43)

      Integer :: l1,m1,lm1
      integer :: ilo
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newsystem(system,nmatp,gsize)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      hlolo(:,:,:,:)=0d0
      Do ilo=1,nlorb (1)
        hlolo(ilo,ilo,1,1)=dble(ilo)*sqrt(4d0*pi)
      End Do

! The answer is a diagonal matrix      
      Do ilo=1,nlorb (1)
        l1=lorbl(ilo,1)
        Do m1=-l1,l1       
          lm1=idxlm(l1,m1)
          hamilton_ref%za(gsize+idxlo(lm1,ilo,1),gsize+idxlo(lm1,ilo,1))=cmplx(ilo,0,8)
        Enddo
      End Do

      Call hmllolon(system,1,1)

      Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(system%hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(system%hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, system%hamilton%za, nmatp, nmatp, tol, 'checking result numbers')
      Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
      Call deletesystem(system)
      Call deletematrix(hamilton_ref)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testHmllolon_SphSym_Serial

!------------------------------------------------------------------------------
! test testHmllolon_SphSymAsym_Serial
!------------------------------------------------------------------------------
! 3rd test, serial
! The purpose is to test whether the radial integrals are handled properly.
! The radial integrals hlolo depend on the local orbital index ilo and the spherical harmonics of the potential.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmllolon_SphSymAsym_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)

      Integer :: l1,m1,lm1,lm2,l3,m3,lm3
      integer :: i,j,ilo,ilo2
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newsystem(system,nmatp,gsize)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      hlolo(:,:,:,:)=0d0
      Do ilo=1,nlorb (1)
        Do ilo2=1,nlorb (1)
          Do lm1=1,lmmaxvr
            hlolo(ilo,ilo2,lm1,1)=cos(dble(ilo))*cos(dble(ilo2))*sin(dble(lm1))
          End Do
        End Do
      End Do

! Preparing the answer... 
      Do ilo=1,nlorb(1)
        l1=lorbl(ilo,1)
        Do m1=-l1,l1
          lm1=idxlm(l1,m1)
          i=idxlo(lm1,ilo,1)
          Do ilo2=1,nlorb(1)
            l3=lorbl(ilo2,1)
            Do m3=-l3,l3
              lm3=idxlm(l3,m3) 
              j=idxlo(lm3,ilo2,1)
              do lm2=1,lmmaxvr
                if (i.le.j) then
                  hamilton_ref%za(gsize+i,gsize+j)=hamilton_ref%za(gsize+i,gsize+j)+gntyry(lm1,lm2,lm3)*cos(dble(ilo))*cos(dble(ilo2))*sin(dble(lm2))
                endif
              enddo
            End Do
          End Do
        End Do
      End Do

      Call hmllolon(system,1,1)

      Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(system%hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(system%hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, system%hamilton%za, nmatp, nmatp, tol, 'checking result numbers')
      Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
      Call deletesystem(system)
      Call deletematrix(hamilton_ref)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testHmllolon_SphSymAsym_Serial


!------------------------------------------------------------------------------
! test testHmllolon_Gnt_1Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 1 proc
! The purpose is to test whether gaunt coefficients are handled properly.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmllolon_Gnt_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=43)
      
      Integer :: l1,m1,l2,m2,l3,m3
      integer :: i,j,ilo,ilo2
      complex(8)               :: test
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref

! Externals
      Complex(8), External :: gauntyry

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
         Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
         Call initGntyry

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

         hlolo(:,:,:,:)=1d0

         i=0
         Do ilo= 1,nlorb (1)
         l1 = lorbl (ilo, 1)
         Do m1 = - l1, l1
            i=i+1
            j=0
            Do ilo2= 1,nlorb (1)
               l3 = lorbl (ilo2, 1)
               Do m3 = - l3, l3
               j=j+1
               if (i.le.j) then
                  test=cmplx(0,0,8)
                  Do l2 = 0, input%groundstate%lmaxvr
                     Do m2 = - l2, l2
                     test = test + gauntyry (l1, l2, l3, m1, m2, m3)
                     End Do
                  End Do
                  hamilton_ref%za(i+gsize,j+gsize)=test
               endif
               End Do
            End Do
         End Do
         End Do

         Call hmllolon(system,1,1)

         Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(system%hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(system%hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, system%hamilton%za, nmatp, nmatp, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(hamilton_ref)
! deallocation of global variables   
         Call freeGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmllolon_Gnt_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmllolon_SphSym_1Proc
!------------------------------------------------------------------------------
! 2nd test, MPI with 1 proc
! The purpose is to test whether the radial integrals are handled properly.
! The radial integrals hlolo depend on the local orbital index ilo.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmllolon_SphSym_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=43)

      Integer :: l1,m1,lm1
      integer :: ilo
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref

! Externals
      Complex(8), External :: gauntyry

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
         Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
         Call initGntyry

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

         hlolo(:,:,:,:)=0d0
         Do ilo=1,nlorb (1)
            hlolo(ilo,ilo,1,1)=dble(ilo)*sqrt(4d0*pi)
         End Do

! The answer is a diagonal matrix      
         Do ilo=1,nlorb (1)
            l1=lorbl(ilo,1)
            Do m1=-l1,l1       
               lm1=idxlm(l1,m1)
               hamilton_ref%za(gsize+idxlo(lm1,ilo,1),gsize+idxlo(lm1,ilo,1))=cmplx(ilo,0,8)
            Enddo
         End Do

         Call hmllolon(system,1,1)

         Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(system%hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(system%hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, system%hamilton%za, nmatp, nmatp, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(hamilton_ref)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmllolon_SphSym_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmllolon_SphSymAsym_1Proc
!------------------------------------------------------------------------------
! 3rd test, MPI with 1 proc
! The purpose is to test whether the radial integrals are handled properly.
! The radial integrals hlolo depend on the local orbital index ilo and the spherical harmonics of the potential.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmllolon_SphSymAsym_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)

      Integer :: l1,m1,lm1,lm2,l3,m3,lm3
      integer :: i,j,ilo,ilo2
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref

! Externals
      Complex(8), External :: gauntyry

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
         Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
         Call initGntyry

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

         hlolo(:,:,:,:)=0d0
         Do ilo=1,nlorb (1)
            Do ilo2=1,nlorb (1)
               Do lm1=1,lmmaxvr
                  hlolo(ilo,ilo2,lm1,1)=cos(dble(ilo))*cos(dble(ilo2))*sin(dble(lm1))
               End Do
            End Do
         End Do

! Preparing the answer... 
         Do ilo=1,nlorb(1)
            l1=lorbl(ilo,1)
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               i=idxlo(lm1,ilo,1)
               Do ilo2=1,nlorb(1)
                  l3=lorbl(ilo2,1)
                  Do m3=-l3,l3
                  lm3=idxlm(l3,m3) 
                  j=idxlo(lm3,ilo2,1)
                  do lm2=1,lmmaxvr
                     if (i.le.j) then
                        hamilton_ref%za(gsize+i,gsize+j)=hamilton_ref%za(gsize+i,gsize+j)+gntyry(lm1,lm2,lm3)*cos(dble(ilo))*cos(dble(ilo2))*sin(dble(lm2))
                     endif
                  enddo
                  End Do
               End Do
            End Do
         End Do

         Call hmllolon(system,1,1)

         Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(system%hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(system%hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, system%hamilton%za, nmatp, nmatp, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(hamilton_ref)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmllolon_SphSymAsym_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmllolon_Gnt_4Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 4 procs
! The purpose is to test whether gaunt coefficients are handled properly.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmllolon_Gnt_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=43)
      
      Integer :: l1,m1,l2,m2,l3,m3
      integer :: i,j,ilo,ilo2
      complex(8)               :: test
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! Externals
      Complex(8), External :: gauntyry

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, ncols_loc

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

! initialisation of global variables
         Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
         Call initGntyry

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

         hlolo(:,:,:,:)=1d0

         hamilton_ref_global = cmplx(0,0,8)
         i=0
         Do ilo= 1,nlorb (1)
            l1 = lorbl (ilo, 1)
            Do m1 = - l1, l1
               i=i+1
               j=0
               Do ilo2= 1,nlorb (1)
                  l3 = lorbl (ilo2, 1)
                  Do m3 = - l3, l3
                  j=j+1
                  if (i.le.j) then
                     test=cmplx(0,0,8)
                     Do l2 = 0, input%groundstate%lmaxvr
                        Do m2 = - l2, l2
                        test = test + gauntyry (l1, l2, l3, m1, m2, m3)
                        End Do
                     End Do
                     hamilton_ref_global(i+gsize,j+gsize)=test
                  endif
                  End Do
               End Do
            End Do
         End Do
         Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal)

         Select Case (MPIglobal%myprocrow)
            Case (0)
               nrows_loc = 22
            Case (1)
               nrows_loc = 21
         End Select
         Select Case (MPIglobal%myproccol)
            Case (0)
               ncols_loc = 22
            Case (1)
               ncols_loc = 21
         End Select

         Call hmllolon(system,1,1)

         Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(system%hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(system%hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, system%hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(hamilton_ref)
! deallocation of global variables   
         Call freeGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmllolon_Gnt_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmllolon_SphSym_4Proc
!------------------------------------------------------------------------------
! 2nd test, MPI with 4 procs
! The purpose is to test whether the radial integrals are handled properly.
! The radial integrals hlolo depend on the local orbital index ilo.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmllolon_SphSym_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=43)

      Integer :: l1,m1,lm1
      integer :: ilo
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! Externals
      Complex(8), External :: gauntyry

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, ncols_loc

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

! initialisation of global variables
         Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
         Call initGntyry

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

         hlolo(:,:,:,:)=0d0
         Do ilo=1,nlorb (1)
            hlolo(ilo,ilo,1,1)=dble(ilo)*sqrt(4d0*pi)
         End Do

! The answer is a diagonal matrix      
         hamilton_ref_global = cmplx(0,0,8)
         Do ilo=1,nlorb (1)
            l1=lorbl(ilo,1)
            Do m1=-l1,l1       
               lm1=idxlm(l1,m1)
               hamilton_ref_global(gsize+idxlo(lm1,ilo,1),gsize+idxlo(lm1,ilo,1))=cmplx(ilo,0,8)
            Enddo
         End Do
         Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal)

         Select Case (MPIglobal%myprocrow)
            Case (0)
               nrows_loc = 22
            Case (1)
               nrows_loc = 21
         End Select
         Select Case (MPIglobal%myproccol)
            Case (0)
               ncols_loc = 22
            Case (1)
               ncols_loc = 21
         End Select

         Call hmllolon(system,1,1)

         Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(system%hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(system%hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, system%hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(hamilton_ref)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmllolon_SphSym_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmllolon_SphSymAsym_4Proc
!------------------------------------------------------------------------------
! 3rd test, MPI with 4 procs
! The purpose is to test whether the radial integrals are handled properly.
! The radial integrals hlolo depend on the local orbital index ilo and the spherical harmonics of the potential.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmllolon_SphSymAsym_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)

      Integer :: l1,m1,lm1,lm2,l3,m3,lm3
      integer :: i,j,ilo,ilo2
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! Externals
      Complex(8), External :: gauntyry

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, ncols_loc

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

! initialisation of global variables
         Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
         Call initGntyry

         Call newsystem(system,nmatp,gsize)
         Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

         hlolo(:,:,:,:)=0d0
         Do ilo=1,nlorb (1)
            Do ilo2=1,nlorb (1)
               Do lm1=1,lmmaxvr
                  hlolo(ilo,ilo2,lm1,1)=cos(dble(ilo))*cos(dble(ilo2))*sin(dble(lm1))
               End Do
            End Do
         End Do

! Preparing the answer... 
         hamilton_ref_global = cmplx(0,0,8)
         Do ilo=1,nlorb(1)
            l1=lorbl(ilo,1)
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               i=idxlo(lm1,ilo,1)
               Do ilo2=1,nlorb(1)
                  l3=lorbl(ilo2,1)
                  Do m3=-l3,l3
                     lm3=idxlm(l3,m3) 
                     j=idxlo(lm3,ilo2,1)
                     do lm2=1,lmmaxvr
                        if (i.le.j) then
                           hamilton_ref_global(gsize+i,gsize+j)=hamilton_ref_global(gsize+i,gsize+j)+gntyry(lm1,lm2,lm3)*cos(dble(ilo))*cos(dble(ilo2))*sin(dble(lm2))
                        endif
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal)

         Select Case (MPIglobal%myprocrow)
            Case (0)
               nrows_loc = 12
            Case (1)
               nrows_loc = 11
         End Select
         Select Case (MPIglobal%myproccol)
            Case (0)
               ncols_loc = 12
            Case (1)
               ncols_loc = 11
         End Select

         Call hmllolon(system,1,1)

         Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(system%hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(system%hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, system%hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(hamilton_ref)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmllolon_SphSymAsym_4Proc
#endif


end module modHmllolon_test
