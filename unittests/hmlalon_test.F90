! !MODULE:  modHmlalon_test
! !DESCRIPTION:
! Modules with test and helper routines for testing hmlalon
!
! !REVISION HISTORY:
!   Created: March 2013 (G. Huhs - BSC, A. Gulans)
!
module modHmlalon_test
    Use fruit
    Use test_helpers
    Use modinput,        Only: input
    Use mod_muffin_tin,  Only: idxlm, lmmaxapw,lmmaxvr,lmmaxmat
    Use mod_APW_LO,      Only: apword, apwordmax, nlomax, lorbl, nlorb, lolmmax
    Use mod_Gkvector,    Only: ngkmax
    Use mod_atoms,       Only: natmtot
    Use mod_eigensystem, Only: gntyry, hloa, idxlo
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix
#ifdef MPI
    use modmpi
#endif

    Double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains


    Subroutine testcaseHmlalonSerial
      Implicit None

      Call set_test_name ('Gaunt coefficients')
      Call testHmlalon_Gnt_Serial
      Call set_test_name ('Spherically symmetric contributions, hloa')
      Call testHmlalon_SphSymLO_Serial
      Call set_test_name ('Spherically symmetric contributions, apwalm')
      Call testHmlalon_SphSymAPW_Serial
      Call set_test_name ('Spherically symmetric and asymmetric contributions, apwalm')
      Call testHmlalon_SphSymAsymAPW_Serial
      Call set_test_name ('Spherically symmetric and asymmetric contributions, hloa')
      Call testHmlalon_SphSymAsymLO_Serial

    End Subroutine testcaseHmlalonSerial


#ifdef MPI
    Subroutine testcaseHmlalon1Proc
      Implicit None

      Call set_test_name ('Gaunt coefficients')
      Call testHmlalon_Gnt_1Proc
      Call set_test_name ('Spherically symmetric contributions, hloa')
      Call testHmlalon_SphSymLO_1Proc
      Call set_test_name ('Spherically symmetric contributions, apwalm')
      Call testHmlalon_SphSymAPW_1Proc
      Call set_test_name ('Spherically symmetric and asymmetric contributions, apwalm')
      Call testHmlalon_SphSymAsymAPW_1Proc
      Call set_test_name ('Spherically symmetric and asymmetric contributions, hloa')
      Call testHmlalon_SphSymAsymLO_1Proc

    End Subroutine testcaseHmlalon1Proc
#endif


#ifdef MPI
    Subroutine testcaseHmlalon4Proc
      Implicit None

      Call set_test_name ('Gaunt coefficients')
      Call testHmlalon_Gnt_4Proc
      Call set_test_name ('Spherically symmetric contributions, hloa')
      Call testHmlalon_SphSymLO_4Proc
      Call set_test_name ('Spherically symmetric contributions, apwalm')
      Call testHmlalon_SphSymAPW_4Proc
      Call set_test_name ('Spherically symmetric and asymmetric contributions, apwalm')
      Call testHmlalon_SphSymAsymAPW_4Proc
      Call set_test_name ('Spherically symmetric and asymmetric contributions, hloa')
      Call testHmlalon_SphSymAsymLO_4Proc
    End Subroutine testcaseHmlalon4Proc
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
      apwordmax=1
      apword(:,:)=1

      nlomax=10

      If (allocated(hloa)) Deallocate (hloa)
      Allocate (hloa(nlomax, apwordmax, 0:input%groundstate%lmaxmat, &
     & lmmaxvr, natmtot))
      
      lolmmax=25
      nlorb(1)=4
      lorbl(1:4,1)=(/1,0,1,2/)
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

      Deallocate(hloa)
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
! test testHmlalon_SpherSymmAsymGnt_Serial
!------------------------------------------------------------------------------
! 1st test, serial
! The purpose is to test whether gaunt coefficients are handled properly.
! The matching coefficients apwalm contain a constant.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmlalon_Gnt_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)
      
      Integer :: l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,g1
      integer :: i,ilo
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8)               :: test
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      hloa(:,:,:,:,:)=1d0

      Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      Do l1=0,lmaxapw
         Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
               apwalm(g1,:,lm1,1)=cmplx(1d0,0,8)
            Enddo
         Enddo
      Enddo

      i=0
      Do ilo= 1,nlorb (1)
         l1 = lorbl (ilo, 1)
         Do m1 = - l1, l1
            i=i+1
            lm1 = idxlm (l1, m1)
            test=cmplx(0,0,8)
            Do l2 = 0, input%groundstate%lmaxvr
               Do m2 = - l2, l2
               lm2 = idxlm (l2, m2)
                  Do l3 = 0, input%groundstate%lmaxmat
                     Do m3 = - l3, l3
                        lm3 = idxlm (l3, m3)
                        test = test + gauntyry (l1, l2, l3, m1, m2, m3)
                     End Do
                  End Do
               End Do
            End Do
            hamilton_ref%za(1:gsize,i+gsize)=conjg(test)
         End Do
      End Do

      Call hmlalon(hamilton,1,1,gsize,apwalm)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmlalon_Gnt_Serial

!------------------------------------------------------------------------------
! test testHmlalon_SphSymLO_Serial
!------------------------------------------------------------------------------
! 2nd test, serial
! The purpose is to test whether whether hloa is handled properly.
! This test involves only the spherically symmetric component, meaning only a part of hloa is tested.
! The matching coefficients apwalm depend on g1.
! The radial integrals (hloa) depend on the LO index and angular momentum of APW.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmlalon_SphSymLO_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)

      Integer :: l1,m1,lm1,lm2,l3,g1
      integer :: ilo
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      hloa(:,:,:,:,:)=1d8
! The condition for the spherically symmetric evaluation
      hloa(:,1,:,2:lmmaxvr,1)=0d0

      Do ilo= 1,nlorb (1)
         l1=lorbl(ilo,1)
         Do l3 = 0, input%groundstate%lmaxmat
            hloa(ilo,:,l3,1,1)=dble(ilo)*sqrt(dble(l1+1))/sqrt(dble(l3+1))*sqrt(4d0*pi)
         End Do
      End Do

      Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
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
            hamilton_ref%za(g1,gsize+idxlo(lm1,ilo,1):gsize+idxlo(lm2,ilo,1))=g1*ilo
         enddo
      enddo

      Call hmlalon(hamilton,1,1,gsize,apwalm)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testHmlalon_SphSymLO_Serial



!------------------------------------------------------------------------------
! test testHmlalon_SphSymAPW_Serial
!------------------------------------------------------------------------------
! 3rd test, serial
! The purpose is to test whether whether apwalm is handled properly.
! This test involves only the spherically symmetric component, meaning only a part of hloa is tested.
! The matching coefficients apwalm depend on g1, angular momentum of APW and LO index.
! The radial integrals (hloa) depend on angular momentum of APW.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
! IMPORTANT NOTE: This test depends on initialisation done in a specific way.
! Make sure that lorbl contains either (/0/), (/0,1/), (/0,1,2/), (/0,1,2,3/) etc.

    Subroutine testHmlalon_SphSymAPW_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=43)

      Integer :: l1,m1,lm1,l3,g1,g2
      integer :: ilo,last,i,l,m,lm
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

      nlorb(1)=4
      lorbl(1:4,1)=(/0,1,2,3/)
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

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      hloa(:,:,:,:,:)=1d8
      hloa(:,:,:,2:lmmaxvr,:)=0d0
      Do ilo= 1,nlorb (1)
         l1=lorbl(ilo,1)
         Do l3 = 0, input%groundstate%lmaxmat
            hloa(ilo,:,l3,1,1)=1d0/sqrt(dble(l3+1))*sqrt(4d0*pi)
         End Do
      End Do

      Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
      Do l1=0,lmaxapw
         Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
               ilo=l1**2+l1+m1+1
               apwalm(g1,:,lm1,1)=cmplx(g1,0,8)*cmplx(sqrt(dble(l1+1)),0,8)*dble(ilo)
            Enddo
         Enddo
      Enddo

      last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
      do g1=1,gsize
         do g2=1,last
            hamilton_ref%za(g1,gsize+g2)=g1*g2
         enddo
      enddo

      Call hmlalon(hamilton,1,1,gsize,apwalm)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testHmlalon_SphSymAPW_Serial



!------------------------------------------------------------------------------
! test testHmlalon_SphSymAsymAPW_Serial
!------------------------------------------------------------------------------
! 4th test, serial
! The purpose is to test whether whether apwalm is handled properly.
! This test involves both the spherically symmetric and asymmetric components.
! This test is a more general version of the 3rd test, but it is also more difficult to follow if something goes wrong here.
! The matching coefficients apwalm depend on depend on g1, the orbital momentum l3 and magnetic moment m3.
! The radial integrals (hloa) are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmlalon_SphSymAsymAPW_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)

      Integer :: l1,m1,lm1,lm2,g1,g2
      integer :: ilo,last
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      hloa(:,:,:,:,:)=1d8
      hloa(:,1,:,:,1)=1d0

      Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
      Do l1=0,lmaxapw
         Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
               apwalm(g1,:,lm1,1)=cmplx(g1,0,8)*cmplx(cos(dble(lm1)),sin(dble(lm1)),8)
            Enddo
         Enddo
      Enddo

      last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
      do ilo=1,nlorb(1)
         l1=lorbl(ilo,1)
         do m1=-l1,l1
            lm1=idxlm(l1,m1)
            g2=idxlo(lm1,ilo,1)
            do lm2=1,lmmaxmat
               hamilton_ref%za(1:gsize,gsize+g2)=conjg(sum(gntyry(lm1,1:lmmaxvr,lm2))*cmplx(cos(dble(lm2)),sin(dble(lm2)),8))+hamilton_ref%za(1:gsize,gsize+g2)
            enddo
         enddo
      enddo
      do g1=1,gsize
         hamilton_ref%za(g1,gsize+1:gsize+last)=hamilton_ref%za(g1,gsize+1:gsize+last)*dble(g1)
      enddo

      Call hmlalon(hamilton,1,1,gsize,apwalm)
      
      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testHmlalon_SphSymAsymAPW_Serial


!------------------------------------------------------------------------------
! test testHmlalon_SphSymAsymLO_Serial
!------------------------------------------------------------------------------
! 5th test, serial
! The purpose is to test whether whether hloa is handled properly.
! This test involves both the spherically symmetric and asymmetric components.
! This test is a more general version of the 2nd test, but it is also more difficult to follow if something goes wrong here.
! The matching coefficients apwalm depend on depend on g1.
! The radial integrals (hloa) depend on the LO index (ilo), the orbital momentum of APW (l3) and the lm couple of the potential (lm2).
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmlalon_SphSymAsymLO_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)

      Integer :: l1,m1,lm1,l2,lm2,l3,lm3,g1,g2
      integer :: ilo,last
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      hloa(:,:,:,:,:)=1d8
      hloa(:,1,:,:,1)=1d0

      Do ilo= 1,nlorb (1)
         l1=lorbl(ilo,1)
         Do l3 = 0, input%groundstate%lmaxmat
            do lm2=1,lmmaxvr
               hloa(ilo,:,l3,lm2,1)=dble(ilo)*cos(dble(l3))*sin(dble(lm2))
            enddo
         End Do
      End Do

      Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
      Do g1=1,gsize
         apwalm(g1,1,:,1)=cmplx(g1,0,8)
      Enddo

      last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
      do ilo=1,nlorb(1)
         l1=lorbl(ilo,1)
         do m1=-l1,l1
            lm1=idxlm(l1,m1)
            g2=idxlo(lm1,ilo,1)
            do l2=0,lmaxmat
               do lm3=1,lmmaxvr
                  hamilton_ref%za(1:gsize,gsize+g2)=conjg(sum(gntyry(lm1,lm3,idxlm(l2,-l2):idxlm(l2,l2))))*dble(ilo)*cos(dble(l2))*sin(dble(lm3))+hamilton_ref%za(1:gsize,gsize+g2)
               enddo
            enddo
         enddo
      enddo
      do g1=1,gsize
        hamilton_ref%za(g1,gsize+1:gsize+last)=hamilton_ref%za(g1,gsize+1:gsize+last)*dble(g1)
      enddo

      Call hmlalon(hamilton,1,1,gsize,apwalm)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testHmlalon_SphSymAsymLO_Serial


!------------------------------------------------------------------------------
! test testHmlalon_SpherSymmAsymGnt_1Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 1 proc
! The purpose is to test whether gaunt coefficients are handled properly.
! The matching coefficients apwalm contain a constant.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlalon_Gnt_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)
      
      Integer :: l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,g1
      Integer :: i,ilo
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Complex(8)               :: test
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Integer, Dimension(nmatp) :: hamilton_loc_cols

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

         Call newmatrix(hamilton,nmatp)
         Call newmatrix(hamilton_ref,nmatp)
         hamilton_loc_cols = (/(i,i=1,nmatp)/)
! initialisation is finished

         hloa(:,:,:,:,:)=1d0

         Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1=1,gsize
                  apwalm(g1,:,lm1,1)=cmplx(1d0,0,8)
               Enddo
            Enddo
         Enddo

         i=0
         Do ilo= 1,nlorb (1)
            l1 = lorbl (ilo, 1)
            Do m1 = - l1, l1
               i=i+1
               lm1 = idxlm (l1, m1)
               test=cmplx(0,0,8)
               Do l2 = 0, input%groundstate%lmaxvr
                  Do m2 = - l2, l2
                     lm2 = idxlm (l2, m2)
                     Do l3 = 0, input%groundstate%lmaxmat
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           test = test + gauntyry (l1, l2, l3, m1, m2, m3)
                        End Do
                     End Do
                  End Do
               End Do
               hamilton_ref%za(1:gsize,i+gsize)=conjg(test)
            End Do
         End Do

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize,gsize,hamilton_loc_cols)

         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm)
! deallocation of global variables   
         Call freeGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_Gnt_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlalon_SphSymLO_1Proc
!------------------------------------------------------------------------------
! 2nd test, MPI with 1 proc
! The purpose is to test whether whether hloa is handled properly.
! This test involves only the spherically symmetric component, meaning only a part of hloa is tested.
! The matching coefficients apwalm depend on g1.
! The radial integrals (hloa) depend on the LO index and angular momentum of APW.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlalon_SphSymLO_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)

      Integer :: l1,m1,lm1,lm2,l3,g1,i
      integer :: ilo
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Integer, Dimension(nmatp) :: hamilton_loc_cols

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

         Call newmatrix(hamilton,nmatp)
         Call newmatrix(hamilton_ref,nmatp)
         hamilton_loc_cols = (/(i,i=1,nmatp)/)
! initialisation is finished

         hloa(:,:,:,:,:)=1d8
! The condition for the spherically symmetric evaluation
         hloa(:,1,:,2:lmmaxvr,1)=0d0

         Do ilo= 1,nlorb (1)
            l1=lorbl(ilo,1)
            Do l3 = 0, input%groundstate%lmaxmat
               hloa(ilo,:,l3,1,1)=dble(ilo)*sqrt(dble(l1+1))/sqrt(dble(l3+1))*sqrt(4d0*pi)
            End Do
         End Do

         Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
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
               hamilton_ref%za(g1,gsize+idxlo(lm1,ilo,1):gsize+idxlo(lm2,ilo,1))=g1*ilo
            enddo
         enddo

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize,gsize,hamilton_loc_cols)

         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_SphSymLO_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlalon_SphSymAPW_1Proc
!------------------------------------------------------------------------------
! 3rd test, MPI with 1 proc
! The purpose is to test whether whether apwalm is handled properly.
! This test involves only the spherically symmetric component, meaning only a part of hloa is tested.
! The matching coefficients apwalm depend on g1, angular momentum of APW and LO index.
! The radial integrals (hloa) depend on angular momentum of APW.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
! IMPORTANT NOTE: This test depends on initialisation done in a specific way.
! Make sure that lorbl contains either (/0/), (/0,1/), (/0,1,2/), (/0,1,2,3/) etc.
#ifdef MPI
    Subroutine testHmlalon_SphSymAPW_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=43)
      Integer :: l1,m1,lm1,l3,g1,g2
      integer :: ilo,last,i,l,m,lm
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Integer, Dimension(nmatp) :: hamilton_loc_cols

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

         nlorb(1)=4
         lorbl(1:4,1)=(/0,1,2,3/)
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

! allocate and generate complex Gaunt coefficient array
         Call initGntyry

         Call newmatrix(hamilton,nmatp)
         Call newmatrix(hamilton_ref,nmatp)
         hamilton_loc_cols = (/(i,i=1,nmatp)/)
! initialisation is finished

         hloa(:,:,:,:,:)=1d8
         hloa(:,:,:,2:lmmaxvr,:)=0d0
         Do ilo= 1,nlorb (1)
            l1=lorbl(ilo,1)
            Do l3 = 0, input%groundstate%lmaxmat
               hloa(ilo,:,l3,1,1)=1d0/sqrt(dble(l3+1))*sqrt(4d0*pi)
            End Do
         End Do

         Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1=1,gsize
                  ilo=l1**2+l1+m1+1
                  apwalm(g1,:,lm1,1)=cmplx(g1,0,8)*cmplx(sqrt(dble(l1+1)),0,8)*dble(ilo)
               Enddo
            Enddo
         Enddo

         last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
         do g1=1,gsize
            do g2=1,last
               hamilton_ref%za(g1,gsize+g2)=g1*g2
            enddo
         enddo

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize,gsize,hamilton_loc_cols)

         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_SphSymAPW_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlalon_SphSymAsymAPW_1Proc
!------------------------------------------------------------------------------
! 4th test, MPI with 1 proc
! The purpose is to test whether whether apwalm is handled properly.
! This test involves both the spherically symmetric and asymmetric components.
! This test is a more general version of the 3rd test, but it is also more difficult to follow if something goes wrong here.
! The matching coefficients apwalm depend on depend on g1, the orbital momentum l3 and magnetic moment m3.
! The radial integrals (hloa) are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlalon_SphSymAsymAPW_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)

      Integer :: l1,m1,lm1,lm2,g1,g2,i
      integer :: ilo,last
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Integer, Dimension(nmatp) :: hamilton_loc_cols

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

         Call newmatrix(hamilton,nmatp)
         Call newmatrix(hamilton_ref,nmatp)
         hamilton_loc_cols = (/(i,i=1,nmatp)/)
! initialisation is finished

         hloa(:,:,:,:,:)=1d8
         hloa(:,1,:,:,1)=1d0

         Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1=1,gsize
                  apwalm(g1,:,lm1,1)=cmplx(g1,0,8)*cmplx(cos(dble(lm1)),sin(dble(lm1)),8)
               Enddo
            Enddo
         Enddo

         last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
         do ilo=1,nlorb(1)
            l1=lorbl(ilo,1)
            do m1=-l1,l1
               lm1=idxlm(l1,m1)
               g2=idxlo(lm1,ilo,1)
               do lm2=1,lmmaxmat
                  hamilton_ref%za(1:gsize,gsize+g2)=conjg(sum(gntyry(lm1,1:lmmaxvr,lm2))*cmplx(cos(dble(lm2)),sin(dble(lm2)),8))+hamilton_ref%za(1:gsize,gsize+g2)
               enddo
            enddo
         enddo
         do g1=1,gsize
            hamilton_ref%za(g1,gsize+1:gsize+last)=hamilton_ref%za(g1,gsize+1:gsize+last)*dble(g1)
         enddo

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize,gsize,hamilton_loc_cols)
         
         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_SphSymAsymAPW_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlalon_SphSymAsymLO_4Proc
!------------------------------------------------------------------------------
! 5th test, MPI with 1 proc
! The purpose is to test whether whether hloa is handled properly.
! This test involves both the spherically symmetric and asymmetric components.
! This test is a more general version of the 2nd test, but it is also more difficult to follow if something goes wrong here.
! The matching coefficients apwalm depend on depend on g1.
! The radial integrals (hloa) depend on the LO index (ilo), the orbital momentum of APW (l3) and the lm couple of the potential (lm2).
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlalon_SphSymAsymLO_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)

      Integer :: l1,m1,lm1,l2,lm2,l3,lm3,g1,g2,i
      integer :: ilo,last
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Integer, Dimension(nmatp) :: hamilton_loc_cols

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

         Call newmatrix(hamilton,nmatp)
         Call newmatrix(hamilton_ref,nmatp)
         hamilton_loc_cols = (/(i,i=1,nmatp)/)
! initialisation is finished

         hloa(:,:,:,:,:)=1d8
         hloa(:,1,:,:,1)=1d0

         Do ilo= 1,nlorb (1)
            l1=lorbl(ilo,1)
            Do l3 = 0, input%groundstate%lmaxmat
               do lm2=1,lmmaxvr
                  hloa(ilo,:,l3,lm2,1)=dble(ilo)*cos(dble(l3))*sin(dble(lm2))
               enddo
            End Do
         End Do

         Allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))
         Do g1=1,gsize
         apwalm(g1,1,:,1)=cmplx(g1,0,8)
         Enddo

         last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
         do ilo=1,nlorb(1)
            l1=lorbl(ilo,1)
            do m1=-l1,l1
               lm1=idxlm(l1,m1)
               g2=idxlo(lm1,ilo,1)
               do l2=0,lmaxmat
                  do lm3=1,lmmaxvr
                  hamilton_ref%za(1:gsize,gsize+g2)=conjg(sum(gntyry(lm1,lm3,idxlm(l2,-l2):idxlm(l2,l2))))*dble(ilo)*cos(dble(l2))*sin(dble(lm3))+hamilton_ref%za(1:gsize,gsize+g2)
                  enddo
               enddo
            enddo
         enddo
         do g1=1,gsize
         hamilton_ref%za(g1,gsize+1:gsize+last)=hamilton_ref%za(g1,gsize+1:gsize+last)*dble(g1)
         enddo

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize,gsize,hamilton_loc_cols)

         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_SphSymAsymLO_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlalon_SpherSymmAsymGnt_4Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 4 procs
! The purpose is to test whether gaunt coefficients are handled properly.
! The matching coefficients apwalm contain a constant.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlalon_Gnt_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23)
      
      Integer :: l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,g1_loc
      integer :: i,ilo
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8)               :: test
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! Externals
      Integer,    External :: NUMROC
      Complex(8), External :: gauntyry

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, gsize_ncols_loc, ncols_loc
      Integer :: gsize_nrows_loc, nmat_ncols_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc2glob, hamilton_loc_cols

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

         Call newmatrix(hamilton, nmatp, DISTRIBUTE_2D)
         Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_2D)

! init datastructures for splitting apwalm
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         gsize_ncols_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         nmat_ncols_loc  = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         Allocate(apwalm1_loc2glob(gsize_nrows_loc), hamilton_loc_cols(nmat_ncols_loc))
         Call getLocalIndices(gsize, nmatp, apwalm1_loc2glob, hamilton_loc_cols, MPIglobal)
! initialisation is finished

         hloa(:,:,:,:,:)=1d0

         hamilton_ref_global(:,:)=cmplx(0,0,8)

         Allocate(apwalm (gsize_nrows_loc, apwordmax, lmmaxapw, natmtot))      
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1_loc=1,gsize_nrows_loc !splitting of apwalm along first dimension
                                           ! according to the processor rows
                  apwalm(g1_loc,:,lm1,1)=cmplx(1d0,0,8)
               Enddo
            Enddo
         Enddo

         hamilton_ref_global = Cmplx(0,0,8)
         i=0
         Do ilo= 1,nlorb (1)
            l1 = lorbl (ilo, 1)
            Do m1 = - l1, l1
               i=i+1
               lm1 = idxlm (l1, m1)
               test=cmplx(0,0,8)
               Do l2 = 0, input%groundstate%lmaxvr
                  Do m2 = - l2, l2
                     lm2 = idxlm (l2, m2)
                     Do l3 = 0, input%groundstate%lmaxmat
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           test = test + gauntyry (l1, l2, l3, m1, m2, m3)
                        End Do
                     End Do
                  End Do
               End Do
               hamilton_ref_global(1:gsize,i+gsize)=conjg(test)
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

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize_nrows_loc,gsize_ncols_loc,hamilton_loc_cols)

         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm, hamilton_loc_cols, apwalm1_loc2glob)
! deallocation of global variables   
         Call freeGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_Gnt_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmlalon_SphSymLO_4Proc
!------------------------------------------------------------------------------
! 2nd test, MPI with 4 procs
! The purpose is to test whether whether hloa is handled properly.
! This test involves only the spherically symmetric component, meaning only a part of hloa is tested.
! The matching coefficients apwalm depend on g1.
! The radial integrals (hloa) depend on the LO index and angular momentum of APW.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlalon_SphSymLO_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23) 

      Integer :: l1,m1,lm1,lm2,l3,g1,g1_loc
      integer :: ilo
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! Externals
      Integer,    External :: NUMROC
      Complex(8), External :: gauntyry

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, gsize_ncols_loc, ncols_loc
      Integer :: gsize_nrows_loc, nmat_ncols_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc2glob, hamilton_loc_cols

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

         Call newmatrix(hamilton, nmatp, DISTRIBUTE_2D)
         Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_2D)

! init datastructures for splitting apwalm
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         gsize_ncols_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         nmat_ncols_loc  = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         Allocate(apwalm1_loc2glob(gsize_nrows_loc), hamilton_loc_cols(nmat_ncols_loc))
         Call getLocalIndices(gsize, nmatp, apwalm1_loc2glob, hamilton_loc_cols, MPIglobal)
! initialisation is finished

         hloa(:,:,:,:,:)=1d8
! The condition for the spherically symmetric evaluation
         hloa(:,1,:,2:lmmaxvr,1)=0d0

         Do ilo= 1,nlorb (1)
            l1=lorbl(ilo,1)
            Do l3 = 0, input%groundstate%lmaxmat
               hloa(ilo,:,l3,1,1)=dble(ilo)*sqrt(dble(l1+1))/sqrt(dble(l3+1))*sqrt(4d0*pi)
            End Do
         End Do

         Allocate(apwalm (gsize_nrows_loc, apwordmax, lmmaxapw, natmtot))
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

         hamilton_ref_global = Cmplx(0,0,8)
         do g1=1,gsize
            do ilo=1,nlorb (1)
               l1=lorbl(ilo,1)
               lm1=idxlm(l1,-l1)
               lm2=idxlm(l1,l1)
               hamilton_ref_global(g1,gsize+idxlo(lm1,ilo,1):gsize+idxlo(lm2,ilo,1))=g1*ilo
            enddo
         enddo
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

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize_nrows_loc,gsize_ncols_loc,hamilton_loc_cols)

         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm, hamilton_loc_cols, apwalm1_loc2glob)
! deallocation of global variables   
         Call freeGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_SphSymLO_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmlalon_SphSymAPW_4Proc
!------------------------------------------------------------------------------
! 3rd test, MPI with 4 procs
! The purpose is to test whether whether apwalm is handled properly.
! This test involves only the spherically symmetric component, meaning only a part of hloa is tested.
! The matching coefficients apwalm depend on g1, angular momentum of APW and LO index.
! The radial integrals (hloa) depend on angular momentum of APW.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
! IMPORTANT NOTE: This test depends on initialisation done in a specific way.
! Make sure that lorbl contains either (/0/), (/0,1/), (/0,1,2/), (/0,1,2,3/) etc.
#ifdef MPI
    Subroutine testHmlalon_SphSymAPW_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=43)

      Integer :: l1,m1,lm1,l3,g1,g1_loc,g2
      integer :: ilo, last,i,l,m,lm
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! Externals
      Integer,    External :: NUMROC
      Complex(8), External :: gauntyry

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, gsize_ncols_loc, ncols_loc
      Integer :: gsize_nrows_loc, nmat_ncols_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc2glob, hamilton_loc_cols

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

! initialisation of global variables
         Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

         nlorb(1)=4
         lorbl(1:4,1)=(/0,1,2,3/)
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

! allocate and generate complex Gaunt coefficient array
         Call initGntyry

         Call newmatrix(hamilton, nmatp, DISTRIBUTE_2D)
         Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_2D)

! init datastructures for splitting apwalm
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         gsize_ncols_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         nmat_ncols_loc  = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         Allocate(apwalm1_loc2glob(gsize_nrows_loc), hamilton_loc_cols(nmat_ncols_loc))
         Call getLocalIndices(gsize, nmatp, apwalm1_loc2glob, hamilton_loc_cols, MPIglobal)
! initialisation is finished

         hloa(:,:,:,:,:)=1d8
         hloa(:,:,:,2:lmmaxvr,:)=0d0
         Do ilo= 1,nlorb (1)
            l1=lorbl(ilo,1)
            Do l3 = 0, input%groundstate%lmaxmat
               hloa(ilo,:,l3,1,1)=1d0/sqrt(dble(l3+1))*sqrt(4d0*pi)
            End Do
         End Do

         Allocate(apwalm(gsize_nrows_loc, apwordmax, lmmaxapw, natmtot))
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1_loc=1,gsize_nrows_loc !splitting of apwalm along first dimension
                                          ! according to the processor rows
                  g1=apwalm1_loc2glob(g1_loc) 
                  ilo=l1**2+l1+m1+1
                  apwalm(g1_loc,:,lm1,1)=cmplx(g1,0,8)*cmplx(sqrt(dble(l1+1)),0,8)*dble(ilo)
               Enddo
            Enddo
         Enddo

         hamilton_ref_global = Cmplx(0,0,8)
         last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
         do g1=1,gsize
            do g2=1,last
               hamilton_ref_global(g1,gsize+g2)=g1*g2
            enddo
         enddo
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

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize_nrows_loc,gsize_ncols_loc,hamilton_loc_cols)

         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm, hamilton_loc_cols, apwalm1_loc2glob)
! deallocation of global variables   
         Call freeGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_SphSymAPW_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmlalon_SphSymAsymAPW_4Proc
!------------------------------------------------------------------------------
! 4th test, MPI with 4 procs
! The purpose is to test whether whether apwalm is handled properly.
! This test involves both the spherically symmetric and asymmetric components.
! This test is a more general version of the 3rd test, but it is also more difficult to follow if something goes wrong here.
! The matching coefficients apwalm depend on depend on g1, the orbital momentum l3 and magnetic moment m3.
! The radial integrals (hloa) are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlalon_SphSymAsymAPW_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23) 

      Integer :: l1,m1,lm1,lm2,g1,g1_loc,g2
      integer :: ilo, last
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! Externals
      Integer,    External :: NUMROC
      Complex(8), External :: gauntyry

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, gsize_ncols_loc, ncols_loc
      Integer :: gsize_nrows_loc, nmat_ncols_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc2glob, hamilton_loc_cols

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

         Call newmatrix(hamilton, nmatp, DISTRIBUTE_2D)
         Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_2D)

! init datastructures for splitting apwalm
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         gsize_ncols_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         nmat_ncols_loc  = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         Allocate(apwalm1_loc2glob(gsize_nrows_loc), hamilton_loc_cols(nmat_ncols_loc))
         Call getLocalIndices(gsize, nmatp, apwalm1_loc2glob, hamilton_loc_cols, MPIglobal)
! initialisation is finished

         hloa(:,:,:,:,:)=1d8
         hloa(:,1,:,:,1)=1d0

         Allocate(apwalm(gsize_nrows_loc, apwordmax, lmmaxapw, natmtot))
         Do l1=0,lmaxapw
            Do m1=-l1,l1
               lm1=idxlm(l1,m1)
               Do g1_loc=1,gsize_nrows_loc !splitting of apwalm along first dimension
                                          ! according to the processor rows
                  g1=apwalm1_loc2glob(g1_loc) 
                  apwalm(g1_loc,:,lm1,1)=cmplx(g1,0,8)*cmplx(cos(dble(lm1)),sin(dble(lm1)),8)
               Enddo
            Enddo
         Enddo

         hamilton_ref_global = Cmplx(0,0,8)
         last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
         do ilo=1,nlorb(1)
            l1=lorbl(ilo,1)
            do m1=-l1,l1
               lm1=idxlm(l1,m1)
               g2=idxlo(lm1,ilo,1)
               do lm2=1,lmmaxmat
                  hamilton_ref_global(1:gsize,gsize+g2)=conjg(sum(gntyry(lm1,1:lmmaxvr,lm2))*cmplx(cos(dble(lm2)),sin(dble(lm2)),8))+hamilton_ref_global(1:gsize,gsize+g2)
               enddo
            enddo
         enddo
         do g1=1,gsize
            hamilton_ref_global(g1,gsize+1:gsize+last)=hamilton_ref_global(g1,gsize+1:gsize+last)*dble(g1)
         enddo
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

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize_nrows_loc,gsize_ncols_loc,hamilton_loc_cols)

         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm, hamilton_loc_cols, apwalm1_loc2glob)
! deallocation of global variables   
         Call freeGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_SphSymAsymAPW_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmlalon_SphSymAsymLO_4Proc
!------------------------------------------------------------------------------
! 5th test, MPI with 4 procs
! The purpose is to test whether whether hloa is handled properly.
! This test involves both the spherically symmetric and asymmetric components.
! This test is a more general version of the 2nd test, but it is also more difficult to follow if something goes wrong here.
! The matching coefficients apwalm depend on depend on g1.
! The radial integrals (hloa) depend on the LO index (ilo), the orbital momentum of APW (l3) and the lm couple of the potential (lm2).
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlalon_SphSymAsymLO_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=23) 

      Integer :: l1,l2,m1,lm1,lm2,lm3,l3,g1,g1_loc,g2
      integer :: ilo, last
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! Externals
      Integer,    External :: NUMROC
      Complex(8), External :: gauntyry

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, gsize_ncols_loc, ncols_loc
      Integer :: gsize_nrows_loc, nmat_ncols_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc2glob, hamilton_loc_cols

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

         Call newmatrix(hamilton, nmatp, DISTRIBUTE_2D)
         Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_2D)

! init datastructures for splitting apwalm
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         gsize_ncols_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         nmat_ncols_loc  = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         Allocate(apwalm1_loc2glob(gsize_nrows_loc), hamilton_loc_cols(nmat_ncols_loc))
         Call getLocalIndices(gsize, nmatp, apwalm1_loc2glob, hamilton_loc_cols, MPIglobal)
! initialisation is finished

         hloa(:,:,:,:,:)=1d8
         hloa(:,1,:,:,1)=1d0

         Do ilo= 1,nlorb (1)
            l1=lorbl(ilo,1)
            Do l3 = 0, input%groundstate%lmaxmat
               do lm2=1,lmmaxvr
                  hloa(ilo,:,l3,lm2,1)=dble(ilo)*cos(dble(l3))*sin(dble(lm2))
               enddo
            End Do
         End Do

         Allocate(apwalm(gsize_nrows_loc, apwordmax, lmmaxapw, natmtot))
         Do g1_loc=1,gsize_nrows_loc !splitting of apwalm along first dimension
                                    ! according to the processor rows
            g1=apwalm1_loc2glob(g1_loc) 
            apwalm(g1_loc,1,:,1)=cmplx(g1,0,8)
         Enddo

         hamilton_ref_global = Cmplx(0,0,8)
         last=idxlo(idxlm(lorbl(nlorb (1),1),lorbl(nlorb (1),1)),nlorb (1),1)
         do ilo=1,nlorb(1)
            l1=lorbl(ilo,1)
            do m1=-l1,l1
               lm1=idxlm(l1,m1)
               g2=idxlo(lm1,ilo,1)
               do l2=0,lmaxmat
                  do lm3=1,lmmaxvr
                     hamilton_ref_global(1:gsize,gsize+g2)=conjg(sum(gntyry(lm1,lm3,idxlm(l2,-l2):idxlm(l2,l2))))*dble(ilo)*cos(dble(l2))*sin(dble(lm3))+hamilton_ref_global(1:gsize,gsize+g2)
                  enddo
               enddo
            enddo
         enddo
         do g1=1,gsize
            hamilton_ref_global(g1,gsize+1:gsize+last)=hamilton_ref_global(g1,gsize+1:gsize+last)*dble(g1)
         enddo
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

         Call hmlalon(hamilton,1,1,gsize,apwalm,gsize_nrows_loc,gsize_ncols_loc,hamilton_loc_cols)

         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(apwalm, hamilton_loc_cols, apwalm1_loc2glob)
! deallocation of global variables   
         Call freeGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
   End Subroutine testHmlalon_SphSymAsymLO_4Proc
#endif


end module modHmlalon_test
