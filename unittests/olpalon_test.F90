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
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix
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

    End Subroutine testcaseOlpalon1Proc
#endif


#ifdef MPI
    Subroutine testcaseOlpalon4Proc
      Implicit None

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
      Type (HermitianMatrix)   :: overlap,overlap_ref

! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,gsize)

      Call newmatrix(overlap,nmatp)
      Call newmatrix(overlap_ref,nmatp)
! initialisation is finished

! The condition for the spherically symmetric evaluation
      oalo(:,:,:)=1d0

      overlap%za(:,:)=cmplx(0,0,8)
      overlap_ref%za(:,:)=cmplx(0,0,8)

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

      Call olpalon(overlap,1,1,gsize,apwalm)

      Call assert_equals(nmatp, overlap%size, 'checking result rank')
      Call assert_equals(nmatp, size(overlap%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(overlap%za,2), 'checking result size cols')
      Call assert_equals(overlap_ref%za, overlap%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(overlap)
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
      Type (HermitianMatrix)   :: overlap,overlap_ref

! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,gsize)

      Call newmatrix(overlap,nmatp)
      Call newmatrix(overlap_ref,nmatp)
! initialisation is finished

! The condition for the spherically symmetric evaluation
      oalo(:,:,:)=1d8
      oalo(:,:,1)=0d0

      overlap%za(:,:)=cmplx(0,0,8)
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

      Call olpalon(overlap,1,1,gsize,apwalm)

      Call assert_equals(nmatp, overlap%size, 'checking result rank')
      Call assert_equals(nmatp, size(overlap%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(overlap%za,2), 'checking result size cols')
      Call assert_equals(overlap_ref%za, overlap%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(overlap)
      Call deletematrix(overlap_ref)
      Deallocate(apwalm)
! deallocation of global variables 
      Call freeGlobals
    End Subroutine testOlpalon_APWLO_Serial

end module modOlpalon_test
