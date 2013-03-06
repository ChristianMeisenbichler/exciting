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

    End Subroutine testcaseOlpaan1Proc
#endif


#ifdef MPI
    Subroutine testcaseOlpaan4Proc
      Implicit None

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

! ! initialisation of global variables
      Call InitOverlapGlobals(lmaxmat,lmaxapw,gsize)

      Call newmatrix(overlap,nmatp)
      Call newmatrix(overlap_ref,nmatp)
! initialisation is finished


      overlap%za(:,:)=cmplx(0,0,8)
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
    End Subroutine testOlpaan_Serial

end module modOlpaan_test
