! !MODULE:  modOlplolon_test
! !DESCRIPTION:
! Modules with test and helper routines for testing olplolon
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
module modOlplolon_test
    Use fruit
    Use test_helpers
    Use modinput,        Only: input
    Use mod_muffin_tin,  Only: idxlm, lmmaxapw,lmmaxmat
    Use mod_APW_LO,      Only: nlomax, lorbl, nlorb, lolmmax
    Use mod_Gkvector,    Only: ngkmax
    Use mod_atoms,       Only: natmtot
    Use mod_eigensystem, Only: ololo, idxlo
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix
#ifdef MPI
    use modmpi
#endif

    Double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains


    Subroutine testcaseOlplolonSerial
      Implicit None

      Call set_test_name ('Radial integrals for sparse overlap matrix, ololo')
      Call testOlplolon_LO_Serial(.true.)
      Call set_test_name ('Radial integrals for dense overlap matrix, ololo')
      Call testOlplolon_LO_Serial(.false.)



    End Subroutine testcaseOlplolonSerial


#ifdef MPI
    Subroutine testcaseOlplolon1Proc
      Implicit None

    End Subroutine testcaseOlplolon1Proc
#endif


#ifdef MPI
    Subroutine testcaseOlplolon4Proc
      Implicit None

    End Subroutine testcaseOlplolon4Proc
#endif


! initialisation of global variables
    Subroutine initGlobals(lmaxmat,lmaxapw,gsize,sparse)
      Implicit None

      Integer, Intent(in) :: lmaxmat,lmaxapw,gsize
      Logical, Intent(in) :: sparse

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

      nlomax=10

      If (allocated(ololo)) Deallocate (ololo)
      Allocate (ololo(nlomax, nlomax, natmtot))
      
      lolmmax=25
      if (sparse) then
        nlorb(1)=7
        lorbl(1:7,1)=(/1,0,1,0,0,2,0/)
      else
        nlorb(1)=10
        lorbl(1:10,1)=(/0,0,1,0,0,0,0,0,0,1/)
      endif
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

      Deallocate(ololo)
      Deallocate(idxlm,idxlo,input%groundstate)

    End Subroutine freeGlobals


!------------------------------------------------------------------------------
! test testHmllolon_SphSym_Serial
!------------------------------------------------------------------------------
! 2nd test, serial
! The purpose is to test whether the radial integrals are handled properly.
! The radial integrals hlolo depend on the local orbital index ilo.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testOlplolon_LO_Serial(sparse)

      Implicit None
      Logical, Intent(in) :: sparse
! Size of the tests
      Integer lmaxmat,lmaxapw,gsize,nmatp
      parameter (lmaxmat=5,lmaxapw=10,gsize=9,nmatp=43)

      Integer :: l1,m1,lm1,l2,m2,lm2
      integer :: ilo,ilo2
      Type (HermitianMatrix)   :: overlap,overlap_ref

! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,gsize,sparse)

      Call newmatrix(overlap,nmatp)
      Call newmatrix(overlap_ref,nmatp)
! initialisation is finished

      ololo(:,:,:)=0d0
      Do ilo=1,nlorb (1)
        Do ilo2=1,nlorb (1)
          ololo(ilo,ilo2,1)=dble(ilo*ilo2)
        EndDo
      End Do

      overlap%za(:,:)=cmplx(0,0,8)
      overlap_ref%za(:,:)=cmplx(0,0,8)

! The answer is a diagonal matrix      
      Do ilo=1,nlorb (1)
        l1=lorbl(ilo,1)
        Do m1=-l1,l1       
          lm1=idxlm(l1,m1)
          Do ilo2=1,nlorb (1)
            l2=lorbl(ilo2,1)
            Do m2=-l2,l2
              lm2=idxlm(l2,m2)
              if ((lm1.eq.lm2).and.(idxlo(lm2,ilo2,1).le.idxlo(lm1,ilo,1))) then
                 overlap_ref%za(gsize+idxlo(lm2,ilo2,1),gsize+idxlo(lm1,ilo,1))=cmplx(ilo*ilo2,0,8)
              endif
            EndDo
          EndDo
        EndDo
      End Do

      Call olplolon(overlap,1,1,gsize)

      Call assert_equals(nmatp, overlap%size, 'checking result rank')
      Call assert_equals(nmatp, size(overlap%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(overlap%za,2), 'checking result size cols')
      Call assert_equals(overlap_ref%za, overlap%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(overlap)
      Call deletematrix(overlap_ref)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testOlplolon_LO_Serial



end module modOlplolon_test
