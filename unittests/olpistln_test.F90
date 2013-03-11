! !MODULE:  modOlpistln_test
! !DESCRIPTION:
! Modules with test and helper routines for testing olpistln
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
module modOlpistln_test
    Use fruit
    Use test_helpers
    Use mod_Gkvector,    Only: ngkmax
    Use mod_Gvector,     Only: ivg,ngvec,cfunig,ngrtot,ivgig
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix
#ifdef MPI
    use modmpi
#endif

    Double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains


    Subroutine testcaseOlpistlnSerial
      Implicit None

      Call set_test_name ('Theta function / cfunig')
      Call testOlpistln_Theta_Serial

    End Subroutine testcaseOlpistlnSerial


#ifdef MPI
    Subroutine testcaseOlpistln1Proc
      Implicit None

    End Subroutine testcaseOlpistln1Proc
#endif


#ifdef MPI
    Subroutine testcaseOlpistln4Proc
      Implicit None

    End Subroutine testcaseOlpistln4Proc
#endif


! initialisation of global variables
    Subroutine initGlobals(maxg,gsize)
      Implicit None
      Integer, Intent(in) :: maxg,gsize
      Integer order1,order2
      Integer i,j,k

      ngkmax=gsize
      ngrtot=(2*maxg+1)**3
      if (allocated(ivgig)) deallocate(ivgig)
      allocate(ivgig(-maxg:maxg,-maxg:maxg,-maxg:maxg))
    
      if (allocated(ivg)) deallocate(ivg)
      allocate(ivg(3,ngrtot))

      order1=(2*maxg+1)
      order2=(2*maxg+1)**2
      Do k=-maxg,maxg
        Do j=-maxg,maxg
          Do i=-maxg,maxg
            ivgig(i,j,k)=1+maxg+i+order1*(maxg+j)+order2*(maxg+k)
            ivg(1,1+maxg+i+order1*(maxg+j)+order2*(maxg+k))=i
            ivg(2,1+maxg+i+order1*(maxg+j)+order2*(maxg+k))=j
            ivg(3,1+maxg+i+order1*(maxg+j)+order2*(maxg+k))=k
          EndDo
        EndDo
      EndDo

      ngvec=ngrtot

      if (allocated(cfunig)) deallocate(cfunig)
      allocate(cfunig(ngrtot))
      cfunig(:)=0d0

    End Subroutine initGlobals


! deallocation of global variables       
    Subroutine freeGlobals
      Implicit None
      deallocate(ivg,ivgig,cfunig)
    End Subroutine freeGlobals

!------------------------------------------------------------------------------
! test testOlpistln_Theta_Serial
!------------------------------------------------------------------------------
! 1st test, serial
! The purpose is to test whether the theta-function (cfunig) is handled properly.
! The theta-function (cfunig) depends on the plane-wave index.
    Subroutine testOlpistln_Theta_Serial

      Implicit None
! Size of the tests
      Integer gsize,nmatp,maxg,maxgk
      parameter (maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10)

      Integer :: igpig (gsize)
      Integer i,j,k
      Type (HermitianMatrix)   ::  overlap, overlap_ref

! initialisation of global variables
      Call initGlobals(maxg,gsize)

! initialisation of hmlistln parameters
      j=0
      Do i=1,ngrtot
        If ((ivg(1,i)**2.le.maxgk**2).and.(ivg(2,i)**2.le.maxgk**2).and.(ivg(3,i)**2.le.maxgk**2)) then
          j=j+1
          igpig(j)=i
        End If
      EndDo

! allocate and generate complex Gaunt coefficient array
      Call newmatrix( overlap,nmatp)
      Call newmatrix( overlap_ref,nmatp)
! initialisation is finished

      do i=1,ngrtot
        cfunig(i)=dble(i)
      enddo

! preparation of the correct answer
      do i=1,gsize
        do k=1,i
           overlap_ref%za(k,i)=ivgig(ivg(1,igpig(k))-ivg(1,igpig(i)),ivg(2,igpig(k))-ivg(2,igpig(i)),ivg(3,igpig(k))-ivg(3,igpig(i)))
        enddo
      enddo

      Call olpistln( overlap,gsize,igpig)

      Call assert_equals(gsize, j, 'Is the test itself set up properly?')
      Call assert_equals(nmatp,  overlap%size, 'checking result rank')
      Call assert_equals(nmatp, size( overlap%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size( overlap%za,2), 'checking result size cols')
      Call assert_equals( overlap_ref%za,  overlap%za, nmatp, nmatp, tol, 'checking result numbers')

! deallocation of global variables   
      Call freeGlobals
    End Subroutine testOlpistln_Theta_Serial

end module modOlpistln_test
