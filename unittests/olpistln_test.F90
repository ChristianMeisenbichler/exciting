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
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix,evsystem,newsystem,deletesystem
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

      Call set_test_name ('Theta function / cfunig')
      Call testOlpistln_Theta_1Proc

    End Subroutine testcaseOlpistln1Proc
#endif


#ifdef MPI
    Subroutine testcaseOlpistln4Proc
      Implicit None

      Call set_test_name ('Theta function / cfunig')
      Call testOlpistln_Theta_4Proc

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
      Type(evsystem)          :: system
      Type(HermitianMatrix)   :: overlap_ref

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
      Call newsystem(system, nmatp, gsize)
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

      Call olpistln(system,igpig)

      Call assert_equals(gsize, j, 'Is the test itself set up properly?')
      Call assert_equals(nmatp,  system%overlap%size, 'checking result rank')
      Call assert_equals(nmatp, size( system%overlap%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size( system%overlap%za,2), 'checking result size cols')
      Call assert_equals( overlap_ref%za,  system%overlap%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletesystem(system)
      Call deletematrix(overlap_ref)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testOlpistln_Theta_Serial


!------------------------------------------------------------------------------
! test testOlpistln_Theta_1Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 1 proc
! The purpose is to test whether the theta-function (cfunig) is handled properly.
! The theta-function (cfunig) depends on the plane-wave index.
#ifdef MPI
    Subroutine testOlpistln_Theta_1Proc

      Implicit None
! Size of the tests
      Integer gsize,nmatp,maxg,maxgk
      parameter (maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10)

      Integer :: igpig (gsize)
      Integer i,j,k
      Type(evsystem)          :: system
      Type(HermitianMatrix)   :: overlap_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIkpt_2D%comm, MPIkpt_2D%context, ierror_t)
      MPIkpt_2D%blocksize = 2

      If (MPIkpt_2D%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIkpt_2D)

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
         Call newsystem(system, nmatp, gsize)
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

         Call olpistln(system,igpig)

         Call assert_equals(gsize, j, 'Is the test itself set up properly?')
         Call assert_equals(nmatp,  system%overlap%size, 'checking result rank')
         Call assert_equals(nmatp, size( system%overlap%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size( system%overlap%za,2), 'checking result size cols')
         Call assert_equals( overlap_ref%za,  system%overlap%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
         Call deletesystem(system)
         Call deletematrix(overlap_ref)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIkpt_2D%comm, MPIkpt_2D%context, ierror_t)
      End If
    End Subroutine testOlpistln_Theta_1Proc
#endif


!------------------------------------------------------------------------------
! test testOlpistln_Theta_4Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 4 procs
! The purpose is to test whether the theta-function (cfunig) is handled properly.
! The theta-function (cfunig) depends on the plane-wave index.
#ifdef MPI
    Subroutine testOlpistln_Theta_4Proc

      Implicit None
! Size of the tests
      Integer gsize,nmatp,maxg,maxgk
      parameter (maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10) ! -> nmatp=37

      Integer :: igpig (gsize)
      Integer i,j,k
      Type(evsystem)          :: system
      Type(HermitianMatrix)   :: overlap_ref
      Complex (8), Dimension(nmatp,nmatp) :: overlap_ref_global

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, ncols_loc

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIkpt_2D%comm, MPIkpt_2D%context, ierror_t)
      MPIkpt_2D%blocksize = 2

      If (MPIkpt_2D%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIkpt_2D)

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
         Call newsystem(system, nmatp, gsize)
         Call newmatrix( overlap_ref,nmatp)
! initialisation is finished

         do i=1,ngrtot
            cfunig(i)=dble(i)
         enddo

! preparation of the correct answer
         overlap_ref_global = cmplx(0,0,8)
         do i=1,gsize
            do k=1,i
               overlap_ref_global(k,i)=ivgig(ivg(1,igpig(k))-ivg(1,igpig(i)),ivg(2,igpig(k))-ivg(2,igpig(i)),ivg(3,igpig(k))-ivg(3,igpig(i)))
            enddo
         enddo
         Call getBlockDistributedLoc(overlap_ref_global, overlap_ref%za, MPIkpt_2D)

         Select Case (MPIkpt_2D%myprocrow)
            Case (0)
               nrows_loc = 19
            Case (1)
               nrows_loc = 18
         End Select
         Select Case (MPIkpt_2D%myproccol)
            Case (0)
               ncols_loc = 19
            Case (1)
               ncols_loc = 18
         End Select

         Call olpistln(system,igpig)

         Call assert_equals(gsize, j, 'Is the test itself set up properly?')
         Call assert_equals(nmatp, system%overlap%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(system%overlap%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(system%overlap%za,2), 'checking result size cols')
         Call assert_equals(overlap_ref%za, system%overlap%za, nrows_loc, ncols_loc, tol, 'checking result numbers')
         Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking hamilton=0')

! finalisation
         Call deletesystem(system)
         Call deletematrix(overlap_ref)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIkpt_2D%comm, MPIkpt_2D%context, ierror_t)
      End If
    End Subroutine testOlpistln_Theta_4Proc
#endif

end module modOlpistln_test
