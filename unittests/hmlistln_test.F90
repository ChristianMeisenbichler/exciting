! !MODULE:  modHmlistln_test
! !DESCRIPTION:
! Modules with test and helper routines for testing hmlistln
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
module modHmlistln_test
    Use fruit
    Use test_helpers
    Use mod_Gkvector,    Only: ngkmax
    Use mod_Gvector,     Only: ivg,ngvec,cfunig,ngrtot,ivgig
    Use mod_potential_and_density, Only : veffig
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix
#ifdef MPI
    use modmpi
#endif

    Double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains


    Subroutine testcaseHmlistlnSerial
      Implicit None

! the matrix dimesion is equal to gsize+10 or equivalently to (2*maxgk+1)**3+10
      Call set_test_name ('Potential / veffig with matrix dimension 11')
      Call testHmlistln_Potential_Serial(0)
      Call set_test_name ('Potential / veffig with matrix dimension 37')
      Call testHmlistln_Potential_Serial(1)
      Call set_test_name ('Potential / veffig with matrix dimension 135')
      Call testHmlistln_Potential_Serial(2)
      Call set_test_name ('Potential / veffig with matrix dimension 353')
      Call testHmlistln_Potential_Serial(3)
      Call set_test_name ('Theta function / cfunig')
      Call testHmlistln_Theta_Serial
      Call set_test_name ('G-vector lengths / vgpc')
      Call testHmlistln_vgpc_Serial

    End Subroutine testcaseHmlistlnSerial


#ifdef MPI
    Subroutine testcaseHmlistln1Proc
      Implicit None

      Call set_test_name ('Potential / veffig with matrix dimension 11')
      Call testHmlistln_Potential_1Proc(0)
      Call set_test_name ('Potential / veffig with matrix dimension 37')
      Call testHmlistln_Potential_1Proc(1)
      Call set_test_name ('Potential / veffig with matrix dimension 135')
      Call testHmlistln_Potential_1Proc(2)
      Call set_test_name ('Potential / veffig with matrix dimension 353')
      Call testHmlistln_Potential_1Proc(3)
      Call set_test_name ('Theta function / cfunig')
      Call testHmlistln_Theta_1Proc
      Call set_test_name ('G-vector lengths / vgpc')
      Call testHmlistln_vgpc_1Proc

    End Subroutine testcaseHmlistln1Proc
#endif


#ifdef MPI
    Subroutine testcaseHmlistln4Proc
      Implicit None

      Call set_test_name ('Potential / veffig with matrix dimension 11')
      Call testHmlistln_Potential_4Proc(0)
      Call set_test_name ('Potential / veffig with matrix dimension 37')
      Call testHmlistln_Potential_4Proc(1)
      Call set_test_name ('Potential / veffig with matrix dimension 135')
      Call testHmlistln_Potential_4Proc(2)
      Call set_test_name ('Potential / veffig with matrix dimension 353')
      Call testHmlistln_Potential_4Proc(3)
      Call set_test_name ('Theta function / cfunig')
      Call testHmlistln_Theta_4Proc
      Call set_test_name ('G-vector lengths / vgpc')
      Call testHmlistln_vgpc_4Proc

    End Subroutine testcaseHmlistln4Proc
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

      if (allocated(veffig)) deallocate(veffig)
      allocate(veffig(ngrtot))
      veffig(:)=0d0
    End Subroutine initGlobals


! deallocation of global variables       
    Subroutine freeGlobals
      Implicit None
      deallocate(ivg,ivgig,cfunig,veffig)
    End Subroutine freeGlobals



!------------------------------------------------------------------------------
! test testHmlistln_Potential_Serial
!------------------------------------------------------------------------------
! 1st test, serial
! The purpose is to test whether the potential (veffig) is handled properly.
! The potential depends on the plane-wave index.
! The kinetic energy term (calculated using cfunig and vgpc) is zero.
    Subroutine testHmlistln_Potential_Serial(maxgk)

      Implicit None
! Size of the tests
! The fruit library happily handles maxgk=0,1,2,3. 
! Above that - hamilton gets too large apparently.
      Integer, Intent(in) :: maxgk
      Integer gsize,nmatp
      Integer maxg
      parameter (maxg=7)

      Integer :: igpig ((2*maxgk+1)**3)
      Real(8) :: vgpc(3,(2*maxgk+1)**3)
      Integer i,j,k
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

      gsize=(2*maxgk+1)**3
      nmatp=gsize+10
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
      vgpc(:,:)=0d0
       
      Do i=1,ngrtot
         veffig(i)=dble(i)
      EndDo      

! initialisation is finished

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)

! preparation of the correct answer
      do i=1,gsize
        do k=1,i
          hamilton_ref%za(k,i)=ivgig(ivg(1,igpig(k))-ivg(1,igpig(i)),ivg(2,igpig(k))-ivg(2,igpig(i)),ivg(3,igpig(k))-ivg(3,igpig(i)))
        enddo
      enddo

      Call hmlistln(hamilton,gsize,igpig,vgpc)

      Call assert_equals(gsize, j, 'Is the test itself set up properly?')
      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')
      
! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmlistln_Potential_Serial


!------------------------------------------------------------------------------
! test testHmlistln_Theta_Serial
!------------------------------------------------------------------------------
! 2st test, serial
! The purpose is to test whether the theta-function (cfunig) is handled properly.
! The potential is set to zero.
! The theta-function (cfunig) depends on index.
! The plane-wave lengths (vgpc) are constant.
    Subroutine testHmlistln_Theta_Serial

      Implicit None
! Size of the tests
      Integer gsize,nmatp,maxg,maxgk
      parameter (maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10)

      Integer :: igpig (gsize)
      Real(8) :: vgpc(3,gsize)
      Integer i,j,k
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

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
      vgpc(:,:)=sqrt(2d0/3d0)

! allocate and generate complex Gaunt coefficient array
      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      do i=1,ngrtot
        cfunig(i)=dble(i)
      enddo

! preparation of the correct answer
      do i=1,gsize
        do k=1,i
          hamilton_ref%za(k,i)=ivgig(ivg(1,igpig(k))-ivg(1,igpig(i)),ivg(2,igpig(k))-ivg(2,igpig(i)),ivg(3,igpig(k))-ivg(3,igpig(i)))
        enddo
      enddo

      Call hmlistln(hamilton,gsize,igpig,vgpc)

      Call assert_equals(gsize, j, 'Is the test itself set up properly?')
      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testHmlistln_Theta_Serial


!------------------------------------------------------------------------------
! test testHmlistln_vgpc_Serial
!------------------------------------------------------------------------------
! 3rd test, serial
! The purpose is to test whether the plane-wave lengths (vgpc) are handled properly.
! The potential is set to zero.
! The theta-function (cfunig) is constant.
! The plane-wave lengths (vgpc) depends on the plane-wave index.
    Subroutine testHmlistln_vgpc_Serial

      Implicit None
! Size of the tests
      Integer gsize,nmatp,maxg,maxgk
      parameter (maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10)

      Integer :: igpig (gsize)
      Real(8) :: vgpc(3,gsize)
      Integer i,j,k
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

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
      do i=1,gsize
        vgpc(:,i)=sqrt(1d0/7d0)*(/1d0,2d0,3d0/)*dble(i)
      enddo

! allocate and generate complex Gaunt coefficient array
      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      do i=1,ngrtot
        cfunig(i)=1d0
      enddo

! preparation of the correct answer
      do i=1,gsize
        do k=1,i
          hamilton_ref%za(k,i)=k*i
        enddo
      enddo

      Call hmlistln(hamilton,gsize,igpig,vgpc)

      Call assert_equals(gsize, j, 'Is the test itself set up properly?')
      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
! deallocation of global variables   
      Call freeGlobals
    End Subroutine testHmlistln_vgpc_Serial


!------------------------------------------------------------------------------
! test testHmlistln_Potential_1Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 1 proc
! The purpose is to test whether the potential (veffig) is handled properly.
! The potential depends on the plane-wave index.
! The kinetic energy term (calculated using cfunig and vgpc) is zero.
#ifdef MPI
    Subroutine testHmlistln_Potential_1Proc(maxgk)

      Implicit None
! Size of the tests
! The fruit library happily handles maxgk=0,1,2,3. 
! Above that - hamilton gets too large apparently.
      Integer, Intent(in) :: maxgk
      Integer gsize,nmatp
      Integer maxg
      parameter (maxg=7)

      Integer :: igpig ((2*maxgk+1)**3)
      Real(8) :: vgpc(3,(2*maxgk+1)**3)
      Integer i,j,k
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer, Dimension(:), Allocatable :: hamilton_loc_rows, hamilton_loc_cols

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

         gsize=(2*maxgk+1)**3
         nmatp=gsize+10
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
         vgpc(:,:)=0d0
         
         Do i=1,ngrtot
            veffig(i)=dble(i)
         EndDo      

! initialisation is finished

         Call newmatrix(hamilton,nmatp)
         Call newmatrix(hamilton_ref,nmatp)
         Allocate(hamilton_loc_rows(nmatp), hamilton_loc_cols(nmatp))
         hamilton_loc_rows = (/(i,i=1,nmatp)/)
         hamilton_loc_cols = (/(i,i=1,nmatp)/)

! preparation of the correct answer
         do i=1,gsize
            do k=1,i
               hamilton_ref%za(k,i)=ivgig(ivg(1,igpig(k))-ivg(1,igpig(i)),ivg(2,igpig(k))-ivg(2,igpig(i)),ivg(3,igpig(k))-ivg(3,igpig(i)))
            enddo
         enddo

         Call hmlistln(hamilton,gsize,igpig,vgpc,gsize,gsize,hamilton_loc_rows,hamilton_loc_cols)

         Call assert_equals(gsize, j, 'Is the test itself set up properly?')
         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')
      
! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(hamilton_loc_rows, hamilton_loc_cols)
! deallocation of global variables   
         Call freeGlobals    
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmlistln_Potential_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlistln_Theta_1Proc
!------------------------------------------------------------------------------
! 2st test, MPI with 1 proc
! The purpose is to test whether the theta-function (cfunig) is handled properly.
! The potential is set to zero.
! The theta-function (cfunig) depends on index.
! The plane-wave lengths (vgpc) are constant.
#ifdef MPI
    Subroutine testHmlistln_Theta_1Proc

      Implicit None
! Size of the tests
      Integer gsize,nmatp,maxg,maxgk
      parameter (maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10)

      Integer :: igpig (gsize)
      Real(8) :: vgpc(3,gsize)
      Integer i,j,k
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer, Dimension(:), Allocatable :: hamilton_loc_rows, hamilton_loc_cols

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

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
         vgpc(:,:)=sqrt(2d0/3d0)

         Call newmatrix(hamilton,nmatp)
         Call newmatrix(hamilton_ref,nmatp)
         Allocate(hamilton_loc_rows(nmatp), hamilton_loc_cols(nmatp))
         hamilton_loc_rows = (/(i,i=1,nmatp)/)
         hamilton_loc_cols = (/(i,i=1,nmatp)/)
! initialisation is finished

         do i=1,ngrtot
            cfunig(i)=dble(i)
         enddo

! preparation of the correct answer
         do i=1,gsize
            do k=1,i
               hamilton_ref%za(k,i)=ivgig(ivg(1,igpig(k))-ivg(1,igpig(i)),ivg(2,igpig(k))-ivg(2,igpig(i)),ivg(3,igpig(k))-ivg(3,igpig(i)))
            enddo
         enddo

         Call hmlistln(hamilton,gsize,igpig,vgpc,gsize,gsize,hamilton_loc_rows,hamilton_loc_cols)

         Call assert_equals(gsize, j, 'Is the test itself set up properly?')
         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(hamilton_loc_rows, hamilton_loc_cols)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmlistln_Theta_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlistln_vgpc_1Proc
!------------------------------------------------------------------------------
! 3rd test, MPI with 1 proc
! The purpose is to test whether the plane-wave lengths (vgpc) are handled properly.
! The potential is set to zero.
! The theta-function (cfunig) is constant.
! The plane-wave lengths (vgpc) depends on the plane-wave index.
#ifdef MPI
    Subroutine testHmlistln_vgpc_1Proc

      Implicit None
! Size of the tests
      Integer gsize,nmatp,maxg,maxgk
      parameter (maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10)

      Integer :: igpig (gsize)
      Real(8) :: vgpc(3,gsize)
      Integer i,j,k
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer, Dimension(:), Allocatable :: hamilton_loc_rows, hamilton_loc_cols

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

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
         do i=1,gsize
            vgpc(:,i)=sqrt(1d0/7d0)*(/1d0,2d0,3d0/)*dble(i)
         enddo

! allocate and generate complex Gaunt coefficient array
         Call newmatrix(hamilton,nmatp)
         Call newmatrix(hamilton_ref,nmatp)
         Allocate(hamilton_loc_rows(nmatp), hamilton_loc_cols(nmatp))
         hamilton_loc_rows = (/(i,i=1,nmatp)/)
         hamilton_loc_cols = (/(i,i=1,nmatp)/)
! initialisation is finished

         do i=1,ngrtot
         cfunig(i)=1d0
         enddo

! preparation of the correct answer
         do i=1,gsize
            do k=1,i
               hamilton_ref%za(k,i)=k*i
            enddo
         enddo

         Call hmlistln(hamilton,gsize,igpig,vgpc,gsize,gsize,hamilton_loc_rows,hamilton_loc_cols)

         Call assert_equals(gsize, j, 'Is the test itself set up properly?')
         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(hamilton_loc_rows, hamilton_loc_cols)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmlistln_vgpc_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlistln_Potential_4Proc
!------------------------------------------------------------------------------
! 1st test, MPI with 4 procs
! The purpose is to test whether the potential (veffig) is handled properly.
! The potential depends on the plane-wave index.
! The kinetic energy term (calculated using cfunig and vgpc) is zero.
#ifdef MPI
    Subroutine testHmlistln_Potential_4Proc(maxgk)

      Implicit None
! Size of the tests
! The fruit library happily handles maxgk=0,1,2,3. 
! Above that - hamilton gets too large apparently.
      Integer, Intent(in) :: maxgk
      Integer gsize,nmatp
      Integer maxg
      parameter (maxg=7)

      Integer :: igpig ((2*maxgk+1)**3)
      Real(8) :: vgpc(3,(2*maxgk+1)**3)
      Integer i,j,k
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(:,:), Allocatable :: hamilton_ref_global

! Externals
      Integer,    External :: NUMROC

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, ncols_loc, gsize_ncols_loc, gsize_nrows_loc
      Integer, Dimension(:), Allocatable :: hamilton_loc_rows, hamilton_loc_cols

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

         gsize=(2*maxgk+1)**3
         nmatp=gsize+10
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
         vgpc(:,:)=0d0
         
         Do i=1,ngrtot
            veffig(i)=dble(i)
         EndDo      

! initialisation is finished

         Call newmatrix(hamilton,nmatp, DISTRIBUTE_2D)
         Call newmatrix(hamilton_ref,nmatp, DISTRIBUTE_2D)
         Allocate(hamilton_ref_global(nmatp,nmatp))

         nrows_loc = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         ncols_loc = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         gsize_ncols_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         Allocate(hamilton_loc_rows(nrows_loc), hamilton_loc_cols(ncols_loc))
         Call getLocalIndices(nmatp, nmatp, hamilton_loc_rows, hamilton_loc_cols, MPIglobal)

! preparation of the correct answer
         hamilton_ref_global = Cmplx(0,0,8)
         do i=1,gsize
            do k=1,i
               hamilton_ref_global(k,i)=ivgig(ivg(1,igpig(k))-ivg(1,igpig(i)),ivg(2,igpig(k))-ivg(2,igpig(i)),ivg(3,igpig(k))-ivg(3,igpig(i)))
            enddo
         enddo
         Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal)

         Call hmlistln(hamilton,gsize,igpig,vgpc,gsize_nrows_loc,gsize_ncols_loc,hamilton_loc_rows,hamilton_loc_cols)

         Call assert_equals(gsize, j, 'Is the test itself set up properly?')
         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')
      
! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(hamilton_ref_global,hamilton_loc_rows,hamilton_loc_cols)  
! deallocation of global variables   
         Call freeGlobals  
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmlistln_Potential_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmlistln_Theta_4Proc
!------------------------------------------------------------------------------
! 2st test, MPI with 4 procs
! The purpose is to test whether the theta-function (cfunig) is handled properly.
! The potential is set to zero.
! The theta-function (cfunig) depends on index.
! The plane-wave lengths (vgpc) are constant.
#ifdef MPI
    Subroutine testHmlistln_Theta_4Proc

      Implicit None
! Size of the tests
      Integer gsize,nmatp,maxg,maxgk
      parameter (maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10)

      Integer :: igpig (gsize)
      Real(8) :: vgpc(3,gsize)
      Integer i,j,k
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(:,:), Allocatable :: hamilton_ref_global

! Externals
      Integer,    External :: NUMROC

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, ncols_loc, gsize_ncols_loc, gsize_nrows_loc
      Integer, Dimension(:), Allocatable :: hamilton_loc_rows, hamilton_loc_cols

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

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
         vgpc(:,:)=sqrt(2d0/3d0)

         Call newmatrix(hamilton,nmatp, DISTRIBUTE_2D)
         Call newmatrix(hamilton_ref,nmatp, DISTRIBUTE_2D)
         Allocate(hamilton_ref_global(nmatp,nmatp))

         nrows_loc = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         ncols_loc = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         gsize_ncols_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         Allocate(hamilton_loc_rows(nrows_loc), hamilton_loc_cols(ncols_loc))
         Call getLocalIndices(nmatp, nmatp, hamilton_loc_rows, hamilton_loc_cols, MPIglobal)
! initialisation is finished

         do i=1,ngrtot
            cfunig(i)=dble(i)
         enddo

! preparation of the correct answer
         hamilton_ref_global = Cmplx(0,0,8)
         do i=1,gsize
            do k=1,i
               hamilton_ref_global(k,i)=ivgig(ivg(1,igpig(k))-ivg(1,igpig(i)),ivg(2,igpig(k))-ivg(2,igpig(i)),ivg(3,igpig(k))-ivg(3,igpig(i)))
            enddo
         enddo
         Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal)

         Call hmlistln(hamilton,gsize,igpig,vgpc,gsize_nrows_loc,gsize_ncols_loc,hamilton_loc_rows,hamilton_loc_cols)

         Call assert_equals(gsize, j, 'Is the test itself set up properly?')
         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(hamilton_ref_global, hamilton_loc_rows, hamilton_loc_cols)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmlistln_Theta_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmlistln_vgpc_4Proc
!------------------------------------------------------------------------------
! 3rd test, MPI with 4 procs
! The purpose is to test whether the plane-wave lengths (vgpc) are handled properly.
! The potential is set to zero.
! The theta-function (cfunig) is constant.
! The plane-wave lengths (vgpc) depends on the plane-wave index.
#ifdef MPI
    Subroutine testHmlistln_vgpc_4Proc

      Implicit None
! Size of the tests
      Integer gsize,nmatp,maxg,maxgk
      parameter (maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10)

      Integer :: igpig (gsize)
      Real(8) :: vgpc(3,gsize)
      Integer i,j,k
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(:,:), Allocatable :: hamilton_ref_global

! Externals
      Integer,    External :: NUMROC

! MPI related variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      Integer :: nrows_loc, ncols_loc, gsize_ncols_loc, gsize_nrows_loc
      Integer, Dimension(:), Allocatable :: hamilton_loc_rows, hamilton_loc_cols

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
         Call getBlacsGridInfo(MPIglobal)

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
         do i=1,gsize
            vgpc(:,i)=sqrt(1d0/7d0)*(/1d0,2d0,3d0/)*dble(i)
         enddo

! allocate and generate complex Gaunt coefficient array
         Call newmatrix(hamilton,nmatp, DISTRIBUTE_2D)
         Call newmatrix(hamilton_ref,nmatp, DISTRIBUTE_2D)
         Allocate(hamilton_ref_global(nmatp,nmatp))

         nrows_loc = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         ncols_loc = NUMROC(nmatp, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         gsize_nrows_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
         gsize_ncols_loc = NUMROC(gsize, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
         Allocate(hamilton_loc_rows(nrows_loc), hamilton_loc_cols(ncols_loc))
         Call getLocalIndices(nmatp, nmatp, hamilton_loc_rows, hamilton_loc_cols, MPIglobal)
! initialisation is finished

         do i=1,ngrtot
         cfunig(i)=1d0
         enddo

! preparation of the correct answer
         hamilton_ref_global = Cmplx(0,0,8)
         do i=1,gsize
            do k=1,i
               hamilton_ref_global(k,i)=k*i
            enddo
         enddo
         Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal)

         Call hmlistln(hamilton,gsize,igpig,vgpc,gsize_nrows_loc,gsize_ncols_loc,hamilton_loc_rows,hamilton_loc_cols)

         Call assert_equals(gsize, j, 'Is the test itself set up properly?')
         Call assert_equals(nmatp, hamilton%size, 'checking result rank')
         Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
         Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
         Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
         Call deletematrix(hamilton)
         Call deletematrix(hamilton_ref)
         Deallocate(hamilton_ref_global, hamilton_loc_rows, hamilton_loc_cols)
! deallocation of global variables   
         Call freeGlobals
! freeing proc grid
         Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      End If
    End Subroutine testHmlistln_vgpc_4Proc
#endif

end module modHmlistln_test
