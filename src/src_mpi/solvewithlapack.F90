    Subroutine solvewithlapack(system,nstfv,evecfv,evalfv)
      use mod_timing
      use modinput
      use modfvsystem
      use mod_eigensystem, only:nmatmax
      type(evsystem)::system
      integer::nstfv

      Real (8), Intent (Out) :: evalfv (nstfv)
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
      !local
      Integer :: is, ia, i, m, np, info ,nmatp
      Real (8) :: vl, vu
      Real (8) :: ts0, ts1
      ! allocatable arrays
      Integer, Allocatable :: iwork (:)
      Integer, Allocatable :: ifail (:)
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: v (:)
      Complex (8), Allocatable :: work (:)
      Call timesec (ts0)

         vl = 0.d0
         vu = 0.d0
         ! LAPACK 3.0 call
         !nmatmax
         nmatp=system%hamilton%size
         Allocate (iwork(5*nmatp))
         Allocate (ifail(nmatp))
         Allocate (w(nmatp))
         Allocate (rwork(7*nmatp))
         Allocate (v(1))
         Allocate (work(2*nmatp))
         !Call zhpgvx (1, 'V', 'I', 'U', nmatp, system%hamilton%zap, &
         !system%overlap%zap, vl, vu, 1, nstfv, &
         !input%groundstate%solver%evaltol, m, w, evecfv, nmatmax, work, &
         !rwork, iwork, ifail, info)
         call ZHEGVX(1, 'V', 'I', 'U', nmatp, system%hamilton%za, nmatp, system%overlap%za, nmatp,&
              vl, vu, 1, nstfv, input%groundstate%solver%evaltol, &
              m, w, evecfv, nmatmax, work, 2*nmatp, rwork, iwork, ifail, info )
         evalfv (1:nstfv) = w (1:nstfv)
         !
         !
         !
         If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(seceqnfv): diagonalisation failed")')
            Write (*, '(" ZHPGVX returned INFO = ", I8)') info
            If (info .Gt. nmatp) Then
               i = info - nmatp
               Write (*, '(" The leading minor of the overlap matrix of or&
                    &der ", I8)'                ) i
               Write (*, '("  is not positive definite")')
               Write (*, '(" Order of overlap matrix : ", I8)') nmatp
               Write (*,*)
            End If
            Stop
         End If
!         Call timesec (ts1)
         !$OMP CRITICAL
!         timefv = timefv + ts1 - ts0
         !$OMP END CRITICAL
         Call deleteystem (system)
         Deallocate (iwork, ifail, w, rwork, v, work)
         Call timesec (ts1)
         timefv = timefv + ts1 - ts0
    end subroutine solvewithlapack
