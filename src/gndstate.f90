! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gndstate
! !INTERFACE:
!
!
Subroutine gndstate
  ! !USES:
      Use modinput
      Use modmain
      Use modmpi
      Use scl_xml_out_Module
      use mpi_allgatherv_ifc_module

!
  ! !DESCRIPTION:
  !   Computes the self-consistent Kohn-Sham ground-state. General information is
  !   written to the file {\tt INFO.OUT}. First- and second-variational
  !   eigenvalues, eigenvectors and occupancies are written to the unformatted
  !   files {\tt EVALFV.OUT}, {\tt EVALSV.OUT}, {\tt EVECFV.OUT}, {\tt EVECSV.OUT}
  !   and {\tt OCCSV.OUT}.
  !
  ! !REVISION HISTORY:
  !   Created October 2002 (JKD)
  !EOP
  !BOC
      Implicit None
  ! local variables
      Logical :: exist
      Integer :: ik, is, ia, idm
      Integer :: n, nwork
      Real (8) :: timetot, et, fm
  ! allocatable arrays
      Real (8), Allocatable :: v (:)
!
      Real (8), Allocatable :: evalfv (:, :)
      Complex (8), Allocatable :: evecfv (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Logical :: force_converged, tibs
  ! time measurements
      Real(8) :: ts0,ts1,tsg0,tsg1,tin1,tin0



!! TIME - Initialisation segment
    Call timesec (tsg0)
    Call timesec (ts0)
     
  ! require forces for structural optimisation
      If ((task .Eq. 2) .Or. (task .Eq. 3)) input%groundstate%tforce = &
     & .True.
  ! initialise global variables
      Call timesec (tin0)
      Call init0
      Call timesec (tin1)
      time_init0=tin1-tin0
      Call timesec (tin0)
      Call init1
      Call timesec (tin1)
      time_init1=tin1-tin0
  ! initialize reference density
      If (allocated(rhomtref)) deallocate (rhomtref)
      Allocate (rhomtref(lmmaxvr, nrmtmax, natmtot))
      If (allocated(rhoirref)) deallocate (rhoir)
      Allocate (rhoirref(ngrtot))
  ! initialise OEP variables if required
      If (input%groundstate%xctypenumber .Lt. 0) Call init2
      If (MPIglobal%rank .Eq. 0) Then
     ! write the real and reciprocal lattice vectors to file
         Call writelat
     ! write interatomic distances to file
         Call writeiad (.False.)
     ! write symmetry matrices to file
         Call writesym
#ifdef XS
     ! write realtion to inverse symmetries
         Call writesymi
     ! write advanced information on symmetry group
         Call writesym2
     ! write out symmetrization matrix for rank 2 tensors
         Call writesymt2
#endif
     ! output the k-point set to file
         Call writekpts
     ! write lattice vectors and atomic positions to file
         Call writegeom (.False.)
         Call writegeometryxml (.False.)
     ! open INFO.OUT file
         Open (60, File='INFO'//trim(filext), Action='WRITE', Form='FOR&
        &MATTED')
     ! open TOTENERGY.OUT
         Open (61, File='TOTENERGY'//trim(filext), Action='WRITE', &
        & Form='FORMATTED')
     ! open FERMIDOS.OUT
         Open (62, File='FERMIDOS'//trim(filext), Action='WRITE', &
        & Form='FORMATTED')
     ! open MOMENT.OUT if required
         If (associated(input%groundstate%spin)) open (63, file='MOMENT&
        &'//trim(filext), action='WRITE', form='FORMATTED')
     ! open FORCEMAX.OUT if required
         If (input%groundstate%tforce) open (64, file='FORCEMAX'//&
        & trim(filext), action='WRITE', form='FORMATTED')
     ! open RMSDVEFF.OUT
         Open (65, File='RMSDVEFF'//trim(filext), Action='WRITE', &
        & Form='FORMATTED')
     ! open DTOTENERGY.OUT
         open(66,file='DTOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
     ! open DFORCEMAX.OUT
         if (input%groundstate%tforce) open(67,file='DFORCEMAX'//trim(filext),action='WRITE',form='FORMATTED')
     ! open CHGDIST.OUT
         open(68,file='CHGDIST'//trim(filext),action='WRITE',form='FORMATTED')
     ! open PCHARGE.OUT
         open(69,file='PCHARGE'//trim(filext),action='WRITE',form='FORMATTED')
     ! write out general information to INFO.OUT
         Call writeinfo (60)

     ! initialise or read the charge density and potentials from file
      End If
      iscl = 0
      If (MPIglobal%rank .Eq. 0) write (60,*)
      If ((task .Eq. 1) .Or. (task .Eq. 3)) Then
         Call readstate
         If (MPIglobal%rank .Eq. 0) write (60, '("Potential read in from STATE.OU&
        &T")')
      Else If (task .Eq. 200) Then
         Call phveff
         If (MPIglobal%rank .Eq. 0) write (60, '("Supercell potential constructed&
        & from STATE.OUT")')
      Else
         Call timesec(tin0)
         Call rhoinit
         Call timesec(tin1)
         time_density_init=tin1-tin0
  ! store density to reference
         rhoirref(:)=rhoir(:)
         rhomtref(:,:,:)=rhomt(:,:,:)
         Call timesec(tin0)
         Call poteff
         Call genveffig
         Call timesec(tin1)
         time_pot_init=tin1-tin0
         If (MPIglobal%rank .Eq. 0) write (60, '("Density and potential initialis&
        &ed from atomic data")')
      End If
      If (MPIglobal%rank .Eq. 0) Call flushifc (60)
  ! size of mixing vector
      n = lmmaxvr * nrmtmax * natmtot + ngrtot
      If (associated(input%groundstate%spin)) n = n * (1+ndmag)
      If (ldapu .Ne. 0) n = n + 2 * lmmaxlu * lmmaxlu * nspinor * &
     & nspinor * natmtot
  ! allocate mixing arrays
      Allocate (v(n))
  ! set stop flag
      tstop = .False.

      Call timesec (ts1)
      timeinit = timeinit+ts1-ts0
!! TIME - End of initialisation segment


10    Continue

!! TIME - Mixer segment
      Call timesec (ts0)
  !call mixing array allocation functions by setting
      nwork = - 1
  !and call interface
      If (MPIglobal%rank .Eq. 0) Call mixerifc (input%groundstate%mixernumber, n, &
     & v, currentconvergence, nwork)
      et = 0.d0
      fm = 0.d0
      Call timesec (ts1)
      timemixer = ts1-ts0+timemixer
!! TIME - End of mixer segment

!! TIME - First IO segment
      Call timesec (ts0)
  ! set last iteration flag
      tlast = .False.
  ! delete any existing eigenvector files
      If ((MPIglobal%rank .Eq. 0) .And. ((task .Eq. 0) .Or. (task &
     & .Eq. 2))) Call delevec
  ! begin the self-consistent loop
      If (MPIglobal%rank .Eq. 0) Then
         Write (60,*)
         Write (60, '("+------------------------------+")')
         Write (60, '("| Self-consistent loop started |")')
         Write (60, '("+------------------------------+")')
      End If
      Do iscl = 1, input%groundstate%maxscl
         If (MPIglobal%rank .Eq. 0) Then
            Write (60,*)
            Write (60, '("+-------------------------+")')
            Write (60, '("| Iteration number : ", I4, " |")') iscl
            Write (60, '("+-------------------------+")')
         End If
         If (iscl .Ge. input%groundstate%maxscl) Then
            If (MPIglobal%rank .Eq. 0) Then
               Write (60,*)
               Write (60, '("Reached self-consistent loops maximum")')
               Write (100,*)
               Write (100, '("Warning(gndstate): Reached self-consistent loops maximum")')
            End If
            tlast = .True.
         End If
         If (MPIglobal%rank .Eq. 0) Call flushifc (60)
         Call timesec (ts1)
         timeio=ts1-ts0+timeio
!! TIME - End of first IO segment

!! TIME - Muffin-tin segment
         Call timesec (ts0)
     ! generate the core wavefunctions and densities
         Call gencore
         Select Case (trim(input%groundstate%findlinentype))
         Case ('simple')
         Case ('advanced')
            If (MPIglobal%rank .Eq. 0) Then
               Write (60,*)
               Write (60, '("Using advanced method for search of linear&
              &ization energies")')
            End If
         End Select
     ! find the new linearisation energies
         Call linengy
     ! write out the linearisation energies
         if (MPIglobal%rank .eq. 0) Call writelinen
     ! generate the APW radial functions
         Call genapwfr
     ! generate the local-orbital radial functions
         Call genlofr
     ! compute the overlap radial integrals
         Call olprad
     ! compute the Hamiltonian radial integrals
         Call hmlrad
     ! zero partial charges
         chgpart(:,:,:)=0.d0

         Call timesec (ts1)
         timemt=ts1-ts0+timemt
!! TIME - End of muffin-tin segment


!! TIME - Second IO segment       
          Call timesec (ts0)
#ifdef MPI
         Call MPI_barrier (MPI_COMM_WORLD, MPIglobal%ierr)
         If (MPIglobal%rank .Eq. 0) Call delevec ()
#endif
#ifdef MPISEC
         splittfile = .True.
         Do ik = firstk (MPIglobal%rank), lastk (MPIglobal%rank)
!
#endif
#ifdef NEVERDEFINED
         End Do
#endif
#ifndef MPISEC
         splittfile = .False.
     ! begin parallel loop over k-points
#ifdef KSMP
     !$OMP PARALLEL DEFAULT(SHARED) &
     !$OMP PRIVATE(evalfv,evecfv,evecsv)
     !$OMP DO
#endif
         Do ik = 1, nkpt
#endif
        ! every thread should allocate its own arrays
            Allocate (evalfv(nstfv, nspnfv))
            Allocate (evecfv(nmatmax, nstfv, nspnfv))
            Allocate (evecsv(nstsv, nstsv))
!! TIME - seceqn does not belong to IO
            Call timesec(ts1)
            timeio=ts1-ts0+timeio            
        ! solve the first- and second-variational secular equations
            Call seceqn (ik, evalfv, evecfv, evecsv)
            Call timesec(ts0)
        ! write the eigenvalues/vectors to file
            Call putevalfv (ik, evalfv)
            Call putevalsv (ik, evalsv(:, ik))
            Call putevecfv (ik, evecfv)
            Call putevecsv (ik, evecsv)
        ! calculate partial charges
            call genpchgs(ik,evecfv,evecsv)
            Deallocate (evalfv, evecfv, evecsv)
         End Do
#ifdef KSMP
     !$OMP END DO
     !$OMP END PARALLEL
#endif
#ifdef MPISEC
         call mpi_allgatherv_ifc(nkpt,nstsv,rbuf=evalsv)
         Call MPI_barrier (MPI_COMM_WORLD, MPIglobal%ierr)
#endif
     ! find the occupation numbers and Fermi energy
         Call occupy
         If (MPIglobal%rank .Eq. 0) Then
        ! write out the eigenvalues and occupation numbers
            Call writeeval
        ! write the Fermi energy to file
            Call writefermi
         End If
     ! set the charge density and magnetisation to zero
         rhomt (:, :, :) = 0.d0
         rhoir (:) = 0.d0
         If (associated(input%groundstate%spin)) Then
            magmt (:, :, :, :) = 0.d0
            magir (:, :) = 0.d0
         End If
#ifdef MPIRHO
         Do ik = firstk (MPIglobal%rank), lastk (MPIglobal%rank)
        !write the occupancies to file
            Call putoccsv (ik, occsv(:, ik))
         End Do
         Do ik = firstk (MPIglobal%rank), lastk (MPIglobal%rank)
#endif
#ifndef MPIRHO
            If (MPIglobal%rank .Eq. 0) Then
               Do ik = 1, nkpt
              !write the occupancies to file
                  Call putoccsv (ik, occsv(:, ik))
               End Do
            End If
#ifdef KSMP
 ! begin parallel loop over k-points
 !$OMP PARALLEL DEFAULT(SHARED) &
 !$OMP PRIVATE(evecfv,evecsv)
 !$OMP DO
#endif
            Do ik = 1, nkpt
#endif
               Allocate (evecfv(nmatmax, nstfv, nspnfv))
               Allocate (evecsv(nstsv, nstsv))
           ! get the eigenvectors from file
               Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
               Call getevecsv (vkl(:, ik), evecsv)

!! TIME - rhovalk does not belong to IO
               Call timesec(ts1)
               timeio=ts1-ts0+timeio
           ! add to the density and magnetisation
               Call rhovalk (ik, evecfv, evecsv)
               Deallocate (evecfv, evecsv)
               Call timesec(ts0)
            End Do

#ifndef MPIRHO
#ifdef KSMP
        !$OMP END DO
        !$OMP END PARALLEL
#endif
#endif
#ifdef MPIRHO
            Call mpisumrhoandmag
#endif
#ifdef MPI
            If (input%groundstate%xctypenumber .Lt. 0) Call &
           & mpiresumeevecfiles ()
#endif
            If (MPIglobal%rank .Eq. 0) Then
        ! write out partial charges
               call writepchgs(69,input%groundstate%lmaxvr)
               call flushifc(69)
            end if

            Call timesec(ts1)
            timeio=ts1-ts0+timeio
!! TIME - End of second IO segment

!! TIME - potential
           Call timesec(ts0)
        ! symmetrise the density
            Call symrf (input%groundstate%lradstep, rhomt, rhoir)
        ! symmetrise the magnetisation
            If (associated(input%groundstate%spin)) Call symrvf &
           & (input%groundstate%lradstep, magmt, magir)
        ! convert the density from a coarse to a fine radial mesh
            Call rfmtctof (rhomt)
        ! convert the magnetisation from a coarse to a fine radial mesh
            Do idm = 1, ndmag
               Call rfmtctof (magmt(:, :, :, idm))
            End Do
        ! add the core density to the total density
            Call addrhocr
        ! calculate the charges
            Call charge
        ! calculate the moments
            If (associated(input%groundstate%spin)) Call moment
        ! normalise the density
            Call rhonorm
        ! LDA+U
            If (ldapu .Ne. 0) Then
           ! generate the LDA+U density matrix
               Call gendmatlu
           ! generate the LDA+U potential matrix
               Call genvmatlu
           ! write the LDA+U matrices to file
               if (MPIglobal%rank .eq. 0) Call writeldapu
            End If
        ! generate charge distance
            call chgdist
        ! store density to reference
            rhoirref(:)=rhoir(:)
            rhomtref(:,:,:)=rhomt(:,:,:)
        ! compute the effective potential
            Call poteff
        ! pack interstitial and muffin-tin effective potential and field into one array
            Call packeff (.True., n, v)
        ! mix in the old potential and field with the new
            If (MPIglobal%rank .Eq. 0) Call mixerifc &
           & (input%groundstate%mixernumber, n, v, currentconvergence, &
           & nwork)
#ifdef MPI
         Call MPI_bcast (v(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPIglobal%ierr)
#endif
 ! unpack potential and field
            Call packeff (.False., n, v)
        ! add the fixed spin moment effect field
            If (getfixspinnumber() .Ne. 0) Call fsmfield
        ! Fourier transform effective potential to G-space
            Call genveffig
        ! reduce the external magnetic fields if required
            If (associated(input%groundstate%spin)) Then
               If (input%groundstate%spin%reducebf .Lt. 1.d0) Then
                  input%groundstate%spin%bfieldc(:) = &
                 & input%groundstate%spin%bfieldc(:) * &
                 & input%groundstate%spin%reducebf
                  Do is = 1, nspecies
                     Do ia = 1, natoms (is)
                        input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) = &
                       & input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) * input%groundstate%spin%reducebf
                     End Do
                  End Do
               End If
            End If
        ! compute the energy components
            Call energy
           Call timesec(ts1)
           timepot=ts1-ts0+timepot
!! TIME - End of potential

!! TIME - Third IO segment
        ! compute the forces (without IBS corrections)
            Call timesec(ts0)
            If (input%groundstate%tforce) Then
               tibs=input%groundstate%tfibs
               input%groundstate%tfibs=.false.
!! forces do no belong to IO
               Call timesec(ts1)
               timeio=ts1-ts0+timeio
               call force
               Call timesec(ts0)
               
               input%groundstate%tfibs=tibs
               If (MPIglobal%rank .Eq. 0) Then
                  
       ! output forces to INFO.OUT
                  Call writeforce (60)
                  ! write maximum force magnitude to FORCEMAX.OUT
                  Write (64, '(G18.10,"   (without IBS correction)")') forcemax
                  Call flushifc (64)
               end if
            end if
            If (MPIglobal%rank .Eq. 0) Then
        ! output energy components
               Call writeengy (60)
               Write (60,*)
               Write (60, '("Density of states at Fermi energy : ", G18&
              &.10)') fermidos
               Write (60, '(" (states/Hartree/unit cell)")')
           ! write total energy to TOTENERGY.OUT and flush
               Write (61, '(G22.12)') engytot
               Call flushifc (61)
           ! write DOS at Fermi energy to FERMIDOS.OUT and flush
               Write (62, '(G18.10)') fermidos
               Call flushifc (62)
           ! output charges and moments
               Call writechg (60)
           ! write total moment to MOMENT.OUT and flush
               If (associated(input%groundstate%spin)) Then
                  Write (63, '(3G18.10)') momtot (1:ndmag)
                  Call flushifc (63)
               End If
           ! output effective fields for fixed spin moment calculations
               If (getfixspinnumber() .Ne. 0) Call writefsm (60)
           ! check for WRITE file
               Inquire (File='WRITE', Exist=Exist)
               If (exist) Then
                  Write (60,*)
                  Write (60, '("WRITE file exists - writing STATE.OUT")&
                 &')
                  Call writestate
                  Open (50, File='WRITE')
                  Close (50, Status='DELETE')
               End If
           ! write STATE.OUT file if required
               If (input%groundstate%nwrite .Ge. 1) Then
                  If (Mod(iscl, input%groundstate%nwrite) .Eq. 0) Then
                     Call writestate
                     Write (60,*)
                     Write (60, '("Wrote STATE.OUT")')
                  End If
               End If
               ! update convergence criteria
               deltae=abs(et-engytot)
               dforcemax=abs(fm-forcemax)
               Call scl_iter_xmlout ()
               If (associated(input%groundstate%spin)) Call &
              & scl_xml_write_moments ()
               Call scl_xml_out_write ()
            End If
            Call timesec(ts1)
            timeio=ts1-ts0+timeio
!! TIME - End of third IO segment

 ! exit self-consistent loop if last iteration is complete
            If (tlast) Go To 20

!! TIME - Fourth IO segment
            Call timesec(ts0)
            If (MPIglobal%rank .Eq. 0) Then
    ! check for convergence
               If (iscl .Ge. 2) Then
                  Write (60,*)
                  Write (60, '("RMS change in effective potential (targ&
                    &et) : ", G18.10, " (", G18.10, ")")') &
                    & currentconvergence, input%groundstate%epspot
                  write(60,'("Absolute change in total energy (target)   : ",G18.10," (",&
                  &G18.10,")")') deltae, input%groundstate%epsengy
                  if (input%groundstate%tforce) then
                    write(60,'("Absolute change in |max. force| (target)   : ",G18.10," (",&
                    &G18.10,")")') dforcemax, input%groundstate%epsforce
                  end if
                  write(60,'("Charge distance (target)                   : ",G18.10," (",&
                  &G18.10,")")') chgdst, input%groundstate%epschg
                  write(66,'(G18.10)') deltae
                  call flushifc(66)
                  if (input%groundstate%tforce) then
                    write(67,'(G18.10)') dforcemax
                    call flushifc(67)
                  end if
                  write(68,'(G18.10)') chgdst
                  call flushifc(68)
                  Write (65, '(G18.10)') currentconvergence
                  Call flushifc (65)
                  If ((currentconvergence .Lt. input%groundstate%epspot).and. &
                 & (deltae .lt. input%groundstate%epsengy).and. &
                 & (chgdst .lt. input%groundstate%epschg).and. &
                 & (dforcemax .lt. input%groundstate%epsforce)) Then
                     Write (60,*)
                     Write (60, '("Convergence targets achieved")')
                     tlast = .True.
                  End If
               End If
               et = engytot
               fm = forcemax
               If (input%groundstate%xctypenumber .Lt. 0) Then
                  Write (60,*)
                  Write (60, '("Magnitude of OEP residual : ", G18.10)') resoep
               End If
           ! check for STOP file
               Inquire (File='STOP', Exist=Exist)
               If (exist) Then
                  Write (60,*)
                  Write (60, '("STOP file exists - stopping self-consis&
                 &tent loop")')
                  tstop = .True.
                  tlast = .True.
                  Open (50, File='STOP')
                  Close (50, Status='DELETE')
               End If
           ! output the current total time
               timetot = timeinit + timemat + timefv + timesv + timerho &
              & + timepot + timefor+timeio+timemt+timemixer
               Write (60,*)
               Write (60, '("Wall time (seconds) : ", F12.2)') timetot
           ! end the self-consistent loop
            End If
#ifdef MPI
            Call MPI_bcast (tstop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPIglobal%ierr)
            Call MPI_bcast (tlast, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPIglobal%ierr)
#endif
            Call timesec(ts1)
            timeio=ts1-ts0+timeio
!! TIME - End of fourth IO segment
         End Do
20       Continue

!! TIME - Fifth IO segment
         Call timesec(ts0)
         If (MPIglobal%rank .Eq. 0) Then
            Write (60,*)
            Write (60, '("+------------------------------+")')
            Write (60, '("| Self-consistent loop stopped |")')
            Write (60, '("+------------------------------+")')
 ! write density and potentials to file only if maxscl > 1
            If (input%groundstate%maxscl .Gt. 1) Then
               Call writestate
               Write (60,*)
               Write (60, '("Wrote STATE.OUT")')
            End If
         end if
 !-----------------------!
 !     compute forces    !
 !-----------------------!

         If (( .Not. tstop) .And. (input%groundstate%tforce)) Then
            Call force
            If (MPIglobal%rank .Eq. 0) Then
       ! output forces to INFO.OUT
                  Call writeforce (60)


              ! write maximum force magnitude to FORCEMAX.OUT
                  Write (64, '(G18.10)') forcemax
                  Call flushifc (64)
            End If
         End If
     !---------------------------------------!
     !     perform structural relaxation     !
     !---------------------------------------!
         If (( .Not. tstop) .And. ((task .Eq. 2) .Or. (task .Eq. 3))) &
        & Then
            If (MPIglobal%rank .Eq. 0) Then
               Write (60,*)
               Write (60, '("Maximum force magnitude (target) : ", G18.&
              &10, " (", G18.10, ")")') forcemax, &
              & input%structureoptimization%epsforce
               Call flushifc (60)
            End If
 ! check force convergence
            force_converged = .False.
            If (forcemax .Le. input%structureoptimization%epsforce) Then
               If (MPIglobal%rank .Eq. 0) Then
                  Write (60,*)
                  Write (60, '("Force convergence target achieved")')
               End If
               force_converged = .True.
            End If
            If (force_converged) Go To 30
 ! update the atomic positions if forces are not converged
            Call updatpos
            If (MPIglobal%rank .Eq. 0) Then
               Write (60,*)
               Write (60, '("+--------------------------+")')
               Write (60, '("| Updated atomic positions |")')
               Write (60, '("+--------------------------+")')
               Do is = 1, nspecies
                  Write (60,*)
                  Write (60, '("Species : ", I4, " (", A, ")")') is, trim (input%structure%speciesarray(is)%species%chemicalSymbol)
                  Write (60, '(" atomic positions (lattice) :")')
                  Do ia = 1, natoms (is)
                     Write (60, '(I4, " : ", 3F14.8)') ia, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
                  End Do
               End Do
           ! add blank line to TOTENERGY.OUT, FERMIDOS.OUT, MOMENT.OUT and RMSDVEFF.OUT
               Write (61,*)
               Write (62,*)
               If (associated(input%groundstate%spin)) write (63,*)
               Write (65,*)
           ! add blank line to DTOTENERGY.OUT, DFORCEMAX.OUT, CHGDIST.OUT and PCHARGE.OUT
               Write (66,*)
               if (input%groundstate%tforce) Write (67,*)
               Write (68,*)
               Write (69,*)
            End If
           ! begin new self-consistent loop with updated positions
!! TIME - stop the IO timer if we return to the SCF loop
            Call timesec(ts1)
            timeio=ts1-ts0+timeio
            Go To 10
         End If
30       Continue
     ! output timing information
         If (MPIglobal%rank .Eq. 0) Then
 ! close the TOTENERGY.OUT file
            Close (61)
 ! close the FERMIDOS.OUT file
            Close (62)
 ! close the MOMENT.OUT file
            If (associated(input%groundstate%spin)) close (63)
 ! close the FORCEMAX.OUT file
            If (input%groundstate%tforce) close (64)
 ! close the RMSDVEFF.OUT file
            Close (65)
 ! close the DTOTENERGY.OUT file
            close(66)
 ! close the DFORCEMAX.OUT file
            if (input%groundstate%tforce) close(67)
 ! close the CHGDIST.OUT file
            close(68)
 ! close the PCHARGE.OUT file
            close(69)
            Call structure_xmlout ()
            Call scl_xml_setGndstateStatus ("finished")
            Call scl_xml_out_write ()
         End If
     !set nwork to -2 to tell interface to call the deallocation functions
         If (MPIglobal%rank .Eq. 0) Call mixerifc (input%groundstate%mixernumber, &
        & n, v, currentconvergence,-2)
         Deallocate (v)
         Call mpiresumeevecfiles ()
         Call timesec(ts1)
         timeio=ts1-ts0+timeio
!! TIME - End of fifth IO segment

 ! close the INFO.OUT file
         If (MPIglobal%rank .Eq. 0) then
            Write (60,*)
            Write (60, '("Timings (seconds) :")')
            Write (60, '(" initialisation", T40, ": ", F12.2)') timeinit
            Write (60, '("            - init0", T40,": ", F13.2)') time_init0
            Write (60, '("            - init1", T40,": ", F13.2)') time_init1
            Write (60, '("            - rhoinit", T40,": ", F13.2)') time_density_init
            Write (60, '("            - potential initialisation", T40,": ", F13.2)') time_pot_init
            Write (60, '("            - others", T40,": ", F13.2)') timeinit-time_init0-time_init1-time_density_init-time_pot_init
            Write (60, '(" Hamiltonian and overlap matrix set up", T40,&
           & ": ", F12.2)') timemat
            Write (60, '("            - hmlaan", T40,": ", F13.2)') time_hmlaan
            Write (60, '("            - hmlalon", T40,": ", F13.2)') time_hmlalon
            Write (60, '("            - hmllolon", T40,": ", F13.2)') time_hmllolon
            Write (60, '("            - olpaan", T40,": ", F13.2)') time_olpaan
            Write (60, '("            - olpalon", T40,": ", F13.2)') time_olpalon
            Write (60, '("            - olplolon", T40,": ", F13.2)') time_olplolon
            Write (60, '("            - hmlistln", T40,": ", F13.2)') time_hmlistln
            Write (60, '("            - olpistln", T40,": ", F13.2)') time_olpistln

            Write (60, '(" first-variational secular equation", T40, ":&
           & ", F12.2)') timefv
            If (associated(input%groundstate%spin)) Then
               Write (60, '(" second-variational calculation", T40, ": ", F12.2)') timesv
            End If
            Write (60, '(" charge density calculation", T40, ": ", F12.&
           &2)') timerho
            Write (60, '(" potential calculation", T40, ": ", F12.2)') &
           & timepot
            Write (60, '(" muffin-tin manipulations", T40, ": ", F12.2)') &
           & timemt
            Write (60, '(" APW matching", T40, ": ", F12.2)') &
           & timematch
            Write (60, '(" disk reads/writes", T40, ": ", F12.2)') &
           & timeio
            Write (60, '(" mixing efforts", T40, ": ", F12.2)') &
           & timemixer
            If (input%groundstate%tforce) Then
               Write (60, '(" force calculation", T40, ": ", F12.2)') &
              & timefor
            End If
            timetot = timeinit + timemat + timefv + timesv + timerho + &
           & timepot + timefor+timeio+timemt+timemixer+timematch
            Write (60, '(" sum", T40, ": ", F12.2)') timetot
            Call timesec(tsg1)
            Write (60, '(" total", T40, ": ", F12.2)') tsg1-tsg0
            Write (60,*)
            Write (60, '("+----------------------------+")')
            Write (60, '("| Groundstate module stopped |")')
            Write (60, '("+----------------------------+")')
            close (60)
         endif

         deallocate(rhomtref,rhoirref)
         Return
   End Subroutine gndstate
!EOC
