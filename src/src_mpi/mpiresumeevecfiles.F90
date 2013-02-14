!
!
!
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: mpiresumeevec
! !INTERFACE:
!
!
Subroutine mpiresumeevecfiles
      Use modmain
#ifdef MPI
      Use modmpi
! !DESCRIPTION:
!  Routine reads the lokal EIGVECk1-k2.OUT files that were created to
!  avoid file system inconsistencies. And writes them into the standard
!  EIGVEC.OUT File
!
! !REVISION HISTORY:
!   Created October SEPT 2006 (MULEOBEN)
!   by Cristian Meisenbichler
!EOP
      Implicit None
  ! arguments
      Integer :: ik, proc, nmatmax_, nstfv_, nspnfv_, nstsv_, recl, &
     & token
      Integer :: recvstatus (MPI_STATUS_SIZE)
      Character (256) :: filetag
      Complex (8) :: evecfv (nmatmax, nstfv, nspnfv)
      Complex (8) :: evecsv (nstsv, nstsv)
      Real (8) :: evalfv (nstfv, nspnfv), vkl_ (3), evalsvp (nstsv), &
     & occsvp (nstsv)
      Character (256), External :: outfilenamestring
      If (MPIglobal%procs .Gt. 1) Call MPI_barrier (MPI_COMM_WORLD, MPIglobal%ierr)
      If (splittfile .And. (MPIglobal%procs .Gt. 1) .and. (nkpt .gt.MPIglobal%rank)) Then
! start a receive in order to pass around a token from rank 0 to max
         If (MPIglobal%rank .Ne. 0) Call mpi_recv (token, 1, MPI_INTEGER, MPIglobal%rank-1, &
        & 1, MPI_COMM_WORLD, recvstatus, MPIglobal%ierr)
!
         filetag = 'EVECFV'
         Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_, &
        & evecfv
         Open (71, File=trim(filetag)//trim(filext), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
         proc = MPIglobal%rank
         Open (77, File=outfilenamestring(filetag, firstk(proc)), &
        & Action='READ', Form='UNFORMATTED', Access='DIRECT', &
        & Recl=Recl)
         Do ik = firstk (proc), lastk (proc)
            Read (77, Rec=ik-firstk(procofk(ik))+1) vkl_, nmatmax_, &
           & nstfv_, nspnfv_, evecfv
            Write (71, Rec=ik) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
         End Do
         Close (77, Status='DELETE')
         Close (71)
!
         filetag = 'EVECSV'
         Inquire (IoLength=Recl) vkl_, nstsv_, evecsv
         Open (71, File=trim(filetag)//trim(filext), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
         proc = MPIglobal%rank
         Open (77, File=outfilenamestring(filetag, firstk(proc)), &
        & Action='READ', Form='UNFORMATTED', Access='DIRECT', &
        & Recl=Recl)
         Do ik = firstk (proc), lastk (proc)
            Read (77, Rec=ik-firstk(procofk(ik))+1) vkl_, nstsv_, &
           & evecsv
            Write (71, Rec=ik) vkl_, nstsv_, evecsv
         End Do
         Close (77, Status='DELETE')
         Close (71)
!
         filetag = 'EVALFV'
         Inquire (IoLength=Recl) vkl_, nstfv_, nspnfv_, evalfv
         Open (71, File=trim(filetag)//trim(filext), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
         proc = MPIglobal%rank
         Open (77, File=outfilenamestring(filetag, firstk(proc)), &
        & Action='READ', Form='UNFORMATTED', Access='DIRECT', &
        & Recl=Recl)
         Do ik = firstk (proc), lastk (proc)
            Read (77, Rec=ik-firstk(procofk(ik))+1) vkl_, nstfv_, &
           & nspnfv_, evalfv
            Write (71, Rec=ik) vkl_, nstfv_, nspnfv_, evalfv
         End Do
         Close (77, Status='DELETE')
         Close (71)
!
         filetag = 'EVALSV'
         Inquire (IoLength=Recl) vkl_, nstsv_, evalsvp
         Open (71, File=trim(filetag)//trim(filext), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
         proc = MPIglobal%rank
         Open (77, File=outfilenamestring(filetag, firstk(proc)), &
        & Action='READ', Form='UNFORMATTED', Access='DIRECT', &
        & Recl=Recl)
         Do ik = firstk (proc), lastk (proc)
            Read (77, Rec=ik-firstk(procofk(ik))+1) vkl_, nstsv_, &
           & evalsvp
            Write (71, Rec=ik) vkl_, nstsv_, evalsvp
         End Do
         Close (77, Status='DELETE')
         Close (71)
!
         filetag = 'OCCSV'
         Inquire (IoLength=Recl) vkl_, nstsv_, occsvp
         Open (71, File=trim(filetag)//trim(filext), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
         proc = MPIglobal%rank
         Open (77, File=outfilenamestring(filetag, firstk(proc)), &
        & Action='READ', Form='UNFORMATTED', Access='DIRECT', &
        & Recl=Recl)
         Do ik = firstk (proc), lastk (proc)
            Read (77, Rec=ik-firstk(procofk(ik))+1) vkl_, nstsv_, &
           & occsvp
            Write (71, Rec=ik) vkl_, nstsv_, occsvp
         End Do
         Close (77, Status='DELETE')
         Close (71)
!
         write(*,*) 'cheers'
         Call SYSTEM ("sync")
         If (MPIglobal%rank .Eq. 0) Then
            Write (60,*) "resumed split files"
            Call flushifc (60)
         End If
         !if I am not the last process pass on the token
         If (MPIglobal%rank .Ne. (MPIglobal%procs-1) .and. (MPIglobal%rank .ne. (nkpt-1))) then
         Call mpi_send (token, 1, MPI_INTEGER, &
        & MPIglobal%rank+1, 1, MPI_COMM_WORLD, MPIglobal%ierr)
         endif
      End If
      call barrier
      Call SYSTEM ("sync")
      splittfile = .False.

#endif
!
      Return
End Subroutine mpiresumeevecfiles
