!
!
!
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genparidxran (typ, n)
      Use modmain
      Use modmpi
      Use modxs
      Implicit None
  ! arguments
      Character (1), Intent (In) :: typ
      Integer, Intent (In) :: n
  ! local variables
      Integer :: np
  ! check if number of processors is greater than set
      If (MPIglobal%procs .Gt. n) Then
         Write (*,*)
         Write (*, '("Error(genparidxran): number of processors exceeds&
        & size of set")')
         Write (*, '(" parallelization type : ", a)') typ
         Write (*, '(" size of set	       : ", i6)') n
         Write (*, '(" number of processors : ", i6)') MPIglobal%procs
         Write (*,*)
         Call terminate
      End If
  ! default values
      wpari = 1
      wparf = nwdf
      qpari = 1
      qparf = nqpt
      kpari = 1
      kparf = nkpt
  ! number of (k,kp) pairs
      np = nkpt * (nkpt+1) / 2
      ppari = 1
      pparf = np
      Select Case (typ)
      Case ('w')
         wpari = firstofset (MPIglobal%rank, n)
         wparf = lastofset (MPIglobal%rank, n)
      Case ('q')
         qpari = firstofset (MPIglobal%rank, n)
         qparf = lastofset (MPIglobal%rank, n)
      Case ('k')
         kpari = firstofset (MPIglobal%rank, n)
         kparf = lastofset (MPIglobal%rank, n)
      Case ('p')
         ppari = firstofset (MPIglobal%rank, n)
         pparf = lastofset (MPIglobal%rank, n)
      Case Default
         Write (*,*)
         Write (*, '("Error(genparidxran): unknown parallelization type&
        &: ", a)') typ
         Write (*,*)
         Call terminate
      End Select
      partype = typ
End Subroutine genparidxran
