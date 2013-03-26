!
! service functions to partition k points
! still kept around but should be replaced by generic functions
! firstofset nofset lastofset ....
Module setpatrtitioning

contains
      Function nofk (process)
         Use modmain, Only: nkpt
         use modmpi_types, only: MPIglobal
         Integer :: nofk
         Integer, Intent (In) :: process
         nofk = nkpt / MPIglobal%procs
         If ((Mod(nkpt, MPIglobal%procs) .Gt. process)) nofk = nofk + 1
      End Function nofk
!
      Function firstk (process)
      use modmpi_types, only: MPIglobal
         Integer :: firstk
         Integer, Intent (In) :: process
         integer::i
         firstk = 1
         Do i = 0, process - 1
            firstk = firstk + nofk (i)
         End Do
      End Function firstk
!
      Function lastk (process)
      use modmpi_types, only: MPIglobal
         Integer :: lastk,i
         Integer, Intent (In) :: process
         lastk = 0
         Do i = 0, process
            lastk = lastk + nofk (i)
         End Do
      End Function lastk
!
      Function procofk (k)
      use modmpi_types, only: MPIglobal
         Integer :: procofk
         Integer, Intent (In) :: k
         Integer :: iproc
         procofk = 0
         Do iproc = 0, MPIglobal%procs - 1
            If (k .Gt. lastk(iproc)) procofk = procofk + 1
         End Do
      End Function procofk
!
!-------------generalized partition!
      Function nofset (process, set)
      use modmpi_types, only: MPIglobal
         Integer :: nofset
         Integer, Intent (In) :: process, set
         nofset = set / MPIglobal%procs
         If ((Mod(set, MPIglobal%procs) .Gt. process)) nofset = nofset + 1
      End Function nofset
!
      Function firstofset (process, set)
      use modmpi_types, only: MPIglobal
         Integer :: firstofset
         Integer, Intent (In) :: process, set
         integer::i
         firstofset = 1
         Do i = 0, process - 1
            firstofset = firstofset + nofset (i, set)
         End Do
      End Function firstofset
!
      Function lastofset (process, set)
      use modmpi_types, only: MPIglobal
         Integer :: lastofset
         Integer, Intent (In) :: process, set
         integer::i
         lastofset = 0
         Do i = 0, process
            lastofset = lastofset + nofset (i, set)
         End Do
      End Function lastofset
!
      Function procofindex (k, set)
      use modmpi_types, only: MPIglobal
         Integer :: procofindex
         Integer, Intent (In) :: k, set
         Integer :: iproc
         procofindex = 0
         Do iproc = 0, MPIglobal%procs - 1
            If (k .Gt. lastofset(iproc, set)) procofindex = procofindex &
           & + 1
         End Do
      End Function procofindex
!
      Function lastproc (row, set)
      use modmpi_types, only: MPIglobal
         Implicit None
         Integer :: lastproc
         Integer, Intent (In) :: row, set
         If (row .Ne. nofset(0, set)) Then
            lastproc = MPIglobal%procs
         Else
            lastproc = modulo (set, MPIglobal%procs)
            If (lastproc .Eq. 0) lastproc = MPIglobal%procs
         End If
         lastproc = lastproc - 1
      End Function lastproc
  end module
