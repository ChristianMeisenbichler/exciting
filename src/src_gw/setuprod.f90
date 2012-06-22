!BOP
!
! !ROUTINE: setuprod
!
! !INTERFACE:
       subroutine setuprod
       
! !DESCRIPTION:
!
! This subroutine calculates the radial product functions
! and their overlap matrix
!
!
! !USES:

      use modinput
      use modmain
      use modgw

! !LOCAL VARIABLES:

      implicit none

      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: ilo1
      integer(4) :: ilo2
      integer(4) :: io1
      integer(4) :: io2
      integer(4) :: ipr1
      integer(4) :: ipr2
      integer(4) :: ir
      integer(4) :: is
      integer(4) :: ist
      integer(4) :: l1
      integer(4) :: l2

      real(8) :: fr(nrmtmax)
      real(8) :: gr(nrmtmax) 
      real(8) :: cf(3,nrmtmax)

! !EXTERNAL ROUTINES:

      external fderiv

! !REVISION HISTORY:
!
! Created 17. May 2006 by RGA
! Revisited 5.05.2011 by DIN 
!
!
!EOP
!BOC
!----------------------------------------------------------------------!
!           Product of core and valence states                         !
!----------------------------------------------------------------------!
!     Set the maximum possible number of product functions
      maxnup=2*(ncmax+lmaxapwmix+nlomax+1)*(lmaxapwmix+nlomax+1)

      write(701,*) 'ncmax,lmaxapw,nlomax:', ncmax,input%groundstate%lmaxapw,nlomax
      write(701,*) 'ncore:', ncore(1)
      write(701,*) 'nlorb:', nlorb(1)
      write(701,*) 'lmaxapwmix:', lmaxapwmix
      write(701,*) 'maxnup:', maxnup

!     Allocate the array for the product functions and initialize
      allocate(uprod(natmtot,maxnup,nrmtmax))
      allocate(eles(natmtot,maxnup,2))
      allocate(nup(natmtot))
      uprod=0.0d0

!     Reset maxnup to zero
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          ipr1=0

!----------------------------------------------------------------------!
!         Product of core and apw functions
!----------------------------------------------------------------------!
          do ist=1,ncore(is)
            do l2=0,lmaxapwmix
              if(mixopt(l2,is,2))then
                ipr1=ipr1+1
                eles(ias,ipr1,1)=spl(ist,is)
                eles(ias,ipr1,2)=l2
                do ir=1,nrmt(is)
                  uprod(ias,ipr1,ir)=apwfr(ir,1,1,l2,ias)*            &
       &                             rwfcr(ir,1,ist,ias)*spr(ir,is)
                enddo ! ir
              endif ! mixopt
            enddo ! l2
            
!----------------------------------------------------------------------!
!           Product of core and local orbital functions
!----------------------------------------------------------------------!
            do ilo2=1,nlorb(is)
              l2=lorbl(ilo2,is)
              if(mixopt(l2,is,2))then
                ipr1=ipr1+1
                eles(ias,ipr1,1)=spl(ist,is)
                eles(ias,ipr1,2)=l2 
                do ir=1,nrmt(is)
                  uprod(ias,ipr1,ir)=lofr(ir,1,ilo2,ias)*             &
       &                             rwfcr(ir,1,ist,ias)*spr(ir,is)
                enddo ! ir
              endif ! mixopt
            enddo ! ilo2
          enddo ! ist
          
!----------------------------------------------------------------------!
!         products between valence/conduction states                   !
!----------------------------------------------------------------------!
          do l1=0,lmaxapwmix
            if(mixopt(l1,is,1))then
              io1=1
              do l2=l1,lmaxapwmix
                if(mixopt(l2,is,2))then
                  io2=1
                  ipr1=ipr1+1
                  eles(ias,ipr1,1)=l1
                  eles(ias,ipr1,2)=l2 
                  do ir=1,nrmt(is)
                    uprod(ias,ipr1,ir)=apwfr(ir,1,io1,l1,ias)*        &
     &                                 apwfr(ir,1,io2,l2,ias)*spr(ir,is)
                  enddo ! ir
                endif ! mixopt2
              enddo ! l2
            endif ! mixopt1
          enddo ! l1

!----------------------------------------------------------------------!
!         product between local orbitals and apw functions
!----------------------------------------------------------------------!
          do ilo1=1,nlorb(is)
            l1=lorbl(ilo1,is)
            if(mixopt(l1,is,1))then
              do l2=0,lmaxapwmix
                if(mixopt(l2,is,2))then
                  io2=1
                  ipr1=ipr1+1
                  eles(ias,ipr1,1)=l1
                  eles(ias,ipr1,2)=l2 
                  do ir=1,nrmt(is)
                    uprod(ias,ipr1,ir)=lofr(ir,1,ilo1,ias)*           &
     &                                 apwfr(ir,1,io2,l2,ias)*spr(ir,is)
                  enddo ! ir
                endif ! mixopt2
              enddo ! l2
            endif ! mixopt1
          enddo ! ilo1

!----------------------------------------------------------------------!         
!        products between local and local orbitals
!----------------------------------------------------------------------!
          do ilo1=1,nlorb(is)
            l1=lorbl(ilo1,is)
            if(mixopt(l1,is,1))then
              do ilo2=1,nlorb(is)
                l2=lorbl(ilo2,is)
                if(mixopt(l2,is,2))then
                  ipr1=ipr1+1
                  eles(ias,ipr1,1)=l1
                  eles(ias,ipr1,2)=l2 
                  do ir=1,nrmt(is)
                    uprod(ias,ipr1,ir)=lofr(ir,1,ilo1,ias)*           &
     &                                 lofr(ir,1,ilo2,ias)*spr(ir,is)
                  enddo ! ir
                endif ! mixopt2
              enddo ! ilo2
            endif ! mixopt1
          enddo ! ilo1

!         reset nup to the exact number of product functions
          nup(ias)=ipr1
          write(701,*) "setuprod: # of uprod for atom",ias,nup(ias)

        enddo ! ia
      enddo !is
!      
!     allocate the overlap matrix and initialize it
!
      allocate(umat(natmtot,maxnup,maxnup))
      umat=0.0d0
!     calculate the overlap matrix of product functions
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          do ipr1=1,nup(ias)
            do ipr2=ipr1,nup(ias)
              do ir=1,nrmt(is)
                fr(ir)=uprod(ias,ipr1,ir)*uprod(ias,ipr2,ir)
              end do
              call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
              umat(ias,ipr1,ipr2)=gr(nrmt(is))
            enddo ! ipr2
          enddo ! ipr1
        enddo ! ia
      enddo !is
      
! DEBUG
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          write(701,*) 
          write(701,*) "###umat for Atom ", ias
          do ipr1=1,nup(ias)
            do ipr2=ipr1,nup(ias)
                write(701,'(2i4,f20.6)') ipr1,ipr2,umat(ias,ipr1,ipr2)

            enddo ! ipr2
          enddo ! ipr1
          write(701,*)        
        enddo ! ia            
      enddo !is
      
      return
      end subroutine setuprod
!EOC
