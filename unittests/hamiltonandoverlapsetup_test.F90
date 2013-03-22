! !MODULE:  modHmlaan_test
! !DESCRIPTION:
! Modules with test and helper routines for testing hmlaan
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
module modHamiltonandoverlapsetup_test
    Use fruit
    Use test_helpers
    Use modinput,        Only: input
    Use mod_muffin_tin,  Only: idxlm, lmmaxapw,lmmaxvr,lmmaxmat, rmt, nrmt
    Use mod_APW_LO,      Only: apword, apwordmax, apwfr, apwdfr
    Use mod_Gkvector,    Only: ngkmax
    Use mod_Gvector,     Only: ivg,ngrtot,ivgig
    Use mod_atoms,       Only: natmtot, idxas
    Use mod_eigensystem, Only: gntyry, haa
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix,evsystem,newsystem,deletesystem
#ifdef MPI
    use modmpi
#endif

    Double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains


    Subroutine testcaseHamiltonandoverlapsetupSerial
      Implicit None

!       Call set_test_name ('contribution of hmlaan')
!       Call testHmlaanContrib_Serial

    End Subroutine testcaseHamiltonandoverlapsetupSerial


#ifdef MPI
    Subroutine testcaseHamiltonandoverlapsetup1Proc
      Implicit None


    End Subroutine testcaseHamiltonandoverlapsetup1Proc
#endif


#ifdef MPI
    Subroutine testcaseHamiltonandoverlapsetup4Proc
      Implicit None


    End Subroutine testcaseHamiltonandoverlapsetup4Proc
#endif


! initialisation of global variables
    Subroutine initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize,maxg)
      Implicit None

      Integer, Intent(in) :: lmaxmat,lmaxapw,lmaxvr,gsize,maxg

      Integer l,m,lm
      Integer order1,order2,i,j,k

      If (.not.associated(input%groundstate)) Allocate(input%groundstate)
      input%groundstate%lmaxmat=lmaxmat
      input%groundstate%lmaxapw=lmaxapw
      input%groundstate%lmaxvr=lmaxvr
      lmmaxapw=(lmaxapw+1)**2
      lmmaxvr=(lmaxvr+1)**2
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
      nrmt(1:natmtot)=10
      apwordmax=1
      apword(:,:)=1

      If (allocated(apwfr)) Deallocate(apwfr)
      Allocate (apwfr(nrmt(1),2,apwordmax,0:lmaxapw,natmtot))

      If (allocated(apwdfr)) Deallocate(apwdfr)
      Allocate (apwdfr(apwordmax,0:lmaxapw,natmtot))

      idxas (:, :) = 1

      ! from *istln
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
       ivgig = 0
       ivg = 0

    End Subroutine initGlobals


! deallocation of global variables       
    Subroutine freeGlobals
      Implicit None

      Deallocate(haa)
      Deallocate(gntyry)
      Deallocate(apwfr,apwdfr)
      Deallocate(idxlm,input%groundstate)
      deallocate(ivg,ivgig)

    End Subroutine freeGlobals



! allocate and generate complex Gaunt coefficient array
! assumes that the global variables are set
! copied from init1.f90
    Subroutine initGntyry
      Implicit None

      Integer l1,m1,lm1,l2,m2,lm2,l3,m3,lm3

      Complex (8) gauntyry
      External gauntyry

      If (allocated(gntyry)) Deallocate (gntyry)
      Allocate (gntyry(lmmaxmat, lmmaxvr, lmmaxapw))
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do l2 = 0, input%groundstate%lmaxvr
               Do m2 = - l2, l2
                  lm2 = idxlm (l2, m2)
                  Do l3 = 0, input%groundstate%lmaxapw
                     Do m3 = - l3, l3
                        lm3 = idxlm (l3, m3)
                        gntyry (lm1, lm2, lm3) = gauntyry (l1, l2, l3, m1, m2, m3)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do


      If (allocated(haa)) Deallocate (haa)
      allocate(haa(apwordmax,0:input%groundstate%lmaxmat,apwordmax, &
     & 0:input%groundstate%lmaxapw, lmmaxvr, natmtot))

    End Subroutine initGntyry


!------------------------------------------------------------------------------
! !TEST: testHmlaan_EkinSurfaceSumL_Serial
!------------------------------------------------------------------------------
! !DESCRIPTION:
! 1st test, serial
! The surface part of the kinetic energy (line 3 in Eq. 24) is under inspection.
! The purpose is to test whether the update of the hamilton and the summation over l is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1 and the orbital momentum l1.
! The radial functions apwfr and apwdfr contain a dependence on the orbital momentum l1.
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
    Subroutine testHmlaanContrib_Serial

      Implicit None
! Size of the tests

      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp,maxg,maxgk
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,maxg=3,maxgk=1,gsize=(2*maxgk+1)**3,nmatp=gsize+10)

      Integer                  :: l1,m1,lm1,g1,g2
      Integer                  :: i,j
      Integer                  :: igpig(gsize)
      Real(8)                  :: vgpc(3,(2*maxgk+1)**3)
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Real(8)                  :: mc
      Type(evsystem)           :: system
      Type (HermitianMatrix)   :: hamilton_ref

! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize,maxg)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newsystem(system,nmatp,gsize)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
      haa(:,:,:,:,:,:)=0d0
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0
      vgpc(:,:)=0d0

      ! this block is from *listln
      j=0
      Do i=1,ngrtot
        If ((ivg(1,i)**2.le.maxgk**2).and.(ivg(2,i)**2.le.maxgk**2).and.(ivg(3,i)**2.le.maxgk**2)) then
          j=j+1
          igpig(j)=i
        End If
      EndDo

      Do l1=0,lmaxapw
        apwfr(nrmt(1),1,:,l1,1)=sqrt(dble(l1+1)) 
        apwdfr(:,l1,1)=sqrt(dble(l1+1))*dble(l1+1)
      Enddo

      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      apwalm(:,:,:,:)=0d0
      Do l1=0,lmaxapw
        Do m1=-l1,l1
          lm1=idxlm(l1,m1)
          Do g1=1,gsize
            mc=dble(g1)/(sqrt(dble(l1+1))*sqrt(dble(2*(lmaxmat+1)*(lmaxmat+2)*(2*l1+1))))
            apwalm(g1,:,lm1,1)=cmplx(mc,mc,8)
          Enddo 
        Enddo
      Enddo

      Do g2=1,gsize
        Do g1=1,gsize      
          hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
        Enddo
      Enddo
#ifndef MPI
Call printglobalmatrix(real(hamilton_ref%za), nmatp, nmatp, 'reference from hmlaan')
#endif

      Call hamiltonandoverlapsetup(system, system%ngp, apwalm, igpig, vgpc)

#ifndef MPI
Call printglobalmatrix(real(system%hamilton%za), nmatp, nmatp, 'system')
#endif
!       Call assert_equals(nmatp, system%hamilton%size, 'checking result rank')
!       Call assert_equals(nmatp, size(system%hamilton%za,1), 'checking result size rows')
!       Call assert_equals(nmatp, size(system%hamilton%za,2), 'checking result size cols')
!       Call assert_equals(hamilton_ref%za, system%hamilton%za, nmatp, nmatp, tol, 'checking result numbers')
!       Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking overlap=0')

! finalisation
      Call deletesystem(system)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmlaanContrib_Serial



end module modHamiltonandoverlapsetup_test
