module modHmlaan_test
    Use fruit
    Use test_helpers
    Use modinput,        Only: input
    Use mod_muffin_tin,  Only: idxlm, lmmaxapw,lmmaxvr,lmmaxmat, rmt, nrmt
    Use mod_APW_LO,      Only: apword, apwordmax, apwfr, apwdfr
    Use mod_Gkvector,    Only: ngkmax
    Use mod_atoms,       Only: natmtot, idxas
    Use mod_eigensystem, Only: gntyry, haa
    Use modfvsystem,     Only: HermitianMatrix,newmatrix,deletematrix
#ifdef MPI
    use modmpi
#endif

    double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains


    subroutine testcaseHmlaanSerial
      Implicit None

      CALL set_test_name ('')
      CALL set_test_name ('surface part of kinetic energy, summation over l')
      CALL testHmlaan_EkinSurfaceSumL_Serial
      CALL set_test_name ('surface part of kinetic energy, summation over m')
      CALL testHmlaan_EkinSurfaceSumM_Serial
      CALL set_test_name ('spherically symmetric contribution, summation over l')
      CALL testHmlaan_SpherSymmSumL_Serial
      CALL set_test_name ('spherically symmetric contribution, summation over m')
      CALL testHmlaan_SpherSymmSumM_Serial
      CALL set_test_name ('spherically symmetric and asymmetric contribution, gaunt coefficients')
      CALL testHmlaan_SpherSymmAsymGnt_Serial
      CALL set_test_name ('spherically symmetric and asymmetric contribution, summation over lm1 and lm3')
      CALL testHmlaan_SpherSymmAsymSumLm1Lm3_Serial

    end subroutine testcaseHmlaanSerial


! initialisation of global variables
    subroutine initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)
      Implicit None

      Integer, Intent(in) :: lmaxmat,lmaxapw,lmaxvr,gsize

      Integer l,m,lm

      If (.not.associated(input%groundstate)) Allocate(input%groundstate)
      input%groundstate%lmaxmat=lmaxmat
      input%groundstate%lmaxapw=lmaxapw
      input%groundstate%lmaxvr=lmaxvr
      lmmaxapw=(lmaxapw+1)**2
      lmmaxvr=(lmaxvr+1)**2
      lmmaxmat=(lmaxmat+1)**2
! copied from init0.f90     
      If (allocated(idxlm)) deallocate (idxlm)
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

      If (allocated(apwfr)) deallocate(apwfr)
      Allocate (apwfr(nrmt(1),2,apwordmax,0:lmaxapw,natmtot))

      If (allocated(apwdfr)) deallocate(apwdfr)
      Allocate (apwdfr(apwordmax,0:lmaxapw,natmtot))

      idxas (:, :) = 1

    end subroutine initGlobals


! deallocation of global variables       
    subroutine freeGlobals
      Implicit None

      Deallocate(haa)
      Deallocate(gntyry)
      Deallocate(apwfr,apwdfr)
      Deallocate(idxlm,input%groundstate)

    end subroutine freeGlobals



! allocate and generate complex Gaunt coefficient array
! assumes that the global variables are set
! copied from init1.f90
    subroutine initGntyry
      Implicit None

      Integer l1,m1,lm1,l2,m2,lm2,l3,m3,lm3

      Complex (8) gauntyry
      External gauntyry

      If (allocated(gntyry)) deallocate (gntyry)
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


      if (allocated(haa)) deallocate (haa)
      allocate(haa(apwordmax,0:input%groundstate%lmaxmat,apwordmax, &
     & 0:input%groundstate%lmaxapw, lmmaxvr, natmtot))

    end subroutine initGntyry


!------------------------------------------------------------------------------
! test testHmlaan_EkinSurfaceSumL_Serial
!------------------------------------------------------------------------------
! 1st test 
! The surface part of the kinetic energy (line 3 in Eq. 24) is under inspection.
! The purpose is to test whether the update of the hamiltonian and the summation over l is done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1 and the orbital momentum l1.
! The radial functions apwfr and apwdfr contain a dependence on the orbital momentum l1.

    subroutine testHmlaan_EkinSurfaceSumL_Serial

      Implicit None
! Size of the tests
      integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Real(8) :: mc
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
     Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call NewMatrix(hamilton,nmatp)
      Call NewMatrix(hamilton_ref,nmatp)
! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0
      do l1=0,lmaxapw
        apwfr(nrmt(1),1,:,l1,1)=sqrt(dble(l1+1)) 
        apwdfr(:,l1,1)=sqrt(dble(l1+1))*dble(l1+1)
      enddo
      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      apwalm(:,:,:,:)=0d0
      do l1=0,lmaxapw
       do m1=-l1,l1
        lm1=idxlm(l1,m1)
        do g1=1,gsize
          mc=dble(g1)/(sqrt(dble(l1+1))*sqrt(dble(2*(lmaxmat+1)*(lmaxmat+2)*(2*l1+1))))
          apwalm(g1,:,lm1,1)=cmplx(mc,mc,8)
        enddo 
       enddo
      enddo

      do g2=1,gsize
       do g1=1,gsize      
         hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
       enddo
      enddo
      Call hmlaan(hamilton,1,1,gsize,apwalm)
      CALL assert_equals(nmatp, hamilton%size, 'checking result rank')
      CALL assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      CALL assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      call DeleteMatrix(hamilton)
      call DeleteMatrix(hamilton_ref)
      deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    end subroutine testHmlaan_EkinSurfaceSumL_Serial



!------------------------------------------------------------------------------
! test testHmlaan_EkinSurfaceSumM_Serial
!------------------------------------------------------------------------------
! 2nd test
! The surface part of the kinetic energy (line 3 in Eq. 24) is under inspection.
! The purpose is to test whether the summation over m is done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial functions apwfr and apwdfr are constant.
    subroutine testHmlaan_EkinSurfaceSumM_Serial

      Implicit None
! Size of the tests
      integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
     Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call NewMatrix(hamilton,nmatp)
      Call NewMatrix(hamilton_ref,nmatp)
! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0
      do l1=0,lmaxapw
        apwfr(nrmt(1),1,:,l1,1)=1d0
        apwdfr(:,l1,1)=1d0
      enddo
      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      do l1=0,lmaxapw
       do m1=-l1,l1
        lm1=idxlm(l1,m1)
        do g1=1,gsize
! The commented line below is another (simpler) option for the test
!          apwalm(g1,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(m1)/dble(2*l1+1)),sin(2d0*pi*dble(m1)/dble(2*l1+1)),8)/sqrt(2d0*dble(2*l1+1))/sqrt(dble((lmaxmat+1)))
           apwalm(g1,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(2d0*dble(lmmaxmat))
        enddo
       enddo
      enddo

      do g2=1,gsize
       do g1=1,gsize
         hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
       enddo
      enddo
      Call hmlaan(hamilton,1,1,gsize,apwalm)
      CALL assert_equals(nmatp, hamilton%size, 'checking result rank')
      CALL assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      CALL assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')     
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      call DeleteMatrix(hamilton)
      call DeleteMatrix(hamilton_ref)
      deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    end subroutine testHmlaan_EkinSurfaceSumM_Serial



!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmSumL_Serial
!------------------------------------------------------------------------------
! 3rd test 
! The spherically symmetric contribution (line 2 in Eq. 24) to the hamiltonian is under inspection.
! The purpose is to test whether the update of the hamiltonian and the summation over l is done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals contain a dependence on the orbital momentum l1.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    subroutine testHmlaan_SpherSymmSumL_Serial

      Implicit None
! Size of the tests
      integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
     Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call NewMatrix(hamilton,nmatp)
      Call NewMatrix(hamilton_ref,nmatp)
! initialisation is finished

      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0

      haa(1,:,1,:,1,1)=1d0
      do l1=0,lmaxmat
        haa(1,l1,1,l1,1,1)=dble(l1+1)
      enddo
      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      do l1=0,lmaxapw
       do m1=-l1,l1
        lm1=idxlm(l1,m1)
        do g1=1,gsize
           apwalm(g1,:,lm1,1)=(4d0*pi)**0.25d0*g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat*(l1+1)))
        enddo
       enddo
      enddo
      do g2=1,gsize
       do g1=1,gsize
         hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
       enddo
      enddo
      Call hmlaan(hamilton,1,1,gsize,apwalm)
      CALL assert_equals(nmatp, hamilton%size, 'checking result rank')
      CALL assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      CALL assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      call DeleteMatrix(hamilton)
      call DeleteMatrix(hamilton_ref)
      deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    end subroutine testHmlaan_SpherSymmSumL_Serial



!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmSumM_Serial
!------------------------------------------------------------------------------
! 4th test
! The spherically symmetric contribution (line 2 in Eq. 24) to the hamiltonian is under inspection.
! The purpose is to test whether the summation over m is done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    subroutine testHmlaan_SpherSymmSumM_Serial

      Implicit None
! Size of the tests
      integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
     Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call NewMatrix(hamilton,nmatp)
      Call NewMatrix(hamilton_ref,nmatp)
! initialisation is finished

      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0

      haa(1,:,1,:,1,1)=1d0

      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      do l1=0,lmaxapw
       do m1=-l1,l1
        lm1=idxlm(l1,m1)
        do g1=1,gsize
           apwalm(g1,:,lm1,1)=(4d0*pi)**0.25d0*g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
        enddo
       enddo
      enddo

      do g2=1,gsize
       do g1=1,gsize
         hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
       enddo
      enddo
      Call hmlaan(hamilton,1,1,gsize,apwalm)
      CALL assert_equals(nmatp, hamilton%size, 'checking result rank')
      CALL assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      CALL assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      call DeleteMatrix(hamilton)
      call DeleteMatrix(hamilton_ref)
      deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    end subroutine testHmlaan_SpherSymmSumM_Serial


!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmAsymGnt_Serial
!------------------------------------------------------------------------------
! 5th test
! The spherically symmetric and asymmetric contributions (lines 2 and 4 in Eq. 24) to the hamiltonian are under inspection.
! The purpose is to test whether gaunt coefficients are handled properly.
! The matching coefficients apwalm contain are constant.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    subroutine testHmlaan_SpherSymmAsymGnt_Serial

      Implicit None
! Size of the tests
      integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      integer l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8) :: test
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

      Complex (8) gauntyry
      External gauntyry

! ! initialisation of global variables
     Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call NewMatrix(hamilton,nmatp)
      Call NewMatrix(hamilton_ref,nmatp)
! initialisation is finished

      haa(:,:,:,:,:,:)=1d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0


      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      do l1=0,lmaxapw
       do m1=-l1,l1
        lm1=idxlm(l1,m1)
        do g1=1,gsize
           apwalm(g1,:,lm1,1)=cmplx(1d0,0,8)
        enddo
       enddo
      enddo

      test=cmplx(0,0,8)
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do l2 = 0, input%groundstate%lmaxvr
               Do m2 = - l2, l2
                  lm2 = idxlm (l2, m2)
                  Do l3 = 0, input%groundstate%lmaxmat
                     Do m3 = - l3, l3
                        lm3 = idxlm (l3, m3)
                        test = test + gauntyry (l1, l2, l3, m1, m2, m3)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do

      do g2=1,gsize
       do g1=1,gsize
         hamilton_ref%za(g1,g2)=test
       enddo
      enddo
      Call hmlaan(hamilton,1,1,gsize,apwalm)
      CALL assert_equals(nmatp, hamilton%size, 'checking result rank')
      CALL assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      CALL assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      call DeleteMatrix(hamilton)
      call DeleteMatrix(hamilton_ref)
      deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    end subroutine testHmlaan_SpherSymmAsymGnt_Serial


!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmAsymSumLm1Lm3_Serial
!------------------------------------------------------------------------------
! 6th test
! The spherically symmetric and asymmetric contribution (lines 2 and 3 in Eq. 24) to the hamiltonian are under inspection.
! The purpose is to test whether the summation over lm1 and lm3 is done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
   subroutine testHmlaan_SpherSymmAsymSumLm1Lm3_Serial

      Implicit None
! Size of the tests
      integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      integer l1,m1,lm1,lm2,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8) :: prefactor
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
     Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call NewMatrix(hamilton,nmatp)
      Call NewMatrix(hamilton_ref,nmatp)
! initialisation is finished

      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0
      
      haa(:,:,:,:,:,:)=cmplx(1,0,8)

      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      apwalm(:,:,:,:)=0d0
      do l1=0,lmaxapw
       do m1=-l1,l1
        lm1=idxlm(l1,m1)
        do g1=1,gsize
           apwalm(g1,1,lm1,1)=cmplx(g1,0,8)*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
        enddo
       enddo
      enddo

      prefactor=cmplx(0,0,8)
      do lm2=1,lmmaxmat
       do lm1=1,lmmaxmat
         prefactor=prefactor+sum(gntyry(lm1,:,lm2))*conjg(apwalm(1,1,lm1,1))*(apwalm(1,1,lm2,1))
       enddo
      enddo

      do g2=1,gsize
       do g1=1,gsize
         hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)*prefactor
       enddo
      enddo
      
      Call hmlaan(hamilton,1,1,gsize,apwalm)
      CALL assert_equals(nmatp, hamilton%size, 'checking result rank')
      CALL assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      CALL assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      call DeleteMatrix(hamilton)
      call DeleteMatrix(hamilton_ref)
      deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    end subroutine testHmlaan_SpherSymmAsymSumLm1Lm3_Serial

end module modHmlaan_test
