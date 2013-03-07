! !MODULE:  modHmlaan_test
! !DESCRIPTION:
! Modules with test and helper routines for testing hmlaan
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
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

    Double precision :: pi
    parameter (pi=3.1415926535897932385d0)

    contains


    Subroutine testcaseHmlaanSerial
      Implicit None

      Call set_test_name ('surface part of kinetic energy, summation over l')
      Call testHmlaan_EkinSurfaceSumL_Serial
      Call set_test_name ('surface part of kinetic energy, summation over m')
      Call testHmlaan_EkinSurfaceSumM_Serial
      Call set_test_name ('spheriCally symmetric contribution, summation over l')
      Call testHmlaan_SpherSymmSumL_Serial
      Call set_test_name ('spheriCally symmetric contribution, summation over m')
      Call testHmlaan_SpherSymmSumM_Serial
      Call set_test_name ('spheriCally symmetric and asymmetric contribution, gaunt coefficients')
      Call testHmlaan_SpherSymmAsymGnt_Serial
      Call set_test_name ('spheriCally symmetric and asymmetric contribution, summation over lm1 and lm3')
      Call testHmlaan_SpherSymmAsymSumLm1Lm3_Serial

    End Subroutine testcaseHmlaanSerial


#ifdef MPI
    Subroutine testcaseHmlaan1Proc
      Implicit None

      Call set_test_name ('surface part of kinetic energy, summation over l')
      Call testHmlaan_EkinSurfaceSumL_1Proc
      Call set_test_name ('surface part of kinetic energy, summation over m')
      Call testHmlaan_EkinSurfaceSumM_1Proc
      Call set_test_name ('spheriCally symmetric contribution, summation over l')
      Call testHmlaan_SpherSymmSumL_1Proc
      Call set_test_name ('spheriCally symmetric contribution, summation over m')
      Call testHmlaan_SpherSymmSumM_1Proc
      Call set_test_name ('spheriCally symmetric and asymmetric contribution, gaunt coefficients')
      Call testHmlaan_SpherSymmAsymGnt_1Proc
      Call set_test_name ('spheriCally symmetric and asymmetric contribution, summation over lm1 and lm3')
      Call testHmlaan_SpherSymmAsymSumLm1Lm3_1Proc

    End Subroutine testcaseHmlaan1Proc
#endif


#ifdef MPI
    Subroutine testcaseHmlaan4Proc
      Implicit None

      Call set_test_name ('surface part of kinetic energy, summation over l')
      Call testHmlaan_EkinSurfaceSumL_4Proc
      Call set_test_name ('surface part of kinetic energy, summation over m')
      Call testHmlaan_EkinSurfaceSumM_4Proc
      Call set_test_name ('spheriCally symmetric contribution, summation over l')
      Call testHmlaan_SpherSymmSumL_4Proc
      Call set_test_name ('spheriCally symmetric contribution, summation over m')
      Call testHmlaan_SpherSymmSumM_4Proc
      Call set_test_name ('spheriCally symmetric and asymmetric contribution, gaunt coefficients')
      Call testHmlaan_SpherSymmAsymGnt_4Proc
      Call set_test_name ('spheriCally symmetric and asymmetric contribution, summation over lm1 and lm3')
      Call testHmlaan_SpherSymmAsymSumLm1Lm3_4Proc

    End Subroutine testcaseHmlaan4Proc
#endif


! initialisation of global variables
    Subroutine initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)
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

    End Subroutine initGlobals


! deallocation of global variables       
    Subroutine freeGlobals
      Implicit None

      Deallocate(haa)
      Deallocate(gntyry)
      Deallocate(apwfr,apwdfr)
      Deallocate(idxlm,input%groundstate)

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
! The purpose is to test whether the update of the hamiltonian and the summation over l is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1 and the orbital momentum l1.
! The radial functions apwfr and apwdfr contain a dependence on the orbital momentum l1.
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
    Subroutine testHmlaan_EkinSurfaceSumL_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Real(8) :: mc
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0
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

      Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmlaan_EkinSurfaceSumL_Serial



!------------------------------------------------------------------------------
! test testHmlaan_EkinSurfaceSumM_Serial
!------------------------------------------------------------------------------
! 2nd test, serial
! The surface part of the kinetic energy (line 3 in Eq. 24) is under inspection.
! The purpose is to test whether the summation over m is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial functions apwfr and apwdfr are constant.
    Subroutine testHmlaan_EkinSurfaceSumM_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
     Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0
      Do l1=0,lmaxapw
        apwfr(nrmt(1),1,:,l1,1)=1d0
        apwdfr(:,l1,1)=1d0
      Enddo
      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      Do l1=0,lmaxapw
        Do m1=-l1,l1
          lm1=idxlm(l1,m1)
          Do g1=1,gsize
! The commented line below is another (simpler) option for the test
!          apwalm(g1,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(m1)/dble(2*l1+1)),sin(2d0*pi*dble(m1)/dble(2*l1+1)),8)/sqrt(2d0*dble(2*l1+1))/sqrt(dble((lmaxmat+1)))
            apwalm(g1,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(2d0*dble(lmmaxmat))
          Enddo
        Enddo
      Enddo

      Do g2=1,gsize
        Do g1=1,gsize
          hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
        Enddo
      Enddo

      Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')     
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmlaan_EkinSurfaceSumM_Serial



!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmSumL_Serial
!------------------------------------------------------------------------------
! 3rd test, serial
! The spheriCally symmetric contribution (line 2 in Eq. 24) to the hamiltonian is under inspection.
! The purpose is to test whether the update of the hamiltonian and the summation over l is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals contain a dependence on the orbital momentum l1.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmlaan_SpherSymmSumL_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
     Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0

      haa(1,:,1,:,1,1)=1d0
      Do l1=0,lmaxmat
        haa(1,l1,1,l1,1,1)=dble(l1+1)
      Enddo
      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      Do l1=0,lmaxapw
        Do m1=-l1,l1
          lm1=idxlm(l1,m1)
          Do g1=1,gsize
            apwalm(g1,:,lm1,1)=(4d0*pi)**0.25d0*g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat*(l1+1)))
          Enddo
        Enddo
      Enddo
      Do g2=1,gsize
        Do g1=1,gsize
          hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
        Enddo
      Enddo

      Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmlaan_SpherSymmSumL_Serial



!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmSumM_Serial
!------------------------------------------------------------------------------
! 4th test, serial
! The spheriCally symmetric contribution (line 2 in Eq. 24) to the hamiltonian is under inspection.
! The purpose is to test whether the summation over m is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmlaan_SpherSymmSumM_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0

      haa(1,:,1,:,1,1)=1d0

      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      Do l1=0,lmaxapw
        Do m1=-l1,l1
          lm1=idxlm(l1,m1)
          Do g1=1,gsize
            apwalm(g1,:,lm1,1)=(4d0*pi)**0.25d0*g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
          Enddo
        Enddo
      Enddo

      Do g2=1,gsize
        Do g1=1,gsize
          hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
        Enddo
      Enddo

      Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmlaan_SpherSymmSumM_Serial


!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmAsymGnt_Serial
!------------------------------------------------------------------------------
! 5th test, serial
! The spheriCally symmetric and asymmetric contributions (lines 2 and 4 in Eq. 24) to the hamiltonian are under inspection.
! The purpose is to test whether gaunt coefficients are handled properly.
! The matching coefficients apwalm contain are constant.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
    Subroutine testHmlaan_SpherSymmAsymGnt_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer :: l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8)               :: test
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! ! initialisation of global variables
     Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      haa(:,:,:,:,:,:)=1d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0

      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      Do l1=0,lmaxapw
        Do m1=-l1,l1
          lm1=idxlm(l1,m1)
          Do g1=1,gsize
            apwalm(g1,:,lm1,1)=cmplx(1d0,0,8)
          Enddo
        Enddo
      Enddo

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

      Do g2=1,gsize
        Do g1=1,gsize
          hamilton_ref%za(g1,g2)=test
        Enddo
      Enddo

      Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmlaan_SpherSymmAsymGnt_Serial


!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmAsymSumLm1Lm3_Serial
!------------------------------------------------------------------------------
! 6th test, serial
! The spheriCally symmetric and asymmetric contribution (lines 2 and 3 in Eq. 24) to the hamiltonian are under inspection.
! The purpose is to test whether the summation over lm1 and lm3 is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
   Subroutine testHmlaan_SpherSymmAsymSumLm1Lm3_Serial

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,lm2,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8) :: prefactor
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! ! initialisation of global variables
      Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
      Call initGntyry

      Call newmatrix(hamilton,nmatp)
      Call newmatrix(hamilton_ref,nmatp)
! initialisation is finished

      haa(:,:,:,:,:,:)=0d0

      hamilton%za(:,:)=cmplx(0,0,8)
      rmt(:)=2d0
      apwfr(:,:,:,:,:)=0d0
      apwdfr(:,:,:)=0d0
      
      haa(:,:,:,:,:,:)=cmplx(1,0,8)

      allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
      apwalm(:,:,:,:)=0d0
      Do l1=0,lmaxapw
        Do m1=-l1,l1
          lm1=idxlm(l1,m1)
            Do g1=1,gsize
              apwalm(g1,1,lm1,1)=cmplx(g1,0,8)*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
            Enddo
        Enddo
      Enddo

      prefactor=cmplx(0,0,8)
      Do lm2=1,lmmaxmat
        Do lm1=1,lmmaxmat
          prefactor=prefactor+sum(gntyry(lm1,:,lm2))*conjg(apwalm(1,1,lm1,1))*(apwalm(1,1,lm2,1))
        Enddo
      Enddo

      Do g2=1,gsize
        Do g1=1,gsize
          hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)*prefactor
        Enddo
      Enddo
      
      Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

      Call assert_equals(nmatp, hamilton%size, 'checking result rank')
      Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
      Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
      Call deletematrix(hamilton)
      Call deletematrix(hamilton_ref)
      Deallocate(apwalm)
! deallocation of global variables   
      Call freeGlobals    
    End Subroutine testHmlaan_SpherSymmAsymSumLm1Lm3_Serial


!------------------------------------------------------------------------------
! !TEST: testHmlaan_EkinSurfaceSumL_1Proc
!------------------------------------------------------------------------------
! !DESCRIPTION:
! 1st test, using MPI with only 1 proc
! The surface part of the kinetic energy (line 3 in Eq. 24) is under inspection.
! The purpose is to test whether the update of the hamiltonian and the summation over l is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1 and the orbital momentum l1.
! The radial functions apwfr and apwdfr contain a dependence on the orbital momentum l1.
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
#ifdef MPI
    Subroutine testHmlaan_EkinSurfaceSumL_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Real(8) :: mc
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
!Integer :: i

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal%blocksize = 2
      MPIglobal_1D%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)
! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0
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
! write (*,*) 'global apwalm'
! Do i=1,gsize
!   write (*,"(i,900f6.1)") i, real(apwalm(i,:,1:10,:))
! End Do


        Do g2=1,gsize
          Do g1=1,gsize      
            hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
          Enddo
        Enddo

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_EkinSurfaceSumL_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlaan_EkinSurfaceSumM_1Proc
!------------------------------------------------------------------------------
! 2nd test, using MPI with only 1 proc
! The surface part of the kinetic energy (line 3 in Eq. 24) is under inspection.
! The purpose is to test whether the summation over m is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial functions apwfr and apwdfr are constant.
#ifdef MPI
    Subroutine testHmlaan_EkinSurfaceSumM_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal%blocksize = 2
      MPIglobal_1D%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)
! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0
        Do l1=0,lmaxapw
          apwfr(nrmt(1),1,:,l1,1)=1d0
          apwdfr(:,l1,1)=1d0
        Enddo
        allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
        Do l1=0,lmaxapw
          Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
! The commented line below is another (simpler) option for the test
!               apwalm(g1,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(m1)/dble(2*l1+1)),sin(2d0*pi*dble(m1)/dble(2*l1+1)),8)/sqrt(2d0*dble(2*l1+1))/sqrt(dble((lmaxmat+1)))
              apwalm(g1,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(2d0*dble(lmmaxmat))
            Enddo
          Enddo
        Enddo

        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
          Enddo
        Enddo

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')     
        Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_EkinSurfaceSumM_1Proc
#endif



!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmSumL_1Proc
!------------------------------------------------------------------------------
! 3rd test, using MPI with only 1 proc
! The spheriCally symmetric contribution (line 2 in Eq. 24) to the hamiltonian is under inspection.
! The purpose is to test whether the update of the hamiltonian and the summation over l is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals contain a dependence on the orbital momentum l1.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlaan_SpherSymmSumL_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal%blocksize = 2
      MPIglobal_1D%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)
! initialisation is finished

        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0

        haa(1,:,1,:,1,1)=1d0
        Do l1=0,lmaxmat
          haa(1,l1,1,l1,1,1)=dble(l1+1)
        Enddo
        allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
        Do l1=0,lmaxapw
          Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
              apwalm(g1,:,lm1,1)=(4d0*pi)**0.25d0*g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat*(l1+1)))
            Enddo
          Enddo
        Enddo
        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
          Enddo
        Enddo

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_SpherSymmSumL_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmSumM_1Proc
!------------------------------------------------------------------------------
! 4th test, using MPI with only 1 proc
! The spheriCally symmetric contribution (line 2 in Eq. 24) to the hamiltonian is under inspection.
! The purpose is to test whether the summation over m is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlaan_SpherSymmSumM_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal%blocksize = 2
      MPIglobal_1D%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)
! initialisation is finished

        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0

        haa(1,:,1,:,1,1)=1d0

        allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
        Do l1=0,lmaxapw
          Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
              apwalm(g1,:,lm1,1)=(4d0*pi)**0.25d0*g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
            Enddo
          Enddo
        Enddo

        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)
          Enddo
        Enddo

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_SpherSymmSumM_1Proc
#endif


!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmAsymGnt_1Proc
!------------------------------------------------------------------------------
! 5th test, using MPI with only 1 proc
! The spheriCally symmetric and asymmetric contributions (lines 2 and 4 in Eq. 24) to the hamiltonian are under inspection.
! The purpose is to test whether gaunt coefficients are handled properly.
! The matching coefficients apwalm contain are constant.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlaan_SpherSymmAsymGnt_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer :: l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8)               :: test
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! Externals
      Complex(8), External :: gauntyry

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal%blocksize = 2
      MPIglobal_1D%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)
! initialisation is finished

        haa(:,:,:,:,:,:)=1d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0


        allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
        Do l1=0,lmaxapw
          Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1=1,gsize
              apwalm(g1,:,lm1,1)=cmplx(1d0,0,8)
            Enddo
          Enddo
        Enddo

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

        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref%za(g1,g2)=test
          Enddo
        Enddo

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_SpherSymmAsymGnt_1Proc
#endif



!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmAsymSumLm1Lm3_1Proc
!------------------------------------------------------------------------------
! 6th test, using MPI with only 1 proc
! The spheriCally symmetric and asymmetric contribution (lines 2 and 3 in Eq. 24) to the hamiltonian are under inspection.
! The purpose is to test whether the summation over lm1 and lm3 is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
   Subroutine testHmlaan_SpherSymmAsymSumLm1Lm3_1Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,lm2,g1,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8) :: prefactor
      Type (HermitianMatrix)   :: hamilton,hamilton_ref

! MPI variables
      Integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal%blocksize = 2
      MPIglobal_1D%blocksize = 2

      If (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)
! initialisation is finished

        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0
        
        haa(:,:,:,:,:,:)=cmplx(1,0,8)

        allocate(apwalm (gsize, apwordmax, lmmaxapw, natmtot))      
        apwalm(:,:,:,:)=0d0
        Do l1=0,lmaxapw
          Do m1=-l1,l1
          lm1=idxlm(l1,m1)
            Do g1=1,gsize
              apwalm(g1,1,lm1,1)=cmplx(g1,0,8)*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
            Enddo
          Enddo
        Enddo

        prefactor=cmplx(0,0,8)
        Do lm2=1,lmmaxmat
          Do lm1=1,lmmaxmat
            prefactor=prefactor+sum(gntyry(lm1,:,lm2))*conjg(apwalm(1,1,lm1,1))*(apwalm(1,1,lm2,1))
          Enddo
        Enddo

        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref%za(g1,g2)=cmplx(g1*g2,0,8)*prefactor
          Enddo
        Enddo
        
        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nmatp, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(nmatp, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nmatp, nmatp, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_SpherSymmAsymSumLm1Lm3_1Proc
#endif


!------------------------------------------------------------------------------
! !TEST: testHmlaan_EkinSurfaceSumL_4Proc
!------------------------------------------------------------------------------
! !DESCRIPTION:
! 1st test, 4 procs
! The surface part of the kinetic energy (line 3 in Eq. 24) is under inspection.
! The purpose is to test whether the update of the hamiltonian and the summation over l is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1 and the orbital momentum l1.
! The radial functions apwfr and apwdfr contain a dependence on the orbital momentum l1.
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC, A. Gulans)
!
#ifdef MPI
    Subroutine testHmlaan_EkinSurfaceSumL_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      Parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer                  :: l1,m1,lm1,g1,g1_idx,g2
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Real(8)                  :: mc
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! MPI related variables
      Integer :: n_procs_test, ierror_t
      Integer :: nrows_loc, ncols_loc
      Integer :: gsize_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc_idx, dummy

! Externals
      Integer, External   :: NUMROC

      n_procs_test = 4
      Call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      If (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

! ! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)

! init datastructures for splitting apwalm
        gsize_loc = NUMROC(gsize, MPIglobal_1D%blocksize, MPIglobal_1D%myproccol, 0, MPIglobal_1D%nproccols)
        allocate(dummy(gsize), apwalm1_loc_idx(gsize_loc))
        Call getLocalIndices(gsize, gsize, dummy, apwalm1_loc_idx, MPIglobal_1D)

! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0
        Do l1=0,lmaxapw
          apwfr(nrmt(1),1,:,l1,1)=sqrt(dble(l1+1)) 
          apwdfr(:,l1,1)=sqrt(dble(l1+1))*dble(l1+1)
        Enddo
        allocate(apwalm(gsize_loc, apwordmax, lmmaxapw, natmtot))      
        apwalm(:,:,:,:)=0d0
        Do l1=0,lmaxapw
          Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1_idx=1,gsize_loc !splitting of apwalm along first dimension!
              g1=apwalm1_loc_idx(g1_idx) 
              mc=dble(g1)/(sqrt(dble(l1+1))*sqrt(dble(2*(lmaxmat+1)*(lmaxmat+2)*(2*l1+1))))
              apwalm(g1_idx,:,lm1,1)=cmplx(mc,mc,8)
            Enddo 
          Enddo
        Enddo

        Do g2=1,gsize
          Do g1=1,gsize      
            hamilton_ref_global(g1,g2)=cmplx(g1*g2,0,8)
          Enddo
        Enddo
        Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 4
          case (1)
            ncols_loc = 3
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        End select

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize_loc)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm, apwalm1_loc_idx, dummy)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_EkinSurfaceSumL_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmlaan_EkinSurfaceSumM_4Proc
!------------------------------------------------------------------------------
! 2nd test, 4 procs
! The surface part of the kinetic energy (line 3 in Eq. 24) is under inspection.
! The purpose is to test whether the summation over m is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial functions apwfr and apwdfr are constant.
#ifdef MPI
    Subroutine testHmlaan_EkinSurfaceSumM_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer                  :: l1,m1,lm1,g1,g1_idx,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! MPI related variables
      Integer :: n_procs_test, ierror_t
      Integer :: nrows_loc, ncols_loc
      Integer :: gsize_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc_idx, dummy

! Externals
      Integer, External   :: NUMROC

      n_procs_test = 4
      Call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      If (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)

! init datastructures for splitting apwalm
        gsize_loc = NUMROC(gsize, MPIglobal_1D%blocksize, MPIglobal_1D%myproccol, 0, MPIglobal_1D%nproccols)
        allocate(dummy(gsize), apwalm1_loc_idx(gsize_loc))
        Call getLocalIndices(gsize, gsize, dummy, apwalm1_loc_idx, MPIglobal_1D)
! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0
        Do l1=0,lmaxapw
          apwfr(nrmt(1),1,:,l1,1)=1d0
          apwdfr(:,l1,1)=1d0
        Enddo
        allocate(apwalm(gsize_loc, apwordmax, lmmaxapw, natmtot))      
        Do l1=0,lmaxapw
          Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1_idx=1,gsize_loc !splitting of apwalm along first dimension!
              g1=apwalm1_loc_idx(g1_idx) 
              apwalm(g1_idx,:,lm1,1)=g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(2d0*dble(lmmaxmat))
            Enddo
          Enddo
        Enddo

        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref_global(g1,g2)=cmplx(g1*g2,0,8)
          Enddo
        Enddo
        Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 4
          case (1)
            ncols_loc = 3
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        End select

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize_loc)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm, apwalm1_loc_idx, dummy)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_EkinSurfaceSumM_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmSumL_4Proc
!------------------------------------------------------------------------------
! 3rd test , 4 procs
! The spheriCally symmetric contribution (line 2 in Eq. 24) to the hamiltonian is under inspection.
! The purpose is to test whether the update of the hamiltonian and the summation over l is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals contain a dependence on the orbital momentum l1.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlaan_SpherSymmSumL_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g1_idx,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! MPI related variables
      Integer :: n_procs_test, ierror_t
      Integer :: nrows_loc, ncols_loc
      Integer :: gsize_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc_idx, dummy

! Externals
      Integer, External   :: NUMROC

      n_procs_test = 4
      Call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      If (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)

! init datastructures for splitting apwalm
        gsize_loc = NUMROC(gsize, MPIglobal_1D%blocksize, MPIglobal_1D%myproccol, 0, MPIglobal_1D%nproccols)
        allocate(dummy(gsize), apwalm1_loc_idx(gsize_loc))
        Call getLocalIndices(gsize, gsize, dummy, apwalm1_loc_idx, MPIglobal_1D)

! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0

        haa(1,:,1,:,1,1)=1d0
        Do l1=0,lmaxmat
          haa(1,l1,1,l1,1,1)=dble(l1+1)
        Enddo
        allocate(apwalm(gsize_loc, apwordmax, lmmaxapw, natmtot))      
        apwalm(:,:,:,:)=0d0
        Do l1=0,lmaxapw
          Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1_idx=1,gsize_loc !splitting of apwalm along first dimension!
              g1=apwalm1_loc_idx(g1_idx) 
              apwalm(g1_idx,:,lm1,1)=(4d0*pi)**0.25d0*g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat*(l1+1)))
            Enddo
          Enddo
        Enddo

        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref_global(g1,g2)=cmplx(g1*g2,0,8)
          Enddo
        Enddo
        Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 4
          case (1)
            ncols_loc = 3
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        End select

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize_loc)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm, apwalm1_loc_idx, dummy)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_SpherSymmSumL_4Proc
#endif


!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmSumM_4Proc
!------------------------------------------------------------------------------
! 4th test, 4 procs
! The spheriCally symmetric contribution (line 2 in Eq. 24) to the hamiltonian is under inspection.
! The purpose is to test whether the summation over m is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlaan_SpherSymmSumM_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer l1,m1,lm1,g1,g1_idx,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! MPI related variables
      Integer :: n_procs_test, ierror_t
      Integer :: nrows_loc, ncols_loc
      Integer :: gsize_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc_idx, dummy

! Externals
      Integer, External   :: NUMROC

      n_procs_test = 4
      Call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      If (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)

! init datastructures for splitting apwalm
        gsize_loc = NUMROC(gsize, MPIglobal_1D%blocksize, MPIglobal_1D%myproccol, 0, MPIglobal_1D%nproccols)
        allocate(dummy(gsize), apwalm1_loc_idx(gsize_loc))
        Call getLocalIndices(gsize, gsize, dummy, apwalm1_loc_idx, MPIglobal_1D)

! initialisation is finished

! The line below makes sure that other terms apart from the surface kinetic energy 
! do not contribute.
        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0

        haa(1,:,1,:,1,1)=1d0

        allocate(apwalm (gsize_loc, apwordmax, lmmaxapw, natmtot))      
        Do l1=0,lmaxapw
          Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1_idx=1,gsize_loc !splitting of apwalm along first dimension!
              g1=apwalm1_loc_idx(g1_idx) 
              apwalm(g1_idx,:,lm1,1)=(4d0*pi)**0.25d0*g1*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
            Enddo
          Enddo
        Enddo

        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref_global(g1,g2)=cmplx(g1*g2,0,8)
          Enddo
        Enddo
        Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 4
          case (1)
            ncols_loc = 3
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        End select

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize_loc)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm, apwalm1_loc_idx, dummy)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_SpherSymmSumM_4Proc
#endif



!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmAsymGnt_4Proc
!------------------------------------------------------------------------------
! 5th test, 4 procs 
! The spheriCally symmetric and asymmetric contributions (lines 2 and 4 in Eq. 24) to the hamiltonian are under inspection.
! The purpose is to test whether gaunt coefficients are handled properly.
! The matching coefficients apwalm contain are constant.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
    Subroutine testHmlaan_SpherSymmAsymGnt_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer ::l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,g1,g1_idx,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8)               :: test
      Type (HermitianMatrix)   :: hamilton,hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! MPI related variables
      Integer :: n_procs_test, ierror_t
      Integer :: nrows_loc, ncols_loc
      Integer :: gsize_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc_idx, dummy

! Externals
      Integer,    External :: NUMROC
      Complex(8), External :: gauntyry

      n_procs_test = 4
      Call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      If (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)

! init datastructures for splitting apwalm
        gsize_loc = NUMROC(gsize, MPIglobal_1D%blocksize, MPIglobal_1D%myproccol, 0, MPIglobal_1D%nproccols)
        allocate(dummy(gsize), apwalm1_loc_idx(gsize_loc))
        Call getLocalIndices(gsize, gsize, dummy, apwalm1_loc_idx, MPIglobal_1D)

! initialisation is finished

        haa(:,:,:,:,:,:)=1d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0

        allocate(apwalm (gsize_loc, apwordmax, lmmaxapw, natmtot))      
        Do l1=0,lmaxapw
          Do m1=-l1,l1
            lm1=idxlm(l1,m1)
            Do g1_idx=1,gsize_loc !splitting of apwalm along first dimension!
              apwalm(g1_idx,:,lm1,1)=cmplx(1d0,0,8)
            Enddo
          Enddo
        Enddo

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

        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref_global(g1,g2)=test
          Enddo
        Enddo
        Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 4
          case (1)
            ncols_loc = 3
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        End select

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize_loc)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm, apwalm1_loc_idx, dummy)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_SpherSymmAsymGnt_4Proc
#endif



!------------------------------------------------------------------------------
! test testHmlaan_SpherSymmAsymSumLm1Lm3_4Proc
!------------------------------------------------------------------------------
! 6th test, 4 procs
! The spheriCally symmetric and asymmetric contribution (lines 2 and 3 in Eq. 24) to the hamiltonian are under inspection.
! The purpose is to test whether the summation over lm1 and lm3 is Done properly.
! The matching coefficients apwalm contain a dependence on the plane-wave index g1, the orbital momentum l1 and the magnetic quantum number m1.
! The radial integrals are constant.
! The table gntyry contains Gaunt coefficients as in a normal exciting run.
#ifdef MPI
   Subroutine testHmlaan_SpherSymmAsymSumLm1Lm3_4Proc

      Implicit None
! Size of the tests
      Integer lmaxmat,lmaxapw,lmaxvr,gsize,nmatp
      parameter (lmaxmat=2,lmaxapw=10,lmaxvr=6,gsize=9,nmatp=11)
      
      Integer :: l1,m1,lm1,lm2,g1,g1_idx,g2
      Complex (8), allocatable :: apwalm (:, :, :, :)
      complex(8)               :: prefactor
      Type (HermitianMatrix)   :: hamilton, hamilton_ref
      Complex (8), Dimension(nmatp,nmatp) :: hamilton_ref_global

! MPI related variables
      Integer :: n_procs_test, ierror_t
      Integer :: nrows_loc, ncols_loc
      Integer :: gsize_loc
      Integer, Dimension(:), Allocatable :: apwalm1_loc_idx, dummy

! Externals
      Integer, External   :: NUMROC

      n_procs_test = 4
      Call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      If (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

! initialisation of global variables
        Call initGlobals(lmaxmat,lmaxapw,lmaxvr,gsize)

! allocate and generate complex Gaunt coefficient array
        Call initGntyry

        Call newmatrix(hamilton, nmatp, DISTRIBUTE_COLS)
        Call newmatrix(hamilton_ref, nmatp, DISTRIBUTE_COLS)

! init datastructures for splitting apwalm
        gsize_loc = NUMROC(gsize, MPIglobal_1D%blocksize, MPIglobal_1D%myproccol, 0, MPIglobal_1D%nproccols)
        allocate(dummy(gsize), apwalm1_loc_idx(gsize_loc))
        Call getLocalIndices(gsize, gsize, dummy, apwalm1_loc_idx, MPIglobal_1D)

! initialisation is finished

        haa(:,:,:,:,:,:)=0d0

        hamilton%za(:,:)=cmplx(0,0,8)
        rmt(:)=2d0
        apwfr(:,:,:,:,:)=0d0
        apwdfr(:,:,:)=0d0
        
        haa(:,:,:,:,:,:)=cmplx(1,0,8)

        allocate(apwalm (gsize_loc, apwordmax, lmmaxapw, natmtot))      
        apwalm(:,:,:,:)=0d0
        Do l1=0,lmaxapw
          Do m1=-l1,l1
          lm1=idxlm(l1,m1)
            Do g1_idx=1,gsize_loc !splitting of apwalm along first dimension!
              g1=apwalm1_loc_idx(g1_idx) 
              apwalm(g1_idx,1,lm1,1)=cmplx(g1,0,8)*cmplx(cos(2d0*pi*dble(lm1)/dble(lmmaxmat)),sin(2d0*pi*dble(lm1)/dble(lmmaxmat)),8)/sqrt(dble(lmmaxmat))
            Enddo
          Enddo
        Enddo

        ! prefactor depends on an element of apwalm which is located at first proc!!
        if (MPIglobal_1D%rank .eq. 0) then
          prefactor=cmplx(0,0,8)
          Do lm2=1,lmmaxmat
            Do lm1=1,lmmaxmat
              prefactor=prefactor+sum(gntyry(lm1,:,lm2))*conjg(apwalm(1,1,lm1,1))*(apwalm(1,1,lm2,1))
            Enddo
          Enddo
        end if
        Call MPI_BCAST(prefactor, 2, MPI_DOUBLE_PRECISION, 0, MPIglobal_1D%comm, ierror_t) 

        Do g2=1,gsize
          Do g1=1,gsize
            hamilton_ref_global(g1,g2)=cmplx(g1*g2,0,8)*prefactor
          Enddo
        Enddo
        Call getBlockDistributedLoc(hamilton_ref_global, hamilton_ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 4
          case (1)
            ncols_loc = 3
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        End select

        Call hmlaan(hamilton,1,1,gsize,apwalm,gsize_loc)

        Call assert_equals(nmatp, hamilton%size, 'checking result rank')
        Call assert_equals(nrows_loc, size(hamilton%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(hamilton%za,2), 'checking result size cols')
        Call assert_equals(hamilton_ref%za, hamilton%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

! finalisation
        Call deletematrix(hamilton)
        Call deletematrix(hamilton_ref)
        Deallocate(apwalm, apwalm1_loc_idx, dummy)
! deallocation of global variables   
        Call freeGlobals    
! freeing proc grid
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      End If
    End Subroutine testHmlaan_SpherSymmAsymSumLm1Lm3_4Proc
#endif


end module modHmlaan_test
