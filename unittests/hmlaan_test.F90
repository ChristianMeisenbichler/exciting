module modHmlaan_test
    use fruit
    use test_helpers
    use modfvsystem
#ifdef MPI
    use modmpi
#endif


    contains

!------------------------------------------------------------------------------
! test testHmlaanSerial
!------------------------------------------------------------------------------
    subroutine testHmlaanSerial
!       Use modinput,        Only: input
!       Use mod_muffin_tin,  Only: idxlm, lmmaxapw, rmt, nrmt
!       Use mod_APW_LO,      Only: apword, apwordmax, apwfr, apwdfr
!       Use mod_Gkvector,    Only: ngkmax
!       Use mod_atoms,       Only: natmtot, idxas
!       Use mod_eigensystem, Only: gntyry, haa
!       Use modfvsystem,     Only: HermitianMatrix


      implicit none


    end subroutine testHmlaanSerial



end module modHmlaan_test
