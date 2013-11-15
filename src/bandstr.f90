! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: bandstr
! !INTERFACE:
!
!
Subroutine bandstr
  ! !USES:
      Use modinput
      Use modmain
      Use FoX_wxml
  ! !DESCRIPTION:
  !   Produces a band structure along the path in reciprocal-space which connects
  !   the vertices in the array {\tt vvlp1d}. The band structure is obtained from
  !   the second-variational eigenvalues and is written to the file {\tt BAND.OUT}
  !   with the Fermi energy set to zero. If required, band structures are plotted
  !   to files {\tt BAND\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species {\tt ss},
  !   which include the band characters for each $l$ component of that atom in
  !   columns 4 onwards. Column 3 contains the sum over $l$ of the characters.
  !   Vertex location lines are written to {\tt BANDLINES.OUT}.
  !
  ! !REVISION HISTORY:
  !   Created June 2003 (JKD)
  !   Modified June 2012 (DIN)
  !EOP
  !BOC
      Implicit None
  ! local variables
      Integer :: lmax, lmmax, l, m, lm
      Integer :: ik, ispn, is, ia, ias, iv, ist
      Real (8) :: emin, emax, sum
      Character (256) :: fname
  ! allocatable arrays
      Real (8), Allocatable :: evalfv (:, :)
      Real (8), Allocatable :: e (:, :)
  ! low precision for band character array saves memory
      Real (4), Allocatable :: bc (:, :, :, :)
      Complex (8), Allocatable :: dmat (:, :, :, :, :)
      Complex (8), Allocatable :: apwalm (:, :, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Character (128) :: buffer, buffer1
      Type (xmlf_t), Save :: xf
!
      integer :: nkpt0, nstsv0, Recl
      real(8),    Allocatable :: vkl0(:,:), ehf(:,:)
      complex(8), allocatable :: e0(:,:), e1(:,:)
      logical :: exist

  ! initialise universal variables
      Call init0
      Call init1

      if (associated(input%groundstate%HartreeFock)) then
        
!----------------------------------------      
! Calculate the HF bandstructure
!----------------------------------------
      
        fname='EVALHF.OUT'
        inquire(File=fname,Exist=exist)
        if (.not.exist) then
          write(*,*)'ERROR(bandstr.f90): File EVALHF.OUT does not exist!'
          stop
        end if
        
        inquire(IoLength=Recl) nkpt0, nstsv0
        open (70, File=fname, Action='READ', Form='UNFORMATTED', &
       &  Access='DIRECT', Recl=Recl)
        read(70, Rec=1) nkpt0, nstsv0
        close(70)
        
        nstsv=min(nstsv,nstsv0)
        
        allocate(vkl0(3,nkpt0))
        allocate(ehf(nstsv0,nkpt0))
        allocate(e0(nkpt0,nstsv))
        allocate(e1(nkpt,nstsv))
        
        inquire(IoLength=Recl) nkpt0, nstsv0, vkl0(:,1), ehf(:,1)
        open (70, File=fname, Action='READ', Form='UNFORMATTED', &
       &  Access='DIRECT', Recl=Recl)
        do ik = 1, nkpt0
          read(70, Rec=ik) nkpt0, nstsv0, vkl0(:,ik), ehf(1:,ik)
        end do ! ik
        close(70)
  
        ! Perform Fourier Interpolation
        do ik = 1, nkpt0
          e0(ik,1:nstsv)=cmplx(ehf(1:nstsv,ik),0.d0,8)
        enddo

        e1(:,:)=zzero
        call fourintp(e0,nkpt0,vkl0,e1,nkpt,vkl,nstsv)
        
        open(88,file='BAND-HF.OUT')
        do ist = 1, nstsv
          do ik = 1, nkpt
            write(88,*) dpp1d(ik), dble(e1(ik,ist))
          end do
          write(88,*)
        end do
        close(88)

        deallocate(vkl0,ehf,e0,e1)
              
      return
      end if 

  ! allocate array for storing the eigenvalues
      Allocate (e(nstsv, nkpt))
  ! maximum angular momentum for band character
      lmax = Min (3, input%groundstate%lmaxapw)
      lmmax = (lmax+1) ** 2
      If (input%properties%bandstructure%character) Then
         Allocate (bc(0:lmax, natmtot, nstsv, nkpt))
      End If
  ! read density and potentials from file
      Call readstate
  ! read Fermi energy from file
      Call readfermi
  ! find the new linearisation energies
      Call linengy
  ! generate the APW radial functions
      Call genapwfr
  ! generate the local-orbital radial functions
      Call genlofr
  ! compute the overlap radial integrals
      Call olprad
  ! compute the Hamiltonian radial integrals
      Call hmlint
      emin = 1.d5
      emax = - 1.d5
  ! begin parallel loop over k-points
#ifdef KSMP
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(evalfv,evecfv,evecsv) &
  !$OMP PRIVATE(dmat,apwalm) &
  !$OMP PRIVATE(ispn,ist,is,ia,ias) &
  !$OMP PRIVATE(l,m,lm,sum)
  !$OMP DO
#endif
      Do ik = 1, nkpt
         Allocate (evalfv(nstfv, nspnfv))
         Allocate (evecfv(nmatmax, nstfv, nspnfv))
         Allocate (evecsv(nstsv, nstsv))
     ! solve the first- and second-variational secular equations
         Call seceqn (ik, evalfv, evecfv, evecsv)
         Do ist = 1, nstsv
        ! subtract the Fermi energy
            e (ist, ik) = evalsv (ist, ik) - efermi
        ! add scissors correction
            If (e(ist, ik) .Gt. 0.d0) e (ist, ik) = e (ist, ik) + &
           & input%properties%bandstructure%scissor
        !$OMP CRITICAL
            emin = Min (emin, e(ist, ik))
            emax = Max (emax, e(ist, ik))
        !$OMP END CRITICAL
         End Do
     ! compute the band characters if required
         If (input%properties%bandstructure%character) Then
            Allocate (dmat(lmmax, lmmax, nspinor, nspinor, nstsv))
            Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, &
           & nspnfv))
        ! find the matching coefficients
            Do ispn = 1, nspnfv
               Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, &
              & ispn, ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, &
              & ispn))
            End Do
        ! average band character over spin and m for all atoms
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
              ! generate the diagonal of the density matrix
                  Call gendmat (.True., .True., 0, lmax, is, ia, ngk(:, &
                 & ik), apwalm, evecfv, evecsv, lmmax, dmat)
                  Do ist = 1, nstsv
                     Do l = 0, lmax
                        sum = 0.d0
                        Do m = - l, l
                           lm = idxlm (l, m)
                           Do ispn = 1, nspinor
                              sum = sum + dble (dmat(lm, lm, ispn, &
                             & ispn, ist))
                           End Do
                        End Do
                        bc (l, ias, ist, ik) = real (sum)
                     End Do
                  End Do
               End Do
            End Do
            Deallocate (dmat, apwalm)
         End If
         Deallocate (evalfv, evecfv, evecsv)
     ! end loop over k-points
      End Do
#ifdef KSMP
  !$OMP END DO
  !$OMP END PARALLEL
#endif
      emax = emax + (emax-emin) * 0.5d0
      emin = emin - (emax-emin) * 0.5d0
  ! output the band structure
      Call xml_OpenFile ("bandstructure.xml", xf, replace=.True., &
     & pretty_print=.True.)
!
Call xml_AddXMLPI(xf,"xml-stylesheet", 'href="'//trim(input%xsltpath)//&
      &'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')
      If ( .Not. input%properties%bandstructure%character) Then
         Open (50, File='BAND.OUT', Action='WRITE', Form='FORMATTED')

         Call xml_NewElement (xf, "bandstructure")
         Call xml_NewElement (xf, "title")
         Call xml_AddCharacters (xf, trim(input%title))
         Call xml_endElement (xf, "title")
         Do ist = 1, nstsv
            Call xml_NewElement (xf, "band")
            Do ik = 1, nkpt
               Write (50, '(2G18.10)') dpp1d (ik), e (ist, ik)
               Call xml_NewElement (xf, "point")
               Write (buffer, '(5G18.10)') dpp1d (ik)
               Call xml_AddAttribute (xf, "distance", &
              & trim(adjustl(buffer)))
               Write (buffer, '(5G18.10)') e (ist, ik)
               Call xml_AddAttribute (xf, "eval", &
              & trim(adjustl(buffer)))
               Call xml_endElement (xf, "point")
            End Do
            Call xml_endElement (xf, "band")
            Write (50, '("     ")')
         End Do
         Close (50)
         Write (*,*)
         Write (*, '("Info(bandstr):")')
         Write (*, '(" band structure plot written to BAND.OUT")')
      Else
         Call xml_NewElement (xf, "bandstructure")
         Call xml_AddAttribute (xf, "character", "true")
         Call xml_NewElement (xf, "title")
         Call xml_AddCharacters (xf, trim(input%title))
         Call xml_endElement (xf, "title")
         Do is = 1, nspecies
            Call xml_NewElement (xf, "species")
            Call xml_AddAttribute (xf, "name", trim(spname(is)))
            Call xml_AddAttribute (xf, "chemicalSymbol", trim(input%structure%speciesarray(is)%species%chemicalSymbol))
            Do ia = 1, natoms (is)
               Call xml_NewElement (xf, "atom")
               Write (buffer, '(5G18.10)') atposc (:, ia, is)
               Call xml_AddAttribute (xf, "coord", &
              & trim(adjustl(buffer)))
               ias = idxas (ia, is)
               Write (fname, '("BAND_S", I2.2, "_A", I4.4, ".OUT")') &
              & is, ia
               Open (50, File=trim(fname), Action='WRITE', Form='FORMAT&
              &TED')
!
               Do ist = 1, nstsv
                  Call xml_NewElement (xf, "band")
                  Do ik = 1, nkpt
                 ! sum band character over l
                     sum = 0.d0
                     Do l = 0, lmax
                        sum = sum + bc (l, ias, ist, ik)
                     End Do
                     Call xml_NewElement (xf, "point")
                     Write (buffer, '(5G18.10)') dpp1d (ik)
                     Call xml_AddAttribute (xf, "distance", &
                    & trim(adjustl(buffer)))
                     Write (buffer, '(5G18.10)') e (ist, ik)
                     Call xml_AddAttribute (xf, "eval", &
                    & trim(adjustl(buffer)))
                     Write (buffer, '(5G18.10)') sum
                     Call xml_AddAttribute (xf, "sum", &
                    & trim(adjustl(buffer)))
                     Do l = 0, lmax
                        Call xml_NewElement (xf, "bc")
                        Write (buffer,*) l
                        Call xml_AddAttribute (xf, "l", &
                       & trim(adjustl(buffer)))
                        Write (buffer, '(5G18.10)') bc (l, ias, ist, &
                       & ik)
                        Call xml_AddAttribute (xf, "character", &
                       & trim(adjustl(buffer)))
                        Call xml_endElement (xf, "bc")
                     End Do
                     Call xml_endElement (xf, "point")
                     Write (50, '(2G18.10, 8F12.6)') dpp1d (ik), e &
                    & (ist, ik), sum, (bc(l, ias, ist, ik), l=0, lmax)
                  End Do
                  Call xml_endElement (xf, "band")
                  Write (50, '("	  ")')
               End Do
               Call xml_endElement (xf, "atom")
               Close (50)
            End Do
            Call xml_endElement (xf, "species")
         End Do
         Write (*,*)
         Write (*, '("Info(bandstr):")')
         Write (*, '(" band structure plot written to BAND_Sss_Aaaaa.OU&
        &T")')
         Write (*, '("	for all species and atoms")')
      End If
      Write (*,*)
      Write (*, '(" Fermi energy is at zero in plot")')
  ! output the vertex location lines
      Open (50, File='BANDLINES.OUT', Action='WRITE', Form='FORMATTED')
      Do iv = 1, nvp1d
         Call xml_NewElement (xf, "vertex")
!
         Write (buffer, '(5G18.10)') dvp1d (iv)
         Call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
         Write (buffer, '(5G18.10)') emax
         Call xml_AddAttribute (xf, "upperboundary", &
        & trim(adjustl(buffer)))
         Write (buffer, '(5G18.10)') emin
         Call xml_AddAttribute (xf, "lowerboundary", &
        & trim(adjustl(buffer)))
         Call xml_AddAttribute (xf, "label", trim(adjustl(input%properties%bandstructure%plot1d%path%pointarray(iv)%point%label)))
         Write (buffer, '(5G18.10)') input%properties%bandstructure%plot1d%path%pointarray(iv)%point%coord
         Call xml_AddAttribute (xf, "coord", trim(adjustl(buffer)))
         Call xml_endElement (xf, "vertex")
         Write (50, '(2G18.10)') dvp1d (iv), emin
         Write (50, '(2G18.10)') dvp1d (iv), emax
         Write (50, '("     ")')
      End Do
      Close (50)
      Write (*,*)
      Write (*, '(" vertex location lines written to BANDLINES.OUT")')
      Write (*,*)
      Deallocate(e)
      If (input%properties%bandstructure%character) deallocate(bc)
      Call xml_close (xf)
      Return
End Subroutine bandstr
!EOC
