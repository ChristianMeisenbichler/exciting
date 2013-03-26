! Copyright (C) 2006-2008 C. Ambrosch-Draxl. C. Meisenbichler S. Sagmeister
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!
! !MODULE:  modmpi
! !DESCRIPTION:
!   MPI variables and interface functions
!   In case of compiled without MPI support it defines the
!   mpi specific variables such that the code behaves exactly as
!   the unmodified scalar version
!
!
! !REVISION HISTORY:
!   Created October 2006 (CHM)
!   Added wrapper routines, 2007-2008 (S. Sagmeister)
!   Added allgatherv interface, August 2010 (S. Sagmeister)
!
!
!
module modmpi

use modmpi_types
use setpatrtitioning

end module
