!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine timesec (ts)
      Use mod_timing
#ifdef MPI
      Use modmpi
#endif
      Implicit None
! arguments
      Real (8), Intent (Out) :: ts
! local variables
      Integer :: count, count_rate, dateA(8),dateB(8)
      real (8) time_cur
      Integer short_year(12),long_year(12),daysA,daysref,yearref
      data short_year /0,31,59,90,120,151,181,212,243,273,304,334/
      data long_year  /0,31,60,91,121,152,182,213,244,274,305,335/
    
      call date_and_time(values=dateA) 
#ifdef MPI
      time_cur = MPI_wtime()
#else
      Call system_clock (count=count, count_rate=count_rate)
      time_cur = dble (count) / dble (count_rate)
#endif
      call date_and_time(values=dateB)

      if ((dateA(3).ne.dateB(3)).and.(time_cur.lt.1d0)) then
        dateA=dateB
      endif


! How many days have elapsed since the reference date? 
      if (mod(date_ref(1),4).eq.0) then
       daysref=long_year(date_ref(2))+date_ref(3)
       yearref=366
      else
       daysref=short_year(date_ref(2))+date_ref(3)
       yearref=365
      endif 
      
      if (mod(dateA(1),4).eq.0) then
       daysA=long_year(dateA(2))+dateA(3)
      else
       daysA=short_year(dateA(2))+dateA(3)
      endif

! Is it the same year? 
! If it is not, we assume the calculation has started during the previous year.
      if (dateA(1).ne.date_ref(1)) then
        daysA=yearref+daysA
      endif

      ts=time_cur-time_ref+24d0*3600d0*dble(daysA-daysref)
      Return
End Subroutine
