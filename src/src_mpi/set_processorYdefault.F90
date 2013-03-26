

      ! taken from SIESTA
      subroutine set_processorYdefault(Nodes,procYdefault)

      ! Finds a sensible default value for the processor grid in the Y
      ! direction to try to get an equal split in X and Y. Currently only
      ! looks for factors of 2, 3 and 5.
      !
      ! Input :
      !
      ! integer Nodes        : total number of processors
      !
      ! Output :
      !
      ! integer procYdefault : default value of Y grid parameter
      !
      ! Written by Julian Gale, November 1999
      !
      ! Passed arguments
            integer :: Nodes, procYdefault
      ! Local variables
            integer :: Nx, Ny, Nrem
            logical :: factor

      ! Initialise values
            Nx = 1
            Ny = 1
            Nrem = Nodes
            factor = .true.

      ! Loop looking for factors
            do while (factor.and.Nrem.gt.1)
              factor = .false.
              if (mod(Nrem,2).eq.0) then
                Nrem = Nrem/2
                factor = .true.
                if (Nx.gt.Ny) then
                  Ny = 2*Ny
                else
                  Nx = 2*Nx
                endif
              endif
              if (mod(Nrem,3).eq.0) then
                Nrem = Nrem/3
                factor = .true.
                if (Nx.gt.Ny) then
                  Ny = 3*Ny
                else
                  Nx = 3*Nx
                endif
              endif
              if (mod(Nrem,5).eq.0) then
                Nrem = Nrem/5
                factor = .true.
                if (Nx.gt.Ny) then
                  Ny = 5*Ny
                else
                  Nx = 5*Nx
                endif
              endif
            enddo

      ! Choose default value as lowest of Nx and Ny
            if (Nx.gt.Ny) then
              procYdefault = Ny
            else
              procYdefault = Nx
            endif

      end subroutine set_processorYdefault
