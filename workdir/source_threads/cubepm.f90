!! cubep3m_dmnu - dm+nu cosmological N-body code
program cubep3m_dmnu
  use omp_lib
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  call mpi_initialize

  sec1 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "starting cubep3m"

  call read_remaining_time
  call variable_initialize
  call coarse_kernel
  call fine_kernel
  call particle_initialize(pid_dm)
  if (injection_step) call particle_initialize(pid_nu)
  call link_list

  if (rank == 0) then
     write(*,*) 
     sec1a = mpi_wtime(ierr)
     write(*,*) 'time taken [hrs] = ',(sec1a-sec1)/3600.
     write(*,*)
     write(*,*) 'starting main loop'
  end if
  
  do 

     call timestep
     call particle_mesh

     if (checkpoint_step) then

        dt_old = 0.0
        call update_position

        call move_grid_back(.true.)
        call link_list

        call checkpoint(kill_step)
        if (halofind_step) call halofind

        if (injection_step) call particle_initialize(pid_nu)

        dt = 0.0

     endif

     if (nts == max_nts .or. final_step .or. a .gt. 1.0) exit
     
  enddo

  call cubepm_fftw(0)
  do ierr=1,cores 
     call cubepm_fftw2('q',ierr)
  enddo

  sec2 = mpi_wtime(ierr)
  if (rank == 0) then
     write(*,*) "stopping cubep3m"
     write(*,*) 'time taken [hrs] = ',(sec2-sec1)/3600.
  end if

  call mpi_finalize(ierr)

end program cubep3m_dmnu
