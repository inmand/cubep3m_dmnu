!! cubep3m_dmnu - dm+nu cosmological N-body code
program cubep3m_dmnu
  use omp_lib
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  call mpi_initialize

  sec1 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STARTING CUBEP3M: ", sec1

  call read_remaining_time
  call variable_initialize
  call coarse_kernel
  call fine_kernel
  call particle_initialize(pid_dm)
  call link_list

  if (rank == 0) write(*,*) 'starting main loop'
  do 
    call timestep
    sec1a = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "TIMESTEP_TIME [hrs] = ", (sec1a - sec1) / 3600.

    call particle_mesh

    !! Determine if it is time to write a checkpoint before being killed
    if (checkpoint_step) then

       !! advance the particles to the end of the current step.
       dt_old = 0.0
       call update_position

       call move_grid_back(.true.)
       call link_list

       call checkpoint(kill_step)
       if (halofind_step) call halofind

       if (a.eq.a_i_nu) call particle_initialize(pid_nu)

       dt = 0.0

    endif

    if (nts == max_nts .or. final_step .or. a .gt. 1.0) exit

 enddo

  call cubepm_fftw(0)
  do ierr=1,cores 
    call cubepm_fftw2('q',ierr)
  enddo

  sec2 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STOPPING CUBEP3M: ", sec2
  if (rank == 0) write(*,*) "ELAPSED CUBEP3M TIME: ", sec2-sec1

  call mpi_finalize(ierr)

end program cubep3m_dmnu
