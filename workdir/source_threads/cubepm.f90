!!  cubep3m - cubical decomposition 2-level particle mesh algorithm with particle-particle interactions
!! Hugh Merz :: merz@cita.utoronto.ca :: 2006 11 02 
program cubep3m
  use omp_lib
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  real(4) :: t_elapsed
  external t_elapsed

  real(8) :: sec1, sec2, sec01, sec02, seconds2wait
  real(8) :: sec1a, sec2a
  logical(kind=4) :: i_continue

  call mpi_initialize
  if (rank == 0) call datestamp

  sec1 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STARTING CUBEP3M: ", sec1

  call read_remaining_time
  call t_start(wc_counter)
  call variable_initialize
  call coarse_kernel
  call fine_kernel
  call particle_initialize
  call link_list

  if (rank == 0) write(*,*) 'starting main loop'
  do 
    call timestep
    sec1a = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "TIMESTEP_TIME [hrs] = ", (sec1a - sec1) / 3600.

    call particle_mesh

    !! Determine if it is time to write a checkpoint before being killed
!    kill_step = .false.
!    sec1a = mpi_wtime(ierr)
!    if (rank == 0) then
!        if ((sec1a - sec1) .ge. kill_time) kill_step = .true.
!    endif
!    call mpi_bcast(kill_step, 1, mpi_logical, 0, mpi_comm_world, ierr)

    if (checkpoint_step.or.projection_step.or.halofind_step.or.kill_step) then

       !! advance the particles to the end of the current step.
       dt_old = 0.0
       call update_position

!!$       force_grid_back = .true.
!!$       call move_grid_back
!!$       force_grid_back = .false.
       call link_list

       if (checkpoint_step .or. kill_step) then
          sec1a = mpi_wtime(ierr)
          if (rank == 0) write(*,*) "STARTING CHECKPOINT: ", sec1a
          call checkpoint(kill_step)
          if (rank == 0) write(*,*) 'finished checkpoint',t_elapsed(wc_counter)
          sec2a = mpi_wtime(ierr)
          if (rank == 0) write(*,*) "STOPPING CHECKPOINT: ", sec2a
          if (rank == 0) write(*,*) "ELAPSED CHECKPOINT TIME: ", sec2a-sec1a
       endif

       if (projection_step.or.halofind_step) then
          call particle_pass
          if (halofind_step) then
             sec1a = mpi_wtime(ierr)
             if (rank == 0) write(*,*) "STARTING HALOFIND: ", sec1a
             call halofind
             if (rank == 0) write(*,*) 'finished halofind',t_elapsed(wc_counter)
             sec2a = mpi_wtime(ierr)
             if (rank == 0) write(*,*) "STOPPING HALOFIND: ", sec2a
             if (rank == 0) write(*,*) "ELAPSED HALOFIND TIME: ", sec2a-sec1a
          endif
          
          if (projection_step) then
             call projection
             if (rank == 0) write(*,*) 'finished projection',t_elapsed(wc_counter)
          endif

          !! Clean up ghost particles
          call delete_particles

       endif

       dt = 0.0

    endif

    if (nts == max_nts .or. final_step .or. a .gt. 1.0) exit

 enddo

#ifdef TIMING
  if (rank==0) then
    print *,'cubep3m finished:' 
    call datestamp
  endif
#endif

  call cubepm_fftw(0)
  do ierr=1,cores 
    call cubepm_fftw2('q',ierr)
  enddo

  sec2 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STOPPING CUBEP3M: ", sec2
  if (rank == 0) write(*,*) "ELAPSED CUBEP3M TIME: ", sec2-sec1

  call mpi_finalize(ierr)

  if (rank == 0) call datestamp

end program cubep3m
