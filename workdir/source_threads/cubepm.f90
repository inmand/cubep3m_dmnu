!! cubep3m_dmnu - dm+nu cosmological N-body code
program cubep3m_dmnu
  use omp_lib
  implicit none
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

  call particle_check

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
     call particle_check

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

contains

  subroutine particle_check
    implicit none
    integer(8) :: n,n_nu,n_dm,n_tot

    n_dm=0;n_nu=0
    do n=1,np_local
       if (PID(n).eq.pid_dm) n_dm=n_dm+1
       if (PID(n).eq.pid_nu) n_nu=n_nu+1
    end do
    if (n_nu+n_dm.ne.np_local) write(*,*) "ERROR IN PARTICLE_CHECK",n_dm,n_nu
    if (rank.eq.0) write(*,*) ''
    call mpi_barrier(mpi_comm_world,ierr)
    write(*,*) 'rank,n_dm',rank,n_dm

    
    if (rank.eq.0) write(*,*) ''
    if (rank.eq.0) write(*,*) 'particle_check'
    call mpi_reduce(n_dm,n_tot,1,mpi_integer8,mpi_sum,0,mpi_comm_world,ierr)
    if (rank.eq.0) write(*,*) 'n_dm=',n_tot
    call mpi_reduce(n_nu,n_tot,1,mpi_integer8,mpi_sum,0,mpi_comm_world,ierr)
    if (rank.eq.0) write(*,*) 'n_nu=',n_tot
    n_dm=np_local
    call mpi_reduce(n_dm,n_tot,1,mpi_integer8,mpi_sum,0,mpi_comm_world,ierr)
    if (rank.eq.0) write(*,*) 'n_tot=',n_tot
    if (rank.eq.0) write(*,*) ''

  end subroutine particle_check

end program cubep3m_dmnu
