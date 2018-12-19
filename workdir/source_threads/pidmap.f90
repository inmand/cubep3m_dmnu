function pidmap(pid) result(id)
  implicit none
# include "cubepm.par"
  integer(kind=pid_int_kind), intent(in) :: pid
  integer(kind=pid_int_kind) :: id

# ifdef UNIQUE_PID
  if (pid.eq.pid_dm) then
     id=pid_dm
  else if (pid.gt.pid_dm) then
     id=pid_nu
  else
     id=pid_nu+1
  end if
#else
  id=pid
#endif

end function pidmap

subroutine particle_check
  implicit none
# include "cubepm.fh"
  integer(8) :: n,n_nu,n_dm,n_tot

  n_dm=0;n_nu=0
  do n=1,np_local
     if (pidmap(PID(n)).eq.pid_dm) n_dm=n_dm+1
     if (pidmap(PID(n)).eq.pid_nu) n_nu=n_nu+1
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
