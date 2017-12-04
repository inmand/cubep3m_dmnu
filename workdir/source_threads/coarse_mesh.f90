!! coarse mesh velocity update
subroutine coarse_mesh
  implicit none
  include 'mpif.h'
#   include "cubepm.fh"

#ifdef DIAG
  integer(4) :: i,j,k
  real(8) :: sumrhoc,sumrhof
#endif

   call coarse_mass

#ifdef DIAG
    sumrhof=0.0
    sumrhoc=0.0
    do k=1,nc_node_dim
      do j=1,nc_node_dim
        do i=1,nc_node_dim
          sumrhof=sumrhof+real(rho_c(i,j,k),kind=8)
        enddo
      enddo
    enddo
    call mpi_reduce(sumrhof,sumrhoc,1,mpi_double_precision, &
                    mpi_sum,0,mpi_comm_world,ierr)
    if (rank == 0) write(*,*) 'sum of rho_c=',sumrhoc
    if (rank == 0) write(*,*) 'rank', rank, 'min/max rho_c =', minval(rho_c), maxval(rho_c)
#endif

   call coarse_force
   call coarse_force_buffer
   call coarse_max_dt
   if (coarse_vel_update) call coarse_velocity

end subroutine coarse_mesh
