!! calculate coarse force on each node 
subroutine coarse_force
  implicit none
# include "cubepm.fh"

  integer(4) :: i,j,k,ii,im

  call system_clock(count=count_i)
  
  call cubepm_fftw(1)

  cmplx_rho_c=slab

  ! first do x direction
  !$omp parallel do default(shared) private(i,j,k,ii,im)
  do k=1,nc_pen+mypadd
     do j=1,nc_node_dim
        do i=1,nc_dim/2
           ii=2*i
           im=ii-1
           slab(im,j,k)=-cmplx_rho_c(ii,j,k)*kern_c(1,i,j,k)
           slab(ii,j,k)=cmplx_rho_c(im,j,k)*kern_c(1,i,j,k)
        enddo
     enddo
  enddo
  !$omp end parallel do

  call cubepm_fftw(-1)
  force_c(1,1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=rho_c

  ! now y direction
  !$omp parallel do default(shared) private(i,j,k,ii,im)
  do k=1,nc_pen+mypadd
     do j=1,nc_node_dim
        do i=1,nc_dim/2
           ii=2*i
           im=ii-1
           slab(im,j,k)=-cmplx_rho_c(ii,j,k)*kern_c(2,i,j,k)
           slab(ii,j,k)=cmplx_rho_c(im,j,k)*kern_c(2,i,j,k)
        enddo
     enddo
  enddo
  !$omp end parallel do

  call cubepm_fftw(-1)
  force_c(2,1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=rho_c

  ! now z direction
  !$omp parallel do default(shared) private(i,j,k,ii,im)
  do k=1,nc_pen+mypadd
     do j=1,nc_node_dim
        do i=1,nc_dim/2
           ii=2*i
           im=ii-1
           slab(im,j,k)=-cmplx_rho_c(ii,j,k)*kern_c(3,i,j,k)
           slab(ii,j,k)=cmplx_rho_c(im,j,k)*kern_c(3,i,j,k)
        enddo
     enddo
  enddo
  !$omp end parallel do

  call cubepm_fftw(-1)
  force_c(3,1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=rho_c

  call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
  call mpi_time_analyze('cm force',real(count_f-count_i)/real(count_r),rank,nodes)
#else
  if (rank==0) write(*,*) 'coarse force finished',real(count_f - count_i)/real(count_r)
#endif

end subroutine coarse_force
