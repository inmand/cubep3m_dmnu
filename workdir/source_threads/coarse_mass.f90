!! calculate coarse mesh density
subroutine coarse_mass
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  integer(4) :: i,j,k,k0,pp,ii,jj,kk
  integer(4), dimension(3) :: i1,i2
  real(4), dimension(3) :: x,dx1,dx2

  integer, parameter :: ijstart = 0
  integer, parameter :: ijstop  = nc_node_dim + 1 
    
  call system_clock(count=count_i)

# ifndef NEUTRINOS
  rho_c = ratio_omega_nu2m
# else
  rho_c= 0.0 !- mass_p * (mesh_scale / 2)**3  
# endif

  do k0 = 0, mesh_scale-1 
     !$omp parallel do schedule(dynamic) default(shared) private(i,j,k,pp)
     do k = k0, nc_node_dim + 1, mesh_scale 
        do j = ijstart, ijstop 
           do i = ijstart, ijstop 
              pp=hoc(i,j,k)
              if (i <= 1 .or. i >= nc_node_dim .or. &
                   j <= 1 .or. j >= nc_node_dim .or. &
                   k <= 1 .or. k >= nc_node_dim) then
                 call coarse_cic_mass_boundry(pp)
              else
                 call coarse_cic_mass(pp)
              endif
           enddo
        enddo
     enddo
     !$omp end parallel do
  enddo

  call system_clock(count=count_f,count_rate=count_r)
# ifdef MPI_TIME
  call mpi_time_analyze('cm  mass',real(count_f-count_i)/real(count_r),rank,nodes)
# else
  if (rank==0) write(*,*) 'coarse mass finished',real(count_f-count_i)/real(count_r)
# endif

end subroutine coarse_mass
