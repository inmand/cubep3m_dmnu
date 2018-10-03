!Preprocessor commands
#undef COARSE_NGP

!! coarse mesh velocity update
subroutine coarse_mesh
  implicit none
# include "cubepm.fh"

# if VERBOSITY>0
  if (rank.eq.0) write(*,*) ':: coarse forces'
# endif

  call coarse_mass
  call coarse_force
  call coarse_force_buffer
  call coarse_max_dt
  call coarse_velocity

end subroutine coarse_mesh

!! calculate coarse mesh density
subroutine coarse_mass
  implicit none
# include "cubepm.fh"

  integer(4) :: i,j,k,k0,pp,ii,jj,kk
  integer(4), dimension(3) :: i1,i2
  real(4), dimension(3) :: x,dx1,dx2

  integer, parameter :: ijstart = 0
  integer, parameter :: ijstop  = nc_node_dim + 1 
    
  rho_c = f_unclustered

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

end subroutine coarse_mass

!! add mass to coarse mesh density
subroutine coarse_cic_mass(pp)
  use omp_lib
  implicit none
# include "cubepm.fh"

  integer(4) :: pp
  integer(4), dimension(3) :: i1,i2
  real(4), dimension(3) :: x,dx1,dx2

  do
     if (pp == 0) exit
     x(:) = (1.0/real(mesh_scale)) * xv(1:3,pp) - 0.5
     i1(:) = floor(x(:)) + 1
     i2(:) = i1(:) + 1
#    ifdef COARSE_NGP
     dx1(:) = 0.0
     dx2(:) = 1.0
#    else
     dx1(:) = i1(:) - x(:)
     dx2(:) = 1.0 - dx1(:)
#    endif

     dx1(1) = mass_p * dx1(1) * mass_p_nudm_fac(PID(pp)) 
     dx2(1) = mass_p * dx2(1) * mass_p_nudm_fac(PID(pp))

     rho_c(i1(1),i1(2),i1(3)) = rho_c(i1(1),i1(2),i1(3)) &
          + dx1(1) * dx1(2) * dx1(3)
     
     rho_c(i2(1),i1(2),i1(3)) = rho_c(i2(1),i1(2),i1(3)) &
          + dx2(1) * dx1(2) * dx1(3)
     
     rho_c(i1(1),i2(2),i1(3)) = rho_c(i1(1),i2(2),i1(3)) &
          + dx1(1) * dx2(2) * dx1(3)

     rho_c(i2(1),i2(2),i1(3)) = rho_c(i2(1),i2(2),i1(3)) &
          + dx2(1) * dx2(2) * dx1(3)

     rho_c(i1(1),i1(2),i2(3)) = rho_c(i1(1),i1(2),i2(3)) &
          + dx1(1) * dx1(2) * dx2(3)
     
     rho_c(i2(1),i1(2),i2(3)) = rho_c(i2(1),i1(2),i2(3)) &
          + dx2(1) * dx1(2) * dx2(3)
     
     rho_c(i1(1),i2(2),i2(3)) = rho_c(i1(1),i2(2),i2(3)) &
          + dx1(1) * dx2(2) * dx2(3)
     
     rho_c(i2(1),i2(2),i2(3)) = rho_c(i2(1),i2(2),i2(3)) &
          + dx2(1) * dx2(2) * dx2(3)
     
     pp = ll(pp)
  enddo
  
end subroutine coarse_cic_mass

!! add mass to coarse mesh density along nodal boundry
subroutine coarse_cic_mass_boundry(pp)
  use omp_lib
  implicit none 
# include "cubepm.fh"

  integer(4) :: pp
  integer(4), dimension(3) :: i1,i2
  real(4), dimension(3) :: x,dx1,dx2

  do
     if (pp == 0) exit
     x(:) = (1.0/real(mesh_scale)) * xv(1:3,pp) - 0.5
     i1(:) = floor(x(:)) + 1
     i2(:) = i1(:) + 1
#    ifdef COARSE_NGP
     dx1(:) = 0.0
     dx2(:) = 1.0
#    else
     dx1(:) = i1(:) - x(:)
     dx2(:) = 1.0 - dx1(:)
#    endif

     dx1(1) = mass_p * dx1(1) * mass_p_nudm_fac(PID(pp)) 
     dx2(1) = mass_p * dx2(1) * mass_p_nudm_fac(PID(pp))

     if (i1(3) >= 1 .and. i1(3) <= nc_node_dim) then
        if (i1(2) >= 1 .and. i1(2) <= nc_node_dim) then
           if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then
              rho_c(i1(1),i1(2),i1(3)) = rho_c(i1(1),i1(2),i1(3)) + &
                   dx1(1) * dx1(2) * dx1(3)
           endif
           if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then
              rho_c(i2(1),i1(2),i1(3)) = rho_c(i2(1),i1(2),i1(3)) + &
                   dx2(1) * dx1(2) * dx1(3)
           endif
        endif
        if (i2(2) >= 1 .and. i2(2) <= nc_node_dim) then
           if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then
              rho_c(i1(1),i2(2),i1(3)) = rho_c(i1(1),i2(2),i1(3)) + &
                   dx1(1) * dx2(2) * dx1(3)
           endif
           if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then
              rho_c(i2(1),i2(2),i1(3)) = rho_c(i2(1),i2(2),i1(3)) + &
                   dx2(1) * dx2(2) * dx1(3)
           endif
        endif
     endif

     if (i2(3) >= 1 .and. i2(3) <= nc_node_dim) then
        if (i1(2) >= 1 .and. i1(2) <= nc_node_dim) then
           if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then
              rho_c(i1(1),i1(2),i2(3)) = rho_c(i1(1),i1(2),i2(3)) + &
                   dx1(1) * dx1(2) * dx2(3)
           endif
           if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then
              rho_c(i2(1),i1(2),i2(3)) = rho_c(i2(1),i1(2),i2(3)) + &
                   dx2(1) * dx1(2) * dx2(3)
           endif
        endif
        if (i2(2) >= 1 .and. i2(2) <= nc_node_dim) then
           if (i1(1) >= 1 .and. i1(1) <= nc_node_dim) then
              rho_c(i1(1),i2(2),i2(3)) = rho_c(i1(1),i2(2),i2(3)) + &
                   dx1(1) * dx2(2) * dx2(3)
           endif
           if (i2(1) >= 1 .and. i2(1) <= nc_node_dim) then
              rho_c(i2(1),i2(2),i2(3)) = rho_c(i2(1),i2(2),i2(3)) + &
                   dx2(1) * dx2(2) * dx2(3)
           endif
        endif
     endif

     pp = ll(pp)
  enddo

end subroutine coarse_cic_mass_boundry

!! calculate coarse force on each node 
subroutine coarse_force
  implicit none
# include "cubepm.fh"

  integer(4) :: i,j,k,ii,im

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

end subroutine coarse_force

!! pass coarse mesh force along boundries to adjacent nodes
subroutine coarse_force_buffer
  implicit none
# include "cubepm.fh"

  integer(4) :: buffer_size
  integer(4) :: tag
  integer(4) :: status(MPI_STATUS_SIZE)

  buffer_size = 3 * (nc_node_dim + 2)**2
  tag=64

  force_c_buffer(:,:,:)=force_c(:,1,:,:)
  !send to node in -x 
  call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
       cart_neighbor(5),tag,cart_neighbor(6), &
       tag,mpi_comm_cart,status,ierr)
  force_c(:,nc_node_dim+1,:,:)=force_c_buffer(:,:,:)
  
  force_c_buffer(:,:,:)=force_c(:,nc_node_dim,:,:)
  !send to node in +x
  call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
       cart_neighbor(6),tag,cart_neighbor(5), &
       tag,mpi_comm_cart,status,ierr)
  force_c(:,0,:,:)=force_c_buffer(:,:,:)
  
  force_c_buffer(:,:,:)=force_c(:,:,1,:)
  !send to node in -y  
  call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
       cart_neighbor(3),tag,cart_neighbor(4), &
       tag,mpi_comm_cart,status,ierr)
  force_c(:,:,nc_node_dim+1,:)=force_c_buffer(:,:,:)
  
  force_c_buffer(:,:,:)=force_c(:,:,nc_node_dim,:)
  !send to node in +y 
  call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
       cart_neighbor(4),tag,cart_neighbor(3), &
       tag,mpi_comm_cart,status,ierr)
  force_c(:,:,0,:)=force_c_buffer(:,:,:)
  
  force_c_buffer(:,:,:)=force_c(:,:,:,1)
  !send to node in -z 
  call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
       cart_neighbor(1),tag,cart_neighbor(2), &
       tag,mpi_comm_cart,status,ierr)
  force_c(:,:,:,nc_node_dim+1)=force_c_buffer(:,:,:)
  
  force_c_buffer(:,:,:)=force_c(:,:,:,nc_node_dim)
  !send to node in +z 
  call mpi_sendrecv_replace(force_c_buffer,buffer_size,mpi_real, &
       cart_neighbor(2),tag,cart_neighbor(1), &
       tag,mpi_comm_cart,status,ierr)
  force_c(:,:,:,0)=force_c_buffer(:,:,:)

end subroutine coarse_force_buffer

!! calculate the maximum dt based on the coarse mesh force
subroutine coarse_max_dt 
  implicit none
# include "cubepm.fh"

  integer(kind=4) :: i,j,k
  real(kind=4) :: force,max_force

  max_force=0.0

  !$omp parallel do default(shared) &
  !$omp private(i,j,k,force) &
  !$omp reduction(MAX:max_force)
  do k=1,nc_node_dim
     do j=1,nc_node_dim
        do i=1,nc_node_dim
           force=sqrt(force_c(1,i,j,k)**2+force_c(2,i,j,k)**2 + &
                force_c(3,i,j,k)**2)
           if (force.gt.max_force) max_force=force
        enddo
     enddo
  enddo
  !$omp end parallel do

  call mpi_reduce(max_force,dt_c_acc,1,mpi_real,mpi_max,0, &
       mpi_comm_world,ierr)
  dt_c_acc=sqrt(real(mesh_scale,kind=4)/(dt_c_acc*a_mid*G))
  call mpi_bcast(dt_c_acc,1,mpi_real,0,mpi_comm_world,ierr)
  
end subroutine coarse_max_dt

!! update coarse mesh velocity
subroutine coarse_velocity
  use omp_lib
  implicit none
# include "cubepm.fh"

  integer(4) :: i,j,k,pp
  integer(4), dimension(3) :: i1,i2
  real(4), dimension(3) :: x,dx1,dx2,dV

  !$omp parallel do default(shared) private(i,j,k,pp,x,i1,i2,dx1,dx2,dV)
  do k=1,nc_node_dim
     do j=1,nc_node_dim
        do i=1,nc_node_dim
           pp=hoc(i,j,k)
           do; if (pp == 0) exit
              x(:) = (1.0/real(mesh_scale)) * xv(1:3,pp) - 0.5
              i1(:) = floor(x(:)) + 1
              i2(:) = i1(:) + 1
#             ifdef COARSE_NGP
              dx1(:) = 0.0
              dx2(:) = 1.0
#             else
              dx1(:) = i1(:) - x(:)
              dx2(:) = 1.0 - dx1(:)
#             endif

              dV = a_mid * G * dt * dx1(1) * dx1(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) + force_c(:,i1(1),i1(2),i1(3)) * dV
              dV = a_mid * G * dt * dx2(1) * dx1(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) + force_c(:,i2(1),i1(2),i1(3)) * dV
              dV = a_mid * G * dt * dx1(1) * dx2(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) + force_c(:,i1(1),i2(2),i1(3)) * dV
              dV = a_mid * G * dt * dx2(1) * dx2(2) * dx1(3)
              xv(4:6,pp) = xv(4:6,pp) + force_c(:,i2(1),i2(2),i1(3)) * dV
              dV = a_mid * G * dt * dx1(1) * dx1(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) + force_c(:,i1(1),i1(2),i2(3)) * dV
              dV = a_mid * G * dt * dx2(1) * dx1(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) + force_c(:,i2(1),i1(2),i2(3)) * dV
              dV = a_mid * G * dt * dx1(1) * dx2(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) + force_c(:,i1(1),i2(2),i2(3)) * dV
              dV = a_mid * G * dt * dx2(1) * dx2(2) * dx2(3)
              xv(4:6,pp) = xv(4:6,pp) + force_c(:,i2(1),i2(2),i2(3)) * dV
              pp = ll(pp)
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do

end subroutine coarse_velocity
