!! update coarse mesh velocity
subroutine coarse_velocity
    use omp_lib
    implicit none

#    include "cubepm.fh"

    integer(4) :: i,j,k,pp
    integer(4), dimension(3) :: i1,i2
    real(4), dimension(3) :: x,dx1,dx2,dV

    call system_clock(count=count_i)

    !$omp parallel do default(shared) private(i,j,k,pp,x,i1,i2,dx1,dx2,dV)
    do k=1,nc_node_dim
      do j=1,nc_node_dim
        do i=1,nc_node_dim
          pp=hoc(i,j,k)
          do; if (pp == 0) exit
            x(:) = (1.0/real(mesh_scale)) * xv(1:3,pp) - 0.5
            i1(:) = floor(x(:)) + 1
            i2(:) = i1(:) + 1
#ifdef COARSE_NGP
            dx1(:) = 0.0
            dx2(:) = 1.0
#else
            dx1(:) = i1(:) - x(:)
            dx2(:) = 1.0 - dx1(:)
#endif
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
#ifdef DEBUG_VEL
            write(*,*) 'pp,aGdt,i1(3),i2(3),dx1(3),dx2(3),force_c(1,8)'
            write(*,*) pp,a_mid*G*dt,i1(:),i2(:),dx1(:),dx2(:), &
                       force_c(1,i1(1):i2(1),i1(2):i2(2),i1(3):i2(3))
#endif
            pp = ll(pp)
          enddo
        enddo
      enddo
    enddo
    !$omp end parallel do

    call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
    call mpi_time_analyze('cm   vel',real(count_f-count_i)/real(count_r),rank,nodes)
#else
    if (rank==0) write(*,*) 'coarse dark matter velocity finished' &
                            ,real(count_f - count_i)/real(count_r)
#endif

  end subroutine coarse_velocity
