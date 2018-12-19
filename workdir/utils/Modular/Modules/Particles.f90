module Particles
  use Parameters
  use Variables
  use mMPI

contains

  subroutine pass_particles(xvp,np,min_x,max_x,buf_x)
    implicit none
    real, dimension(:,:), intent(inout) :: xvp
    integer, intent(inout) :: np
    real, intent(in) :: min_x,max_x,buf_x

    real, dimension(:,:), allocatable :: xvb
    integer :: n,max_np, nb,max_nb,nbp

    integer :: tag=29,status(MPI_STATUS_SIZE)

    if (size(xvp,dim=1).ne.6) then
       write(*,*) 'size(x,dim=1)=',size(xvp,dim=1),6
       call particles_error_stop('[Subroutine - pass_particles] size(xcp,dim=1) != 6')
    end if

    max_np=size(xvp,dim=2)
    if (np.gt.max_np) then
       write(*,*) 'np,max_np=',np,max_np
       call particles_error_stop('[Subroutine - pass_particles] np>max_np')
    end if

    if (minval(xvp(1:3,:np)).lt.min_x) then
       write(*,*) 'minval(xvp(1:3,:np)),min_x=',minval(xvp(1:3,:np)),min_x
       call particles_error_stop('[Subroutine - pass_particles] minval(xvp(1:3,:np)<min_x')
    end if

    if (maxval(xvp(1:3,:np)).gt.max_x) then
       write(*,*) 'maxval(xvp(1:3,:np)),max_x=',maxval(xvp(1:3,:np)),max_x
       call particles_error_stop('[Subroutine - pass_particles] maxval(xvp(1:3,:np)>max_x')
    end if

    nb=max_np-np
    call mpi_allreduce(nb,max_nb,1,mpi_integer,mpi_max,mpi_comm_world,ierr)
    nb=0
    allocate(xvb(6,max_nb))

    !-x
    xvb=0;nb=0
    do n=1,np
       if (xvp(1,n).lt.min_x+buf_x) then
          nb=nb+1
          if (nb.gt.max_nb) then
             write(*,*) 'nb,max_nb=',nb,max_nb
             call particles_error_stop('[Subroutine - pass_particles] nb>max_nb in -x pass')
          end if
          xvb(1:6,nb)=xvp(1:6,n)
       end if
    end do

    call mpi_sendrecv_replace(xvb,size(xvb),mpi_real,cart_neighbor(5),tag,cart_neighbor(6),tag,mpi_comm_cart,status,ierr)
    call mpi_sendrecv_replace(nb,1,mpi_integer,cart_neighbor(5),tag,cart_neighbor(6),tag,mpi_comm_cart,status,ierr)

    do n=1,nb
       np=np+1
       if (np.gt.max_np) then
          write(*,*) 'np,max_np=',np,max_np
          call particles_error_stop('[Subroutine - pass_particles] np>max_np in -x pass')
       end if
       xvp(1:6,np)=xvb(1:6,n)
       xvp(1,np)=xvp(1,np)+(max_x-min_x)
    end do

    !+x
    xvb=0;nb=0
    do n=1,np
       if (xvp(1,n).gt.max_x-buf_x .and. xvp(1,n).lt.max_x) then
          nb=nb+1
          if (nb.gt.max_nb) then
             write(*,*) 'nb,max_nb=',nb,max_nb
             call particles_error_stop('[Subroutine - pass_particles] nb>max_nb in +x pass')
          end if
          xvb(1:6,nb)=xvp(1:6,n)
       end if
    end do

    call mpi_sendrecv_replace(xvb,size(xvb),mpi_real,cart_neighbor(6),tag,cart_neighbor(5),tag,mpi_comm_cart,status,ierr)
    call mpi_sendrecv_replace(nb,1,mpi_integer,cart_neighbor(6),tag,cart_neighbor(5),tag,mpi_comm_cart,status,ierr)          

    do n=1,nb
       np=np+1
       if (np.gt.max_np) then
          write(*,*) 'np,max_np=',np,max_np
          call particles_error_stop('[Subroutine - pass_particles] np>max_np in +x pass')
       end if
       xvp(1:6,np)=xvb(1:6,n)
       xvp(1,np)=xvp(1,np)-(max_x-min_x)
    end do

    !-y
    xvb=0;nb=0
    do n=1,np
       if (xvp(2,n).lt.min_x+buf_x) then
          nb=nb+1
          if (nb.gt.max_nb) then
             write(*,*) 'nb,max_nb=',nb,max_nb
             call particles_error_stop('[Subroutine - pass_particles] nb>max_nb in -y pass')
          end if
          xvb(1:6,nb)=xvp(1:6,n)
       end if
    end do

    call mpi_sendrecv_replace(xvb,size(xvb),mpi_real,cart_neighbor(3),tag,cart_neighbor(4),tag,mpi_comm_cart,status,ierr)
    call mpi_sendrecv_replace(nb,1,mpi_integer,cart_neighbor(3),tag,cart_neighbor(4),tag,mpi_comm_cart,status,ierr)

    do n=1,nb
       np=np+1
       if (np.gt.max_np) then
          write(*,*) 'np,max_np=',np,max_np
          call particles_error_stop('[Subroutine - pass_particles] np>max_np in -y pass')
       end if
       xvp(1:6,np)=xvb(1:6,n)
       xvp(2,np)=xvp(2,np)+(max_x-min_x)
    end do

    !+y
    xvb=0;nb=0
    do n=1,np
       if (xvp(2,n).gt.max_x-buf_x .and. xvp(2,n).lt.max_x) then
          nb=nb+1
          if (nb.gt.max_nb) then
             write(*,*) 'nb,max_nb=',nb,max_nb
             call particles_error_stop('[Subroutine - pass_particles] nb>max_nb in +y pass')
          end if
          xvb(1:6,nb)=xvp(1:6,n)
       end if
    end do

    call mpi_sendrecv_replace(xvb,size(xvb),mpi_real,cart_neighbor(4),tag,cart_neighbor(3),tag,mpi_comm_cart,status,ierr)
    call mpi_sendrecv_replace(nb,1,mpi_integer,cart_neighbor(4),tag,cart_neighbor(3),tag,mpi_comm_cart,status,ierr)          

    do n=1,nb
       np=np+1
       if (np.gt.max_np) then
          write(*,*) 'np,max_np=',np,max_np
          call particles_error_stop('[Subroutine - pass_particles] np>max_np in +y pass')
       end if
       xvp(1:6,np)=xvb(1:6,n)
       xvp(2,np)=xvp(2,np)-(max_x-min_x)
    end do

    !-z
    xvb=0;nb=0
    do n=1,np
       if (xvp(3,n).lt.min_x+buf_x) then
          nb=nb+1
          if (nb.gt.max_nb) then
             write(*,*) 'nb,max_nb=',nb,max_nb
             call particles_error_stop('[Subroutine - pass_particles] nb>max_nb in -z pass')
          end if
          xvb(1:6,nb)=xvp(1:6,n)
       end if
    end do

    call mpi_sendrecv_replace(xvb,size(xvb),mpi_real,cart_neighbor(1),tag,cart_neighbor(2),tag,mpi_comm_cart,status,ierr)
    call mpi_sendrecv_replace(nb,1,mpi_integer,cart_neighbor(1),tag,cart_neighbor(2),tag,mpi_comm_cart,status,ierr)

    do n=1,nb
       np=np+1
       if (np.gt.max_np) then
          write(*,*) 'np,max_np=',np,max_np
          call particles_error_stop('[Subroutine - pass_particles] np>max_np in -z pass')
       end if
       xvp(1:6,np)=xvb(1:6,n)
       xvp(3,np)=xvp(3,np)+(max_x-min_x)
    end do

    !+z
    xvb=0;nb=0
    do n=1,np
       if (xvp(3,n).gt.max_x-buf_x .and. xvp(3,n).lt.max_x) then
          nb=nb+1
          if (nb.gt.max_nb) then
             write(*,*) 'nb,max_nb=',nb,max_nb
             call particles_error_stop('[Subroutine - pass_particles] nb>max_nb in +z pass')
          end if
          xvb(1:6,nb)=xvp(1:6,n)
       end if
    end do

    call mpi_sendrecv_replace(xvb,size(xvb),mpi_real,cart_neighbor(2),tag,cart_neighbor(1),tag,mpi_comm_cart,status,ierr)
    call mpi_sendrecv_replace(nb,1,mpi_integer,cart_neighbor(2),tag,cart_neighbor(1),tag,mpi_comm_cart,status,ierr)          

    do n=1,nb
       np=np+1
       if (np.gt.max_np) then
          write(*,*) 'np,max_np=',np,max_np
          call particles_error_stop('[Subroutine - pass_particles] np>max_np in +z pass')
       end if
       xvp(1:6,np)=xvb(1:6,n)
       xvp(3,np)=xvp(3,np)-(max_x-min_x)
    end do

    deallocate(xvb)

  end subroutine pass_particles

  subroutine delete_particles(xvp,np,min_x,max_x)
    implicit none
    real, dimension(:,:), intent(inout) :: xvp
    integer, intent(inout) :: np
    real, intent(in) :: min_x,max_x
    integer :: n

    n=1
    do 
       if (n.gt.np) exit 
       if (minval(xvp(1:3,n)).lt.min_x .or. maxval(xvp(1:3,n)).gt.max_x) then
          xvp(1:6,n)=xvp(1:6,np)
          np=np-1
          cycle
       end if
       n=n+1
    end do

  end subroutine delete_particles

  subroutine particles_error_stop(expl)
    implicit none
    character(len=*), intent(in) :: expl
    write(*,*) '[MOD - Particles] rank=',rank
    write(*,*) '-->'//expl
    call mpi_abort(mpi_comm_world, ierr, ierr)
  end subroutine particles_error_stop

end module Particles
