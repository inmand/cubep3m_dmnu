module Interpolation
  use Parameters
  use Variables
  use mMPI

  logical, parameter :: interpolation_verbosity=.true.

contains

  subroutine ngp_interpolate(xp,gn,wp)
    implicit none
    real, dimension(:,:), intent(in) :: xp
    real, dimension(:,:,:), intent(inout) :: gn
    real, dimension(:), intent(in), optional :: wp

    integer :: n,np,nc
    integer, dimension(3) :: i1,i2
    real :: w
    real, dimension(3) :: x,dx1,dx2

    real(8), dimension(:,:,:), allocatable :: bg
    real(8), dimension(:,:), allocatable :: bgs
    integer :: bs,tag,status(MPI_STATUS_SIZE)

    real(8) :: ncl,ncg,npl,npg,sgl,sgg

    !Determine number of particles, np, and number of cells, nc
    if (size(xp,dim=1).ne.3) call interpolation_error_stop('[Subroutine - cic_interpolate] size(x,dim=1)!=3')
    np = size(xp,dim=2)
    if (present(wp)) then
       if (size(wp).ne.np) call interpolation_error_stop('[Subroutine - cic_interpolate] size(w)!=np')
    end if

    nc=size(gn,dim=1)
    if (nc.gt.Ncells .or. nc.ne.size(gn,dim=2) .or. nc.ne.size(gn,dim=3)) &
         &call interpolation_error_stop('[Subroutine - cic_interpolate] nc incorrect')

    if (minval(xp).lt.0) call interpolation_error_stop('[Subroutine - cic_interpolate] minval(xp)<0')
    if (maxval(xp).gt.nc) call interpolation_error_stop('[Subroutine - cic_interpolate] maxval(xp)>nc')

    allocate(bg(0:nc+1,0:nc+1,0:nc+1),bgs(0:nc+1,0:nc+1))

    !Interpolate particles to bg array
    bg=0.
    do n=1,np

       x=xp(1:3,n)-0.5
       i1=1+floor(x)
       i2=1+i1
       dx1=0.!i1-x
       dx2=1.-dx1

       if (present(wp)) then
          w=wp(n)
       else
          w=1.
       end if

       bg(i1(1),i1(2),i1(3))=bg(i1(1),i1(2),i1(3))+dx1(1)*dx1(2)*dx1(3)*w
       bg(i2(1),i1(2),i1(3))=bg(i2(1),i1(2),i1(3))+dx2(1)*dx1(2)*dx1(3)*w

       bg(i1(1),i2(2),i1(3))=bg(i1(1),i2(2),i1(3))+dx1(1)*dx2(2)*dx1(3)*w
       bg(i2(1),i2(2),i1(3))=bg(i2(1),i2(2),i1(3))+dx2(1)*dx2(2)*dx1(3)*w

       bg(i1(1),i1(2),i2(3))=bg(i1(1),i1(2),i2(3))+dx1(1)*dx1(2)*dx2(3)*w
       bg(i2(1),i1(2),i2(3))=bg(i2(1),i1(2),i2(3))+dx2(1)*dx1(2)*dx2(3)*w

       bg(i1(1),i2(2),i2(3))=bg(i1(1),i2(2),i2(3))+dx1(1)*dx2(2)*dx2(3)*w
       bg(i2(1),i2(2),i2(3))=bg(i2(1),i2(2),i2(3))+dx2(1)*dx2(2)*dx2(3)*w

    end do

    !Buffer bg array
    !!cart_neighbor(1)=-z
    !!cart_neighbor(2)=+z
    !!cart_neighbor(3)=-y
    !!cart_neighbor(4)=+y
    !!cart_neighbor(5)=-x
    !!cart_neighbor(6)=+x

    bs=size(bgs)
    tag=29

    !-x
    bgs=bg(0,:,:)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(5),tag,cart_neighbor(6),tag,mpi_comm_cart,status,ierr)
    bg(nc,:,:)=bg(nc,:,:)+bgs

    !+x
    bgs=bg(nc+1,:,:)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(6),tag,cart_neighbor(5),tag,mpi_comm_cart,status,ierr)
    bg(1,:,:)=bg(1,:,:)+bgs

    !-y
    bgs=bg(:,0,:)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(3),tag,cart_neighbor(4),tag,mpi_comm_cart,status,ierr)
    bg(:,nc,:)=bg(:,nc,:)+bgs

    !+y
    bgs=bg(:,nc+1,:)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(4),tag,cart_neighbor(3),tag,mpi_comm_cart,status,ierr)
    bg(:,1,:)=bg(:,1,:)+bgs

    !-z
    bgs=bg(:,:,0)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(1),tag,cart_neighbor(2),tag,mpi_comm_cart,status,ierr)
    bg(:,:,nc)=bg(:,:,nc)+bgs
    
    !+z
    bgs=bg(:,:,nc+1)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(2),tag,cart_neighbor(1),tag,mpi_comm_cart,status,ierr)
    bg(:,:,1)=bg(:,:,1)+bgs

    !Normalize
    ncl=(nc*1.d0)**3
    call mpi_allreduce(ncl,ncg,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

    npl=np*1.d0
    call mpi_allreduce(npl,npg,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

    sgl=sum(bg(1:nc,1:nc,1:nc))
    call mpi_allreduce(sgl,sgg,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

    if (interpolation_verbosity .and. rank.eq.0) then
       write(*,*) 'number particles: ',npg
       write(*,*) 'grid sum: ',sgg
       write(*,*) 'number of cells: ',ncg
       write(*,*) 'grid sum/number of cells: ',sgg/ncg
    end if

    bg=bg*ncg/sgg
    
    sgl=sum(bg(1:nc,1:nc,1:nc))
    call mpi_allreduce(sgl,sgg,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

    if (interpolation_verbosity .and. rank.eq.0) then
       write(*,*) 'Normalized'
       write(*,*) 'grid sum: ',sgg
       write(*,*) 'grid sum/number of cells: ',sgg/ncg
    end if

    gn=real(bg(1:nc,1:nc,1:nc),kind=4)

    deallocate(bg,bgs)

  end subroutine ngp_interpolate

  subroutine cic_interpolate(xp,gn,wp)
    implicit none
    real, dimension(:,:), intent(in) :: xp
    real, dimension(:,:,:), intent(inout) :: gn
    real, dimension(:), intent(in), optional :: wp

    integer :: n,np,nc
    integer, dimension(3) :: i1,i2
    real :: w
    real, dimension(3) :: x,dx1,dx2

    real(8), dimension(:,:,:), allocatable :: bg
    real(8), dimension(:,:), allocatable :: bgs
    integer :: bs,tag,status(MPI_STATUS_SIZE)

    real(8) :: ncl,ncg,npl,npg,sgl,sgg

    !Determine number of particles, np, and number of cells, nc
    if (size(xp,dim=1).ne.3) call interpolation_error_stop('[Subroutine - cic_interpolate] size(x,dim=1)!=3')
    np = size(xp,dim=2)
    if (present(wp)) then
       if (size(wp).ne.np) call interpolation_error_stop('[Subroutine - cic_interpolate] size(w)!=np')
    end if

    nc=size(gn,dim=1)
    if (nc.gt.Ncells .or. nc.ne.size(gn,dim=2) .or. nc.ne.size(gn,dim=3)) &
         &call interpolation_error_stop('[Subroutine - cic_interpolate] nc incorrect')

    if (minval(xp).lt.0) call interpolation_error_stop('[Subroutine - cic_interpolate] minval(xp)<0')
    if (maxval(xp).gt.nc) call interpolation_error_stop('[Subroutine - cic_interpolate] maxval(xp)>nc')

    allocate(bg(0:nc+1,0:nc+1,0:nc+1),bgs(0:nc+1,0:nc+1))

    !Interpolate particles to bg array
    bg=0.
    do n=1,np

       x=xp(1:3,n)-0.5
       i1=1+floor(x)
       i2=1+i1
       dx1=i1-x
       dx2=1.-dx1

       if (present(wp)) then
          w=wp(n)
       else
          w=1.
       end if

       bg(i1(1),i1(2),i1(3))=bg(i1(1),i1(2),i1(3))+dx1(1)*dx1(2)*dx1(3)*w
       bg(i2(1),i1(2),i1(3))=bg(i2(1),i1(2),i1(3))+dx2(1)*dx1(2)*dx1(3)*w

       bg(i1(1),i2(2),i1(3))=bg(i1(1),i2(2),i1(3))+dx1(1)*dx2(2)*dx1(3)*w
       bg(i2(1),i2(2),i1(3))=bg(i2(1),i2(2),i1(3))+dx2(1)*dx2(2)*dx1(3)*w

       bg(i1(1),i1(2),i2(3))=bg(i1(1),i1(2),i2(3))+dx1(1)*dx1(2)*dx2(3)*w
       bg(i2(1),i1(2),i2(3))=bg(i2(1),i1(2),i2(3))+dx2(1)*dx1(2)*dx2(3)*w

       bg(i1(1),i2(2),i2(3))=bg(i1(1),i2(2),i2(3))+dx1(1)*dx2(2)*dx2(3)*w
       bg(i2(1),i2(2),i2(3))=bg(i2(1),i2(2),i2(3))+dx2(1)*dx2(2)*dx2(3)*w

    end do

    !Buffer bg array
    !!cart_neighbor(1)=-z
    !!cart_neighbor(2)=+z
    !!cart_neighbor(3)=-y
    !!cart_neighbor(4)=+y
    !!cart_neighbor(5)=-x
    !!cart_neighbor(6)=+x

    bs=size(bgs)
    tag=29

    !-x
    bgs=bg(0,:,:)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(5),tag,cart_neighbor(6),tag,mpi_comm_cart,status,ierr)
    bg(nc,:,:)=bg(nc,:,:)+bgs

    !+x
    bgs=bg(nc+1,:,:)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(6),tag,cart_neighbor(5),tag,mpi_comm_cart,status,ierr)
    bg(1,:,:)=bg(1,:,:)+bgs

    !-y
    bgs=bg(:,0,:)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(3),tag,cart_neighbor(4),tag,mpi_comm_cart,status,ierr)
    bg(:,nc,:)=bg(:,nc,:)+bgs

    !+y
    bgs=bg(:,nc+1,:)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(4),tag,cart_neighbor(3),tag,mpi_comm_cart,status,ierr)
    bg(:,1,:)=bg(:,1,:)+bgs

    !-z
    bgs=bg(:,:,0)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(1),tag,cart_neighbor(2),tag,mpi_comm_cart,status,ierr)
    bg(:,:,nc)=bg(:,:,nc)+bgs
    
    !+z
    bgs=bg(:,:,nc+1)
    call mpi_sendrecv_replace(bgs,bs,mpi_real8,cart_neighbor(2),tag,cart_neighbor(1),tag,mpi_comm_cart,status,ierr)
    bg(:,:,1)=bg(:,:,1)+bgs

    !Normalize
    ncl=(nc*1.d0)**3
    call mpi_allreduce(ncl,ncg,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

    npl=np*1.d0
    call mpi_allreduce(npl,npg,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

    sgl=sum(bg(1:nc,1:nc,1:nc))
    call mpi_allreduce(sgl,sgg,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

    if (interpolation_verbosity .and. rank.eq.0) then
       write(*,*) 'number particles: ',npg
       write(*,*) 'grid sum: ',sgg
       write(*,*) 'number of cells: ',ncg
       write(*,*) 'grid sum/number of cells: ',sgg/ncg
    end if

    bg=bg*ncg/sgg
    
    sgl=sum(bg(1:nc,1:nc,1:nc))
    call mpi_allreduce(sgl,sgg,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

    if (interpolation_verbosity .and. rank.eq.0) then
       write(*,*) 'Normalized'
       write(*,*) 'grid sum: ',sgg
       write(*,*) 'grid sum/number of cells: ',sgg/ncg
    end if

    gn=real(bg(1:nc,1:nc,1:nc),kind=4)

    deallocate(bg,bgs)

  end subroutine cic_interpolate

  subroutine interpolation_error_stop(expl)
    implicit none
    character(len=*), intent(in) :: expl
    write(*,*) '[MOD - Interpolation] rank=',rank
    write(*,*) '-->'//expl
    call mpi_abort(mpi_comm_world, ierr, ierr)
  end subroutine interpolation_error_stop

    subroutine buffer(bg)
    implicit none
    real, dimension(0:Ncells+1,0:Ncells+1,0:Ncells+1), intent(inout) :: bg
    real, dimension(0:Ncells+1,0:Ncells+1) :: bgs
    integer, parameter :: bs=size(bgs),tag=29
    integer :: i,status(MPI_STATUS_SIZE)
    
    do i=1,2
       !-x
       bgs=bg(1,:,:)
       call mpi_sendrecv_replace(bgs,bs,mpi_real,cart_neighbor(5),tag,cart_neighbor(6),tag,mpi_comm_cart,status,ierr)
       bg(Ncells+1,:,:)=bgs

       !+x
       bgs=bg(Ncells,:,:)
       call mpi_sendrecv_replace(bgs,bs,mpi_real,cart_neighbor(6),tag,cart_neighbor(5),tag,mpi_comm_cart,status,ierr)
       bg(0,:,:)=bgs

       !-y
       bgs=bg(:,1,:)
       call mpi_sendrecv_replace(bgs,bs,mpi_real,cart_neighbor(3),tag,cart_neighbor(4),tag,mpi_comm_cart,status,ierr)
       bg(:,Ncells+1,:)=bgs

       !+y
       bgs=bg(:,Ncells,:)
       call mpi_sendrecv_replace(bgs,bs,mpi_real,cart_neighbor(4),tag,cart_neighbor(3),tag,mpi_comm_cart,status,ierr)
       bg(:,0,:)=bgs

       !-z
       bgs=bg(:,:,1)
       call mpi_sendrecv_replace(bgs,bs,mpi_real,cart_neighbor(1),tag,cart_neighbor(2),tag,mpi_comm_cart,status,ierr)
       bg(:,:,Ncells+1)=bgs

       !+z
       bgs=bg(:,:,Ncells+1)
       call mpi_sendrecv_replace(bgs,bs,mpi_real,cart_neighbor(2),tag,cart_neighbor(1),tag,mpi_comm_cart,status,ierr)
       bg(:,:,0)=bgs
    
    end do

  end subroutine buffer
  
end module Interpolation
