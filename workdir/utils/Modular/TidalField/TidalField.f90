#define NGP
#define DECONVOLVE
!#define PAIRTEST
program TidalField
  use Parameters
  use Variables
  use mMPI
  use Interpolation
  use Pencil
  use Deriv3D
  use Class_LinkedList
  use Particles
  implicit none

  !File IO
  !!Checkpoint
  character(len=*), parameter :: fn_cp = '/scratch/dbi208/Tidal_Field/tide3/workdir/input/checkpoints_dmbh'
  integer, parameter :: n_cp_tot = 99 !wc -l fn_cp
  real, dimension(n_cp_tot) :: cp
  integer :: n_cp
  character(len=100) :: cp_s
  !!Output directory
  character(len=*), parameter :: odir = './'

  !Simulation parameters
  integer, parameter :: n_cl = (nodes_dim*Ncells)**3
  integer, parameter :: n_dm = n_cl/4
  integer, parameter :: n_bh = 100
  real, parameter :: f_bh = n_bh*10.**(-5.)
  real, parameter :: Oc=0.26,Ob=0.05,O0=Oc+Ob,fc=Oc/O0,fb=Ob/O0
  real, parameter :: m_dm = (1.-f_bh)*fc*n_cl/n_dm, m_bh=f_bh*fc*n_cl/n_bh

  !Particle arrays
  !!Particle header
  real, dimension(11) :: header
  !!DM particle array
  integer, parameter :: max_np_dm = nint(1.5*n_dm)
  real, dimension(6,max_np_dm) :: xvp_dm
  integer :: np_dm
  !!BH particle array
  integer, parameter :: max_np_bh = n_bh
  real, dimension(6,max_np_bh) :: xvp_bh
  integer(8), dimension(max_np_bh) :: pid_bh,pid_recv
  integer :: np_bh,np_recv

  !Fields
  !!Density fields
  real, dimension(Ncells,Ncells,Ncells) :: delta_dm,delta_bh
  real, dimension(0:Ncells+1,0:Ncells+1,0:Ncells+1) :: delta_m
  !!PM-Tidal field
  real, dimension(0:Ncells+1,0:Ncells+1,0:Ncells+1,3,3) :: tidal_tensor
  real, dimension(3,3,max_np_bh) :: tide_x_pm,tide_recv
  real, dimension(3,3,max_np_bh,n_cp_tot) :: tide_x_pm_z
  real :: trace

  !PP-Tidal Field 
  real, dimension(3,3,max_np_bh) :: tide_x_pp
  real, dimension(3,3,max_np_bh,n_cp_tot) :: tide_x_pp_z

  !Output variables
  real, dimension(max_np_bh) :: pm,pp,p3m

  !!Min/max radius for pp search
  !real, parameter :: eps_pp=0.1, r_ewald=2., buf_pp=6.*r_ewald
  real, parameter :: eps_pp=0.1, r_ewald=64., buf_pp=2.!6.*r_ewald

  !Control variables
  logical, parameter :: correct_kernel = .true.

  !General use
  integer :: n,nn,i,j,k,stat
  integer, dimension(3) :: ix,i1,i2
  real :: r
  real, dimension(3) :: x,dx,dx1,dx2
  character(len=1000) :: fn
  integer :: tag,status(MPI_STATUS_SIZE) 
  
#ifdef PAIRTEST
  real, dimension(3,3,max_np_bh) :: tide_x_tr
  real :: r_pair
  real, parameter :: xs = 0., ys=0.25
#endif

  !!!!!!!!!!!
  !MAIN CODE!
  !!!!!!!!!!!

  call start_mpi
  call setup_pencil
  call omp_set_num_threads(Nomp)

  if (rank.eq.0) write(*,*) 'Program: TidalField'

  !Read in checkpoints
  if (rank.eq.0) then
     write(*,*) 'Reading file: '//fn_cp
     open(unit=11,file=fn_cp)
     do n=1,n_cp_tot
        read(11,*) cp(n)
     end do
     close(11)
     write(*,*) 'Max/Min redshift: ',cp(1),cp(n_cp_tot)
  end if
  call mpi_bcast(cp,size(cp),mpi_real,0,mpi_comm_world,ierr)

#ifdef PAIRTEST
  if (rank.eq.0) then
     write(*,*) 'Performing PAIR TEST instead of TIDAL FIELD'
     open(unit=71,file='pair_pm.txt',status='replace',recl=1000)
     open(unit=72,file='pair_pp.txt',status='replace',recl=1000)
     open(unit=73,file='pair_tr.txt',status='replace',recl=1000)
  end if
#endif

  !Loop over checkpoints
  do n_cp=1,n_cp_tot

     write(cp_s,'(F10.3)') cp(n_cp)
     if (rank.eq.0) write(*,*) 'Redshift: '//trim(adjustl(cp_s))

     !Read in DM particles
     xvp_dm=0.
     fn=dir//'node'//trim(adjustl(rank_s))//'/'//trim(adjustl(cp_s))//'xv'//trim(adjustl(rank_s))//'.dat'
     if (rank.eq.0) write(*,*) 'Reading DM file: '//trim(adjustl(fn))
     open(unit=11,file=trim(adjustl(fn)),status='old',access='stream',iostat=stat)
     read(11) np_dm
     if (np_dm.gt.max_np_dm) write(*,*) 'ERROR: np_dm>max_np_dm',np_dm,max_np_dm,rank
     read(11) header
     do n=1,np_dm
        read(11) xvp_dm(:,n)
     end do
     close(11)

#ifdef PAIRTEST
     np_dm=0
     xvp_dm=0.
#endif

     !Read in BH particles
     xvp_bh=0.
     fn=dir//'node'//trim(adjustl(rank_s))//'/'//trim(adjustl(cp_s))//'xv'//trim(adjustl(rank_s))//'_nu.dat'
     if (rank.eq.0) write(*,*) 'Reading BH file: '//trim(adjustl(fn))
     open(unit=11,file=trim(adjustl(fn)),status='old',access='stream',iostat=stat)
     read(11) np_bh
     if (np_bh.gt.max_np_bh) write(*,*) 'ERROR: np_bh>max_np_bh',np_bh,max_np_bh,rank
     read(11) header
     do n=1,np_bh
        read(11) xvp_bh(:,n)
     end do
     close(11)

     !Read in BH PID
     pid_bh=0
     fn=dir//'node'//trim(adjustl(rank_s))//'/'//trim(adjustl(cp_s))//'pid'//trim(adjustl(rank_s))//'_nu.dat'
     if (rank.eq.0) write(*,*) 'Reading BH pids: '//trim(adjustl(fn))
     open(unit=11,file=trim(adjustl(fn)),status='old',access='stream',iostat=stat)
     read(11) n
     if (n.ne.np_bh) write(*,*) 'ERROR: np_bh_pid != np_bh',n,np_bh,rank
     do n=1,np_bh
        read(11) pid_bh(n)
     end do
     close(11)

#ifdef PAIRTEST
     xvp_bh=0.
     if (rank.eq.0) then
        np_bh=2
        xvp_bh(1,1)=xs+(n_cp-1.)/(n_cp_tot-1.)*(Ncells/2.-xs)
        xvp_bh(2,1)=ys+Ncells/2.
        xvp_bh(3,1)=Ncells/2.
        xvp_bh(1:3,2)=Ncells/2.

        r_pair=sqrt(sum( (xvp_bh(1:3,1)-xvp_bh(1:3,2))**2. ))
        write(*,*) '>r_pair,dx',r_pair,xvp_bh(1,1)-xvp_bh(1,2)
     else
        np_bh=0
     end if
#endif
     
     !Interpolate particles to grids
     delta_dm=0;delta_bh=0;delta_m=0

#ifndef PAIRTEST
     if (rank.eq.0) write(*,*) 'Interpolating DM particles to grid'
#ifdef NGP
     call ngp_interpolate(xvp_dm(1:3,:np_dm),delta_dm)
#else
     call cic_interpolate(xvp_dm(1:3,:np_dm),delta_dm)
#endif
     delta_dm=delta_dm-1.
#endif

     if (rank.eq.0) write(*,*) 'Interpolating BH particles to grid'
#ifdef NGP
     call ngp_interpolate(xvp_bh(1:3,:np_bh),delta_bh)
#else
     call cic_interpolate(xvp_bh(1:3,:np_bh),delta_bh)
#endif
     delta_bh=delta_bh-1.

#ifdef PAIRTEST
     delta_bh=delta_bh/(n_bh/2.)
#endif

     delta_m(1:Ncells,1:Ncells,1:Ncells)=(1.-f_bh)*delta_dm+(f_bh)*delta_bh
     delta_m=delta_m*fc !Not including baryons
     if (correct_kernel) call kernel_l(delta_m(1:Ncells,1:Ncells,1:Ncells))
     call buffer(delta_m)

     !Compute PM-tidal tensor
     if (rank.eq.0) write(*,*) 'Computing PM tidal tensor'
     tidal_tensor=0.
     do j=1,3
        do i=1,3
           if (rank.eq.0) write(*,*) '>>i,j=',i,j
           call gradIgradJinvL(delta_m(1:Ncells,1:Ncells,1:Ncells),tidal_tensor(1:Ncells,1:Ncells,1:Ncells,i,j),i,j)
           call buffer(tidal_tensor(:,:,:,i,j))
        end do
     end do
     
     !Interpolate quantities to PM coordinates
     tide_x_pm=0.
     do n=1,np_bh
        x=xvp_bh(1:3,n)-0.5
        i1=1+floor(x)
        i2=1+i1
#ifdef NGP
        dx1=0.
#else
        dx1=i1-x
#endif
        dx2=1.-dx1
        
        if (minval(i1).lt.0) write(*,*) 'ERROR: minval(i1)<0',minval(i1),x
        if (maxval(i2).gt.Ncells+1) write(*,*) 'ERROR: maxval(i2)>Ncells+1',maxval(i2),x

        do j=1,3
           do i=1,3
           
              tide_x_pm(i,j,n) = tidal_tensor(i1(1),i1(2),i1(3),i,j)*dx1(1)*dx1(2)*dx1(3)&
                   &+ tidal_tensor(i2(1),i1(2),i1(3),i,j)*dx2(1)*dx1(2)*dx1(3)&
                   &+ tidal_tensor(i1(1),i2(2),i1(3),i,j)*dx1(1)*dx2(2)*dx1(3)&
                   &+ tidal_tensor(i2(1),i2(2),i1(3),i,j)*dx2(1)*dx2(2)*dx1(3)&
                   &+ tidal_tensor(i1(1),i1(2),i2(3),i,j)*dx1(1)*dx1(2)*dx2(3)&
                   &+ tidal_tensor(i2(1),i1(2),i2(3),i,j)*dx2(1)*dx1(2)*dx2(3)&
                   &+ tidal_tensor(i1(1),i2(2),i2(3),i,j)*dx1(1)*dx2(2)*dx2(3)&
                   &+ tidal_tensor(i2(1),i2(2),i2(3),i,j)*dx2(1)*dx2(2)*dx2(3)

           end do
        end do

        trace=tide_x_pm(1,1,n)+tide_x_pm(2,2,n)+tide_x_pm(3,3,n)
        tide_x_pm(1,1,n)=tide_x_pm(1,1,n)-trace/3.
        tide_x_pm(2,2,n)=tide_x_pm(2,2,n)-trace/3.
        tide_x_pm(3,3,n)=tide_x_pm(3,3,n)-trace/3.
        
     end do

     !Compute PP-tidal force
     if (rank.eq.0) write(*,*) 'Passing particles'
     call pass_particles(xvp_dm,np_dm,0.,real(Ncells),real(buf_pp))
     call pass_particles(xvp_bh,np_bh,0.,real(Ncells),real(buf_pp))

     !Loop over bh
     tide_x_pp=0.
#ifdef PAIRTEST
     tide_x_tr=0.
#endif
     do n=1,np_bh
        !Make sure on this node
        if ( minval(xvp_bh(1:3,n)).lt.0 .or. maxval(xvp_bh(1:3,n)).gt.Ncells ) cycle

        !Loop over DM particles
        do nn=1,np_dm
           dx=xvp_bh(1:3,n)-xvp_dm(1:3,nn)
           ix=floor(xvp_bh(1:3,n))-floor(xvp_dm(1:3,nn))
           r=sqrt(sum(dx**2))
           !if (r.gt.eps_pp .and. all(abs(ix).le.buf_pp)) then
           if (r.gt.eps_pp .and. r.le.buf_pp) then
              dx=dx/r
              do j=1,3
                 do i=1,3
                    tide_x_pp(i,j,n)=tide_x_pp(i,j,n)-(1./4./pi)*m_dm*(1./r**3)*kernel_s(r,dx,i,j)
                 end do
              end do
           end if
        end do

        !Loop over other BH particles
        do nn=1,np_bh
           if (n.eq.nn) cycle
           dx=xvp_bh(1:3,n)-xvp_bh(1:3,nn)
           ix=floor(xvp_bh(1:3,n))-floor(xvp_bh(1:3,nn))
           r=sqrt(sum(dx**2))
           dx=dx/r
           !if (r.gt.eps_pp .and. all(abs(ix).le.buf_pp)) then
           if (r.gt.eps_pp .and. r.le.buf_pp) then
              do j=1,3
                 do i=1,3
                    tide_x_pp(i,j,n)=tide_x_pp(i,j,n)-(1./4./pi)*m_bh*(1./r**3)*kernel_s(r,dx,i,j)
                 end do
              end do
           end if

#ifdef PAIRTEST
           if ( minval(xvp_bh(1:3,nn)).lt.0 .or. maxval(xvp_bh(1:3,nn)).gt.Ncells ) cycle
           do j=1,3
              do i=1,3
                 tide_x_tr(i,j,n)=tide_x_tr(i,j,n)-(1./4./pi)*m_bh*(1./r**3)*(3.*dx(i)*dx(j)-dij(i,j))
              end do
           end do
#endif

        end do

     end do

     call delete_particles(xvp_dm,np_dm,0.,real(Ncells))
     call delete_particles(xvp_bh,np_bh,0.,real(Ncells))

#ifdef PAIRTEST
     if (rank.eq.0) then
        write(71,*) r_pair,tide_x_pm(:,:,1)
        write(72,*) r_pair,tide_x_pp(:,:,1)
        write(73,*) r_pair,tide_x_tr(:,:,1)
     end if
#else
     !Snychronize with PID
     if (rank.eq.0) then
        write(*,*) 'Synchronizing...'
        tide_x_pm_z(:,:,:,n_cp)=0.
        tide_x_pp_z(:,:,:,n_cp)=0.
        do n=1,np_bh
           tide_x_pm_z(:,:,pid_bh(n),n_cp)=tide_x_pm(:,:,n)
           tide_x_pp_z(:,:,pid_bh(n),n_cp)=tide_x_pp(:,:,n)
        end do
     end if

     do n=1,nodes_dim**3-1
        tag=n
        if (rank.eq.n) then
           call mpi_send(np_bh,1,mpi_integer,0,tag,mpi_comm_world,ierr)
           call mpi_send(pid_bh,size(pid_bh),mpi_integer8,0,tag,mpi_comm_world,ierr)
           call mpi_send(tide_x_pm,size(tide_x_pm),mpi_real,0,tag,mpi_comm_world,ierr)
           call mpi_send(tide_x_pp,size(tide_x_pp),mpi_real,0,tag,mpi_comm_world,ierr)
        else if (rank.eq.0) then

           write(*,*) '>>n',n

           call mpi_recv(np_recv,1,mpi_integer,n,tag,mpi_comm_world,status,ierr)
           call mpi_recv(pid_recv,size(pid_recv),mpi_integer8,n,tag,mpi_comm_world,status,ierr)
           
           !First pm
           call mpi_recv(tide_recv,size(tide_recv),mpi_real,n,tag,mpi_comm_world,status,ierr)
           do nn=1,np_recv
              tide_x_pm_z(:,:,pid_recv(nn),n_cp)=tide_recv(:,:,nn)
           end do

           call mpi_recv(tide_recv,size(tide_recv),mpi_real,n,tag,mpi_comm_world,status,ierr)
           do nn=1,np_recv
              tide_x_pp_z(:,:,pid_recv(nn),n_cp)=tide_recv(:,:,nn)
           end do

        end if

        call mpi_barrier(mpi_comm_world,ierr)

     end do
     call mpi_bcast(tide_x_pm_z(:,:,:,n_cp),size(tide_x_pm_z(:,:,:,n_cp)),mpi_real,0,mpi_comm_world,ierr)
     call mpi_bcast(tide_x_pp_z(:,:,:,n_cp),size(tide_x_pp_z(:,:,:,n_cp)),mpi_real,0,mpi_comm_world,ierr)
#endif

  end do

#ifdef PAIRTEST
  close(71)
  close(72)
  close(73)
#else
  !Output some stuff
  if (rank.eq.0) then
     write(*,*) 'Writing to file'

     open(unit=11,file=odir//'tide_pm.bin',status='replace',access='stream')
     write(11) tide_x_pm_z
     close(11)

     open(unit=11,file=odir//'tide_pp.bin',status='replace',access='stream')
     write(11) tide_x_pp_z
     close(11)

     open(unit=11,file=odir//'tide_p3m.bin',status='replace',access='stream')
     write(11) tide_x_pm_z+tide_x_pp_z
     close(11)

     open(unit=11,file=odir//'tide_z.bin',status='replace',access='stream')
     write(11) cp
     close(11) 

     open(unit=12,file=odir//'tide_pm.txt',recl=10000,status='replace')
     open(unit=13,file=odir//'tide_pp.txt',recl=10000,status='replace')
     open(unit=14,file=odir//'tide_p3m.txt',recl=10000,status='replace')
     do n_cp=1,n_cp_tot

        do n=1,n_bh
           pm(n)=sum(tide_x_pm_z(:,:,n,n_cp)**2.)
           pp(n)=sum(tide_x_pp_z(:,:,n,n_cp)**2.)
           p3m(n)=sum((tide_x_pm_z(:,:,n,n_cp)+tide_x_pp_z(:,:,n,n_cp))**2.)
        end do

        write(12,*) cp(n_cp),pm
        write(13,*) cp(n_cp),pp
        write(14,*) cp(n_cp),p3m

     end do

     close(12)
     close(13)
     close(14)

  end if
#endif

  if (rank.eq.0) write(*,*) 'Finished TidalField'

  call end_mpi

contains

  subroutine kernel_l(grid)
    implicit none

    real, dimension(Ncells,Ncells,Ncells), intent(inout) ::grid

    integer :: n,k,j,i,kg,jg,ig,mg,l,ind,dx,dxy
    real :: kz,ky,kx,kr,kphys,sinc,sincx,sincy,sincz

    real, parameter :: ncr = 1.0*nc
    integer, parameter :: hc = nc/2

    call reset_pencil
    cube=grid
    call cp_fftw(1)

    ind = 0
    dx = fsize(1)
    dxy = dx * fsize(2)

    do k = 1, nc_pen+mypadd
       ind = (k-1)*nc_node_dim*nc/2
       do j = 1, nc_node_dim
          do i = 1, nc, 2
             kg = ind / dxy
             mg = ind - kg * dxy
             jg = mg / dx
             ig = mg - jg * dx
             kg = fstart(3) + kg
             jg = fstart(2) + jg
             ig = 2 * (fstart(1) + ig) - 1
             ind = ind + 1
             if (kg < hc+2) then
                kz=kg-1
             else
                kz=kg-1-nc
             endif
             if (jg < hc+2) then
                ky=jg-1
             else
                ky=jg-1-nc
             endif
             kx = (ig-1)/2

             kr=sqrt(kx**2+ky**2+kz**2)

             slab(i:i+1,j,k)=slab(i:i+1,j,k)*exp(-(2.*pi*kr*r_ewald/ncr)**2)

#ifdef DECONVOLVE
             sincx = merge(sin(pi*kx/ncr)/(pi*kx/ncr),1.0,kx/=0.0)
             sincy = merge(sin(pi*ky/ncr)/(pi*ky/ncr),1.0,ky/=0.0)
             sincz = merge(sin(pi*kz/ncr)/(pi*kz/ncr),1.0,kz/=0.0)
             sinc = (sincx*sincy*sincz)
#ifdef NGP
             slab(i:i+1,j,k)=slab(i:i+1,j,k)/sinc
#else
             slab(i:i+1,j,k)=slab(i:i+1,j,k)/sinc**2
#endif
#endif

          end do
       end do
    end do

    call cp_fftw(-1)
    grid=cube

    return
  end subroutine kernel_l 

  function dij(i,j) result(d)
    implicit none
    integer, intent(in) :: i,j
    integer :: d
    d=merge(1,0,i.eq.j)
  end function dij

  function kernel_s(r,dx,i,j) result(k)
    implicit none
    real, intent(in) :: r
    real, dimension(3), intent(in) :: dx
    integer, intent(in) :: i,j
    real :: k,l

    k=3.*dx(i)*dx(j)-dij(i,j)  
    if (correct_kernel) then
       l=r/r_ewald/2.
       k=k*(erfc(l)+(2.*l/sqrt(pi))*exp(-l**2.))+(4.*l**3/sqrt(pi))*dx(i)*dx(j)*exp(-l**2.)
    end if
  end function kernel_s

end program TidalField
