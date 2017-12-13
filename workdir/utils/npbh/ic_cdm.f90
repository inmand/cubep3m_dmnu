program dist_init 
  use omp_lib
  implicit none

  include 'mpif.h'
#  include "../../parameters"

  integer, parameter  :: num_threads = cores*nested_threads 
  logical, parameter :: generate_seeds=.false.
  logical, parameter :: correct_kernel=.true.

  real, parameter :: ns = n_s
  real, parameter :: As = A_s
  real, parameter :: ko = k_o
  real, parameter :: s8 = sigma_8
  real, parameter :: omegal=omega_l 
  real, parameter :: omegam=1.0-omegal 

  real, parameter :: redshift=z_i
  real, parameter :: scalefactor=1/(1+redshift)

  !! NEUTRINO sims may use this:
  real(4), parameter :: Vphys2sim = 1.0/(300. * sqrt(omega_m) * box * (1. + redshift) / 2. / nc)

  !! Pi
  real, parameter :: pi=3.141592654

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
  integer, parameter :: np=hc 
  real, parameter    :: npr=np

  real, parameter :: knyq=pi*npr/box

  !! internal parallelization parameters
  integer(4), parameter :: nc_node_dim = nc/nodes_dim
  integer(4), parameter :: np_node_dim = np/nodes_dim
  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim

  !! parallelization variables
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local,wc_counter, count_i,count_f,count_r
  logical :: firstfftw
  integer(8) :: plan, iplan

  !! Power spectrum arrays
  real, dimension(3,nk) :: tf    !CLASS
  real, dimension(2,nc) :: pkm,pkn

  !! For pencils decomposition:
  integer(4), parameter   :: nodes_pen = nodes_dim
  integer(4), parameter   :: nc_pen = nc_node_dim / nodes_dim
  integer(4), parameter   :: dim_y = nodes_dim
  integer(4), parameter   :: dim_z = nodes_dim**2
  integer(4) :: pen_dims(2), istart(3), iend(3), isize(3), fstart(3), fend(3), fsize(3), mypadd
  integer(4), dimension(0:nodes_dim-1) :: pen_neighbor_to
  integer(4), dimension(0:nodes_dim-1) :: pen_neighbor_fm

  !! Fourier transform arrays
  real, dimension(nc_node_dim,nc_node_dim,nc_node_dim) :: cube
  real, dimension(nc, nc_node_dim, nc_pen+2) :: slab, slab_work
  real, dimension(nc_node_dim, nc_node_dim, nc_pen, 0:nodes_pen-1)      :: recv_cube
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1,0:nc_node_dim+1) :: phi
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1) :: phi_buf

  !! Particles arrays for subroutine dm
  real, dimension(6,np_node_dim,np_node_dim*(1+bcc),num_threads) :: xvp

  !! Timing variables
  real(8) :: sec1, sec2

  !! Equivalence arrays to save memory
  equivalence (phi,slab_work,recv_cube)
  equivalence (slab,cube)

  !! Common block
  common /rvar/ tf, pkm, pkn, phi_buf
  common / equiv1 / phi
  common / equiv2 / slab
  common / xvar / xvp

  call mpi_initialize

  sec1 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STARTING INITIAL CONDITIONS: ", sec1

  call omp_set_num_threads(num_threads)
  if (rank == 0) print*, 'Note: compile with openmp: num_threads=', num_threads

  if (rank == 0) call writeparams

  call initvar
  call transferfnc
  call noisemap
  call deltafield
  call potentialfield

  call dm(0)
  if (rank == 0) call writepowerspectra
  call veltransfer
  call dm(1)

  call di_fftw(0)

  sec2 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STOPPING INITIAL CONDITIONS: ", sec2
  if (rank == 0) write(*,*) "ELAPSED TIME: ", sec2-sec1
  call mpi_barrier(mpi_comm_world, ierr)

  call mpi_finalize(ierr)

contains

  subroutine mpi_initialize
    implicit none
    
    integer(4) :: i, j, nodes_returned
    integer(4) :: dims(3), ndim
    logical :: periodic(3), reorder
  
    call mpi_init(ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
    if (nodes_returned /= nodes ) then
      write(*,*) 'dist_init compiled for a different number of nodes'
      write(*,*) 'mpirun nodes=',nodes_returned,'dist_init nodes=',nodes 
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    if (mod(nc,nodes_dim**2) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into pencils'
      write(*,*) 'nc=',nc, 'nodes_dim**2=',nodes_dim**2,'mod(nc,nodes_dim**2)=', mod(nc,nodes_dim**2)
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    if (rank==0) then
      write(*,*) 'dist_init running on',nodes,'nodes'
      write(*,*) 'using cubic distribution:',nodes_dim,'nodes per dimension'
      write(*,*) nc,'cells in mesh'
    endif

    slab_coord(3) = rank / nodes_slab
    slab_rank = rank - slab_coord(3) * nodes_slab
    slab_coord(2) = slab_rank / nodes_dim
    slab_coord(1) = slab_rank - slab_coord(2) * nodes_dim
    do j = 0, nodes_dim - 1
        pen_neighbor_to(j) = nodes_slab*slab_coord(3) + slab_coord(2) + j*nodes_dim
        pen_neighbor_fm(j) = nodes_slab*slab_coord(3) + j + nodes_dim*slab_coord(1)
    enddo

    dims(:) = nodes_dim
    periodic(:) = .true.
    reorder = .false.
    ndim = 3

    call mpi_cart_create(mpi_comm_world, ndim,dims, periodic, &
                       reorder, mpi_comm_cart, ierr)
    call mpi_comm_rank(mpi_comm_cart, cart_rank, ierr)
    call mpi_cart_coords(mpi_comm_cart, cart_rank, ndim,  &
                         cart_coords, ierr)

    do i = 0, ndim-1
      call mpi_cart_shift(mpi_comm_cart, i, 1, cart_neighbor(2*(i+1)-1), &
                          cart_neighbor(2*(i+1)), ierr)
    enddo

  end subroutine mpi_initialize

!-------------------------------------------------------------------!

subroutine pack_pencils
    implicit none

    integer(4) :: i,j,k,i0,i1,k1
    integer(4) :: pen_slice,tag,rtag
    integer(8) :: num_elements_i8
    integer(4) :: num_elements 
    integer(4), dimension(2*nodes_dim) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim) :: wait_status
    integer(4) nc_pen_break, breakup
    real(4) :: passGB

    breakup = 1
    num_elements_i8 = int(nc_node_dim,kind=8) * nc_node_dim * nc_pen 
    passGB = 4. * num_elements_i8 / 1024.**3
    if (passGB > 1.) then
        breakup = 2**ceiling(log(passGB)/log(2.))
    endif
    num_elements = num_elements_i8 / breakup

    do k = 1, breakup
        nc_pen_break = nc_pen/breakup*(k-1)
        do j = 0, nodes_dim - 1
            pen_slice = j
            tag  = rank**2
            rtag = pen_neighbor_fm(j)**2
            call mpi_isend(cube(1,1, pen_slice*nc_pen + nc_pen_break + 1), num_elements, &
                           mpi_real, pen_neighbor_to(j), tag, mpi_comm_world, &
                           requests(pen_slice+1),ierr)
            call mpi_irecv(recv_cube(1,1,1+nc_pen_break,pen_slice), &
                           num_elements, mpi_real, pen_neighbor_fm(j),rtag, &
                           mpi_comm_world, requests(pen_slice+1+nodes_dim), &
                           ierr)
        enddo

        call mpi_waitall(2*nodes_dim, requests, wait_status, ierr)

    enddo

    do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        pen_slice = i

        do k = 1, nc_pen
            do j = 1, nc_node_dim
                slab(i0:i1,j,k) = recv_cube(:,j,k,pen_slice)
            enddo
        enddo

    enddo

end subroutine pack_pencils

!-------------------------------------------------------------------!

subroutine unpack_pencils
    implicit none

    integer(4) :: i,j,k,i0,i1,k1
    integer(4) :: pen_slice,num_elements,tag,rtag
    integer(8) :: num_elements_i8
    integer(4), dimension(2*nodes_dim) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim) :: wait_status
    integer(4) nc_pen_break, breakup
    real(4) :: passGB

    do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        pen_slice = i
        do k = 1, nc_pen
            do j = 1, nc_node_dim
                recv_cube(:, j, k, pen_slice) = slab(i0:i1, j, k)
            enddo
        enddo
    enddo

    breakup = 1
    num_elements_i8 = int(nc_node_dim,kind=8) * nc_node_dim * nc_pen 
    passGB = 4. * num_elements_i8 / 1024.**3
    if (passGB > 1.) then
        breakup = 2**ceiling(log(passGB)/log(2.))
    endif
    num_elements = num_elements_i8 / breakup

    do k = 1, breakup
        nc_pen_break = nc_pen/breakup*(k-1)
        do j = 0, nodes_dim - 1
            pen_slice = j
            tag  = rank**2
            rtag = pen_neighbor_to(j)**2
            call mpi_isend(recv_cube(1,1,1+nc_pen_break,pen_slice), num_elements, &
                           mpi_real, pen_neighbor_fm(j), tag, mpi_comm_world, &
                           requests(pen_slice+1),ierr)
            call mpi_irecv(cube(1,1,pen_slice*nc_pen + nc_pen_break + 1), &
                           num_elements, mpi_real, pen_neighbor_to(j),rtag, &
                           mpi_comm_world, requests(pen_slice+1+nodes_dim), &
                           ierr)
        enddo

        call mpi_waitall(2*nodes_dim,requests, wait_status, ierr)

    enddo

end subroutine unpack_pencils

!-------------------------------------------------------------------!

subroutine di_fftw(command)
    use p3dfft
    implicit none

    integer(4) :: i
    integer(4) :: command

    if (firstfftw) then
        pen_dims = (/dim_y,dim_z/)
        call p3dfft_setup(pen_dims, nc, nc, nc, .true.)
        call p3dfft_get_dims(istart, iend, isize, 1, mypadd)
        call p3dfft_get_dims(fstart, fend, fsize, 2)
        firstfftw=.false.
    endif

    if (command /= 0) then
       if (command > 0) then
          call pack_pencils
          call ftran_r2c(slab, slab, "fft")
       else
          call btran_c2r(slab, slab, "tff")
          slab=slab/(real(nc)*real(nc)*real(nc))
          call unpack_pencils
       endif
    else
       call p3dfft_clean
    endif

end subroutine di_fftw

!-------------------------------------------------------------------!

  subroutine writeparams
    implicit none

    real time1,time2
    call cpu_time(time1)

    write(*,*) 'nodes   ', nodes
    write(*,*) 'nc      ', nc
    write(*,*) 'np      ', np
    write(*,*)
    write(*,*) 'n_s      ',ns
    write(*,*) 'sigma_8  ',s8
    write(*,*) 'omega_m  ',omegam
    write(*,*) 'omega_l  ',omegal
    write(*,*)
    write(*,*) 'box      ',box
    write(*,*) 'redshift ',redshift
    write(*,*)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write params'
    return
  end subroutine writeparams

!-------------------------------------------------------------------!

  subroutine writepowerspectra
    implicit none
    integer      :: k
    real         :: kr
    character*180 :: fn

    real time1,time2
    call cpu_time(time1)

    !! Output power spectrum
    !! 1st column is k
    !! 2nd is matter p(k)
    !! 3rd is standard deviation
    !! 6th is noise p(k)
    !! 7th is noise standard deviation
    fn=output_path//'pk.init'

    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,hc+1
       kr=2*pi*(k-1)/box
       write(11,*) kr,pkm(:,k),pkn(:,k)
    enddo
    close(11)

    !! Output cmbfast power spectrum
    fn=output_path//'pk0.init'

    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,2*nc+1
       kr=2*pi*(k-1)/(2*box)
       write(11,*) kr,power(kr,1,2),power(kr,1,3)
    enddo
    close(11)
    
    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write power spectra'
    return
  end subroutine writepowerspectra

!!------------------------------------------------------------------!!

  subroutine transferfnc
    implicit none
    integer :: i,k
    real    :: kr,kmax

    real time1,time2
    call cpu_time(time1)

    if (rank == 0) then
      write(*,*) 'Reading ',fntf
      open(11,file=fntf)
      do k=1,nk
         read(11,*) tf(1,k),tf(2,k),tf(3,k)
      end do
      close(11)

      !Compute power spectrum
      tf(2,:) = As*(tf(1,:)/ko)**(ns-1.)*tf(2,:)**2 !Delta**2
      
    endif

    call mpi_barrier(mpi_comm_world,ierr)
    call mpi_bcast(tf,size(tf),mpi_real,0,mpi_comm_world,ierr) 

    !! tf(1,i) stores k
    !! tf(2,i) stores \Delta^2_c

    call cpu_time(time2)
    time2=time2-time1
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called transfer fnc'
    return
  end subroutine transferfnc

!------------------------------------------------------------------!

  subroutine noisemap
    implicit none
    character(len=6) :: rank_s
    integer      :: i,j,k,seedsize,clock
    integer, dimension(8) :: values

    real         :: x,x1,x2
    character*180 :: fn
    integer(4),allocatable,dimension(:) :: iseed
    integer(4),allocatable,dimension(:) :: iseed_all

    real time1,time2
    call cpu_time(time1)

    write(rank_s,'(i6)') rank
    rank_s=adjustl(rank_s)

    call random_seed() 
#ifdef IBM
    call random_seed(generator=2)
#endif
    call random_seed(size=seedsize)

    allocate(iseed(seedsize))
    allocate(iseed_all(seedsize*nodes))
 
    call system_clock(count=clock)
    iseed = clock + 37 * (/ (i - 1, i = 1, 2) /)
 
    call date_and_time(values=values)
    if(rank==0) write(*,*) values(7:8), iseed

    call random_seed(put=values(7:8)) 
    call random_seed(put=iseed(1:seedsize))

    if (generate_seeds) then
       if (rank==0) then
         write(*,*) 'Generating seeds'
         do j=1,nodes
           do k=1,seedsize
             call random_number(x)
             iseed_all((j-1)*seedsize+k)=int(x*huge(0))
           enddo
         enddo
       endif
       call mpi_scatter(iseed_all,seedsize,mpi_integer,iseed,seedsize,mpi_integer,0,mpi_comm_world,ierr)
    else
       fn=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'seed'//rank_s(1:len_trim(rank_s))//'.init'
       if (rank == 0) print *, 'rank',rank,'Reading ',fn(1:len_trim(fn))
       open(11,file=fn)
       do k=1,seedsize
          read(11,*) i,iseed(k)
       enddo
       close(11)
    endif

    !! Generate random reals between 0 and 1
    call random_seed(put=iseed(1:seedsize))
    call random_number(cube)

    fn=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'seed'//rank_s(1:len_trim(rank_s))//'.init'
    if (rank == 0) print *, 'rank',rank,'Writing ',fn(1:len_trim(fn))
    open(11,file=fn)
    do k=1,seedsize
       write(11,*) k,iseed(k)
    enddo
    close(11)
   
    deallocate(iseed)
    deallocate(iseed_all)

    if (rank==0) then
       write(*,*) 'starting generation of Gaussian random numbers'    
    end if

    !! Generate random Gaussian numbers
    !$omp parallel do default(shared) private(i,j,k,x1,x2)
    do k=1,nc_node_dim
       do j=1,nc_node_dim
          do i=1,nc_node_dim,2
             x1=2*pi*cube(i,j,k)
             if (cube(i+1,j,k)<=0.0) cube(i+1,j,k)=0.0001
             x2=sqrt(-2*log(cube(i+1,j,k)))
             cube(i,j,k)=x2*cos(x1)
             cube(i+1,j,k)=x2*sin(x1)
          enddo
       enddo
    enddo
    !$omp end parallel do

    if (rank==0) then
       write(*,*) 'finished generation of Gaussian random numbers'
    end if

    if (rank==0) then
       write(*,*) 'starting forward FFT'
    end if

    !! Forward FFT white noise field
    call di_fftw(1)
    !! noise is now in slab array

    if (rank==0) then
       write(*,*) 'finished forward FFT'
    end if

    !! Generate noise spectrum
    call powerspectrum(pkn)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called noisemap'
    return
  end subroutine noisemap

!-----------------------------------------------------------------------!

  subroutine deltafield
    implicit none
    integer :: i,j,k,kg,mg,jg,ig,fstat
    real    :: kr,kx,ky,kz
    real    :: powb,powm
    real    :: d,dmin,dmax,dmint,dmaxt
    real*8  :: dsum,dvar,dsumt,dvart
    character(len=4) :: rank_string
    character(len=100) :: check_name

    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0

    !! Interpolate \Delta^2 to grid
    !! Determine Fourier modes \delta(k)
    do k = 1, nc_pen+mypadd
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
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if (kr .eq. 0) then
                   slab(i:i+1,j,k)=0.
!!$                else if ( 2.*pi*kr/box .gt. knyq) then
!!$                   !Don't turn on super-nyquist modes
!!$                   slab(i:i+1,j,k)=0.
                else
                    powm=power(2*pi*kr/box,1,2)/(4*pi*kr**3)
                    slab(i:i+1,j,k)=sqrt(powm*ncr**3)*slab(i:i+1,j,k)
                endif
            enddo
        enddo
    enddo

    !! Calculate matter delta power spectrum

    call powerspectrum(pkm)

    !! Calculate matter delta field statistics

    call di_fftw(-1)

    dmin=0
    dmax=0
    dsum=0
    dvar=0

    !$omp parallel do default(shared) private(i,j,k,d) &
    !$omp& reduction(min:dmin) &
    !$omp& reduction(max:dmax) &
    !$omp& reduction(+:dsum,dvar)
    do k=1,nc_node_dim
       do j=1,nc_node_dim
          do i=1,nc_node_dim
             d=cube(i,j,k)
             dsum=dsum+d
             dvar=dvar+d*d
             dmin=min(dmin,d)
             dmax=max(dmax,d)
          enddo
       enddo
    enddo
    !$omp end parallel do

    call mpi_reduce(dsum,dsumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dvar,dvart,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dmin,dmint,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
    call mpi_reduce(dmax,dmaxt,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

    if (rank==0) then
      dsum=dsumt/real(nc)**3
      dvar=sqrt(dvart/real(ncr)**3)
      write(*,*)
      write(*,*) 'Delta min    ',dmint
      write(*,*) 'Delta max    ',dmaxt
      write(*,*) 'Delta sum ',real(dsum)
      write(*,*) 'Delta var ',real(dvar)
      write(*,*)
    endif

#ifdef write_den 
    if (rank == 0) then
        print *,'Writing density contrast to file'
    endif

    write(rank_string,'(i4)') rank
    rank_string = adjustl(rank_string)

    check_name = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'initdeltafield'// &
        rank_string(1:len_trim(rank_string))//'.bin'

    open(unit=21, file=check_name, status='replace', iostat=fstat, access='stream')

    if (fstat /= 0) then
        write(*,*) 'error opening density file'
        write(*,*) 'rank', rank, 'file:', check_name
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    do k = 1, nc_node_dim
        do j = 1, nc_node_dim

            write(21) cube(:, j, k)

        enddo
    enddo

#endif

    call di_fftw(1)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called delta field'
    return
  end subroutine deltafield

!!------------------------------------------------------------------!!

  subroutine potentialfield
    implicit none

    integer :: i,j,k,ioerr
    integer :: im,ip,ig,jm,jp,jg,km,kp,kg,mg
    real    :: r,x,y,z
    real    :: kr,ksq,kx,ky,kz
    real    :: phi8,phi8tot
    character*180 :: fn
!! for more than 1+1e6 processes this will have to be increased
    character(len=6) :: rank_s

    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

    if (correct_kernel) then
      !! write delta to disk so we can build kernel
      if (rank == 0) write(*,*) 'Caching Delta on disk'
      write(rank_s,'(i6)') rank
      rank_s=adjustl(rank_s)
      fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'delta'//rank_s(1:len_trim(rank_s))

      open(11,file=fn,status='replace',iostat=ioerr, access='stream')

      if (ioerr /= 0) then
        print *,'error opening Delta cache file:',fn
        stop
      endif
      do k=1,nc_pen+2
        write(11) slab(:,:,k)
      enddo
      close(11)

    else
      slab_work=slab
    endif

    !! Construct uncorrected potential kernel in Fourier space

    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0

    do k = 1, nc_pen+mypadd
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
                kz=2*sin(pi*kz/ncr)
                if (jg < hc+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc
                endif
                ky=2*sin(pi*ky/ncr)
                kx = (ig-1)/2
                kx=2*sin(pi*kx/ncr)
                ksq=kx**2+ky**2+kz**2
                ip=i+1
                if (ksq .eq. 0) then
                    slab(i:ip,j,k)=0
                else
                    slab(i,j,k)=-4*pi/ksq
                    slab(ip,j,k)=0
                endif
            enddo
        enddo
    enddo

    if (correct_kernel) then

       !! Inverse FFT potential kernel
       call di_fftw(-1) 

       phi8=0.0

       if (nc_node_dim < 9) then
         print *,'warning: mesh too small in potential kernel correction'
         call mpi_abort(mpi_comm_world,ierr,ierr)
       endif

       if (cart_coords(1) == 0 .and. cart_coords(2) == 0 .and. &
           cart_coords(3) == 0) phi8=cube(9,1,1)+cube(1,9,1)+cube(1,1,9)

       if (cart_coords(3) == nodes_dim-1 .and. cart_coords(2) == 0 .and. &
           cart_coords(1) == 0) phi8=phi8+cube(nc_node_dim-7,1,1)

       if (cart_coords(3) == 0 .and. cart_coords(2) == nodes_dim-1 .and. &
           cart_coords(1) == 0) phi8=phi8+cube(1,nc_node_dim-7,1)

       if (cart_coords(3) == 0 .and. cart_coords(2) == 0 .and. &
           cart_coords(1) == nodes_dim -1) phi8=phi8+cube(1,1,nc_node_dim-7)

       call mpi_reduce(phi8,phi8tot,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       if (rank == 0) phi8=phi8tot/6.0      
       call mpi_bcast(phi8,1,mpi_real,0,mpi_comm_world,ierr)

 
       !! Construct Ewald potential kernel in real space
       !$omp parallel do default(shared) private(i,j,k) &
       !$omp& private(r,x,y,z,ig,jg,kg)
       do k=1,nc_node_dim
          kg=k+nc_node_dim*cart_coords(1)
          if (kg .lt. hc+2) then
             z=kg-1
          else
             z=kg-1-nc
          endif
          do j=1,nc_node_dim
             jg=j+nc_node_dim*cart_coords(2)
             if (jg .lt. hc+2) then
                y=jg-1
             else
                y=jg-1-nc
             endif
             do i=1,nc_node_dim
                ig=i+nc_node_dim*cart_coords(3)
                if (ig .lt. hc+2) then
                   x=ig-1
                else
                   x=ig-1-nc
                endif
                r=sqrt(x**2+y**2+z**2)
                if (r .gt. 8) then
                   cube(i,j,k)=cube(i,j,k)-(phi8+1/8.)
                else
                   if (r .eq. 0) then
                      cube(i,j,k)=-2.5
                   else
                      cube(i,j,k)=-1/r
                   endif
                endif
             enddo
          enddo
       enddo
       !$omp end parallel do

       !! Forward FFT potential kernel
       call di_fftw(1)

      open(11,file=fn,status='old',iostat=ioerr,access='stream')
      if (ioerr /= 0) then
        print *,'error opening Delta cache file:',fn
        stop
      endif
      do k=1,nc_pen+2
        read(11) slab_work(:,:,k)
      enddo
      close(11)

    endif

    !! Complex multiply density field with potential kernel
    !$omp parallel do default(shared) private(i,j,k)
    do k=1,nc_pen+2
       do j=1,nc_node_dim
          do i=1,nc,2
             slab(i:i+1,j,k)=slab_work(i:i+1,j,k)*slab(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    !! Inverse FFT potential field
    call di_fftw(-1)

    !! put cube in phi
    phi=0.0
    phi(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=cube
    call mesh_buffer

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called potential field'
    return
  end subroutine potentialfield

!!------------------------------------------------------------------!!

subroutine veltransfer
    !
    ! Replacing Fourier density field delta(k) with delta(k)*T_velocity(k)/T_density(k) 
    ! where T_velocity(k) is the neutrino velocity transfer function. 
    ! Store results in phi and pass along to dm subroutine with COMMAND == 1.
    !    

    implicit none

    integer :: ii, k, j, i, l, kg, jg, ig, mg, ioerr
    real(4) :: kz, ky, kx, kr, interpVTF,interpMTF, kphys, w1, w2
    integer :: dx, dxy, ind
    real, dimension(2,nk) :: vtf

    character*180 :: fn
    character(len=6) :: rank_s

    real time1,time2
    call cpu_time(time1)

    ! Read density and velocity transfer function
    if (rank == 0) then
        write(*,*) 'Reading ',fntf
        open(11,file=fntf)
        do k=1,nk
           read(11,*) tf(1,k),tf(2,k),tf(3,k)
        end do
        close(11)
        !! Put vtf into simulation units
        tf(3,:) = tf(3,:)*Vphys2sim
    endif

    call mpi_bcast(tf, size(tf), mpi_real, 0, mpi_comm_world, ierr)

    call di_fftw(1)

    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
    
    !$omp parallel default(shared) private(k, j, i, kg, mg, jg, ig, ind, kz, ky, kx, kphys, kr, interpMTF, interpVTF, w1, w2, ii)
    !$omp do schedule(dynamic)
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

                kphys = 2*pi*sqrt(kx**2+ky**2+kz**2)/box
                kx = 2*sin(pi*kx/ncr)
                ky = 2*sin(pi*ky/ncr)
                kz = 2*sin(pi*kz/ncr)

                kr = sqrt(kx**2+ky**2+kz**2)

                !! Interpolate 
                if (kphys <= tf(1,1)) then
                    interpMTF = tf(2,1)
                    interpVTF = tf(3,nk)
                else if (kphys >= tf(1, nk)) then
                    interpMTF = tf(2,nk)
                    interpVTF = tf(3,nk)
                else
                    do ii = 1, nk-1
                        if (kphys <= tf(1,ii+1)) then
                            w1 = tf(1,ii+1) - kphys
                            w2 = kphys - tf(1,ii)
                            interpMTF = (w1*tf(2,ii)+w2*tf(2,ii+1))/(tf(1,ii+1)-tf(1,ii))
                            interpVTF = (w1*tf(3,ii)+w2*tf(3,ii+1))/(tf(1,ii+1)-tf(1,ii))
                           exit
                        endif
                    enddo
                endif

                !! Multiply slab by interpVTF/interpMTF*kr**2/kr
                slab(i:i+1, j, k) = slab(i:i+1, j, k) * interpVTF/interpMTF * kr * (-1.0) !* (-2.0) * pi / ncr
                if (kr .EQ. 0) slab(i:i+1,j,k) = 0.0
            enddo
        enddo
    enddo
    !$omp end do 
    !$omp end parallel 

    call di_fftw(-1)

    phi=0.0
    phi(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=cube
    call mesh_buffer

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called veltransfer'
    return

end subroutine veltransfer

!!------------------------------------------------------------------!!
    !! Dark matter data
    !! xvp(1:3) = x,y,z in grid units from 0 to nc
    !! xvp(4:6) = vx,vy,vz in km/s
  subroutine dm(COMMAND)
    implicit none
    integer :: i,j,k,ioerr,sc
    integer(8) :: pos_in,pos_out 
    integer :: i1,j1,k1,lb,ub
    real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt
    real*8  :: dsum,dvar,dsumt,dvart
    real, dimension(3) :: dis,x
    real*8, dimension(3) :: xav
    character*250 :: fn
    character(len=6) :: rank_s !! for more than 1+1e6 processes this will have to be increased
    character(len=500) :: z_s
    integer(4), dimension(11) :: header_garbage
    integer(4) :: thread
    integer(4), parameter :: xsize1 = 4*3*np_node_dim**2*(bcc+1) 
    integer(4), parameter :: xsize2 = 2*xsize1

    integer :: COMMAND

    real time1,time2
    call cpu_time(time1)

    write(z_s,'(f10.3)') redshift
    z_s=adjustl(z_s)
 
    write(rank_s,'(i6)') rank
    rank_s=adjustl(rank_s)

    !! Open input file
    if (COMMAND == 1) then
       fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'delta'//rank_s(1:len_trim(rank_s))
       open(unit=31,file=fn,status='old',iostat=ioerr,access='stream')
       if (ioerr /= 0) then
          print *,'error opening Delta cache file:',fn
          stop
       endif
    endif

    !! Open output file
    if (COMMAND == 0) then !! Temp file
       if (rank == 0) write(*,*) 'Caching xvp on disk'
       fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'delta'//rank_s(1:len_trim(rank_s))
       open(unit=11,file=fn,status='replace',iostat=ioerr,access='stream')
       if (ioerr /= 0) then
          print *,'error opening xvp cache file:',fn
          stop
       endif
    else !! xv file

       fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'

       open(unit=11,file=fn,status='replace',iostat=ioerr,access='stream')
       if (ioerr /= 0) then
          print *,'error opening:',fn
          stop
       endif
    endif

    !! Write the header in checkpoint format with zeros for all variables (except np_local)
    np_local=np_node_dim**3*(1+bcc)
    header_garbage(:) = 0

    if (COMMAND==1) then

       write(11) np_local,scalefactor,0.0,-3./sqrt(scalefactor),0,1000.,1000.,1000.,1,1,1,1.

    endif

    !! Displace particles
    !! Finite-difference potential to get displacement field
    !$omp parallel default(shared) private(k,k1,j,j1,i,i1,pos_in,pos_out,thread)
    thread = 1
    thread = omp_get_thread_num() + 1
    if (COMMAND == 0) then 
        !! First time we call dm We only need the positions
        !$omp do schedule(dynamic)
        do k=1,np_node_dim
            pos_out = 1 + xsize1*(k-1)
            do j=1,np_node_dim
                do i=1,np_node_dim
                    do sc=0,bcc

                       k1=(nc/np)*(k-1)+1+sc
                       j1=(nc/np)*(j-1)+1+sc
                       i1=(nc/np)*(i-1)+1+sc

                       xvp(1,i,j+np_node_dim*sc,thread)=(phi(i1-1,j1,k1)-phi(i1+1,j1,k1))/2./(4.*pi)+(i1-0.5)
                       xvp(2,i,j+np_node_dim*sc,thread)=(phi(i1,j1-1,k1)-phi(i1,j1+1,k1))/2./(4.*pi)+(j1-0.5)
                       xvp(3,i,j+np_node_dim*sc,thread)=(phi(i1,j1,k1-1)-phi(i1,j1,k1+1))/2./(4.*pi)+(k1-0.5)

                 end do

                enddo
            enddo
            write(unit=11,pos=pos_out) xvp(1:3,:,:,thread) !! save in temp file 
        enddo
        !$omp end do
    else
        !! Second time we call dm we want velocity
        !$omp do schedule(dynamic)
        do k=1,np_node_dim
            pos_out = 1 + sizeof(np_local) + sizeof(header_garbage) + xsize2*(k-1)
            pos_in = 1 + xsize1*(k-1)
            read(unit=31,pos=pos_in) xvp(1:3,:,:,thread) !! catch from temp file
            do j=1,np_node_dim
                do i=1,np_node_dim
                    do sc=0,bcc

                       k1=(nc/np)*(k-1)+1+sc
                       j1=(nc/np)*(j-1)+1+sc
                       i1=(nc/np)*(i-1)+1+sc

                       xvp(4,i,j+np_node_dim*sc,thread)=(phi(i1-1,j1,k1)-phi(i1+1,j1,k1))/2./(4.*pi)
                       xvp(5,i,j+np_node_dim*sc,thread)=(phi(i1,j1-1,k1)-phi(i1,j1+1,k1))/2./(4.*pi)
                       xvp(6,i,j+np_node_dim*sc,thread)=(phi(i1,j1,k1-1)-phi(i1,j1,k1+1))/2./(4.*pi)

                    end do

                enddo
            enddo
            write(unit=11,pos=pos_out) xvp(:,:,:,thread)
        enddo
        !$omp end do
    endif
    !$omp end parallel

    !! Close I/O files
    close(11)
    if (command == 1) close(31)
    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return
  end subroutine dm
!!------------------------------------------------------------------!!

!!--------------------------------------------------------------!!

  subroutine powerspectrum(pk)
    implicit none
    real, dimension(2,nc)       :: pk

    integer :: i,j,k,kg,mg,jg,ig
    integer :: k1,k2
    real    :: kr,kx,ky,kz,w1,w2,pow
    real, dimension(3,nc,nc_pen+2) :: pkt
    real, dimension(3,nc) :: pktsum
    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0

    pkt=0.0
    pktsum=0.0

    !! Compute power spectrum
    !! Cannot thread because of ind index
    do k = 1, nc_pen+mypadd 
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
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if (kr .ne. 0) then
                    k1=ceiling(kr)
                    k2=k1+1
                    w1=k1-kr
                    w2=1-w1
                    pow=sum((slab(i:i+1,j,k)/ncr**3)**2)
                    pkt(1,k1,k)=pkt(1,k1,k)+w1*pow
                    pkt(2,k1,k)=pkt(2,k1,k)+w1*pow**2
                    pkt(3,k1,k)=pkt(3,k1,k)+w1
                    pkt(1,k2,k)=pkt(1,k2,k)+w2*pow
                    pkt(2,k2,k)=pkt(2,k2,k)+w2*pow**2
                    pkt(3,k2,k)=pkt(3,k2,k)+w2
                endif
            enddo
        enddo
    enddo

    !! Merge power spectrum from threads
    do k=2,nc_pen+2
       pkt(:,:,1)=pkt(:,:,1)+pkt(:,:,k)
    enddo

    !! Reduce to rank 0
    call mpi_reduce(pkt(:,:,1),pktsum,3*nc,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

    !! Divide by weights
    !! pk(1,k) stores pk(k)
    !! pk(2,k) stores standard deviation
    if (rank == 0) then
      do k=1,nc
        if (pktsum(3,k) .eq. 0) then
          pk(:,k)=0
        else
          pk(:,k)=pktsum(1:2,k)/pktsum(3,k)
          pk(2,k)=sqrt(abs((pk(2,k)-pk(1,k)**2)/(pktsum(3,k)-1)))
          pk(1:2,k)=4.*pi*(real(k)-1.)**3*pk(1:2,k)
       endif
      enddo
    endif

    call mpi_bcast(pk,2*nc,mpi_real,0,mpi_comm_world,ierr)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called power spectrum'
    return
  end subroutine powerspectrum

!!------------------------------------------------------------------!!

  subroutine initvar
    implicit none
    integer :: k

    real time1,time2
    call cpu_time(time1)

    !! Initialize large arrays in parallel
    !$omp parallel default(shared) private(k)
    !$omp do
    do k=0,nc_node_dim+1
       phi(:,:,k)=0
    enddo
    !$omp end do
    !$omp do

    do k=1,nc_pen+2

       slab(:,:,k)=0

    enddo
    !$omp end do
    !$omp do
    do k=1,nk
       tf(:,k)=0
    enddo
    !$omp end do 
    !$omp do
    do k=1,nc
       pkm(:,k)=0
    enddo
    !$omp end do 
    !$omp do
    do k=1,nc
       pkn(:,k)=0
    enddo
    !$omp end do 
    !$omp end parallel
    
    !! Initialize fftw so that it generates the plans!
    firstfftw=.true.

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar

!!------------------------------------------------------------------!!

subroutine mesh_buffer
  implicit none

  integer(4) :: buffer_size
  integer(4) :: tag
  integer(4) :: status(MPI_STATUS_SIZE)

  buffer_size = (nc_node_dim + 2)**2

  if (rank==0) write(*,*) 'buffer_size =', buffer_size


  tag=64

!! send to node in -x

      phi_buf(:,:)=phi(1,:,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(5),tag,cart_neighbor(6), &
                              tag,mpi_comm_cart,status,ierr)

      phi(nc_node_dim+1,:,:)=phi(nc_node_dim+1,:,:)+phi_buf(:,:)

!! send to node in +x
   
      phi_buf(:,:)=phi(nc_node_dim,:,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(6),tag,cart_neighbor(5), &
                              tag,mpi_comm_cart,status,ierr)

      phi(0,:,:)=phi(0,:,:)+phi_buf(:,:)

!! send to node in -y

      phi_buf(:,:)=phi(:,1,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(3),tag,cart_neighbor(4), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,nc_node_dim+1,:)=phi(:,nc_node_dim+1,:)+phi_buf(:,:)

!! send to node in +y

      phi_buf(:,:)=phi(:,nc_node_dim,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(4),tag,cart_neighbor(3), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,0,:)=phi(:,0,:)+phi_buf(:,:)

!! send to node in -z
    
      phi_buf(:,:)=phi(:,:,1)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(1),tag,cart_neighbor(2), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,:,nc_node_dim+1)=phi(:,:,nc_node_dim+1)+phi_buf(:,:)

!! send to node in +z

      phi_buf(:,:)=phi(:,:,nc_node_dim)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(2),tag,cart_neighbor(1), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,:,0)=phi(:,:,0)+phi_buf(:,:)

  end subroutine mesh_buffer

!!------------------------------------------------------------------!!

  function power(kr,ix,iy)
   implicit none
    real    :: kr
    integer :: ix,iy

    integer :: i,i1,i2
    real    :: x,y,x1,x2,y1,y2
    real    :: power

    i1=1
    i2=nk
    do while (i2-i1 .gt. 1)
       i=(i1+i2)/2
       if (kr .gt. tf(ix,i)) then
          i1=i
       else
          i2=i
       endif
    enddo

    x1=log(tf(ix,i1))
    y1=log(tf(iy,i1))
    x2=log(tf(ix,i2))
    y2=log(tf(iy,i2))
    x =log(kr)
    y =y1+(y2-y1)*(x-x1)/(x2-x1)
    power=exp(y)
    
    return
  end function power

!!------------------------------------------------------------------!!
  function tophat(x)
    implicit none
    real :: x,tophat

    if (x .ne. 0) then
       tophat=3*(sin(x)-cos(x)*x)/x**3
    else
       tophat=1
    endif

    return
  end function tophat

end program dist_init 
