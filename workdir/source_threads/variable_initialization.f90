!! initialize variables
subroutine variable_initialize
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  integer(4) :: i,fstat

  ierr=-1
  fstat=0

  !! initial scalefactor, see cubepm.par
  a=a_i
  tau=-3.0/sqrt(a_i)

  !! set these variables prohibitively high such that they don't 
  !! bugger up the first half-step update 
  !! (1st sweep starts with a half-step position update)

  dt_f_acc=1000.0
  dt_c_acc=1000.0

  dt_pp_acc=1000.0
#ifdef PP_EXT
  dt_pp_ext_acc=1000.0
#endif

  dt_vmax=1000.0

  !! zero everything else 
  dt_old=0.0
  dt=0.0
  da=0.0
  nts=0
  t=0.0
  np_local=0
  mass_p=0.0
  mass_p_nudm_fac=0.0
    
  rho_f=0.0
  kern_f=0.0
  force_f=0.0
  kern_c=0.0
  slab=0.0
  send_buf=0.0
  recv_buf=0.0
  xv=0.0
  ll=0
  hoc=0
  np_buf=0
  final_step=.false.
  kill_step=.false.
  injection_step= (a_i.eq.a_i_nu) .and. (n_bh.gt.0)
  shake_offset=0.0

  PID=0
  send_buf_PID=0.0
  recv_buf_PID=0.0
  
  if (rank == 0) then

     !! print cosmology
     write(*,*) 'cosmology: '
     write(*,*) '>>omega_m: ',omega_m
     write(*,*) '>>omega_r: ',omega_r
     write(*,*) '>>omega_l: ',omega_l
     write(*,*) '>>omega_k: ',1.-omega_m-omega_r-omega_l

     !! read in when to store checkpoints
     open(11,file=checkpoints,status='old',iostat=fstat)
     if (fstat /= 0) then
        write(*,*) 'error opening checkpoint list file'
        write(*,*) 'rank',rank,'file:',checkpoints
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     do num_checkpoints=1,max_input
        read(unit=11,err=51,end=41,fmt='(f10.4)') z_checkpoint(num_checkpoints)
     enddo
41   num_checkpoints=num_checkpoints-1
51   close(11)

     if (z_checkpoint(1).eq.z_i) then
        z_checkpoint(1:num_checkpoints-1) = z_checkpoint(2:num_checkpoints)
        num_checkpoints=num_checkpoints-1
     end if

     if (num_checkpoints.eq.1) then
        write(*,*) 'problem reading checkpoints '
     endif

     if (num_checkpoints.eq.max_input) then
        write(*,*) 'too many checkpoints to store > ',max_input
        call mpi_abort(mpi_comm_world,i,ierr)
     endif

     if (z_checkpoint(1).gt.z_i) then
        write(*,*) 'z_initial less than first checkpoint, exiting',z_i,z_checkpoint(1)
        call mpi_abort(mpi_comm_world,i,ierr)
     endif

     do i=1,num_checkpoints
        a_checkpoint(i)=1.0/(1.0+z_checkpoint(i))
     enddo

     write(*,*)
     write(*,*) 'Starting simulation at:'
     write(*,*) 'z      a'
     write(*,'(f10.3,2x,f10.3)') z_i,a_i
     if (num_checkpoints > 0) then
        write(*,*)
        write(*,*) 'Checkpointing performed at:'
        do i=1,num_checkpoints
           write(*,'(f10.3,2x,f10.3)') z_checkpoint(i),a_checkpoint(i)
        enddo
     else
        write(*,*) 'no checkpoints to be stored'
     endif

     !Usually just do halofinds at checkpoints
     num_halofinds=num_checkpoints
     z_halofind=z_checkpoint
     a_halofind=1.0/(1.0+z_halofind)

     if (num_halofinds > 0) then
        write(*,*)
        write(*,*) 'Halo catalogs generated at:'
        write(*,*) 'z        a'
        do i=1,num_halofinds
           write(*,'(f10.3,2x,f10.3)') z_halofind(i),a_halofind(i)
        enddo
     else
        a_halofind(1)=100.0
        write(*,*) 'no halo catalogs to be stored'
     endif

  endif

  cur_halofind=1
  cur_checkpoint=1

  ! Initialize first time through fftw flag
  firstfftw=.true.
  firstfftw_nest=.true.
  firstfftw2=.true.

  ! NEW TRICK BY JHD TO REMOVE MEMORY RACING CONDITION ON THREADED PLAN
  ! CREATION. 
  do i = 1,cores 
     call cubepm_fftw2('o',i)
  enddo

  ! Initialize halo finding arrays
  call initialize_halofind

  call omp_set_num_threads(cores*nested_threads)
  call omp_set_nested(.true.)
  if (rank == 0) then
     write(*,*)
     write(*,*) 'omp setup'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'finished variable init'
  end if

end subroutine variable_initialize
