!write checkpoints to disk
subroutine checkpoint(dokill)
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  logical, intent(in) :: dokill

  character (len=max_path) :: ofile,ofile2
  character (len=6) :: rank_s
  character (len=7) :: z_s  

  integer(kind=4) :: i,j,fstat,blocksize,num_writes,nplow,nphigh
  integer(kind=4) :: cur_proj,cur_halo
  real(kind=4) :: z_write
#ifdef NEUTRINOS
  integer(4) :: np_dm, np_nu, ind_check1, ind_check2
  character (len=max_path) :: ofile_nu, ofile2_nu
  integer, parameter :: pdm = 1
#endif

  !! label files with the same z as in the checkpoints file
  if (rank == 0) then
     if (dokill) then
        if (rank == 0) z_write = 1./a - 1.
     else
        z_write=z_checkpoint(cur_checkpoint)
     end if
  end if
  call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

  !! most linux systems choke when writing more than 2GB of data
  !! in one write statement, so break up into blocks < 2GB 
  !    blocksize=(2047*1024*1024)/24
  ! reduced to 32MB chunks because of intel compiler
  blocksize=(32*1024*1024)/24
  num_writes=np_local/blocksize+1

  ! Checkpoint file names
  write(rank_s,'(i6)') rank
  rank_s=adjustl(rank_s)

  write(z_s,'(f7.3)') z_write
  z_s=adjustl(z_s)

  !! Dark matter file
  ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
       rank_s(1:len_trim(rank_s))//'.dat'

#ifdef NEUTRINOS
  !! Neutrino file
  ofile_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
       rank_s(1:len_trim(rank_s))//'_nu.dat'
#endif

  ! Increment checkpoint counter so restart works on next checkpoint
  if (.not. dokill) then
     cur_checkpoint=cur_checkpoint+1
     cur_proj = cur_projection
     cur_halo = cur_halofind
     if (projection_step) cur_proj = cur_proj + 1
     if (halofind_step) cur_halo = cur_halo + 1
  end if
  
  ! Open checkpoint files
  
  !! Dark matter 
  open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")
  if (fstat /= 0) then
     write(*,*) 'error opening checkpoint file for write'
     write(*,*) 'rank',rank,'file:',ofile
     call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

#ifdef NEUTRINOS
  !! Neutrinos 
  open(unit=22, file=ofile_nu, status="replace", iostat=fstat, access="stream") 
  if (fstat /= 0) then
     write(*,*) 'error opening checkpoint file for write'
     write(*,*) 'rank',rank,'file:',ofile_nu
     call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

  ! Write file headers
  !! Determine how many dark matter and neutrino particles this rank has
  np_dm = count(PID(1:np_local) == pdm)
  np_nu = np_local - np_dm 
  if (rank == 0) write(*,*) "checkpoint np_dm, np_nu = ", np_dm, np_nu
  write(12) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p
  write(22) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p
#else
  write(12) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_proj,cur_halo,mass_p
#endif

#ifdef NEUTRINOS
  ! Write data for neutrino simulation
  ind_check1 = 0
  ind_check2 = 0
  
  do i=1,num_writes
     nplow=(i-1)*blocksize+1
     nphigh=min(i*blocksize,np_local)
     do j=nplow,nphigh
        if (PID(j) == 1) then
           write(12) xv(1:3,j) - shake_offset
           write(12) xv(4:6,j)
           ind_check1 = ind_check1 + 1
        else
           write(22) xv(1:3,j) - shake_offset
           write(22) xv(4:6,j)
           ind_check2 = ind_check2 + 1
        endif
     enddo
  enddo

  close(12)
  close(22)

  !! Consistency check
  if (ind_check1 .ne. np_dm .or. ind_check2 .ne. np_nu) then
     write(*,*) "Dark Matter checkpoint error: ind_checks ", ind_check1, np_dm, ind_check2, np_nu
     call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

#else
  ! Write data for non-neutrino simulation

  do j = 1, np_local
     xv(1:3,j) = xv(1:3,j) - shake_offset
  enddo

  do i=1,num_writes
     nplow=(i-1)*blocksize+1
     nphigh=min(i*blocksize,np_local)
     write(12) xv(:,nplow:nphigh)
  enddo

  close(12)

  do j = 1, np_local
     xv(1:3,j) = xv(1:3,j) + shake_offset
  enddo

#endif
  
  if (dokill) then
     if (rank.eq.0) then
        write(*,*) 'checkpoint kill completed'
        write(*,*) 'exiting code'
     end if
     call mpi_barrier(mpi_comm_world,ierr)
     call mpi_abort(mpi_comm_world,ierr,ierr)
  else
     ! Print info to screen
     if (rank==0) then
        print*, 'current steps recorded in xv file:'
        print*, 'cur_checkpoint =', cur_checkpoint
        print*, 'cur_projection =', cur_proj
        print*, 'cur_halofind   =', cur_halo

        write(*,*) "BPOS: ", minval(xv(1:3,1:np_local)), maxval(xv(1:3,1:np_local))
        write(*,*) "BVEL: ", minval(abs(xv(4:6,1:np_local))), maxval(abs(xv(4:6,1:np_local)))
     endif
  end if

  checkpoint_step=.false.

end subroutine checkpoint

! Subroutine that reads the first line of killtime_path which is to contain
! the amount of time (in seconds) remaining in the job. Will checkpoint with
! kill_remaining (parameter in cubepm.par) seconds remaining.
subroutine read_remaining_time
  implicit none
  include 'mpif.h'
# include "cubepm.fh"
  
  integer(4) :: time_lefti
  real(4) :: time_left
  real(4) :: time_left_default = 48.*3600. !! Maximum walltime of GPC is taken if problems exist

  if (rank == 0) then !! Only master rank will read this file                                                                                          

     open(31, file=killtime_path, status="old", iostat=ierr)
     if (ierr /= 0) then !! If file does not exist then assume default time remaining
        time_left = time_left_default
        write(*,*) "WARNING: Could not open file ", killtime_path
        write(*,*) "Assuming the time left (in hours) is ", time_left/3600.
     else !! Read the file                                                                                                                            
        read(unit=31,fmt='(I20)') time_lefti
        close(31)
        time_left = 1.*time_lefti
     endif

     !! Kill the job in this amount of time from now.                                                                                                 
     kill_time = time_left - kill_remaining
     write(*,*) "Killing job ", kill_time/3600., " hours from now"

     !! Consistency check                                                                                                                             
     if (kill_time <= 0.) then
        kill_time = time_left_default
        write(*,*) "WARNING: kill_time did not make sense. Taking kill_time (in hours): ", kill_time/3600.
        write(*,*) "You should check these: ", time_left, time_left_default, kill_remaining
     endif

  endif

  !! Broadcast to everyone                                                                                                                             
  call mpi_bcast(kill_time, 1, mpi_real, 0, mpi_comm_world, ierr)
  
end subroutine read_remaining_time
