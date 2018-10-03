!write checkpoints to disk
subroutine checkpoint(dokill)
  implicit none
# include "cubepm.fh"

  logical, intent(in) :: dokill

  logical :: checkpoint_nu

  character (len=max_path) :: ofile
  character (len=6) :: rank_s
  character (len=100) :: z_s  

  integer(kind=4) :: i,j,fstat,blocksize,num_writes,nplow,nphigh
  integer(kind=4) :: cur_halo
  real(kind=4) :: z_write

  integer(4) :: np_dm, np_nu, ind_check1, ind_check2
  character (len=max_path) :: ofile_nu, ofile_nu_pid

  if (rank.eq.0) write(*,*) 'starting checkpoint'
  if (rank.eq.0 .and. kill_step) write(*,*) 'checkpoint kill'

  !Whether to open neutrino file or not
  checkpoint_nu = a .gt. a_i_nu
  np_dm=0
  do i=1,np_local
     if (pidmap(PID(i)).eq.pid_dm) np_dm=np_dm+1
  end do
  np_nu = np_local - np_dm 
  if (rank == 0) write(*,*) "checkpoint np_dm, np_nu = ", np_dm, np_nu

  !! label files with the same z as in the checkpoints file
  if (rank == 0) then
     if (dokill) then
        if (rank == 0) z_write = 1./a - 1.
     else
        z_write=z_checkpoint(cur_checkpoint)
     end if
  end if
  call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)

  ! Increment checkpoint counter so restart works on next checkpoint
  cur_halo = cur_halofind
  if (.not. dokill) then
     cur_checkpoint=cur_checkpoint+1
     !cur_halo = cur_halofind
     if (halofind_step) cur_halo = cur_halo + 1
  end if

  !! most linux systems choke when writing more than 2GB of data
  !! in one write statement, so break up into blocks < 2GB 
  !    blocksize=(2047*1024*1024)/24
  ! reduced to 32MB chunks because of intel compiler
  blocksize=(32*1024*1024)/24
  num_writes=np_local/blocksize+1

  ! Checkpoint file names
  write(rank_s,'(i6)') rank
  rank_s=adjustl(rank_s)

  write(z_s,'(f10.3)') z_write
  z_s=adjustl(z_s)

  ! Open checkpoint files
  !! Dark matter file
  ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
       rank_s(1:len_trim(rank_s))//'.dat'
  open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")
  if (fstat /= 0) then
     write(*,*) 'error opening checkpoint file for write'
     write(*,*) 'rank',rank,'file:',ofile
     call mpi_abort(mpi_comm_world,ierr,ierr)
  endif
  !Write header
  write(12) np_dm,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,1,cur_halo,mass_p

  if ( checkpoint_nu ) then
     !! Neutrino file
     ofile_nu=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
          rank_s(1:len_trim(rank_s))//'_nu.dat'
     open(unit=22, file=ofile_nu, status="replace", iostat=fstat, access="stream") 
     if (fstat /= 0) then
        write(*,*) 'error opening checkpoint file for write'
        write(*,*) 'rank',rank,'file:',ofile_nu
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     ! Write header
     write(22) np_nu,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,1,cur_halo,mass_p

#    ifdef UNIQUE_PID     
     ofile_nu_pid=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'pid'// &
          rank_s(1:len_trim(rank_s))//'_nu.dat'
     open(unit=23, file=ofile_nu_pid, status="replace", iostat=fstat, access="stream")
     if (fstat /= 0) then
        write(*,*) 'error opening checkpoint pid file for write'
        write(*,*) 'rank',rank,'file:',ofile_nu_pid
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     write(23) np_nu
#    endif
  end if

  ! Write data for neutrino simulation
  ind_check1 = 0
  ind_check2 = 0
  
  do i=1,num_writes
     nplow=(i-1)*blocksize+1
     nphigh=min(i*blocksize,np_local)
     do j=nplow,nphigh
        if (pidmap(PID(j)) == pid_dm) then
           write(12) xv(1:3,j) - shake_offset
           write(12) xv(4:6,j)
           ind_check1 = ind_check1 + 1
        else
           !Make sure we don't try to write to a non-opened file
           if (.not. checkpoint_nu) then
              if (rank.eq.0) write(*,*) 'ERROR IN CHECKPOINT! FOUND AN UNWANTED NEUTRINO'
              call mpi_abort(mpi_comm_world,ierr,ierr)
           end if
           write(22) xv(1:3,j) - shake_offset
           write(22) xv(4:6,j)
#          ifdef UNIQUE_PID
           write(23) PID(j)
#          endif
           ind_check2 = ind_check2 + 1
        endif
     enddo
  enddo

  close(12)
  if (checkpoint_nu) then
     close(22)
#    ifdef UNIQUE_PID
     close(23)
#    endif
  end if
  !! Consistency check
  if (ind_check1 .ne. np_dm .or. ind_check2 .ne. np_nu) then
     write(*,*) "Dark Matter checkpoint error: ind_checks ", ind_check1, np_dm, ind_check2, np_nu
     call mpi_abort(mpi_comm_world,ierr,ierr)
  endif
  
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
        print*, 'current steps recorded in xv file:',nts
        print*, 'cur_checkpoint =', cur_checkpoint
        print*, 'cur_halofind   =', cur_halo
        sec2a = mpi_wtime(ierr)
        write(*,*) 'stopping checkpoint',sec2a
        write(*,*) 'checkpoint time: ',sec2a-sec1a
     endif
  end if

  checkpoint_step=.false.

end subroutine checkpoint

! Subroutine that reads the first line of killtime_path which is to contain
! the amount of time (in seconds) remaining in the job. Will checkpoint with
! kill_remaining (parameter in cubepm.par) seconds remaining.
subroutine read_remaining_time
  implicit none
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
     write(*,*) "Killing CUBEP3M ", kill_time/3600., " hours from now"

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
