!! initialize particle list
subroutine particle_initialize
  use omp_lib
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  real(4) :: rnum,z_write,dummy
  integer(4) :: i,j,k,pp,fstat,blocksize,num_writes,nplow,nphigh
  integer*8 :: np_total,npl8
  character(len=max_path) :: ofile
  character(len=6) :: rank_s
  character(len=7) :: z_s, z_s2
  integer(4) :: np_nu

#ifdef NEUTRINOS
  integer(4) :: np_dm
  integer(8) :: np_total_nu

  !! Check that the PID flag is also defined
#ifndef PID_FLAG
  write(*,*) "ERROR: Using Neutrinos but PID is not enabled !!"
  call mpi_abort(mpi_comm_world,ierr,ierr)
#endif
#endif

  print*,'particle_initialize'
  fstat=0

  np_local=(nf_physical_node_dim/2)**3

  !! set flag in cubepm.par to select desired initial conditions
  if (restart_ic) then

     ! ---------------------------------------------------------------------------------------------
     ! Read in checkpoint
     ! ---------------------------------------------------------------------------------------------
     if (rank == 0) z_write = z_checkpoint(restart_checkpoint)
     call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
     write(z_s,'(f7.3)') z_write
     z_s=adjustl(z_s)

     write(rank_s,'(i6)') rank
     rank_s=adjustl(rank_s)

     ! Read dark matter particles

     ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
          rank_s(1:len_trim(rank_s))//'.dat'
     if (rank==0) print*, 'opening dark matter checkpoint file:',ofile

     open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
     
     if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif

     read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
          cur_projection,cur_halofind,mass_p
     
     if (rank == 0) print *,'restarting simulation from z=',z_checkpoint(cur_checkpoint-1)

     if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif

     !reduced to 32MB chunks because of intel compiler
     blocksize=(32*1024*1024)/24
     num_writes=np_local/blocksize+1

     do i = 1,num_writes
        nplow=(i-1)*blocksize+1
        nphigh=min(i*blocksize,np_local)
        do j=nplow,nphigh
           read(21) xv(:,j)
        enddo
     enddo
     close(21)

      ! Read neutrino particles 
#    ifdef NEUTRINOS
     ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
          rank_s(1:len_trim(rank_s))//'_nu.dat'
     if (rank==0) print*, 'opening neutrino checkpoint file:',ofile

     open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
     if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif

     !! Only store the local number of particles. All other info already read from dark matter checkpoint file.
     read(21) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
     
     if (rank == 0) print *,'neutrinos restarting simulation from z=',z_checkpoint(cur_checkpoint-1)

     if (np_local+np_nu > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif

     !reduced to 32MB chunks because of intel compiler
     blocksize=(32*1024*1024)/24
     num_writes=np_nu/blocksize+1
     do i = 1,num_writes
        nplow=(i-1)*blocksize+1 + np_local
        nphigh=min(i*blocksize,np_nu) + np_local
        do j=nplow,nphigh
           read(21) xv(:,j)
        enddo
     enddo
     close(21)
#    endif

#ifdef CHECKPOINT_KILL

     ! ---------------------------------------------------------------------------------------------
     ! Read in checkpoint kill
     ! ---------------------------------------------------------------------------------------------

  elseif (restart_kill) then

     write(rank_s,'(i6)') rank
     rank_s=adjustl(rank_s)

     ! Read dark matter particles
     ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'xvres'// &
          rank_s(1:len_trim(rank_s))//'.dat'

     open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
     if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif

     read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
          cur_projection,cur_halofind,mass_p

     if (rank == 0) print *,'restarting simulation from z=', reskill_prefix
     if (rank == 0) print *,'current checkpoint, proj and halo entries are:', cur_checkpoint, &
          cur_projection,cur_halofind
     if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif

     !reduced to 32MB chunks because of intel compiler
     blocksize=(32*1024*1024)/24
     num_writes=np_local/blocksize+1
     
     do i = 1,num_writes
        nplow=(i-1)*blocksize+1
        nphigh=min(i*blocksize,np_local)
        do j=nplow,nphigh
           read(21) xv(:,j)
        enddo
     enddo
     close(21)

#ifdef NEUTRINOS
     ! Read neutrino particles 
     ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//reskill_prefix//'xvres'// &
          rank_s(1:len_trim(rank_s))//'_nu.dat'

     open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
     if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif

     !! Only store the local number of particles. All other info already read from dark matter checkpoint file.
     read(21) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy

     if (rank == 0) print *,'neutrinos restarting simulation from z=',reskill_prefix

     if (np_local+np_nu > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
        
     !reduced to 32MB chunks because of intel compiler
     blocksize=(32*1024*1024)/24
     num_writes=np_nu/blocksize+1
     
     do i = 1,num_writes
        nplow=(i-1)*blocksize+1 + np_local
        nphigh=min(i*blocksize,np_nu) + np_local
        do j=nplow,nphigh
           read(21) xv(:,j)
        enddo
     enddo

     close(21)

#endif

  else

     ! ---------------------------------------------------------------------------------------------
     ! Read in ICs
     ! ---------------------------------------------------------------------------------------------
       
     write(z_s,'(f7.3)') z_i
     z_s=adjustl(z_s)

     write(rank_s,'(i6)') rank
     rank_s=adjustl(rank_s)

     ! Read dark matter particles
     ofile=ic_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'
     if (rank==0) print *,'opening particle list:',ofile(1:len_trim(ofile))
       
     open(unit=20, file=ofile, status="old", iostat=fstat, access="stream")
     if (fstat /= 0) then
        write(*,*) 'error opening initial conditions'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     read(20) np_local,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
     if (np_local > max_np) then
        write(*,*) 'too many particles to store'
        write(*,*) 'rank',rank,'np_local',np_local,'max_np',max_np
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     read(20) xv(:,:np_local)
     close(20)
       
  endif

  !! calculate total number of particles and particle mass
  npl8=int(np_local,kind=8)
  call mpi_reduce(npl8,np_total,1,MPI_INTEGER8, &
       mpi_sum,0,mpi_comm_world,ierr)
  if (rank == 0) write(*,*) 'number of particles =', np_total
  call mpi_bcast(np_total,1,MPI_INTEGER8,0,mpi_comm_world,ierr)
  
  if (.not.restart_ic) then
     mass_p = real(nf_physical_dim)**3 / real(np_total)
  endif

  if (rank == 0) write(*,*) 'particle mass=', mass_p
  if (rank == 0) write(*,*) 'total dark matter mass =', mass_p * np_total

#ifdef NEUTRINOS
  ! Assign 1 byte integer PIDs
  do i = 1, np_local
     PID(i) = 1 !! All dark matter will have a PID of 1
  enddo
  do i = np_local+1, np_local+np_nu
     PID(i) = 2 !! All neutrinos will have a PID of 2 
  enddo

  ! Print some stats to screen
  
  npl8=int(np_nu,kind=8)
  call mpi_reduce(npl8,np_total_nu,1,MPI_INTEGER8, mpi_sum,0,mpi_comm_world,ierr)
  !! Append np_nu to np_local (must be done after PIDs and after mass calculation)
  np_dm = np_local
  np_local = np_local + np_nu
  if (rank == 0) then
     write(*,*) "np_dm = ", np_dm
     write(*,*) "np_nu = ", np_nu
     write(*,*) "np_local = ", np_local
     write(*,*) "np_dm_total = ", np_total
     write(*,*) "np_nu_total = ", np_total_nu
  endif
#endif

end subroutine particle_initialize
