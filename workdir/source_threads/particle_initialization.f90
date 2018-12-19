!! initialize particle list
subroutine particle_initialize(ispec)
  use omp_lib
  implicit none
# include "cubepm.fh"

  integer(4), intent(in) :: ispec

  real(4) :: rnum,z_write,dummy
  integer(4) :: i,j,k,pp,fstat,blocksize,num_writes,nplow,nphigh
  integer*8 :: np_total,npl8
  character(len=max_path) :: ofile
  character(len=6) :: rank_s
  character(len=100) :: z_s, z_s2
  integer(4) :: np_nu

  integer, dimension(nodes_dim**3) :: np_nu_rank

  integer(4) :: np_dm
  integer(8) :: np_total_nu

  if (rank.eq.0) then
     write(*,*)
     print*,'particle_initialize'
  end if

  fstat=0

  write(rank_s,'(i6)') rank
  rank_s=adjustl(rank_s)

  !! set flag in cubepm.par to select desired initial conditions
  if (restart_ic) then
     !Checkpoint
     if (rank == 0) z_write = z_checkpoint(restart_checkpoint)
     call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
     write(z_s,'(f10.3)') z_write
     z_s=adjustl(z_s)
  elseif (restart_kill) then
     !Checkpoint kill
     z_s=reskill_prefix
     z_s=adjustl(z_s)
  elseif (ispec .ne. pid_dm) then
     !Neutrino IC
     if (rank==0) z_write = z_i_nu!z_checkpoint(cur_checkpoint-1) !checkpoint comes first!
     call mpi_bcast(z_write,1,mpi_real,0,mpi_comm_world,ierr)
     write(z_s,'(f10.3)') z_write
     z_s=adjustl(z_s)
  else
     !Dark matter IC
     write(z_s,'(f10.3)') z_i
     z_s=adjustl(z_s)
  endif
  
  np_dm=0
  np_nu=0

  !Amount of unclustered matter
  f_unclustered = 1.d0-omega_c/omega_m

  if (ispec.eq.pid_dm) then

     !Read dark matter particles
     if (rank.eq.0) write(*,*) 'reading in dm particles'

     ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
          rank_s(1:len_trim(rank_s))//'.dat'

     if (rank==0) print*, 'opening dark matter file:',trim(adjustl(ofile))
     open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
     if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif

     if (restart_ic .or. restart_kill) then
        read(21) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
             dummy,cur_halofind,mass_p
        if (rank == 0) print *,'restarting simulation from z=',1./a-1.
     else
        read(21) np_local,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
        if (rank == 0) print *,'starting simulation from z=',z_i
     endif
     
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

     !! calculate total number of particles and particle mass
     npl8=int(np_local,kind=8)
     call mpi_reduce(npl8,np_total,1,MPI_INTEGER8, &
          mpi_sum,0,mpi_comm_world,ierr)
     if (rank == 0) write(*,*) 'number of particles =', np_total
     call mpi_bcast(np_total,1,MPI_INTEGER8,0,mpi_comm_world,ierr)
  
     if (.not.restart_ic .and. .not.restart_kill) then
        !Main code does not care about this parameter, but halofinder does
        mass_p = real(nf_physical_dim,kind=8)**3 / real(np_total)
     endif

     if (rank == 0) write(*,*) 'particle mass=', mass_p
     if (rank == 0) write(*,*) 'total dark matter mass =', mass_p * np_total

     ! Assign 1 byte integer PIDs
     do i = 1, np_local
        PID(i) = pid_dm
     enddo

     !Set appropriate mass 
     mass_p_nudm_fac(pid_dm) = (real(nf_physical_dim,kind=8)**3/np_total/mass_p) * omega_c/omega_m

  else

     !Read neutrino particles 
     if (rank.eq.0) write(*,*) 'reading in nu particles'

     ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'// &
          rank_s(1:len_trim(rank_s))//'_nu.dat'
     if (rank==0) print*, 'opening neutrino file:',trim(adjustl(ofile))

     open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
     if (fstat /= 0) then
        write(*,*) 'error opening checkpoint'
        write(*,*) 'rank',rank,'file:',ofile
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif

     !! Only store the local number of particles. All other info already read from dark matter checkpoint file.
     read(21) np_nu,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
     
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

#    ifdef UNIQUE_PID
     if (restart_ic .or. restart_kill) then
        !Read in PID file
        ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'pid'// &
            rank_s(1:len_trim(rank_s))//'_nu.dat'
        open(unit=21, file=ofile, status="old", iostat=fstat, access="stream")
        if (fstat /= 0) then
           write(*,*) 'error opening checkpoint pid'
           write(*,*) 'rank',rank,'file:',ofile
           call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
        read(21) i
        if (i.ne. np_nu) then
           write(*,*) 'error in pid file: np_nu does not match'
           write(*,*) 'rank',rank,'expected/received np_nu',np_nu,i
           call mpi_abort(mpi_comm_world,ierr,ierr)
        end if
        do i=np_local+1,np_local+np_nu
           read(21) PID(i)
        end do
        close(21)
     else
        !Assign unique 8 byte integer PIDs, assumes each node starts with np_nu particles, write to file
        ofile=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'pid'// &
             rank_s(1:len_trim(rank_s))//'_nu.dat'
        open(unit=23, file=ofile, status="replace", iostat=fstat, access="stream")
        if (fstat /= 0) then
           write(*,*) 'error opening checkpoint pid file for write'
           write(*,*) 'rank',rank,'file:',ofile
           call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
        write(23) np_nu

        !Find how many bh in previous nodes
        call mpi_allgather(np_nu,1,mpi_integer,np_nu_rank,1,mpi_integer,mpi_comm_world,ierr)
        np_nu_rank(1)=sum(np_nu_rank(1:rank+1))-np_nu_rank(rank+1)
        do i = np_local+1, np_local+np_nu
           PID(i) = int(i-np_local-1+pid_nu,kind=pid_int_kind)+int(np_nu_rank(1),kind=pid_int_kind)
           write(23) PID(i)
        end do

        close(23)

     end if
#    else
     ! Assign identical 1 byte integer PIDs
     do i = np_local+1, np_local+np_nu
        PID(i) = pid_nu
     enddo
#    endif

     !! calculate total number of particles 
     npl8=int(np_local,kind=8)
     call mpi_reduce(npl8,np_total,1,MPI_INTEGER8, &
          mpi_sum,0,mpi_comm_world,ierr)
     if (rank == 0) write(*,*) 'number of particles =', np_total
     call mpi_bcast(np_total,1,MPI_INTEGER8,0,mpi_comm_world,ierr)

     npl8=int(np_nu,kind=8)
     call mpi_reduce(npl8,np_total_nu,1,MPI_INTEGER8, &
          mpi_sum,0,mpi_comm_world,ierr)
     if (rank == 0) write(*,*) 'number of particles =', np_total_nu
     call mpi_bcast(np_total_nu,1,MPI_INTEGER8,0,mpi_comm_world,ierr)     

     !Set appropriate mass
     mass_p_nudm_fac(pid_nu) = (real(nf_physical_dim,kind=8)**3/np_total_nu/mass_p) * omega_nu/omega_m

     !Correct unclustered amount
     f_unclustered = 1.d0-omega_c/omega_m-omega_nu/omega_m     

  end if

#ifdef TIDAL_FIELD
  call initialize_tidal_field
#endif

  !Print some stats to screen
  np_dm = np_local
  np_local = np_local + np_nu
  if (rank == 0) then
     write(*,*) "np_dm = ", np_dm
     write(*,*) "np_nu = ", np_nu
     write(*,*) "np_local = ", np_local
     write(*,*) "np_dm_total = ", np_total
     write(*,*) "np_nu_total = ", np_total_nu

     write(*,*) "mass_p = ",mass_p
     write(*,*) "mass_p_nudm_fac: dm = ",mass_p_nudm_fac(pid_dm)
     write(*,*) "mass_p_nudm_fac: nu = ",mass_p_nudm_fac(pid_nu)

     write(*,*) 'Dark matter mass on grid',1.d0*np_total*mass_p*mass_p_nudm_fac(pid_dm)/real(nf_physical_dim,kind=8)**3
     write(*,*) 'Neutrino mass on grid',1.d0*np_total_nu*mass_p*mass_p_nudm_fac(pid_nu)/real(nf_physical_dim,kind=8)**3
     write(*,*) 'Unclustered mass on grid',1.d0*f_unclustered
     write(*,*) 'Total=1: ',real(1.d0*f_unclustered+&
          &(1.d0*np_total*mass_p_nudm_fac(pid_dm)+1.d0*np_total_nu*mass_p_nudm_fac(pid_nu))&
          &*mass_p/real(nf_physical_dim,kind=8)**3)
  endif

  if (rank == 0) write(*,*) 'finished initializing particles'

end subroutine particle_initialize
