!! pass particles to adjacent nodes
subroutine particle_pass
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  real(4), parameter :: rnf_buf = nf_buf
  integer(4) :: np_max
  integer(4) :: np_buf_max,np_buf_max_dir
  integer(4) :: np_local0
  
  integer(4) :: i,j,k,pp,tag
  integer(4) :: nppx,nppy,nppz,npmx,npmy,npmz
  integer(4), dimension(mpi_status_size) :: status,sstatus,rstatus
  integer(4) :: srequest,rrequest,sierr,rierr
  integer(4) :: ikill, ikill_loc

  np_buf_max_dir = 0

  call mpi_barrier(mpi_comm_world,ierr) !necessary?

  !! Keep track of np_local at the start in case we need to checkpoint_kill
  np_local0 = np_local

  !! If this variable becomes 1 then we must checkpoint kill
  ikill_loc = 0

  ! pass +x
  tag=11
  np_buf = 0
  do k=hoc_nc_l,hoc_nc_h
     do j=hoc_nc_l,hoc_nc_h
        do i=hoc_nc_h - hoc_pass_depth, hoc_nc_h 
           pp = hoc(i,j,k)
           do while (pp /= 0)
              if (xv(1,pp) >= nf_physical_node_dim - rnf_buf) then
                 np_buf = np_buf + 1
                 send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
                 send_buf_PID(np_buf)=PID(pp)
              endif
              pp = ll(pp)
           enddo
        enddo
     enddo
  enddo

  !! Check to see if we need to checkpoint kill
  call check_buf_num(np_local0,'+x pass')

  nppx = np_buf
  
  call mpi_reduce(nppx,np_buf_max_dir,1,mpi_integer, &
       mpi_max,0,mpi_comm_world,ierr)
  if(rank==0) then
     if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
  endif

  call mpi_sendrecv_replace(nppx,1,mpi_integer,cart_neighbor(6), &
       tag,cart_neighbor(5),tag,mpi_comm_world, &
       status,ierr)

  !! Check to see if we need to checkpoint kill
  call check_par_num(np_local,nppx,np_local0,'+x pass')

  call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(6), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf,nppx*6,mpi_real,cart_neighbor(5), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(6), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf_PID,nppx,MPI_integer1,cart_neighbor(5), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  do i=1,nppx
     xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
     xv(1,np_local+i)=max(xv(1,np_local+i)-nf_physical_node_dim,-rnf_buf)
     PID(np_local+i)=recv_buf_PID(i)
  enddo

  np_local=np_local+nppx

  ! pass -x
  np_buf = 0
  do k=hoc_nc_l,hoc_nc_h
     do j=hoc_nc_l,hoc_nc_h
        do i=hoc_nc_l,hoc_nc_l + hoc_pass_depth
           pp = hoc(i,j,k)
           do while (pp /= 0)
              if (xv(1,pp) < rnf_buf) then
                 np_buf = np_buf + 1
                 send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
                 send_buf_PID(np_buf)=PID(pp)
              endif
              pp = ll(pp)
           enddo
        enddo
     enddo
  enddo

  !! Check to see if we need to checkpoint kill
  call check_buf_num(np_local0,'-x pass')

  npmx = np_buf
  
  call mpi_reduce(npmx,np_buf_max_dir,1,mpi_integer, &
       mpi_max,0,mpi_comm_world,ierr)
  if(rank==0) then
     if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
  endif

  call mpi_sendrecv_replace(npmx,1,mpi_integer,cart_neighbor(5), &
       tag,cart_neighbor(6),tag,mpi_comm_world, &
       status,ierr)

  !! Check to see if we need to checkpoint kill
  call check_par_num(np_local,npmx,np_local0,'-x pass')

  call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(5), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf,npmx*6,mpi_real,cart_neighbor(6), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(5), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf_PID,npmx,MPI_integer1,cart_neighbor(6), &
       tag,mpi_comm_world,rrequest,rierr)
    
  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)


  do i=1,npmx
     xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
     PID(np_local+i)=recv_buf_PID(i)
     if (abs(xv(1,np_local+i)).lt.eps) then
        if (xv(1,np_local+i) < 0.0) then
           xv(1,np_local+i)=-eps
        else
           xv(1,np_local+i)=eps
        endif
     endif
     xv(1,np_local+i)=min(xv(1,np_local+i)+real(nf_physical_node_dim,4), &
          nf_physical_node_dim+rnf_buf-eps)
  enddo
  
  np_local=np_local+npmx

  ! add additional particles to linked list
  !! should add/subtract offsets here!

  pp=np_local-npmx-nppx+1
  do
     if (pp > np_local) exit
     i=floor(xv(1,pp)/mesh_scale)+1
     j=floor(xv(2,pp)/mesh_scale)+1
     k=floor(xv(3,pp)/mesh_scale)+1
     if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
          j < hoc_nc_l .or. j > hoc_nc_h .or. &
          k < hoc_nc_l .or. k > hoc_nc_h) then
        write(*,*) 'x-pass particle out of link list range - deleted'
        write(*,*) rank,xv(:,pp),pp,i,j,k,hoc_nc_l,hoc_nc_h
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     ll(pp)=hoc(i,j,k)
     hoc(i,j,k)=pp
     pp=pp+1
  enddo

  ! pass -y
  np_buf = 0
  do k=hoc_nc_l,hoc_nc_h
     do j=hoc_nc_l,hoc_nc_l + hoc_pass_depth
        do i=hoc_nc_l,hoc_nc_h
           pp = hoc(i,j,k)
           do while (pp /= 0)
              if (xv(2,pp) < rnf_buf) then
                 np_buf = np_buf + 1
                 send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
                 send_buf_PID(np_buf)=PID(pp)
              endif
              pp = ll(pp)
           enddo
        enddo
     enddo
  enddo

  !! Check to see if we need to checkpoint kill
  call check_buf_num(np_local0,'-y pass')

  npmy = np_buf

  call mpi_reduce(npmy,np_buf_max_dir,1,mpi_integer, &
       mpi_max,0,mpi_comm_world,ierr)
  if(rank==0) then
     if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
  endif

  call mpi_sendrecv_replace(npmy,1,mpi_integer,cart_neighbor(3), &
       tag,cart_neighbor(4),tag,mpi_comm_world, &
       status,ierr)

  !! Check to see if we need to checkpoint kill
  call check_par_num(np_local,npmy,np_local0,'-y pass')

  call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(3), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf,npmy*6,mpi_real,cart_neighbor(4), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(3), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf_PID,npmy,MPI_integer1,cart_neighbor(4), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  do i=1,npmy
     xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
     PID(np_local+i)=recv_buf_PID(i)
     if (abs(xv(2,np_local+i)).lt.eps) then
        if (xv(2,np_local+i) < 0.0) then
           xv(2,np_local+i)=-eps
        else
           xv(2,np_local+i)=eps
        endif
     endif
     xv(2,np_local+i)=min(xv(2,np_local+i)+real(nf_physical_node_dim,4), &
          nf_physical_node_dim+rnf_buf-eps)
  enddo
  np_local=np_local+npmy


  ! pass +y
  tag=11
  np_buf = 0
  do k=hoc_nc_l,hoc_nc_h
     do j=hoc_nc_h - hoc_pass_depth, hoc_nc_h 
        do i=hoc_nc_l,hoc_nc_h
           pp = hoc(i,j,k)
           do while (pp /= 0)
              if (xv(2,pp) >= nf_physical_node_dim - rnf_buf) then
                 np_buf = np_buf + 1
                 send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
                 send_buf_PID(np_buf)=PID(pp)
              endif
              pp = ll(pp)
           enddo
        enddo
     enddo
  enddo

  !! Check to see if we need to checkpoint kill
  call check_buf_num(np_local0,'+y pass')

  nppy = np_buf

  call mpi_reduce(nppy,np_buf_max_dir,1,mpi_integer, &
       mpi_max,0,mpi_comm_world,ierr)
  if(rank==0) then
     if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
  endif

  call mpi_sendrecv_replace(nppy,1,mpi_integer,cart_neighbor(4), &
       tag,cart_neighbor(3),tag,mpi_comm_world, &
       status,ierr)

  !! Check to see if we need to checkpoint kill
  call check_par_num(np_local,nppy,np_local0,'+y pass')

  call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(4), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf,nppy*6,mpi_real,cart_neighbor(3), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(4), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf_PID,nppy,MPI_integer1,cart_neighbor(3), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  do i=1,nppy
     xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
     xv(2,np_local+i)=max(xv(2,np_local+i)-nf_physical_node_dim,-rnf_buf)
     PID(np_local+i)=recv_buf_PID(i)
  enddo
  np_local=np_local+nppy

  ! add additional particles to linked list
  pp=np_local-npmy-nppy+1
  do
     if (pp > np_local) exit
     i=floor(xv(1,pp)/mesh_scale)+1
     j=floor(xv(2,pp)/mesh_scale)+1
     k=floor(xv(3,pp)/mesh_scale)+1
     if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
          j < hoc_nc_l .or. j > hoc_nc_h .or. &
          k < hoc_nc_l .or. k > hoc_nc_h) then
        write(*,*) 'y-pass particle out of link list range - deleted'
        write(*,*) rank,xv(:,pp),pp,i,j,k,hoc_nc_l,hoc_nc_h
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     ll(pp)=hoc(i,j,k)
     hoc(i,j,k)=pp
     pp=pp+1
  enddo

! pass +z
  tag=11
  np_buf = 0
  do k=hoc_nc_h - hoc_pass_depth, hoc_nc_h 
     do j=hoc_nc_l, hoc_nc_h
        do i=hoc_nc_l, hoc_nc_h
           pp = hoc(i,j,k)
           do while (pp /= 0)
              if (xv(3,pp) >= nf_physical_node_dim - rnf_buf) then
                 np_buf = np_buf + 1
                 send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
                 send_buf_PID(np_buf)=PID(pp)
              endif
              pp = ll(pp)
           enddo
        enddo
     enddo
  enddo

  !! Check to see if we need to checkpoint kill
  call check_buf_num(np_local0,'+z pass')

  nppz = np_buf

  call mpi_reduce(nppz,np_buf_max_dir,1,mpi_integer, &
       mpi_max,0,mpi_comm_world,ierr)
  if(rank==0) then
     if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
  endif

  call mpi_sendrecv_replace(nppz,1,mpi_integer,cart_neighbor(2), &
       tag,cart_neighbor(1),tag,mpi_comm_world, &
       status,ierr)

  !! Check to see if we need to checkpoint kill
  call check_par_num(np_local,nppz,np_local0,'+z pass')

  call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(2), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf,nppz*6,mpi_real,cart_neighbor(1), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(2), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf_PID,nppz,MPI_integer1,cart_neighbor(1), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  do i=1,nppz
     xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
     xv(3,np_local+i)=max(xv(3,np_local+i)-nf_physical_node_dim,-rnf_buf)
     PID(np_local+i)=recv_buf_PID(i)
  enddo

  np_local=np_local+nppz

  ! pass -z
  np_buf = 0
  do k=hoc_nc_l,hoc_nc_l + hoc_pass_depth
     do j=hoc_nc_l,hoc_nc_h
        do i=hoc_nc_l,hoc_nc_h
           pp = hoc(i,j,k)
           do while (pp /= 0)
              if (xv(3,pp) < rnf_buf) then
                 np_buf = np_buf + 1
                 send_buf((np_buf-1)*6+1:np_buf*6)=xv(:,pp)
                 send_buf_PID(np_buf)=PID(pp)
              endif
              pp = ll(pp)
           enddo
        enddo
     enddo
  enddo

  !! Check to see if we need to checkpoint kill
  call check_buf_num(np_local0,'-z pass')

  npmz = np_buf

  call mpi_reduce(npmz,np_buf_max_dir,1,mpi_integer, &
       mpi_max,0,mpi_comm_world,ierr)
  if(rank==0) then
     if(np_buf_max_dir > np_buf_max) np_buf_max = np_buf_max_dir
  endif

  call mpi_sendrecv_replace(npmz,1,mpi_integer,cart_neighbor(1), &
       tag,cart_neighbor(2),tag,mpi_comm_world, &
       status,ierr)

  !! Check to see if we need to checkpoint kill
  call check_par_num(np_local,npmz,np_local0,'-z pass')

  call mpi_isend(send_buf,np_buf*6,mpi_real,cart_neighbor(1), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf,npmz*6,mpi_real,cart_neighbor(2), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_isend(send_buf_PID,np_buf,MPI_integer1,cart_neighbor(1), &
       tag,mpi_comm_world,srequest,sierr)
  call mpi_irecv(recv_buf_PID,npmz,MPI_integer1,cart_neighbor(2), &
       tag,mpi_comm_world,rrequest,rierr)

  call mpi_wait(srequest,sstatus,sierr)
  call mpi_wait(rrequest,rstatus,rierr)

  do i=1,npmz
     xv(:,np_local+i)=recv_buf((i-1)*6+1:(i-1)*6+6)
     PID(np_local+i)=recv_buf_PID(i)
     if (abs(xv(3,np_local+i)).lt.eps) then
        if (xv(3,np_local+i) < 0.0) then
           xv(3,np_local+i)=-eps
        else
           xv(3,np_local+i)=eps
        endif
     endif
     xv(3,np_local+i)=min(xv(3,np_local+i)+nf_physical_node_dim, &
          nf_physical_node_dim+rnf_buf-eps)
  enddo

  np_local=np_local+npmz

  ! add additional particles to linked list
  pp=np_local-npmz-nppz+1
  do
     if (pp > np_local) exit
     i=floor(xv(1,pp)/mesh_scale)+1
     j=floor(xv(2,pp)/mesh_scale)+1
     k=floor(xv(3,pp)/mesh_scale)+1
     if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
          j < hoc_nc_l .or. j > hoc_nc_h .or. &
          k < hoc_nc_l .or. k > hoc_nc_h) then
        write(*,*) 'z-pass particle out of link list range - deleted'
        write(*,*) rank,xv(:,pp),pp,i,j,k,hoc_nc_l,hoc_nc_h
        call mpi_abort(mpi_comm_world,ierr,ierr)
     endif
     ll(pp)=hoc(i,j,k)
     hoc(i,j,k)=pp
     pp=pp+1
  enddo
  
  call mpi_reduce(np_local,np_max,1,mpi_integer, &
       mpi_max,0,mpi_comm_world,ierr)
  min_den_buf = real(np_max)*density_buffer/real(max_np)
!!$  if(rank==0) write(*,*) '*************** Density_buffer Analysis *************'
!!$  if(rank==0) write(*,*) '*** max np allowed             = ' , max_np, '    ***'
!!$  if(rank==0) write(*,*) '*** max np_local (with ghosts) = ' , np_max, '    ***'
!!$  if(rank==0) write(*,*) '*** min density_buffer allowed = ', min_den_buf, ' ***'
!!$  
!!$  if(rank==0) write(*,*) '*************** SendRecv Analysis *******************'
!!$  if(rank==0) write(*,*) '*** max np_buf                 = ' , np_buf_max, '    ***'
!!$  if(rank==0) write(*,*) '*** max allowed                = ' , max_buf/6 , '    ***'
!!$  if(rank==0) write(*,*) '*****************************************************'

end subroutine particle_pass

subroutine check_buf_num(npl0,astr)
  implicit none
  include 'mpif.h'
# include "cubepm.fh"
  character(len=*), intent(in) :: astr
  integer(4), intent(in) :: npl0
  integer(4) :: ikill, ikill_loc
  ikill_loc = 0

  !! Check to see if we need to checkpoint kill
  if (np_buf*6 > max_buf) then
     write(*,*) 'rank:',rank,'not enough buffer space in pass',np_buf*6,max_buf
     ikill_loc = 1 
  endif
  call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  if (ikill > 0) then
     if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting: "//astr! (+x pass) ... "
     !! Reset np_local to its starting point so that we don't write duplicates. 
     np_local = npl0 !np_local0
     call checkpoint(.true.)
     call mpi_barrier(mpi_comm_world,ierr)
     call mpi_abort(mpi_comm_world,ierr,ierr)
  endif

end subroutine check_buf_num

subroutine check_par_num(npl,npb,npl0,astr)
  implicit none
  include 'mpif.h'
# include "cubepm.fh"
  character(len=*), intent(in) :: astr
  integer(4), intent(in) :: npl,npb,npl0
  integer(4) :: ikill, ikill_loc,np_check
  ikill_loc = 0

  np_check=npl+npb
  if (np_check > max_np) then
     write(*,*) 'rank:',rank,'exceeded max_np in pass',np_check,max_np
     write(*,*) 'rank:',rank,'particles in node/buffer',npl,npb
     write(*,*) 'rank:',rank,'dm particles in buffer',count(send_buf_PID(:npb).eq.pid_dm)
     ikill_loc = 1
  endif
  call mpi_allreduce(ikill_loc, ikill, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  if (ikill > 0) then 
     if (rank == 0) write(*,*) "Calling checkpoint_kill and then aborting: "//astr
     !! Reset np_local to its starting point so that we don't write duplicates. 
     np_local = npl0
     call checkpoint(.true.)
     call mpi_barrier(mpi_comm_world,ierr)
     call mpi_abort(mpi_comm_world,ierr,ierr)
  endif
end subroutine check_par_num
