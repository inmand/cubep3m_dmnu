#ifdef TIDAL_FIELD
subroutine compute_tidal_field
  use omp_lib
  implicit none
# include "cubepm.fh"

  integer :: n,nn,i,j
  real :: r
  real, dimension(3) :: dx

  integer :: tag,status(MPI_STATUS_SIZE) 

  if (rank.eq.0) write(*,*) 'Computing tidal field'

  !Store scalefactor
  tide_a(nts)=a_mid

  !Set to zero
  tide_bh(:,:,:,nts)=0.

  !Loop over all particles and compute tidal field
  !$omp parallel do default(none) shared(nts,xv,PID,np_local,tide_bh,mass_p,mass_p_nudm_fac) private(n,nn,dx,r,i,j)
  do n=1,np_local

     !Select only bh
     if ( pidmap(PID(n)) .eq. pid_dm ) cycle

     !Select only bh on this node
     if ( minval(xv(1:3,n)).lt.0 .or. maxval(xv(1:3,n)).gt.nf_physical_node_dim ) cycle
        
     !Find all nearby particles
     do nn=1,np_local

        !Don't self count (in case rsoft=0)
        if (n.eq.nn) cycle
           
        dx=xv(1:3,n)-xv(1:3,nn)
        r=sqrt(sum(dx**2))
        if (r.le.tide_rmin .or. r.gt.tide_rmax) cycle
                   
        !Compute Tij
        dx=dx/r 
        do j=1,3
           do i=1,3
              tide_bh(i,j,PID(n),nts) = tide_bh(i,j,PID(n),nts) - &
                   &(1./4./pi)*(1./r**3)*(3.*dx(i)*dx(j)-merge(1,0,i.eq.j))*mass_p*mass_p_nudm_fac(pidmap(PID(nn)))
           end do
        end do

     end do

  end do
  !$omp end parallel do
  
  !Synchronize tidal field
  if (rank.eq.0) write(*,*) 'Synchronizing tidal field'
  do n=1,nodes_dim**3-1
     tag=n
     tide_recv=0.
     if (rank.eq.n) then
        !Send to head node
        call mpi_send(tide_bh(:,:,:,nts),size(tide_bh(:,:,:,nts)),mpi_real,0,tag,mpi_comm_world,ierr)
     else if (rank.eq.0) then
        !Receive on head node
        call mpi_recv(tide_recv,size(tide_recv),mpi_real,n,tag,mpi_comm_world,status,ierr)
        !Transfer to tide_bh array
        do nn=1,n_bh
           if (sum(tide_recv(:,:,nn)).ne.0) then
              !Double check not double counted
              if (sum(tide_bh(:,:,nn,nts)).ne.0) then
                 write(*,*) 'Error: tide_bh double counted bh - rank,n,nn',rank,n,nn
                 write(*,*) '>nn,tide_recv',tide_recv(:,:,nn)
                 write(*,*) '>nn,tide_bh',tide_bh(:,:,nn,nts)
                 call mpi_abort(mpi_comm_world,ierr,ierr)
              end if
              !Add to array
              tide_bh(:,:,nn,nts)=tide_recv(:,:,nn)
           end if
        end do
     end if
     call mpi_barrier(mpi_comm_world,ierr)
  end do
  !Pass to all other nodes
  call mpi_bcast(tide_bh(:,:,:,nts),size(tide_bh(:,:,:,nts)),mpi_real,0,mpi_comm_world,ierr)

  if (rank.eq.0) write(*,*) 'Finished tidal field'

end subroutine compute_tidal_field

subroutine checkpoint_tidal_field
  implicit none
# include "cubepm.fh"

  integer :: n

  open(unit=11,file=output_path//'tide_bh.bin',access='stream',status='replace')
  write(11) tide_bh
  close(11)

  open(unit=11,file=output_path//'tide_a.bin',access='stream',status='replace')
  write(11) tide_a
  close(11)

  open(unit=11,file=output_path//'tide_a.txt',status='replace')
  do n=1,nts
     write(11,*) tide_a(n)
  end do
  close(11)

end subroutine checkpoint_tidal_field

subroutine initialize_tidal_field
  implicit none
# include "cubepm.fh"

  tide_bh=0.; tide_a=0.
  if (restart_ic .or. restart_kill) then

     if (rank.eq.0) then
        open(unit=11,file=output_path//'tide_bh.bin',access='stream',status='replace')
        write(11) tide_bh
        close(11)
        
        open(unit=11,file=output_path//'tide_a.bin',access='stream',status='replace')
        write(11) tide_a
        close(11)
     end if

     call mpi_bcast(tide_bh,size(tide_bh),mpi_real,0,mpi_comm_world,ierr)
     call mpi_bcast(tide_a,size(tide_a),mpi_real,0,mpi_comm_world,ierr)
     
  end if

end subroutine initialize_tidal_field
#endif
