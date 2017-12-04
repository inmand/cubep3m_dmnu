!! generate linked list for particles based on their positions
!! within the coarse mesh
subroutine link_list
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  integer(4) :: i,j,k,pp,pc
  integer(4) :: omp_get_thread_num,omp_get_num_threads
  external omp_get_thread_num,omp_get_num_threads

  call system_clock(count=count_i)

  hoc(:,:,:)=0
  np_buf=0

  pp=1
  do
     if (pp > np_local) exit
     i=floor(xv(1,pp)/mesh_scale)+1
     j=floor(xv(2,pp)/mesh_scale)+1
     k=floor(xv(3,pp)/mesh_scale)+1
     if (i < hoc_nc_l .or. i > hoc_nc_h .or. &
          j < hoc_nc_l .or. j > hoc_nc_h .or. &
          k < hoc_nc_l .or. k > hoc_nc_h) then
        print*, 'link_list: particle moved out of buffer!',xv(:,pp)
        print*, 'check timestep & update_position!'
        call mpi_abort(mpi_comm_world,ierr,ierr)
        cycle
     else
        ll(pp)=hoc(i,j,k)
        hoc(i,j,k)=pp
     endif
     pp=pp+1
  enddo
  
#ifdef DIAG
  if (np_buf.gt.0) write(*,*) rank,'deleted',np_buf,'particles in ll'
#endif

  call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
  call mpi_time_analyze('linklist',real(count_f-count_i)/real(count_r),rank,nodes)
#else
  if (rank==0) write(*,*) 'link list finished',real(count_f-count_i)/real(count_r)
#endif

end subroutine link_list
