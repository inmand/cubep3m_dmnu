!! update particle positions
subroutine update_position
  use omp_lib
  implicit none
# ifdef DISP_MESH 
  include 'mpif.h'
# endif
# include "cubepm.fh"

  integer(4) :: i,j

# ifdef DISP_MESH 
  real(4), dimension(3) :: offset

  if (rank==0) then
       
     call random_number(offset)

# ifdef NEUTRINOS
     offset=(offset-0.5)*mesh_scale ! no shake offset
# else
     offset=(offset-0.5)*mesh_scale*4.0  - shake_offset
# endif
     shake_offset=shake_offset+offset
     print*,'current shake offset:',shake_offset
  
  endif

  call mpi_bcast(offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(shake_offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  
# endif

  call system_clock(count=count_i)
    
  !$omp parallel do default(shared) private(i)
  do i=1,np_local
#ifdef DISP_MESH 
     xv(1:3,i)=xv(1:3,i)+xv(4:6,i)*0.5*(dt + dt_old)+offset(:)
#else
     xv(1:3,i)=xv(1:3,i)+xv(4:6,i)*0.5*(dt + dt_old)
#endif
  enddo
  !$omp end parallel do

  call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
  call mpi_time_analyze('pos updt',real(count_f-count_i)/real(count_r),rank,nodes)
#else
  if (rank==0) write(*,*) 'position update finished',real(count_f-count_i)/real(count_r)
#endif

end subroutine update_position


!!$       !This will always use the same random number at each time step. 
!!$       !It surely introduces a bias, but is good for testing code.        
!!$       call random_seed
!!$       call random_seed(size=seedsize)
!!$
!!$       allocate(iseed(seedsize))
!!$       !allocate(old(seedsize))
!!$
!!$
!!$       seedfile = ic_path//'node0/seed0.init'  !or any other seed files available
!!$       open(11,file=seedfile)
!!$       write(*,*) 'opened ',seedfile
!!$       do i = 1,seedsize 
!!$          read(11,*) j,iseed(i)
!!$       enddo
!!$       close(11)
!!$
!!$       call random_seed(put=iseed(1:seedsize))
