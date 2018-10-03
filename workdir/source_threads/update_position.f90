!! update particle positions
subroutine update_position
  use omp_lib
  implicit none
# include "cubepm.fh"

  integer(4) :: i,j

  real(4), dimension(3) :: offset

# if VERBOSITY>0
  if (rank.eq.0) write(*,*) ':: update position'
# endif

  if (rank==0) then
       
     call random_number(offset)

     if (a.gt.a_i_nu) then
        !offset=(offset-0.5)*mesh_scale ! no shake offset
        offset=(offset-0.5)*mesh_scale*offset_boost ! no shake offset
     else
        offset=(offset-0.5)*mesh_scale*offset_boost  - shake_offset
     end if

     shake_offset=shake_offset+offset
  
  endif

  call mpi_bcast(offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(shake_offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  
  !$omp parallel do default(shared) private(i)
  do i=1,np_local
     xv(1:3,i)=xv(1:3,i)+xv(4:6,i)*0.5*(dt + dt_old)+offset(:)
  enddo
  !$omp end parallel do

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
