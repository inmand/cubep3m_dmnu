program init_pbh
  implicit none
#include "../../parameters"

  real :: r
  integer :: p,d,stat

  character(len=1000) :: str

  !Seed stuff
  integer(4) seedsize
  integer(4), allocatable :: iseed(:)
  
  call random_seed(size=seedsize)
  seedsize=max(seedsize,12)
  allocate(iseed(seedsize))
    open(11,file=scratch_path//'node0/seed0.init',status='old',access='stream')
  read(11) iseed
  close(11)
  call random_seed(put=iseed)

  write(str,'(f10.3)') z_i
  str=scratch_path//'/node0/'//trim(adjustl(str))//'xv0_nu.dat'
  write(*,*) 'Writing to file: '//trim(adjustl(str))
  open(unit=11,file=trim(adjustl(str)),status='replace',iostat=stat,access='stream')

  write(11) npbh,1.0/(1.0+z_i),0.0,-3./sqrt(1.0/(1.0+z_i)),0,1000.,1000.,1000.,1,1,1,1. !Header garbage

  do p=1,Npbh

     !Positions
     do d=1,3
        call random_number(r)
        write(11) r*nc
     end do

     !Velocities
     do d=1,3
        write(11) 0.0
     end do

  end do

  close(11)

end program init_pbh
