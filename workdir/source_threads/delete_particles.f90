!! delete particles which lie outside of the physical region
subroutine delete_particles
  implicit none
  include 'mpif.h'
# include "cubepm.fh"

  integer(4) :: pp,i

  pp=1
  do
42   if (pp > np_local) exit
     if (xv(1,pp) >= nf_physical_node_dim .or. xv(1,pp) < 0.0 .or. &
          xv(2,pp) >= nf_physical_node_dim .or. xv(2,pp) < 0.0 .or. &
          xv(3,pp) >= nf_physical_node_dim .or. xv(3,pp) < 0.0) then                    
        xv(:,pp)=xv(:,np_local)
        PID(pp)=PID(np_local)
        np_local=np_local-1
        goto 42
     endif
     pp=pp+1
  enddo
        
end subroutine delete_particles
