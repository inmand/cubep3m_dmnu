!-*- mode: f90-*-
module Class_LinkedList
  implicit none
  private
  
  public :: LinkedList
  public :: setup_linked_list, destroy_linked_list
        
  type LinkedList
     integer, dimension(:,:,:), allocatable :: hoc
     integer, dimension(:), allocatable :: ll
  end type LinkedList

contains

  subroutine setup_linked_list(tll, x)
    implicit none
    type(LinkedList), intent(inout) :: tll
    real, dimension(:,:), intent(in) :: x
    integer :: p, n, maxc, minc
    integer, dimension(3) :: nx
    
    n=size(x,dim=2)
    if (size(x,dim=1).ne.3) then
       write(*,*) 'ERROR in setup_linked_list: x dimension not 3',size(x,dim=1)
       ERROR STOP
    end if
    minc=1+floor(minval(x))
    maxc=1+floor(maxval(x))

    allocate(tll%ll(n), tll%hoc(minc:maxc,minc:maxc,minc:maxc))

    !$omp workshare
    tll%ll = 0
    tll%hoc = 0
    !$omp end workshare

    do p=1,n
       nx(:) = 1+floor(x(:,p))
       tll%ll(p) = tll%hoc(nx(1),nx(2),nx(3))
       tll%hoc(nx(1),nx(2),nx(3)) = p
    end do

  end subroutine setup_linked_list

  subroutine destroy_linked_list(tll)
    implicit none 
    type(LinkedList), intent(inout) :: tll
    deallocate(tll%ll, tll%hoc)
  end subroutine destroy_linked_list

end module Class_LinkedList
