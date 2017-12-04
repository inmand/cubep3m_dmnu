!! calculate coarse mesh density
  subroutine coarse_mass
   implicit none

    include 'mpif.h'
#   include "cubepm.fh"

   integer(4) :: i,j,k,k0,pp,ii,jj,kk
   integer(4), dimension(3) :: i1,i2
   real(4), dimension(3) :: x,dx1,dx2

#ifdef COARSEPROJ
    character (len=max_path) :: ofile
    character (len=6) :: step_string
    character (len=6) :: rank_s
    integer :: fstat
    real(8) :: cpsum_local, cpsum
    integer, parameter :: ijstart = -1
    integer, parameter :: ijstop  = nc_node_dim + 2
#else
    integer, parameter :: ijstart = 0
    integer, parameter :: ijstop  = nc_node_dim + 1 
#endif

   call system_clock(count=count_i)

#ifndef NEUTRINOS
   rho_c = ratio_omega_nu2m
#else
   rho_c= 0.0 !- mass_p * (mesh_scale / 2)**3  
#endif

#ifdef COARSEPROJ
    !! Determine if want to do projection here
    doCoarseProj = .false.
    if (mod(nts,writeCoarseProjEverySteps) == 0 .and. a >= writeCoarseProjAboveA) then 
        doCoarseProj = .true.
        !! Initialize grids
        crhoproj = 0. ; crhoprojsum = 0.        
#ifdef NEUTRINOS
        crhoproj_nu = 0. ; crhoprojsum_nu = 0.
#endif
        !! Completely remove z offset while binning x, y to nearest coarse mesh cell
        call mpi_bcast(shake_offset, 3, mpi_real, 0, mpi_comm_world, ierr)
        soffcproj(1) = (shake_offset(1) - mesh_scale*int(shake_offset(1)/mesh_scale))/mesh_scale 
        soffcproj(2) = (shake_offset(2) - mesh_scale*int(shake_offset(2)/mesh_scale))/mesh_scale
        soffcproj(3) = shake_offset(3)/mesh_scale

        if (rank == 2) write(*,*) "shake_offset = ", shake_offset
        if (rank == 2) write(*,*) "soffcproj = ", soffcproj

    endif
#endif 

    do k0 = 0, mesh_scale-1 
        !$omp parallel do schedule(dynamic) default(shared) private(i,j,k,pp)
        do k = k0, nc_node_dim + 1, mesh_scale 
            do j = ijstart, ijstop 
                do i = ijstart, ijstop 
                    pp=hoc(i,j,k)
                    if (i <= 1 .or. i >= nc_node_dim .or. &
                        j <= 1 .or. j >= nc_node_dim .or. &
                        k <= 1 .or. k >= nc_node_dim) then
                      call coarse_cic_mass_boundry(pp)
                    else
                        call coarse_cic_mass(pp)
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    enddo

#ifdef COARSEPROJ
    if (doCoarseProj) then
        !! Sum up grid as consistency check
        cpsum_local = sum(crhoproj(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim))
        call mpi_reduce(cpsum_local, cpsum, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr) 
        if (rank == 0) write(*,*) "Dark matter coarsened projection sum = ", cpsum
#ifdef NEUTRINOS 
        cpsum_local = sum(crhoproj_nu(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim))
        call mpi_reduce(cpsum_local, cpsum, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
        if (rank == 0) write(*,*) "Neutrino coarsened projection sum = ", cpsum
#endif
        !! Write data using only those on the top xy slab of the box
        if (cart_coords(1) == 0) then
            !! Rank and time step strings
            write(rank_s,'(i6)') rank
            rank_s = adjustl(rank_s)
            write(step_string, "(i6.6)") nts
            !! Sum up grids into projections
            do k = nc_coarseproj_start, nc_coarseproj_stop 
                do j = 1, nc_node_dim
                    do i = 1, nc_node_dim
                        crhoprojsum(i, j) = crhoprojsum(i,j) + crhoproj(i,j,k) 
#ifdef NEUTRINOS
                        crhoprojsum_nu(i,j) = crhoprojsum_nu(i,j) + crhoproj_nu(i,j,k)
#endif
                    enddo
                enddo
            enddo
            !! Write dark matter projection
            ofile=output_path//'coarseproj/node'//rank_s(1:len_trim(rank_s))//'/cproj'//rank_s(1:len_trim(rank_s))//'_'//step_string//'.dat'
            open(unit=12, file=ofile, status="replace", iostat=fstat, access="stream")
            if (fstat /= 0) then
                write(*,*) 'error opening coarse projection file for write'
                write(*,*) 'rank',rank,'file:',ofile
                call mpi_abort(mpi_comm_world,ierr,ierr)
            endif
            write(12) shake_offset
            write(12) soffcproj
            write(12) crhoprojsum 
            close(12)
#ifdef NEUTRINOS
            !! Write neutrino projection
            ofile=output_path//'coarseproj/node'//rank_s(1:len_trim(rank_s))//'/cproj'//rank_s(1:len_trim(rank_s))//'_'//step_string//'_nu.dat'
            open(unit=13, file=ofile, status="replace", iostat=fstat, access="stream")
            if (fstat /= 0) then
                write(*,*) 'error opening coarse projection file for write'
                write(*,*) 'rank',rank,'file:',ofile
                call mpi_abort(mpi_comm_world,ierr,ierr)
            endif
            write(13) shake_offset
            write(13) soffcproj
            write(13) crhoprojsum_nu
            close(13)
#endif
        endif !! cart_coord
        if (rank == 0) write(*,*) "Coarse Projection written"
        call mpi_barrier(mpi_comm_world, ierr)
    endif !! doCoarseProj
#endif

  call system_clock(count=count_f,count_rate=count_r)
#ifdef MPI_TIME
  call mpi_time_analyze('cm  mass',real(count_f-count_i)/real(count_r),rank,nodes)
#else
  if (rank==0) write(*,*) 'coarse mass finished',real(count_f-count_i)/real(count_r)
#endif


  end subroutine coarse_mass
