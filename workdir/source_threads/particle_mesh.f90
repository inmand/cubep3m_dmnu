!! main particle mesh subroutine
subroutine particle_mesh
  implicit none
# include "cubepm.fh"

  integer(4) :: i, j, k, cur_tile, thread
  integer(4), dimension(3) :: tile
  real(4), dimension(3) :: offset
  integer(4) :: k0, i3
  integer(4), dimension(3) :: cic_l, cic_h 
  integer(4), dimension(3, cores) :: cic_fine_l, cic_fine_h
  real(4) :: force_mag
  real(4) :: f_force_max_node
  real(4) :: pp_force_max_node
  real(4) :: pp_ext_force_max_node
  integer(4) :: omp_get_thread_num
  external omp_get_thread_num
  real(4), dimension(nested_threads) :: pp_ext_sum

  !! Variables that appear in the nested parts (give them all _n suffix) 
  integer(4) :: i_n, j_n, k_n, thread_n, n_pairs_n 
  real(4), dimension(3) :: x_n, offset_n, dx1_n, dx2_n
  integer(4) :: pp_n, pp1_n, pp2_n
  integer(4), dimension(3) :: i1_n, i2_n
  integer(4) :: ii_n, jj_n, kk_n, im_n, jm_n, km_n, ip_n, jp_n, kp_n
  integer(4) :: ip_min_n, ip_max_n, jp_min_n, jp_max_n, kp_min_n, kp_max_n
  integer(4) :: ipl_n(mesh_scale, mesh_scale, mesh_scale)
  real(4), dimension(3) :: sep_n, force_pp_n
  real(4) :: pp_force_mag_n, force_mag_n, rmag_n, dVc_n, pp_ext_sum_k_n
  real(4) :: fpp1_n, fpp2_n !! Only used for neutrino simulations

# if VERBOSITY>0
  if (rank.eq.0) write(*,*) 'particle mesh'
# endif

  !! start of particle mesh.  All particles are within (1:nc_node_dim]
  call update_position

  !! particles must not have advanced past hoc_nc_l:hoc_nc_h
  call link_list
  call mpi_barrier(mpi_comm_world, ierr)
  call particle_pass
  call mpi_barrier(mpi_comm_world, ierr)

# if VERBOSITY>0
  if (rank.eq.0) write(*,*) ':: fine + pp forces'
# endif

  !$omp  parallel num_threads(cores) default(shared) &
  !$omp& private(cur_tile, i, j, k, k0, tile, thread, i3, cic_l, cic_h, offset, force_mag, pp_ext_sum)
  thread = 1
  thread = omp_get_thread_num() + 1
  f_mesh_mass(thread)=0.0
  f_force_max(thread, :)=0.0
  pp_force_max(thread, :)=0.0
  pp_ext_force_max(thread)=0.
  !$omp do schedule(dynamic) 
  do cur_tile=1,tiles_node
     tile(3) = (cur_tile-1) / (tiles_node_dim * tiles_node_dim)
     j = cur_tile - tile(3) * tiles_node_dim * tiles_node_dim
     tile(2) = (j-1) /  tiles_node_dim
     j = j - tile(2) * tiles_node_dim
     tile(1) = j - 1
     
     !fine_mass
     !! normalize fine mesh density
     rho_f(:,:,:,thread)= f_unclustered

     !! calculate coarse mesh offsets
     cic_l(:) = nc_tile_dim * tile(:) + 2 - nc_buf
     cic_h(:) = nc_tile_dim * (tile(:) + 1) + nc_buf - 1

     !! calculate fine mesh density for tile
     do k0 = 0, mesh_scale-1 
        !$omp parallel num_threads(nested_threads) default(shared) private(i_n, j_n, k_n, pp_n, offset_n, x_n, i1_n) 
        !$omp do schedule(dynamic)
        do k_n = cic_l(3) + k0, cic_h(3), mesh_scale
           do j_n = cic_l(2), cic_h(2)
              do i_n = cic_l(1), cic_h(1)
                 pp_n = hoc(i_n, j_n, k_n)
                 offset_n(:)= - tile(:) * nf_physical_tile_dim + nf_buf

                 do; if (pp_n == 0) exit
                    x_n = xv(1:3, pp_n) + offset_n(:)
                    i1_n(:) = floor(x_n(:)) + 1
                    rho_f(i1_n(1),i1_n(2),i1_n(3),thread) = rho_f(i1_n(1),i1_n(2),i1_n(3),thread)+mass_p*mass_p_nudm_fac(pidmap(PID(pp_n)))
                    pp_n = ll(pp_n)
                 enddo !! pp_n

              enddo !! i_n
           enddo !! j_n
        enddo !! k_n
        !$omp end do
        !$omp end parallel
     enddo !! k0

     !fine force
     !! transform and calculate fine mesh force 
     call cubepm_fftw2('f',thread)
     cmplx_rho_f(:,:,:,thread)=rho_f(:,:,:,thread)

     do i3 = 1, 3
        !$omp parallel num_threads(nested_threads) default(shared) private(k_n, j_n, i_n, ii_n, im_n)
        !$omp do
        do k_n = 1, nf_tile
           do j_n = 1, nf_tile
              do i_n = 1, nf_tile/2+1
                 ii_n = 2*i_n
                 im_n = ii_n - 1
                 rho_f(im_n, j_n, k_n, thread) = -cmplx_rho_f(ii_n, j_n, k_n, thread) * kern_f(i3, i_n, j_n, k_n)
                 rho_f(ii_n, j_n, k_n, thread) =  cmplx_rho_f(im_n, j_n, k_n, thread) * kern_f(i3, i_n, j_n, k_n)
              enddo
           enddo
        enddo
        !$omp end do 
        !$omp end parallel

        call cubepm_fftw2('b',thread)

        !$omp parallel num_threads(nested_threads) default(shared) private(i_n, j_n, k_n)
        !$omp do
        do k_n = nf_buf-1, nf_tile-nf_buf+1
           do j_n = nf_buf-1, nf_tile-nf_buf+1
              do i_n = nf_buf-1, nf_tile-nf_buf+1
                 force_f(i3, i_n, j_n, k_n, thread) = rho_f(i_n, j_n, k_n, thread)
              enddo
           enddo
        enddo
        !$omp end do 
        !$omp end parallel

     enddo
     
     !! calculate max fine mesh force   
     force_mag=0.0
     do k=nf_buf-1,nf_tile-nf_buf+1
        do j=nf_buf-1,nf_tile-nf_buf+1
           do i=nf_buf-1,nf_tile-nf_buf+1
              force_mag=sqrt(force_f(1,i,j,k,thread)**2 &
                   +force_f(2,i,j,k,thread)**2+force_f(3,i,j,k,thread)**2)
              if (force_mag > maxval(f_force_max(thread, :))) f_force_max(thread, 1)=force_mag
           enddo
        enddo
     enddo
     
     !fine velocity
     offset(:) = real(nf_buf) - tile(:) * nf_physical_tile_dim
     do k0 = 0, mesh_scale-1 
        !$omp  parallel num_threads(nested_threads) default(shared) &
        !$omp& private(i_n, j_n, k_n, pp_n, x_n, i1_n, im_n, jm_n, km_n, i2_n, force_mag_n, dx1_n, dx2_n, dVc_n, &
        !$omp&         ipl_n, ip_n, jp_n, kp_n, pp1_n, pp2_n, sep_n, rmag_n, force_pp_n, pp_force_mag_n, thread_n, &
        !$omp&         fpp1_n, fpp2_n)
        thread_n = 1
        thread_n = omp_get_thread_num() + 1
        !$omp do
        do k_n = k0 + tile(3) * nc_tile_dim + 1, (tile(3) + 1) * nc_tile_dim, mesh_scale 
           do j_n = tile(2) * nc_tile_dim + 1, (tile(2) + 1) * nc_tile_dim
              do i_n = tile(1) * nc_tile_dim + 1, (tile(1) + 1) * nc_tile_dim
                 pp_n = hoc(i_n, j_n, k_n)
                 ipl_n = 0

                 do; if (pp_n == 0) exit
                    x_n = xv(1:3, pp_n) + offset(:)
                    i1_n(:) = floor(x_n(:)) + 1

                    do im_n = 1, 3
                       i2_n(im_n) = mod(i1_n(im_n)-1,mesh_scale) + 1
                    enddo
                    ipl_n(i2_n(1),i2_n(2),i2_n(3)) = ipl_n(i2_n(1),i2_n(2),i2_n(3))+1
                    if (ipl_n(i2_n(1),i2_n(2),i2_n(3))>max_llf) then
                       print *,'exceeded max_llf',max_llf,i1_n,i2_n,ipl_n
                       stop
                    endif
                    llf(ipl_n(i2_n(1),i2_n(2),i2_n(3)),i2_n(1),i2_n(2),i2_n(3),thread,thread_n)=pp_n

                    xv(4:6,pp_n)=xv(4:6,pp_n)+force_f(1:3,i1_n(1),i1_n(2),i1_n(3),thread) * a_mid * G * dt
                    
                    pp_n = ll(pp_n)
                 enddo
                 
                 do km_n = 1, mesh_scale
                    do jm_n = 1, mesh_scale
                       do im_n = 1, mesh_scale
                          
                          pp_force_accum(:, :ipl_n(im_n,jm_n,km_n), thread, thread_n) = 0.
                          do ip_n = 1, ipl_n(im_n,jm_n,km_n) - 1
                             pp1_n = llf(ip_n, im_n, jm_n, km_n, thread, thread_n)
                             fpp1_n = mass_p_nudm_fac(pidmap(PID(pp1_n))) 

                             do jp_n = ip_n+1, ipl_n(im_n,jm_n,km_n)
                                pp2_n = llf(jp_n, im_n, jm_n, km_n, thread, thread_n)
                                fpp2_n = mass_p_nudm_fac(pidmap(PID(pp2_n))) 

#                               ifdef MASSLESS_PARTICLES
                                !Don't include massless-massless interactions
                                if (fpp1_n.eq.0 .and. fpp2_n.eq.0) cycle
#                               endif                              

                                sep_n = xv(:3,pp1_n) - xv(:3,pp2_n)                      
                                rmag_n = sqrt(sep_n(1)*sep_n(1) + sep_n(2)*sep_n(2) + sep_n(3)*sep_n(3))
                                if (rmag_n > rsoft) then

                                   force_pp_n = mass_p*(sep_n/(rmag_n*pp_bias)**3)  !mass_p divides out below
                                   
                                   pp_force_accum(:, ip_n, thread, thread_n) = pp_force_accum(:, ip_n, thread, thread_n) - force_pp_n*fpp2_n
                                   pp_force_accum(:, jp_n, thread, thread_n) = pp_force_accum(:, jp_n, thread, thread_n) + force_pp_n*fpp1_n

                                   xv(4:,pp1_n) = xv(4:,pp1_n) - force_pp_n*a_mid*G*dt*fpp2_n
                                   xv(4:,pp2_n) = xv(4:,pp2_n) + force_pp_n*a_mid*G*dt*fpp1_n

                                endif
                             enddo
                          enddo

                          do ip_n=1, ipl_n(im_n,jm_n,km_n)
                             pp_force_mag_n = sqrt(&
                                  &pp_force_accum(1,ip_n,thread,thread_n)**2 + &
                                  &pp_force_accum(2,ip_n,thread,thread_n)**2 + &
                                  pp_force_accum(3,ip_n,thread,thread_n)**2 )
                             if (pp_force_mag_n > pp_force_max(thread, thread_n)) pp_force_max(thread, thread_n)=pp_force_mag_n
                          enddo

                       enddo !im_n
                    enddo !jm_n
                 enddo !km_n

              enddo !i_n
           enddo !j_n
        enddo !k_n
        !$omp end do
        !$omp end parallel
     enddo !k0

#ifdef PP_EXT

     !fine link list
     hoc_fine(:,:,:,thread)=0
 
     ! Get physical coordinate of boundaries of tiles in fine mesh units
     cic_fine_l(:,thread) = tile(:)*nf_physical_tile_dim + 1
     cic_fine_h(:,thread) = (tile(:) + 1)*nf_physical_tile_dim
     
     ! Include pp_range
     cic_fine_l(:,thread) = cic_fine_l(:,thread) - pp_range 
     cic_fine_h(:,thread) = cic_fine_h(:,thread) + pp_range 
      
     !$omp parallel num_threads(nested_threads) default(shared) private(i_n, j_n, k_n, ii_n, jj_n, kk_n, pp_n)
     !$omp do schedule(dynamic)
     do k_n = cic_l(3), cic_h(3)
        do j_n = cic_l(2), cic_h(2)
           do i_n = cic_l(1), cic_h(1)
              
              pp_n = hoc(i_n, j_n, k_n)
              
              do; if (pp_n == 0) exit

                 ii_n = floor(xv(1,pp_n))+1
                 jj_n = floor(xv(2,pp_n))+1
                 kk_n = floor(xv(3,pp_n))+1
                 
                 if (ii_n >= cic_fine_l(1,thread) .and. ii_n <= cic_fine_h(1,thread) .and. &
                      jj_n >= cic_fine_l(2,thread) .and. jj_n <= cic_fine_h(2,thread) .and. &
                      kk_n >= cic_fine_l(3,thread) .and. kk_n <= cic_fine_h(3,thread)) then
                       
                    ll_fine(pp_n,thread) = hoc_fine(ii_n-cic_fine_l(1,thread)+1, jj_n-cic_fine_l(2,thread)+1, kk_n-cic_fine_l(3,thread)+1,thread)
                    hoc_fine(ii_n-cic_fine_l(1,thread)+1, jj_n-cic_fine_l(2,thread)+1, kk_n-cic_fine_l(3,thread)+1,thread)=pp_n
                    
                 endif

                 pp_n = ll(pp_n)
              enddo
              
           enddo
        enddo
     enddo
     !$omp end do
     !$omp end parallel

     !$omp parallel num_threads(nested_threads) default(shared) private(k_n)
     !$omp do
     do k_n = 1, max_np 
        pp_ext_force_accum(:,k_n,thread) = 0 
     enddo
     !$omp end do
     !$omp end parallel

     if(pp_range .ne. 0)then

        do k0 = 0, pp_range
           !$omp  parallel num_threads(nested_threads) default(shared) & 
           !$omp& private(i_n, j_n, k_n, pp1_n, kp_min_n, kp_max_n, kp_n, jp_min_n, jp_max_n, fpp1_n, fpp2_n, &
           !$omp&         jp_n, ip_min_n, ip_max_n, ip_n, pp2_n, n_pairs_n, sep_n, rmag_n, force_pp_n) 
           !$omp do
           do k_n=1+k0,nf_physical_tile_dim+pp_range,pp_range+1 ! We never loop towards smaller z
              do j_n=1,nf_physical_tile_dim+2*pp_range
                 do i_n=1,nf_physical_tile_dim+2*pp_range
                    
                    pp1_n = hoc_fine(i_n, j_n, k_n, thread) 
                    if(pp1_n == 0) cycle 
 
                    kp_min_n = k_n
                    kp_max_n = k_n + pp_range
                       
                    do kp_n = kp_min_n, kp_max_n

                       if(kp_n == k_n) then
                          jp_min_n = j_n
                       else
                          jp_min_n = j_n - pp_range
                          if(jp_min_n <= 0) jp_min_n = 1
                       endif
                       jp_max_n = j_n + pp_range
                       if(jp_max_n > nf_physical_tile_dim+2*pp_range) jp_max_n = nf_physical_tile_dim+2*pp_range
                       
                       do jp_n = jp_min_n, jp_max_n 

                          if((kp_n == k_n .and. jp_n == j_n))then
                             ip_min_n = i_n+1
                          else
                             ip_min_n = i_n - pp_range
                             if(ip_min_n <= 0) ip_min_n = 1
                          endif
                          ip_max_n = i_n + pp_range
                          if(ip_max_n > nf_physical_tile_dim+2*pp_range) ip_max_n = nf_physical_tile_dim+2*pp_range
                        
                          do ip_n = ip_min_n, ip_max_n 
                           
                             pp2_n = hoc_fine(ip_n, jp_n, kp_n, thread) 
                             if(pp2_n == 0) cycle                                                     

                             n_pairs_n = 0
                                
                             do; if(pp1_n == 0)exit
                                fpp1_n = mass_p_nudm_fac(pidmap(PID(pp1_n)))
                                
                                do ; if(pp2_n == 0)exit
                                   fpp2_n = mass_p_nudm_fac(pidmap(PID(pp2_n)))

#                                  ifdef MASSLESS_PARTICLES
                                   !Don't include massless-massless interactions
                                   if (fpp1_n.eq.0 .and. fpp2_n.eq.0) then
                                      pp2_n = ll_fine(pp2_n, thread)
                                      cycle
                                   end if
#                                  endif
                                      
                                   n_pairs_n = n_pairs_n + 1
                           
                                   sep_n = xv(:3,pp1_n) - xv(:3,pp2_n)
                                   rmag_n = sqrt(sep_n(1)*sep_n(1) + sep_n(2)*sep_n(2) + sep_n(3)*sep_n(3))

                                   if (rmag_n > rsoft) then
                                      
                                      if(rmag_n>real(nf_cutoff)+sqrt(3.0))then
                                         force_pp_n = mass_p*(sep_n/(rmag_n*pp_bias)**3)
                                      else
                                         force_pp_n = mass_p*(sep_n/(rmag_n*pp_bias)**3)*(1 - (7.0/4.0)*(rmag_n*pp_bias/(nf_cutoff))**3 + &
                                              (3.0/4.0)*(rmag_n*pp_bias/(nf_cutoff))**5)  !mass_p divides out below
                                      endif

                                      pp_ext_force_accum(:, pp1_n, thread) = pp_ext_force_accum(:, pp1_n, thread) - force_pp_n*fpp2_n
                                      pp_ext_force_accum(:, pp2_n, thread) = pp_ext_force_accum(:, pp2_n, thread) + force_pp_n*fpp1_n
                              
                                      if (pp_ext_force_flag) then
                                 
                                         ! Update only particles in physical space
                                         if((pp_range<i_n).and.(i_n<=nf_physical_tile_dim+pp_range) .and.&
                                              (pp_range<j_n).and.(j_n<=nf_physical_tile_dim+pp_range).and.&
                                              (pp_range<k_n).and.(k_n<=nf_physical_tile_dim+pp_range)) then 
                                               
                                            xv(4:,pp1_n) = xv(4:,pp1_n) - force_pp_n*a_mid*G*dt*fpp2_n

                                         endif
                                            
                                         if((pp_range<ip_n).and.(ip_n<=nf_physical_tile_dim+pp_range) .and.&
                                              (pp_range<jp_n).and.(jp_n<=nf_physical_tile_dim+pp_range).and.&
                                              (pp_range<kp_n).and.(kp_n<=nf_physical_tile_dim+pp_range)) then 

                                            xv(4:,pp2_n)=xv(4:,pp2_n) + force_pp_n*a_mid*G*dt*fpp1_n

                                         endif
                             
                                      endif !! pp_ext_force_flag
                                      
                                   endif !! (rmag_n > rsoft)
                           
                                   pp2_n = ll_fine(pp2_n, thread)                           
                           
                                enddo !! pp2_n

                                pp2_n = hoc_fine(ip_n, jp_n, kp_n, thread)
                                pp1_n = ll_fine(pp1_n, thread)

                             enddo !! pp1_n

                             ! Restore pp1_n value for the next pp2_n iteration
                             pp1_n = hoc_fine(i_n, j_n, k_n, thread)
                        
                          enddo !! ip_n                    
                       enddo !! jp_n
                    enddo !! kp_n
                  
                 enddo !! i_n
              enddo !! j_n
           enddo !! k_n
           !$omp end do
           !$omp end parallel
        enddo !! k0

     endif !! (pp_range .ne. 0)

     pp_ext_sum(:) = 0.
     !$omp parallel num_threads(nested_threads) default(shared) private(k_n, pp_ext_sum_k_n, thread_n)
     thread_n = 1
     thread_n = omp_get_thread_num() + 1
     !$omp do
     do k_n = 1, max_np
        pp_ext_sum_k_n = pp_ext_force_accum(1,k_n,thread)**2 + pp_ext_force_accum(2,k_n,thread)**2 + pp_ext_force_accum(3,k_n,thread)**2
        if (pp_ext_sum_k_n > pp_ext_sum(thread_n)) pp_ext_sum(thread_n) = pp_ext_sum_k_n
     enddo
     !$omp end do
     !$omp end parallel

     pp_ext_force_max(thread) = max(pp_ext_force_max(thread), sqrt(maxval(pp_ext_sum(:))))

#endif
      
  enddo
  !$omp end do
  !$omp end parallel
  
  call mpi_barrier(mpi_comm_world, ierr) 

  !! calculate maximum dt from fine mesh force
  f_force_max_node=maxval(f_force_max)
  call mpi_reduce(f_force_max_node,dt_f_acc,1,mpi_real,mpi_max,0, &
       mpi_comm_world,ierr)
  if (rank == 0) then
     dt_f_acc=1.0/sqrt(max(0.0001,dt_f_acc)*a_mid*G)
  endif
  call mpi_bcast(dt_f_acc,1,mpi_real,0,mpi_comm_world,ierr)

  !! calculate maximum dt from particle-particle force
  pp_force_max_node=maxval(pp_force_max)
  call mpi_reduce(pp_force_max_node,dt_pp_acc,1,mpi_real,mpi_max,0, &
       mpi_comm_world,ierr)
  if (rank == 0) then
     dt_pp_acc=sqrt(dt_pp_scale*rsoft)/max(sqrt(dt_pp_acc*a_mid*G),1e-3)
  endif
  call mpi_bcast(dt_pp_acc,1,mpi_real,0,mpi_comm_world,ierr)

#ifdef PP_EXT
  pp_ext_force_max_node=maxval(pp_ext_force_max)
  call mpi_reduce(pp_ext_force_max_node,dt_pp_ext_acc,1,mpi_real,mpi_max,0, &
       mpi_comm_world,ierr)
  if (rank == 0) then
     dt_pp_ext_acc=sqrt(dt_pp_scale*rsoft)/max(sqrt(dt_pp_ext_acc*a_mid*G),1e-3)
  endif
  call mpi_bcast(dt_pp_ext_acc,1,mpi_real,0,mpi_comm_world,ierr)
#endif
     
  call coarse_mesh

#ifdef TIDAL_FIELD
  call compute_tidal_field
#endif

  !! delete all particles outside (1:nc_node_dim]
  call delete_particles

  !! Prevent drift from going too far
  call move_grid_back(.false.)

end subroutine particle_mesh
