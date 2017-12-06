!! increment timestep 
subroutine timestep
  implicit none 
  include 'mpif.h'
# include "cubepm.fh"

  real(4) :: ra,da_1,da_2,dt_e,am
  integer(4) :: n

  real(4) :: vmax, vmax_local
  real(4) :: Dx

  nts = nts + 1
  if (nts /= 1) dt_old = dt

  !! Compute maximum timestep allowed by the maximum velocity
  vmax_local = 0.
  do n = 1, np_local
     vmax_local = max(vmax_local, maxval(abs(xv(4:6,n))))
  enddo
  call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)

#ifdef NEUTRINOS
  Dx=real(nf_buf)-0.5*mesh_scale ! 5.5 coarse grids
#else
  Dx = real(nf_buf)-4.*mesh_scale ! 2.0 coarse grids
#endif

  !! fbuf is the fraction of the buffer that the fastest particle is allowed to move. 
  !! As long as the maximum velocity increases no more than zeta = 2/fbuf - 1 compared 
  !! to the previous time step then this method will work.
  dt_vmax = fbuf * Dx / vmax

  if (rank == 0) write(*,*) 'vmax and maximum timestep from vmax=',vmax,dt_vmax
  
  if (rank == 0) then
     
     if (cosmo) then

        dt_e=dt_max

        n=0
        do
           n=n+1
           call expansion(a,dt_e,da_1,da_2)
           da=da_1+da_2
           ra=da/(a+da)
           if (ra.gt.ra_max) then
              dt_e=dt_e*(ra_max/ra)
           else
              exit
           endif
           if (n .gt. 10) exit
        enddo

        ! take the minimum of all the limits 

#ifdef PPINT
#ifdef PP_EXT
        dt = min(dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_pp_ext_acc,dt_c_acc)
#else
        dt = min(dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_c_acc)
#endif
#else
        dt = min(dt_e,dt_f_acc,dt_vmax,dt_c_acc)
#endif

        dt = dt * dt_scale

        call expansion(a,dt,da_1,da_2)

        da=da_1+da_2 

        ! Check to see if we are checkpointing / projecting / halofinding this step 

        checkpoint_step=.false.
        projection_step=.false.
        halofind_step=.false.
        
        am=min(a_checkpoint(cur_checkpoint),a_projection(cur_projection),a_halofind(cur_halofind))
        
        if (a_checkpoint(cur_checkpoint)==am) then
           
          if (a+da > a_checkpoint(cur_checkpoint)) then
            checkpoint_step=.true.
            dt=dt*(a_checkpoint(cur_checkpoint)-a)/da
            call expansion(a,dt,da_1,da_2)
            if (cur_checkpoint == num_checkpoints) final_step=.true.
            if (a_projection(cur_projection) == am .and. cur_projection <= num_projections) projection_step=.true. 
            if (a_halofind(cur_halofind) == am .and. cur_halofind <= num_halofinds) halofind_step=.true.
          endif
          
        elseif (a_projection(cur_projection)==am .and. cur_projection <= num_projections) then 
           
          if (a+da > a_projection(cur_projection)) then
            projection_step=.true.
            dt=dt*(a_projection(cur_projection)-a)/da
            call expansion(a,dt,da_1,da_2)
            if (a_halofind(cur_halofind) == am .and. cur_halofind <= num_halofinds) halofind_step=.true.
          endif

        elseif (a_halofind(cur_halofind) == am  .and. cur_halofind <= num_halofinds) then
           
          if (a+da > a_halofind(cur_halofind)) then
            halofind_step=.true.
            dt=dt*(a_halofind(cur_halofind)-a)/da
            call expansion(a,dt,da_1,da_2)
          endif
          
        endif
        
        !! Check to see whether we should perform checkpoint kill
        kill_step = .false.
        sec1a = mpi_wtime(ierr)
        if (rank == 0) then
           if ((sec1a - sec1) .ge. kill_time) kill_step = .true.
        endif
        call mpi_bcast(kill_step, 1, mpi_logical, 0, mpi_comm_world, ierr)

        !! Calculate timestep parameters to be used

        dt_gas=dt/4      
        da=da_1+da_2 
        ra=da/(a+da)
        a_mid=a+(da/2) !da_1

        write(*,*)
        write(*,*) 'Sweep number: ',nts
        write(*,*) 'Tau         : ',tau,tau+dt
        write(*,*) 'Redshift    : ',1.0/a-1.0,1.0/(a+da)-1.0
        write(*,*) 'Scale factor: ',a,a_mid,a+da
        write(*,*) 'Expansion   : ',ra

#ifdef PPINT
#ifdef PP_EXT
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_pp_ext_acc,dt_c_acc
#else
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_vmax,dt_pp_acc,dt_c_acc
#endif
#else
        write(*,*) 'Time step   : ',dt,dt_e,dt_f_acc,dt_vmax,dt_c_acc
#endif

        tau=tau+dt
        t=t+dt
        a=a+da
        
     else ! not cosmo
   
        a = 1.0
        a_mid = a
        da = 0.0
#ifdef PPINT
#ifdef PP_EXT
        dt = min(1.0,dt_f_acc,dt_vmax,dt_pp_acc,dt_pp_ext_acc,dt_c_acc)
#else
        dt = min(1.0,dt_f_acc,dt_vmax,dt_pp_acc,dt_c_acc)
#endif
#else
        dt = min(1.0,dt_f_acc,dt_vmax,dt_c_acc)
#endif

        t = t + dt
        if (rank == 0) write(*,*) 'nts=',nts,'t=',t,'dt=',dt,dt_f_acc,dt_vmax,dt_pp_acc, dt_pp_ext_acc,dt_c_acc,dt_max_v,0.1/sqrt(G*mass_p/cur_sep**2)
  
     endif

  endif

  ! broadcast timestep variables 

  call mpi_bcast(a,1,mpi_real,0,mpi_comm_world,ierr)
  call mpi_bcast(a_mid,1,mpi_real,0,mpi_comm_world,ierr)
  call mpi_bcast(dt,1,mpi_real,0,mpi_comm_world,ierr)
  call mpi_bcast(dt_gas,1,mpi_real,0,mpi_comm_world,ierr)
  call mpi_bcast(checkpoint_step,1,mpi_logical,0,mpi_comm_world,ierr)
  call mpi_bcast(projection_step,1,mpi_logical,0,mpi_comm_world,ierr)
  call mpi_bcast(halofind_step,1,mpi_logical,0,mpi_comm_world,ierr)
  call mpi_bcast(final_step,1,mpi_logical,0,mpi_comm_world,ierr)

  if (rank == 0) write(*,*) 'finished particle mesh'

end subroutine timestep

!! Expansion subroutine :: Hy Trac -- trac@cita.utoronto.ca
!! Added Equation of State for Dark Energy :: Pat McDonald -- pmcdonal@cita.utoronto.ca
subroutine expansion(a0,dt0,da1,da2)
  implicit none    
# include "cubepm.par"

  real(4) :: a0,dt0,dt_x,da1,da2
  real(8) :: a_x,adot,addot,atdot,arkm,a3rlm,omHsq
  real(8), parameter :: e = 2.718281828459046
  
  !! Expand Friedman equation to third order and integrate
  dt_x=dt0/2
  a_x=a0
  omHsq=4.0/9.0
  a3rlm=a_x**(-3*wde)*omega_l/omega_m
  arkm=a_x*(1.0-omega_m-omega_l)/omega_m

  adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
  addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
  atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)

  da1=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

  a_x=a0+da1
  omHsq=4.0/9.0
  a3rlm=a_x**(-3*wde)*omega_l/omega_m
  arkm=a_x*(1.0-omega_m-omega_l)/omega_m

  adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
  addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
  atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)
  
  da2=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

end subroutine expansion
