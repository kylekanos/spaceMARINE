program spacemarine

   use boxlib
   use multifab_module
   use bl_IO_module
   use ml_layout_module
   use problem_module
   use write_plotfile_module
   use advance_module
   use define_bc_module
   use make_new_grids_module
   use regrid_module
   use solvers, only: solver_init
   use source_module, only: source_init
   
   implicit none
   
   ! stuff you can set with the inputs file (otherwise use default values below)
   integer    :: max_levs, dim, nsteps, plot_int, n_cell, max_grid_size
   integer    :: amr_buf_width, cluster_minwidth, cluster_blocking_factor
   real(dp_t) :: cluster_min_eff
   integer    :: regrid_int
   integer    :: bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi
   
   ! dummy indices using for reading in inputs file
   integer :: un, farg, narg
   logical :: need_inputs_file, found_inputs_file
   character(len=128) :: inputs_file_name
   
   ! will be allocated with dim components
   integer       , allocatable :: lo(:), hi(:)
   logical       , allocatable :: is_periodic(:)
   real(dp_t)    , allocatable :: prob_lo(:), prob_hi(:)
   
  ! will be allocated with (dim,2) components
   integer       , allocatable :: phys_bc(:,:)
   
  ! will be allocated with max_levs components
   real(dp_t)    , allocatable :: dx(:)
   type(multifab), allocatable :: phi(:)
   
   integer    :: istep,i,n,nl,nlevs,nvars
   logical    :: new_grid
   real(dp_t) :: dt,time,start_time,run_time,run_time_IOproc
   
   type(box)         :: bx
   type(ml_boxarray) :: mba
   type(ml_layout)   :: mla
   type(layout), allocatable :: la_array(:)
   
   type(bc_tower) :: the_bc_tower
   
   namelist /probin/ max_levs, nvars, dim, nsteps, plot_int, n_cell, max_grid_size, amr_buf_width, &
        cluster_minwidth, cluster_blocking_factor, cluster_min_eff, regrid_int, &
        bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
   call boxlib_initialize()
   
   if( parallel_IOProcessor() ) call print_name()

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
   start_time = parallel_wtime()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! default values - will get overwritten by the inputs file
   max_levs      = 4
   dim           = 2
   nvars         = 5        ! always include vz, even if not actually used
   nsteps        = 1000
   plot_int      = 100
   n_cell        = 32
 ! AMR stuff
   max_grid_size = n_cell/2
   amr_buf_width = 4
   cluster_minwidth = 4
   cluster_blocking_factor = 4
   cluster_min_eff = 0.8d0
   regrid_int = 4

  ! allowable options for this example are
  ! -1 = PERIODIC
  ! 12 = OUTLET (dphi/dn=0 at boundary)
  ! 15 = NO_SLIP_WALL (wall with fixed phi=1)
   bc_x_lo       = 15 ! NO_SLIP
   bc_x_hi       = 12 ! OUTLET
   bc_y_lo       = 15 ! NO_SLIP
   bc_y_hi       = 12 ! OUTLET
   bc_z_lo       = -1 ! PERIODIC
   bc_z_hi       = -1 ! PERIODIC

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! read inputs file and overwrite any default values
   narg = command_argument_count()
   need_inputs_file = .true.
   farg = 1
   if( narg >= 1 ) then
      call get_command_argument(farg, value = inputs_file_name)
      inquire(file = inputs_file_name, exist = found_inputs_file )
      if( found_inputs_file ) then
         farg = farg + 1
         un = unit_new()
         open(unit=un, file = inputs_file_name, status = 'old', action = 'read')
         read(unit=un, nml = probin)
         close(unit=un)
         need_inputs_file = .false.
      endif
   endif

  ! now that we have dim, we can allocate these
   allocate(lo(dim),hi(dim))
   allocate(is_periodic(dim))
   allocate(prob_lo(dim),prob_hi(dim))
   allocate(phys_bc(dim,2))

  ! now that we have max_levs, we can allocate these
   allocate(dx(max_levs))
   allocate(phi(max_levs))
   allocate(la_array(max_levs))
   
  ! put all the domain boundary conditions into phys_bc
   phys_bc(1,1) = bc_x_lo
   phys_bc(1,2) = bc_x_hi
   phys_bc(2,1) = bc_y_lo
   phys_bc(2,2) = bc_y_hi
   if(dim == 3) then
      phys_bc(3,1) = bc_z_lo
      phys_bc(3,2) = bc_z_hi
   endif

  ! build an array indicating periodicity in each direction
   is_periodic(:) = .false.
   do i=1,dim
      if(phys_bc(i,1) == -1 .and. phys_bc(i,2) .ne. -1) then
         call bl_error("Invalid BC's - both lo and hi need to be periodic")
      endif
      if(phys_bc(i,2) == -1 .and. phys_bc(i,1) .ne. -1) then
         call bl_error("Invalid BC's - both lo and hi need to be periodic")
      endif
      if(phys_bc(i,1) == -1 .and. phys_bc(i,2) == -1) then
         is_periodic(i) = .true.
      endif
   enddo
 
   call cluster_set_minwidth(cluster_minwidth)
   call cluster_set_blocking_factor(cluster_blocking_factor)
   call cluster_set_min_eff(cluster_min_eff)

  ! tell mba about max_levs and dimensionality of problem
   call ml_boxarray_build_n(mba,max_levs,dim)

  ! tell mba about the ref_ratio between levels
  ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
  ! we use refinement ratio of 2 in every direction between all levels
   do n=2,max_levs
      mba%rr(n-1,:) = 2
   enddo

  ! physical problem is a box on (-1,-1) to (1,1)
   prob_lo(:) =  0.d0
   prob_hi(:) =  1.d0

  ! set grid spacing at each level
  ! the grid spacing is the same in each direction
   dx(1) = (prob_hi(1)-prob_lo(1)) / n_cell
   do n=2,max_levs
      dx(n) = dx(n-1) / mba%rr(n-1,1)
   enddo
   
  ! initialize the variables
   call solver_init(nvars)
   call source_init()

  ! tell the_bc_tower about max_levs, dim, and phys_bc
   call bc_tower_init(the_bc_tower,max_levs,dim,phys_bc)

  ! create a box from (0,0) to (n_cell-1,n_cell-1)
   lo(:) = 0
   hi(:) = n_cell-1
   bx = make_box(lo,hi)

  ! tell mba about the problem domain at each level
   mba%pd(1) = bx
   do n=2,max_levs
      mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
   enddo

  ! initialize the boxarray at level 1 to be one single box
   call boxarray_build_bx(mba%bas(1),bx)

  ! overwrite the boxarray at level 1 to respect max_grid_size
   call boxarray_maxsize(mba%bas(1),max_grid_size)

  ! build the level 1 layout
   call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),is_periodic)

  ! build the level 1 multifab with NVARs component and 2 ghost cell
   call multifab_build(phi(1),la_array(1),nvars,2)

  ! define level 1 of the_bc_tower
   call bc_tower_level_build(the_bc_tower,1,la_array(1))

  ! initialize phi on level 1
   call init_phi_on_level(phi(1),dx(1),prob_lo,nvars,the_bc_tower%bc_tower_array(1))

   nl = 1
   new_grid = .true.

   do while( (nl < max_levs) .and. (new_grid) )
     ! determine whether we need finer grids based on tagging criteria
     ! if so, return new_grid=T and the la_array(nl+1)
      call make_new_grids(new_grid,la_array(nl),la_array(nl+1),phi(nl),dx(nl), &
                          amr_buf_width,mba%rr(nl,1),nl,max_grid_size)
     
      if(new_grid) then
        ! tell mba about the finer level boxarray
         call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

        ! Build the level nl+1 data
         call multifab_build(phi(nl+1),la_array(nl+1),nvars,2)
        
        ! define level nl+1 of the_bc_tower
         call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))
            
        ! initialize phi on level nl+1
         call init_phi_on_level(phi(nl+1),dx(nl+1),prob_lo,nvars,the_bc_tower%bc_tower_array(nl+1))

        ! increment current level counter
         nl = nl+1
      endif
   enddo

  ! the current number of levels in the simulation is nlevs, not necessarily max_levs
   nlevs = nl

  ! destroy phi - we are going to build it again using the new multilevel
  ! layout after we have tested and reconfigured the grids due to proper nesting
   do n=1,nlevs
      call multifab_destroy(phi(n))
   enddo

   if(nlevs >= 3) then
     ! check for proper nesting
      call enforce_proper_nesting(mba,la_array,max_grid_size)
   endif

   do n=1,nlevs
      call destroy(la_array(n))
   enddo

  ! tell mla that there are nlevs levels, not max_levs
   call ml_layout_restricted_build(mla,mba,nlevs,is_periodic)
     
  ! this makes sure the boundary conditions are properly defined everywhere
   do n = 1,nlevs
      call bc_tower_level_build(the_bc_tower,n,mla%la(n))
   enddo

   do n=1,nlevs
      call multifab_build(phi(n),mla%la(n),nvars,2)
   enddo

   call init_phi(mla,phi,dx,prob_lo,nvars,the_bc_tower)

   call destroy(mba)

   istep = 0
   time = 0.d0

  ! choose a time step with a diffusive CFL of 0.9 base on resolution
  ! at max_levs, even if nlevs < max_levs
   dt = 0.9d0*dx(max_levs)**2/(2.d0*dim)

  ! write out plotfile 0
   !call regrid(mla,phi,nlevs,max_levs,dx,the_bc_tower,amr_buf_width,max_grid_size)
   call write_plotfile(mla,phi,istep,dx,nvars,time,prob_lo,prob_hi)
   !stop

   do istep=1,nsteps

     ! regrid
      if( istep > 1 .and. max_levs > 1 ) then
         if( regrid_int > 0 ) then
            if( (mod(istep-1,regrid_int) == 0) ) then
               call regrid(mla,phi,nlevs,max_levs,dx,the_bc_tower,amr_buf_width,max_grid_size)
            endif
         endif
      endif

     ! we only want one processor to write to screen
      if( parallel_IOProcessor() ) then
         print*,'Advancing time step',istep,'with dt=',dt
      endif
     
     ! advance phi & then update dt for the new
      call advance(mla,phi,dx,dt,nvars,the_bc_tower)
      time = time + dt
      call get_new_timestep(mla, phi, dx, dt)

      if(mod(istep,plot_int) == 0 .or. istep == nsteps) then
        ! write out plotfile
         call write_plotfile(mla,phi,istep,dx,nvars,time,prob_lo,prob_hi)
      endif
   enddo

  ! make sure to destroy the multifab or you'll leak memory
   do n=1,nlevs
      call destroy(phi(n))
   enddo
   call destroy(mla)
   call bc_tower_destroy(the_bc_tower)

  deallocate(lo,hi,is_periodic,prob_lo,prob_hi)

  ! deallocate temporary boxarrays and communication mappings
  call layout_flush_copyassoc_cache ()

  ! check for memory that should have been deallocated
  if( parallel_IOProcessor() ) then
     print*, 'MEMORY STATS AT END OF PROGRAM'
     print*, ' '
  endif
  call print(multifab_mem_stats(),    "    multifab")
  call print(fab_mem_stats(),         "         fab")
  call print(boxarray_mem_stats(),    "    boxarray")
  call print(layout_mem_stats(),      "      layout")
  call print(boxassoc_mem_stats(),    "    boxassoc")
  call print(fgassoc_mem_stats(),     "     fgassoc")
  call print(syncassoc_mem_stats(),   "   syncassoc")
  call print(copyassoc_mem_stats(),   "   copyassoc")
  call print(fluxassoc_mem_stats(),   "   fluxassoc")

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  run_time = parallel_wtime() - start_time

  ! collect run_time from each processor and store the maximum
  call parallel_reduce(run_time_IOproc, run_time, MPI_MAX, &
                       proc = parallel_IOProcessorNode())

  if( parallel_IOProcessor() ) then
     print*,"Run time (s) =",run_time_IOproc
  endif

  call boxlib_finalize()
  
contains

  !> print name of code
  subroutine print_name()
    print*,">----------------------------------------------------------------------------<"
    print*,"||   _____                      __  __          _____  _____ _   _ ______   ||"
    print*,"||  / ____|                    |  \/  |   /\   |  __ \|_   _| \ | |  ____|  ||"
    print*,"|| | (___  _ __   __ _  ___ ___| \  / |  /  \  | |__) | | | |  \| | |__     ||"
    print*,"||  \___ \| '_ \ / _` |/ __/ _ \ |\/| | / /\ \ |  _  /  | | | . ` |  __|    ||"
    print*,"||  ____) | |_) | (_| | (_|  __/ |  | |/ ____ \| | \ \ _| |_| |\  | |____   ||"
    print*,"|| |_____/| .__/ \__,_|\___\___|_|  |_/_/    \_\_|  \_\_____|_| \_|______|  ||"
    print*,"||        | |                                                               ||"
    print*,"||        |_|                   By Joshua D Wood, Clemson Univ              ||"
    print*,"||                                 Kyle Kanos                               ||"
    print*,">----------------------------------------------------------------------------<"
    
  end subroutine print_name

end program spacemarine