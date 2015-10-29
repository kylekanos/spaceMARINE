module advance_module

   use multifab_module
   use ml_layout_module
   use define_bc_module
   use bc_module
   use multifab_physbc_module
   use multifab_fill_ghost_module
   use ml_cc_restriction_module
   use problem_module
   use solvers
   use source_module
   use physics_declarations
   
   implicit none
   
   private
   
   public :: advance_grid, get_new_timestep
   integer :: dm

 contains
   !> @brief main advance routine
   subroutine advance_grid(mla, phi, dx, dt, nvars, the_bc_tower)
   
     type(ml_layout), intent(in) :: mla
     type(multifab) , intent(inout) :: phi(:)
     real(dp_t), intent(in) :: dx(:)
     real(dp_t), intent(in) :: dt
     integer, intent(in) :: nvars
     type(bc_tower) , intent(in) :: the_bc_tower
   
     ! local variables
     integer :: i, n, nlevs, ng
     real(dp_t) :: hdt
   
     ! an array of multifabs; one for each direction
     type(multifab) :: psi(mla%nlevel)
   
     dm = mla%dim
     nlevs = mla%nlevel
     ng = nghost(phi(1))
     hdt = 0.5_dp_t*dt
   
     ! build the flux(:,:) multifabs
     do n=1, nlevs
        ! psi(n) has NVARS components & same ghost cells as phi
        call multifab_build(psi(n), mla%la(n), nvars, ng)
        ! copy phi into psi
        call multifab_copy_c(psi(n), 1, phi(n), 1, nvars, ng)
     enddo !- n
   
   ! RK2 step 1:
     call advect(mla, hdt, dx, 1, the_bc_tower, phi, psi)
   ! update sources
     call source_box(mla, hdt, dx, the_bc_tower, psi)
   ! update boundaries
     call update_bcs(mla, psi, the_bc_tower)
   ! RK2 step 2:
     call advect(mla,  dt, dx, 2, the_bc_tower, psi, phi)
   ! update source
     call source_box(mla, hdt, dx, the_bc_tower, phi)
   ! update boundaries
     call update_bcs(mla, phi, the_bc_tower)
   
     ! make sure to destroy the multifab or you'll leak memory
     do n=1, nlevs
        call multifab_destroy(psi(n))
     enddo !- n
   
   end subroutine advance_grid
   
   !> @brief integrate the grid by one dt, level-by-level
   subroutine advect(mla, dt, dxn, step, tower, U, Unew)
     real(dp_t), intent(in) :: dt, dxn(:)
     integer, intent(in) :: step
     type(ml_layout), intent(in) :: mla
     type(bc_tower) , intent(in) :: tower
     type(multifab), intent(in) :: U(:)
     type(multifab), intent(inout) :: Unew(:)
     type(multifab), dimension(mla%nlevel, mla%dim) :: flux
     real(dp_t), pointer, dimension(:,:,:,:) :: pp, fxp, fyp, fzp
     real(dp_t) :: dx, dq
     integer :: lo(mla%dim), hi(mla%dim)
     integer :: n, i, j, k, nf, nlevs, dm, nvar, ng_p, ng_f
     
     dm = mla%dim
     nlevs = mla%nlevel
     nvar = ncomp(U(1))
     
     ! generate flux multifabs
     do n=1, nlevs
        do i=1, dm
           ! variable, layout, # nvars, # ghosts, direction
           call multifab_build_edge(flux(n, i), mla%la(n), nvar, 1, i)
        enddo !- i
     enddo !- n
     
     ng_f = flux(1, 1)%ng
     ng_p = U(1)%ng
     
     ! compute the flux for each level
     do n=1, nlevs
        dx = dxn(n)
        do i=1, nfabs(U(n))
           pp  => dataptr(U(n), i)
           fxp => dataptr(flux(n, 1), i)
           fyp => dataptr(flux(n, 2), i)
           lo  =  lwb(get_box(U(n), i))
           hi  =  upb(get_box(U(n), i))
           
           select case(dm)
             case(2)
                call compute_flux_2d(pp(:,:, 1,:), ng_p, &
                                     fxp(:,:, 1,:), fyp(:,:, 1,:), ng_f, &
                                     lo, hi, dx, step)
             case(3)
                fzp => dataptr(flux(n, 3), i)
                call compute_flux_3d(pp(:,:,:,:), ng_p, &
                                     fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:), ng_f, &
                                     lo, hi, dx)
           end select
        enddo !- i
     enddo !- n
     
     ! match fluxes between levels
     do n=nlevs, 2,-1
        do i=1, dm
           call ml_edge_restriction_c(flux(n-1, i), nvar, flux(n, i), nvar, mla%mba%rr(n-1,:), i, nvar)
        enddo !- i
     enddo !- n
     
     ! update the grid
     do n=1, nlevs
        dq = -dt/dxn(n)
        do nf=1, nfabs(U(n))
           pp  => dataptr(Unew(n), nf)
           fxp => dataptr(flux(n, 1), nf)
           fyp => dataptr(flux(n, 2), nf)
           lo  =  lwb(get_box(U(n), nf))
           hi  =  upb(get_box(U(n), nf))
           select case(dm)
             case(2)
                do j=lo(2), hi(2)
                   do i=lo(1), hi(1)
                      pp(i, j, 1,:) = pp(i, j, 1,:) + dq*(fxp(i, j, 1,:) - fxp(i+1, j, 1,:) + &
                                                          fyp(i, j, 1,:) - fyp(i, j+1, 1,:))
!~                       print '(a,i0,a,i0,a,5(1x,f18.6))',"pp(",i,",",j,":)=",pp(i,j,1,:)
                   enddo !- i
                enddo !- j
             case(3)
                fzp => dataptr(flux(n, 3), nf)
                do k=lo(3), hi(3)
                   do j=lo(2), hi(2)
                      do i=lo(1), hi(1)
                         pp(i, j, k,:) = pp(i, j, k,:) + dq*(fxp(i, j, k,:) - fxp(i+1, j, k,:) + &
                                                             fyp(i, j, k,:) - fyp(i, j+1, k,:) + &
                                                             fzp(i, j, k,:) - fzp(i, j, k+1,:))
                      enddo !- i
                   enddo !- j
                enddo !- k
           end select
        enddo !- i
     enddo !- n
     
     ! destroy multifabs to prevent leaks
     do n=1, nlevs
        do i=1, dm
           call multifab_destroy(flux(n, i))
        enddo !- i
     enddo !- n
   end subroutine advect
   
   !> @brief computes the flux for 2D runs
   subroutine compute_flux_2d(phi, ng_p, fx, fy, ng_f, lo, hi, dx, step)
     integer    :: lo(2), hi(2), ng_p, ng_f
     real(dp_t) :: phi(lo(1)-ng_p:, lo(2)-ng_p:, 1:)
     real(dp_t) ::  fx(lo(1)-ng_f:, lo(2)-ng_f:, 1:)
     real(dp_t) ::  fy(lo(1)-ng_f:, lo(2)-ng_f:, 1:)
     real(dp_t) :: dx
     integer    :: step
   
     ! local variables
     real(dp_t), dimension(:,:,:,:), allocatable :: dv
     real(dp_t), dimension(:,:,:), allocatable :: q
     real(dp_t), dimension(:), allocatable :: pl, pr
     real(dp_t) :: dmm
     integer :: i, j, n, nvar, ierr
     
     nvar = ubound(phi, dim=3)
     allocate(pl(1:nvar), pr(1:nvar))
     allocate( q(1:nvar, lo(1)-2:hi(1)+2, lo(2)-2:hi(2)+2), stat=ierr);     if(ierr /= 0) print *,"  q not allocated"
     allocate( dv(1:nvar, 2, lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1), stat=ierr); if(ierr /= 0) print *," dv not allocated"
     
     
     !$omp parallel
   ! we'll limit the primitive variables
     !$omp do private(i, j) collapse(2)
     do j=lo(2)-2, hi(2)+2
        do i=lo(1)-2, hi(1)+2
           q(:, i, j) = primitive(phi(i, j,:))
        enddo !- i
     enddo !- j
     !$omp end do
     
   ! find differentials (uses monotonized central)
     dv = zero
     !$omp do private(i, j) collapse(2)
     do j=lo(2)-1, hi(2)+1
        do i=lo(1)-1, hi(1)+1
           dv(:, 1, i, j) = slope_lim(q(:, i-1, j), q(:, i, j), q(:, i+1, j))
           dv(:, 2, i, j) = slope_lim(q(:, i, j-1), q(:, i, j), q(:, i, j+1))
        enddo !- i
     enddo !- j
     !$omp end do
     
   ! reconstruct the cell interfaces & compute fluxes in the x-dir
     fx = 0.0_dp_t
     !$omp do private(i, j) collapse(2)
     do j=lo(2)-1,hi(2)+1
        do i=lo(1),hi(1)+1

           pl(1:nvar) = q(1:nvar, i-1, j) + half*dv(1:nvar, 1, i-1, j)
           pr(1:nvar) = q(1:nvar, i  , j) - half*dv(1:nvar, 1, i  , j)
           
           fx(i, j, 1:nvar) = solver(nvar, pl(1:nvar), pr(1:nvar))
        enddo !- i
     enddo !- j
     !$omp end do

   ! reconstruct the cell interfaces & compute fluxes in the y-dir
     fy = 0.0_dp_t
     !$omp do private(i, j) collapse(2)
     do j=lo(2),hi(2)+1
        do i=lo(1)-1,hi(1)+1

           pl(1:nvar) = q(1:nvar, i, j-1) + half*dv(1:nvar, 2, i, j-1)
           pr(1:nvar) = q(1:nvar, i, j  ) - half*dv(1:nvar, 2, i, j  )

           call swapy(pl)
           call swapy(pr)

           fy(i, j, 1:nvar) = solver(nvar,pl(1:nvar), pr(1:nvar))
           call swapy(fy(i,j,:))
        enddo !- i
     enddo !- j
     !$omp end do
          
     !$omp do private(i, j, n) collapse(3)
     do n=1,nvar
        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
           
              if(fx(i,j,n) /= fx(i,j,n)) then
                 print '(3(a,i0),a,i0)',"NAN found in fx(",i,",",j,",",n,") on step ",step
              endif
              if(fy(i,j,n) /= fy(i,j,n)) then
                 print '(3(a,i0),a,i0)',"NAN found in fy(",i,",",j,",",n,") on step ",step
              endif
           enddo !- i
        enddo !- j
     enddo !- n
     !$omp end do
     !$omp end parallel
!~      stop
     
     
!~      print *,"deallocating..."
!~      deallocate(dv);  print *,"  deallocated  dv"
!~      deallocate(dvm); print *,"  deallocated dvm"
!~      deallocate(dvp); print *,"  deallocated dvp"
!~      deallocate(qMx); print *,"  deallocated qMx"
!~      deallocate(qPx); print *,"  deallocated qPx"
!~      deallocate(qPy); print *,"  deallocated qPy"
!~      deallocate(qMy); print *,"  deallocated qMy"
!~      print *,"done deallocating"
     
   end subroutine compute_flux_2d
   
   
   subroutine compute_flux_3d(phi, ng_p, fx, fy, fz, ng_f, lo, hi, dx)
     integer    :: lo(3), hi(3), ng_p, ng_f
     real(dp_t) :: phi(lo(1)-ng_p:, lo(2)-ng_p:, lo(3)-ng_p:, 1:)
     real(dp_t) ::  fx(lo(1)-ng_f:, lo(2)-ng_f:, lo(3)-ng_f:, 1:)
     real(dp_t) ::  fy(lo(1)-ng_f:, lo(2)-ng_f:, lo(3)-ng_f:, 1:)
     real(dp_t) ::  fz(lo(1)-ng_f:, lo(2)-ng_f:, lo(3)-ng_f:, 1:)
     real(dp_t) :: dx
   
!~      ! local variables
!~      real(dp_t), dimension(:,:,:,:), allocatable :: q, qPx, qMx, qMy, qPy, qMz, qPz
!~      real(dp_t), dimension(:), allocatable :: sL, sC, sR
!~      integer :: i, j, k, nvar
!~      
!~      nvar = ubound(phi, dim=3)
!~      allocate(sL(nvar), sC(nvar), sR(nvar))
!~      allocate(qPx(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:nvar),&
!~               qMx(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:nvar),&
!~               qPy(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:nvar),&
!~               qMy(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:nvar),&
!~               qPz(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:nvar),&
!~               qMz(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1, 1:nvar))
!~      
!~      
!~      ! convert conservative to primitive
!~      !$omp parallel     
!~      ! set PWL reconstruction...maybe fix later to PPM?
!~      !$omp do private(i, j, k, sL, sR, sC)
!~      do k=lo(3)-1, hi(3)+1
!~         do j=lo(2)-1, hi(2)+1
!~            do i=lo(1)-1, hi(1)+1
!~               sL(:) = q(i  , j, k,:) - q(i-1, j, k,:)
!~               sR(:) = q(i+1, j, k,:) - q(i  , j, k,:)
!~               sC(:) = minmod(sL(:), sR(:))
!~               qPx(i  , j, k,:) = q(i, j, k,:) - sC(:)*0.5_dp_t
!~               qMx(i+1, j, k,:) = q(i, j, k,:) + sC(:)*0.5_dp_t
!~               
!~               sL(:) = q(i, j  , k,:) - q(i, j-1, k,:)
!~               sR(:) = q(i, j+1, k,:) - q(i, j  , k,:)
!~               sC(:) = minmod(sL(:), sR(:))
!~               qPy(i, j  , k,:) = q(i, j, k,:) - sC(:)*0.5_dp_t
!~               qMy(i, j+1, k,:) = q(i, j, k,:) + sC(:)*0.5_dp_t
!~               
!~               sL(:) = q(i, j, k  ,:) - q(i, j, k-1,:)
!~               sR(:) = q(i, j, k+1,:) - q(i, j, k  ,:)
!~               sC(:) = minmod(sL(:), sR(:))
!~               qPy(i, j, k  ,:) = q(i, j, k,:) - sC(:)*0.5_dp_t
!~               qPy(i, j, k+1,:) = q(i, j, k,:) + sC(:)*0.5_dp_t
!~            enddo !- i
!~         enddo !- j
!~      enddo !- k
!~      !$omp end do
!~      
!~      ! compute fluxes
!~      !$omp do private(i, j, k)
!~      do k=lo(3)-1, hi(3)+1
!~         do j=lo(2)-1, hi(2)+1
!~            do i=lo(1)-1, hi(1)+1
!~               call roe(nvar, 1, qPx(i, j, k,:), qMx(i, j, k,:), fx(i, j, k,:))
!~               call roe(nvar, 2, qPy(i, j, k,:), qMy(i, j, k,:) , fy(i, j, k,:))
!~               call roe(nvar, 3, qPz(i, j, k,:), qMz(i, j, k,:) , fz(i, j, k,:))
!~            enddo !- i
!~         enddo !- j
!~      enddo !- k
!~      !$omp end do
!~      !$omp end parallel
!~      
!~      deallocate(q, qMx, qPx, qPy, qMy, qPz, qMz)
   
   end subroutine compute_flux_3d
   
   !> @brief finds new timestep
   subroutine get_new_timestep(mla, U, dx, dt)
     type(ml_layout), intent(in) :: mla
     type(multifab) , intent(in) :: U(:)
     real(dp_t), intent(in) :: dx(:)
     real(dp_t), intent(inout) :: dt
     real(dp_t), pointer :: q(:,:,:,:)
     real(dp_t) :: max_speed, CFL, cs, a, c2, p, cfx, cfy, cfz
     integer :: i, j, k, nvar, lo(3), hi(3), n, nf
     integer :: dm, nlevs
     
     dm = mla%dim
     nlevs = mla%nlevel
     nvar = ncomp(U(1))
     
     CFL = 0.3d0
     max_speed = 0d0
     
     ! will need MPI calls here to make global dt
     
     do n=1, nlevs
        do nf=1, nfabs(U(n))
           lo = 1; hi = 1
           lo(1:mla%dim) = lwb(get_box(U(n), nf))
           hi(1:mla%dim) = upb(get_box(U(n), nf))
           q => dataptr(U(n), nf)
           
           do k=lo(3), hi(3)
              do j=lo(2), hi(2)
                 do i=lo(1), hi(1)
                    p = press(q(i, j, k,:))
                    c2 = gamma*p/q(i, j, k, 1)
                    cs = sqrt(c2)
                    cfx = cs + abs(q(i, j, k, 2))
                    cfy = cs + abs(q(i, j, k, 3))
                    cfz = cs + abs(q(i, j, k, 4))
                    max_speed = max(max_speed, cfx, cfy, cfz)
                 enddo !- i
              enddo !- j
           enddo !- k
        enddo !- nf
     enddo !- n
     
     dt = CFL * dx(1) / max_speed
     
   end subroutine get_new_timestep

   !> @brief update the boundaries
   subroutine update_bcs(mla, phi, the_bc_tower)
     type(ml_layout), intent(in) :: mla
     type(multifab) , intent(inout) :: phi(:)
     type(bc_tower) , intent(in) :: the_bc_tower

     integer nlevs, n

     nlevs = mla%nlevel
     
     
     if (nlevs == 1) then
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
        call multifab_fill_boundary(phi(nlevs))
   
       ! fill non-periodic domain boundary ghost cells
        call multifab_physbc(phi(nlevs),1,1,1,the_bc_tower%bc_tower_array(nlevs))
     else
       ! the loop over nlevs must count backwards to make sure the finer grids are done first
        do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
           call ml_cc_restriction(phi(n-1),phi(n),mla%mba%rr(n-1,:))
           ! fill level n ghost cells using interpolation from level n-1 data
           ! note that multifab_fill_boundary and multifab_physbc are called for
           ! both levels n-1 and n
           call multifab_fill_ghost_cells(phi(n),phi(n-1),phi(n)%ng,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1), &
                                          the_bc_tower%bc_tower_array(n),1,1,1)
        enddo
     endif
   end subroutine update_bcs
   

end module advance_module

