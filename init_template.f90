module problem_module
   use multifab_module
   use ml_layout_module
   use define_bc_module
   use multifab_physbc_module
   use multifab_fill_ghost_module
   use ml_restriction_module
   use physics_declarations
   
   implicit none
   
   private
   public :: init_phi_on_level, init_phi

 contains
   
   !> @brief Initializes phi per level
   subroutine init_phi_on_level(phi,dx,prob_lo,nvars,the_bc_level)
    ! input/output variables
     type(multifab) , intent(inout) :: phi
     real(dp_t), intent(in) :: dx, prob_lo(phi%dim)
     type(bc_level) , intent(in) :: the_bc_level
     integer, intent(in) :: nvars
    ! local
     integer i,ng,dm
     integer :: lo(phi%dim), hi(phi%dim)

     real(dp_t), pointer :: dp(:,:,:,:)

     ng = phi%ng
     dm = phi%dim

     do i=1,nfabs(phi)
        dp => dataptr(phi,i)
        lo = lwb(get_box(phi,i))
        hi = upb(get_box(phi,i))
        select case(dm)
          case (2)
             call init_phi_2d(dp(:,:,1,:), ng, lo, hi, prob_lo, nvars, dx)
          case (3)
             call init_phi_3d(dp(:,:,:,:), ng, lo, hi, prob_lo, nvars, dx)
        end select
     enddo

     call multifab_fill_boundary(phi)
    
     call multifab_physbc(phi,1,1,1,the_bc_level)

   end subroutine init_phi_on_level
  
   !> @brief initialize all phi's
   subroutine init_phi(mla,phi,dx,prob_lo,nvars,the_bc_tower)
     type(ml_layout), intent(in   ) :: mla
     type(multifab) , intent(inout) :: phi(:)
     real(dp_t), intent(in) :: dx(:), prob_lo(mla%dim)
     type(bc_tower) , intent(in) :: the_bc_tower
     integer, intent(in) :: nvars   
    ! local variables
     integer :: lo(mla%dim), hi(mla%dim)
     integer :: nlevs, dm, ng, i, n
   
     real(dp_t), pointer :: dp(:,:,:,:)
   
     ng = phi(1)%ng
     dm = mla%dim
     nlevs = mla%nlevel
   
     do n=1,nlevs
        do i=1,nfabs(phi(n))
           dp => dataptr(phi(n),i)
           lo = lwb(get_box(phi(n),i))
           hi = upb(get_box(phi(n),i))
           select case(dm)
             case (2)
                call init_phi_2d(dp(:,:,1,:), ng, lo, hi, prob_lo, nvars, dx(n))
             case (3)
                call init_phi_3d(dp(:,:,:,:), ng, lo, hi, prob_lo, nvars, dx(n))
           end select
        enddo
     enddo
   
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
   end subroutine init_phi

   !> @brief set up the 2D problem
   subroutine init_phi_2d(phi, ng, lo, hi, prob_lo, nvars, dx)
     integer    :: lo(2), hi(2), ng, nvars
     real(dp_t) :: phi(lo(1)-ng:,lo(2)-ng:,1:)
     real(dp_t) :: prob_lo(2)
     real(dp_t) :: dx
   
     ! local varables
     integer :: i,j

     !$omp parallel do private(i,j,x,r2)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           
        enddo !- i
     enddo !- j
     !$omp end parallel do
   
     end subroutine init_phi_2d
   
     subroutine init_phi_3d(phi, ng, lo, hi, prob_lo, nvars, dx)
   
     integer    :: lo(3), hi(3), ng, nvars
     real(dp_t) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,1:)
     real(dp_t) :: prob_lo(3)
     real(dp_t) :: dx
   
     ! local varables
     integer    :: i,j,k
   
     !$omp parallel do private(i,j,k,x,r2)
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)

           enddo !- i
        enddo !- j
     enddo !- k
     !$omp end parallel do
   
   end subroutine init_phi_3d
   
   !> returns position vector of domain
   !! @param lo holds LOWER edge of cell, so we adjust by 1/2 dx
   pure function CellPos3D(lo, dx, i, j, k) result(x)
     real(dp_t), intent(in) :: lo(3), dx
     integer, intent(in) :: i, j, k
     real(dp_t) :: x(3)
     x(1) = lo(1) + (real(i) + 0.5_dp_t)*dx
     x(2) = lo(2) + (real(j) + 0.5_dp_t)*dx
     x(3) = lo(3) + (real(k) + 0.5_dp_t)*dx
   end function CellPos3D
   
   !> returns position vector of domain
   pure function CellPos2D(lo, dx, i, j) result(x)
     real(dp_t), intent(in) :: lo(2), dx
     integer, intent(in) :: i, j
     real(dp_t) :: x(2)
     x(1) = lo(1) + (real(i) + 0.5_dp_t)*dx
     x(2) = lo(2) + (real(j) + 0.5_dp_t)*dx
   end function CellPos2D

end module problem_module
