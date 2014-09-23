module advance_module

   use multifab_module
   use ml_layout_module
   use define_bc_module
   use bc_module
   use multifab_physbc_module
   use multifab_fill_ghost_module
   use ml_restriction_module
   use problem_module
   use solvers
   use source_module
   use physics_declarations
   
   implicit none
   
   private
   
   public :: advance, get_new_timestep
   integer :: dm

 contains
   !> @brief main advance routine
   subroutine advance(mla,phi,dx,dt,nvars,the_bc_tower)
   
     type(ml_layout), intent(in) :: mla
     type(multifab) , intent(inout) :: phi(:)
     real(dp_t), intent(in) :: dx(:)
     real(dp_t), intent(in) :: dt
     integer, intent(in) :: nvars
     type(bc_tower) , intent(in) :: the_bc_tower
   
     ! local variables
     integer :: i,n,nlevs,ng
     real(dp_t) :: hdt
   
     ! an array of multifabs; one for each direction
     type(multifab) :: psi(mla%nlevel)
   
     dm = mla%dim
     nlevs = mla%nlevel
     ng = nghost(phi(1))
     hdt = 0.5_dp_t*dt
   
     ! build the flux(:,:) multifabs
     do n=1,nlevs
        ! psi(n) has NVARS components & same ghost cells as phi
        call multifab_build(psi(n),mla%la(n),nvars,ng)
        ! copy phi into psi
        call multifab_copy_c(psi(n),1,phi(n),1,nvars,ng)
     enddo !- n
   
   ! RK2 step 1:
     call advect(mla, hdt, dx, the_bc_tower, phi, psi)
   ! update sources
     call source_box(mla, hdt, dx, the_bc_tower, psi)
   ! RK2 step 2:
     call advect(mla,  dt, dx, the_bc_tower, psi, phi)
   ! update source
     call source_box(mla, hdt, dx, the_bc_tower, phi)
   
     ! make sure to destroy the multifab or you'll leak memory
     do n=1,nlevs
        call multifab_destroy(psi(n))
     enddo !- n
   
   end subroutine advance
   
   !> @brief integrate the grid by one dt, level-by-level
   subroutine advect(mla, dt, dxn, tower, U, Unew)
     real(dp_t), intent(in) :: dt, dxn(:)
     type(ml_layout), intent(in) :: mla
     type(bc_tower) , intent(in) :: tower
     type(multifab), intent(in) :: U(:)
     type(multifab), intent(inout) :: Unew(:)
     type(multifab), dimension(mla%nlevel,mla%dim) :: flux
     real(dp_t), pointer, dimension(:,:,:,:) :: pp, fxp, fyp, fzp
     real(dp_t) :: dx, dq
     integer :: lo(mla%dim), hi(mla%dim)
     integer :: n, i, j, k, nf, nlevs, dm, nvar, ng_p, ng_f
     
     dm = mla%dim
     nlevs = mla%nlevel
     nvar = ncomp(U(1))
     
     ! generate flux multifabs
     do n=1,nlevs
        do i=1,dm
           ! variable, layout, # nvars, # ghosts, direction
           call multifab_build_edge(flux(n,i), mla%la(n), nvar, 1, i)
        enddo !- i
     enddo !- n
     
     ng_f = flux(1,1)%ng
     ng_p = U(1)%ng
     
     ! compute the flux for each level
     do n=1,nlevs
        dx = dxn(n)
        do i=1,nfabs(U(n))
           pp => dataptr(U(n),i)
           fxp => dataptr(flux(n,1),i)
           fyp => dataptr(flux(n,2),i)
           lo = lwb(get_box(U(n),i))
           hi = upb(get_box(U(n),i))
           
           select case(dm)
             case(2)
                call compute_flux_2d(pp(:,:,1,:), ng_p, &
                                     fxp(:,:,1,:), fyp(:,:,1,:), ng_f, &
                                     lo, hi, dx)
             case(3)
                fzp => dataptr(flux(n,3),i)
                call compute_flux_3d(pp(:,:,:,:), ng_p, &
                                     fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:), ng_f, &
                                     lo, hi, dx)
           end select
        enddo !- i
     enddo !- n
     
     ! match fluxes between levels
     do n=nlevs,2,-1
        do i=1,dm
           call ml_edge_restriction_c(flux(n-1,i),nvar,flux(n,i),nvar,mla%mba%rr(n-1,:),i,nvar)
        enddo !- i
     enddo !- n
     
     ! update the grid
     do n=1,nlevs
        dq = -dt/dxn(n)
        do nf=1,nfabs(U(n))
           pp => dataptr(Unew(n),nf)
           fxp => dataptr(flux(n,1),nf)
           fyp => dataptr(flux(n,2),nf)
           lo = lwb(get_box(U(n),nf))
           hi = upb(get_box(U(n),nf))
           select case(dm)
             case(2)
                do j=lo(2),hi(2)
                   do i=lo(1),hi(1)
                      pp(i,j,1,:) = pp(i,j,1,:) + dq*(fxp(i,j,1,:) - fxp(i+1,j,1,:) + &
                                                      fyp(i,j,1,:) - fyp(i,j+1,1,:))
                   enddo !- i
                enddo !- j
             case(3)
                fzp => dataptr(flux(n,3),nf)
                do k=lo(3),hi(3)
                   do j=lo(2),hi(2)
                      do i=lo(1),hi(1)
                         pp(i,j,k,:) = pp(i,j,k,:) + dq*(fxp(i,j,k,:) - fxp(i+1,j,k,:) + &
                                                         fyp(i,j,k,:) - fyp(i,j+1,k,:) + &
                                                         fzp(i,j,k,:) - fzp(i,j,k+1,:))
                      enddo !- i
                   enddo !- j
                enddo !- k
           end select
        enddo !- i
     enddo !- n
     
     ! destroy multifabs to prevent leaks
     do n=1,nlevs
        do i=1,dm
           call multifab_destroy(flux(n,i))
        enddo !- i
     enddo !- n
   end subroutine advect
   
   !> @brief computes the flux for 2D runs
   subroutine compute_flux_2d(phi, ng_p, fx, fy, ng_f, lo, hi, dx)
     integer    :: lo(2), hi(2), ng_p, ng_f
     real(dp_t) :: phi(lo(1)-ng_p:,lo(2)-ng_p:,1:)
     real(dp_t) ::  fx(lo(1)-ng_f:,lo(2)-ng_f:,1:)
     real(dp_t) ::  fy(lo(1)-ng_f:,lo(2)-ng_f:,1:)
     real(dp_t) :: dx
   
     ! local variables
     real(dp_t), dimension(:,:,:,:), allocatable :: dv
     real(dp_t), dimension(:,:,:), allocatable :: q,qPx, qMx, qMy, qPy
     real(dp_t), dimension(:), allocatable :: dvp, dvm
     real(dp_t) :: dmm
     integer :: i,j,nvar,ierr
     
     nvar = ubound(phi,dim=3)
     allocate(dvp(1:nvar), dvm(1:nvar))
     allocate(  q(1:nvar,lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2),stat=ierr); if(ierr /= 0) print *,"  q not allocated"
     allocate(qPx(1:nvar,lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1),stat=ierr); if(ierr /= 0) print *,"qPx not allocated"
     allocate(qMx(1:nvar,lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1),stat=ierr); if(ierr /= 0) print *,"qMx not allocated"
     allocate(qPy(1:nvar,lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1),stat=ierr); if(ierr /= 0) print *,"qPy not allocated"
     allocate(qMy(1:nvar,lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1),stat=ierr); if(ierr /= 0) print *,"qMy not allocated"
     allocate( dv(1:nvar,2,lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1),stat=ierr); if(ierr /= 0) print *," dv not allocated"
     
     
     !$omp parallel
   ! we'll limit the primitive variables
     !$omp do private(i,j) collapse(2)
     do j=lo(2)-2,hi(2)+2
        do i=lo(1)-2,hi(1)+2
           q(:,i,j) = primitive(phi(i,j,:))
        enddo !- i
     enddo !- j
     !$omp end do
     
   ! find differentials (these will be needed later for protection)
     !$omp do private(i,j) collapse(2)
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           dv(:,1,i,j) = q(:,i+1,j) - q(:,i,j)
           dv(:,2,i,j) = q(:,i,j+1) - q(:,i,j)
        enddo !- i
     enddo !- j
     !$omp end do
     
     
   ! 3rd order reconstruction via 3-point stencil (i-1, i, i+1)
     !$omp do private(i,j,dvp,dvm,dmm) collapse(2)
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
         ! x-dir
           dvp(:) = dv(:,1,i  ,j)
           dvm(:) = dv(:,1,i-1,j)
           qPx(:,i,j) = q(:,i,j) + half*dvp(:)*limiter(dvp(:), dvm(:), dx)
           qMx(:,i,j) = q(:,i,j) - half*dvm(:)*limiter(dvp(:), dvm(:), dx)
         ! y-dir
           dvp(:) = dv(:,2,i,j  )
           dvm(:) = dv(:,2,i,j-1)
           qPy(:,i,j) = q(:,i,j) + half*dvp(:)*limiter(dvp(:), dvm(:), dx)
           qMy(:,i,j) = q(:,i,j) - half*dvm(:)*limiter(dvp(:), dvm(:), dx)
        enddo !- i
     enddo !- j
     !$omp end do
     
   ! this limiter is not strictly TVD, so we need to ensure positivity of rho & P
     !$omp do private(i,j,dmm) collapse(2)
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
         ! density test
           if((qMx(1,i,j) < zero) .or. (qPx(1,i,j) < zero)) then
              dmm = minmod(dv(1,1,i,j),dv(1,1,i+1,j))
              qPx(1,i,j) = q(1,i,j) + half*dmm
              qMx(1,i,j) = q(1,i,j) - half*dmm
           endif
           if((qMy(1,i,j) < zero) .or. (qPy(1,i,j) < zero)) then
              dmm = minmod(dv(1,2,i,j),dv(1,2,i,j+1))
              qPy(1,i,j) = q(1,i,j) + half*dmm
              qMy(1,i,j) = q(1,i,j) - half*dmm
           endif
         ! pressure test
           if((qMx(5,i,j) < zero) .or. (qPx(5,i,j) < zero)) then
              dmm = minmod(dv(5,1,i,j),dv(5,1,i+1,j))
              qPx(5,i,j) = q(5,i,j) + half*dmm
              qMx(5,i,j) = q(5,i,j) - half*dmm
           endif
           if((qMy(5,i,j) < zero) .or. (qPy(5,i,j) < zero)) then
              dmm = minmod(dv(5,2,i,j),dv(5,2,i,j+1))
              qPy(5,i,j) = q(5,i,j) + half*dmm
              qMy(5,i,j) = q(5,i,j) - half*dmm
           endif
        enddo !- i
     enddo !- j
     !$omp end do
   
   ! convert primitive interpolants to conservative ones (Riemann solvers expect conservative)
!~      !$omp do private(i,j) collapse(2)
!~      do j=lo(2)-1,hi(2)+1
!~         do i=lo(1)-1,hi(1)+1
!~            qPx(:,i,j) = conservative(qPx(:,i,j))
!~            qMx(:,i,j) = conservative(qMx(:,i,j))
!~            qPy(:,i,j) = conservative(qPy(:,i,j))
!~            qMy(:,i,j) = conservative(qMy(:,i,j))
!~         enddo !- i
!~      enddo !- j
!~      !$omp end do
     
   ! compute fluxes
     !$omp do private(i,j) collapse(2)
     do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
           call roe(nvar, 1, qPx(:,i,j), qMx(:,i,j), fx(i,j,:))
           call roe(nvar, 2, qPy(:,i,j), qMy(:,i,j) ,fy(i,j,:))
        enddo !- i
     enddo !- j
     !$omp end do
     
     
     !$omp do private(i,j) collapse(2)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           if(any(fx(i,j,:)/=fx(i,j,:))) then
              print *,"NAN found in fx at ",i,j
           endif
           if(any(fy(i,j,:)/=fy(i,j,:))) then
              print *,"NAN found in fy at ",i,j
           endif
        enddo !- i
     enddo !- j
     !$omp end do
     !$omp end parallel
     
     !$om
     
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
     real(dp_t) :: phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:, 1:)
     real(dp_t) ::  fx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:, 1:)
     real(dp_t) ::  fy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:, 1:)
     real(dp_t) ::  fz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:, 1:)
     real(dp_t) :: dx
   
     ! local variables
     real(dp_t), dimension(:,:,:,:), allocatable :: q, qPx, qMx, qMy, qPy, qMz, qPz
     real(dp_t), dimension(:), allocatable :: sL, sC, sR
     integer :: i,j,k,nvar
     
     nvar = ubound(phi,dim=3)
     allocate(sL(nvar),sC(nvar),sR(nvar))
     allocate(qPx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:nvar),&
              qMx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:nvar),&
              qPy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:nvar),&
              qMy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:nvar),&
              qPz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:nvar),&
              qMz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:nvar))
     
     
     ! convert conservative to primitive
     !$omp parallel     
     ! set PWL reconstruction...maybe fix later to PPM?
     !$omp do private(i,j,k,sL,sR,sC)
     do k=lo(3)-1,hi(3)+1
        do j=lo(2)-1,hi(2)+1
           do i=lo(1)-1,hi(1)+1
              sL(:) = q(i  ,j,k,:) - q(i-1,j,k,:)
              sR(:) = q(i+1,j,k,:) - q(i  ,j,k,:)
              sC(:) = minmod(sL(:),sR(:))
              qPx(i  ,j,k,:) = q(i,j,k,:) - sC(:)*0.5_dp_t
              qMx(i+1,j,k,:) = q(i,j,k,:) + sC(:)*0.5_dp_t
              
              sL(:) = q(i,j  ,k,:) - q(i,j-1,k,:)
              sR(:) = q(i,j+1,k,:) - q(i,j  ,k,:)
              sC(:) = minmod(sL(:),sR(:))
              qPy(i,j  ,k,:) = q(i,j,k,:) - sC(:)*0.5_dp_t
              qMy(i,j+1,k,:) = q(i,j,k,:) + sC(:)*0.5_dp_t
              
              sL(:) = q(i,j,k  ,:) - q(i,j,k-1,:)
              sR(:) = q(i,j,k+1,:) - q(i,j,k  ,:)
              sC(:) = minmod(sL(:), sR(:))
              qPy(i,j,k  ,:) = q(i,j,k,:) - sC(:)*0.5_dp_t
              qPy(i,j,k+1,:) = q(i,j,k,:) + sC(:)*0.5_dp_t
           enddo !- i
        enddo !- j
     enddo !- k
     !$omp end do
     
     ! compute fluxes
     !$omp do private(i,j,k)
     do k=lo(3)-1,hi(3)+1
        do j=lo(2)-1,hi(2)+1
           do i=lo(1)-1,hi(1)+1
              call roe(nvar, 1, qPx(i,j,k,:), qMx(i,j,k,:), fx(i,j,k,:))
              call roe(nvar, 2, qPy(i,j,k,:), qMy(i,j,k,:) ,fy(i,j,k,:))
              call roe(nvar, 3, qPz(i,j,k,:), qMz(i,j,k,:) ,fz(i,j,k,:))
           enddo !- i
        enddo !- j
     enddo !- k
     !$omp end do
     !$omp end parallel
     
     deallocate(q,qMx,qPx,qPy,qMy,qPz,qMz)
   
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
     
     do n=1,nlevs
        do nf=1,nfabs(U(n))
           lo = 1; hi = 1
           lo(1:mla%dim) = lwb(get_box(U(n),nf))
           hi(1:mla%dim) = upb(get_box(U(n),nf))
           q => dataptr(U(n),nf)
           
           do k=lo(3),hi(3)
              do j=lo(2),hi(2)
                 do i=lo(1),hi(1)
                    p = press(q(i,j,k,:))
                    c2 = gamma*p/q(i,j,k,1)
                    cfx = sqrt(c2) + abs(q(i,j,k,2))
                    cfy = sqrt(c2) + abs(q(i,j,k,3))
                    cfz = sqrt(c2) + abs(q(i,j,k,4))
                    max_speed = max(max_speed, cfx, cfy, cfz)
                 enddo !- i
              enddo !- j
           enddo !- k
        enddo !- nf
     enddo !- n
     
     dt = CFL * dx(1) / max_speed / real(dm,8)
     
   end subroutine get_new_timestep
   

end module advance_module

