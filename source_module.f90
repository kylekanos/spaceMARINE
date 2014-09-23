module source_module
   use multifab_module
   use ml_layout_module
   use define_bc_module
   use bc_module
   use multifab_physbc_module
   use multifab_fill_ghost_module
   use ml_restriction_module
   use physics_declarations
   
   private
   public :: source_box, source_init
   
   type TableDef
      integer :: nPoints=1
      real(dp_t) :: xmin, xmax, dx
      real(dp_t), pointer :: data(:)
      integer :: bc(2)
      logical :: Initialized=.false.
   end type TableDef
   type(TableDef), save :: CoolTable
   logical :: lCooling = .false.
   logical :: lGravity = .false.
   logical :: lSources = .false.
 contains
   !> @brief initializes sources
   subroutine source_init
      integer :: ierr, i
     open(unit=10,file='cooling.tab',iostat=ierr)
     if(ierr/=0) then
        print *,"Error: unable to find cooling.tab, no cooling run enabled"
        lCooling = .false.
        lSources = .false.
        return
     endif
     print *,"table found, reading..."
     lCooling = .true.
   ! read number of points
     read(10,*) CoolTable%nPoints
     allocate(CoolTable%data(CoolTable%nPoints))
   ! read size of table, step between & boundary conditions
     read(10,*) CoolTable%xmin, CoolTable%xmax, CoolTable%dx, CoolTable%bc(1), CoolTable%bc(2)
   ! read in cooling rate
     do i=1,CoolTable%nPoints
        read(10,*) CoolTable%data(i)
     enddo !- i
     close(10)
     
     lSources = lCooling .or. lGravity
   end subroutine source_init

   !> @brief controls sources
   subroutine source_box(mla, dt, dxn, tower, phi)
     real(dp_t), intent(in) :: dt, dxn(:)
     type(ml_layout), intent(in) :: mla
     type(bc_tower) , intent(in) :: tower
     type(multifab), intent(in) :: phi(:)
   ! locals
     real(dp_t), pointer :: q(:,:,:,:)
     real(dp_t), allocatable :: p(:)
     real(dp_t) :: T, S
     integer :: dm, i, j, k, lo(3), hi(3), n, nlevs, nf, nvars
     
     if(.not.lSources) return
     
     nlevs = mla%nlevel
     dm = mla%dim
     nvars = ncomp(phi(1))
     
     allocate(p(nvars))
     
     do n=1,nlevs
        do nf=1,nfabs(phi(n))
           q => dataptr(phi(n),nf)
           lo=1; hi=1
           lo(1:dm) = lwb(get_box(phi(n),i))
           hi(1:dm) = upb(get_box(phi(n),i))
         ! loop over domain for updating sources
           do k=lo(3),hi(3)
              do j=lo(2),hi(2)
                 do i=lo(1),hi(1)
                    p(:) = primitive(q(i,j,k,:))
                    T = Temperature(p)
                    if(lCooling) then
                       call cooling(p, T, S)
                    endif
                    q(i,j,k,5) = q(i,j,k,5)  - dt*S
                 enddo !- i
              enddo !- j
           enddo !- k
        enddo !- nf
     enddo !- n
   end subroutine source_box
   
   !> @brief interpolates cooling function
   subroutine cooling(p, T, S)
     real(dp_t), intent(in) :: p(:), T
     real(dp_t), intent(inout) :: S
     real(dp_t) :: tpos, CoolRate
     integer :: it1, it2
     
     if(T > CoolTable%xmin) then
        if(T < CoolTable%xmax) then
           tpos = (log10(T) - CoolTable%xmin)/CoolTable%dx + 1.0
           it1 = min( CoolTable%nPoints - 1, max(1, int(tpos)))
           it2 = min(CoolTable%nPoints, it1+1)
           CoolRate = CoolTable%data(it1) + (CoolTable%data(it2) - CoolTable%data(it1))*CoolTable%dx
           S = S + p(1)**2*(10.0_dp_t**CoolRate)*ScaleCool
        endif
     endif
   end subroutine cooling
   
   !> @brief computes gravitational source
   subroutine  gravity(p, S)
     real(dp_t), intent(in) :: p(:)
     real(dp_t), intent(inout) :: S
   end subroutine gravity
   
   !> @brief convert primitive variables to physical temperature
   pure real(dp_t) function Temperature(p)
     real(dp_t), intent(in) :: p(:)
      Temperature=TempScale*p(5)/p(1)
   end function Temperature
end module source_module
