module solvers
   use physics_declarations
   implicit none
   integer, parameter :: dp = kind(1d0)
   integer, dimension(:,:), allocatable, save :: ivar
   
   private
   
   public :: solver_init, minmod, limiter, hlle, roe, press
contains
   !> @brief initialize quick thing here
   subroutine solver_init(nvar)
     integer, intent(in) :: nvar
     integer :: i
     
     allocate(ivar(nvar,3))
     ivar(:,1) = [(i,i=1,nvar)]
     ivar(:,2) = ivar(:,1); ivar(1:5,2) = [1,3,4,2,5]
     ivar(:,3) = ivar(:,1); ivar(1:5,3) = [1,4,2,3,5]
   end subroutine solver_init
   
   !> @brief minmod limiter
   elemental real(dp) function minmod(a,b)
     real(dp), intent(in) :: a,b
     minmod = (sign(1._dp,a)+sign(1._dp,b))*min(abs(a),abs(b))/2._dp
   end function minmod
   
   !> @brief 3rd order limiter function
   elemental real(dp) function limiter(dvp, dvm, dx)
     real(dp), intent(in) :: dvp, dvm, dx
     real(dp) :: r, a, b, c, q, th, eta, psi, eps
     
     r = 0.1d0; eps = 1d-12
     th = dvm/(dvp + 1d-16)
     q = (2d0 + th)*0.3333333333333333333d0
     
     a = min(1.5d0, 2d0*th)
     a = min(q, a)
     b = max(-half*th, a)
     c = min(q, b)
     psi = max(zero, c)
     
     eta = (dvm**2 + dvp**2)/(r*dx)**2
     if(eta <= one - eps) then
        limiter = q
     else if(eta >= one + eps) then
        limiter = psi
     else
        limiter = half*((one - (eta - one)/eps)*q + (one + (eta - one)/eps)*psi)
     endif
   end function limiter
   
   pure function press(q) result (p)
     real(dp), intent(in) :: q(:)
     real(dp) :: p
     p = (gamma-1d0)*(q(5) - half*sum(q(2:4)**2/q(1)))
   end function press
   
   !> @brief roe solver
   subroutine roe(nvar, ic, ucl, ucr, f)
     integer, intent(in) :: ic, nvar
     real(dp), dimension(nvar), intent(in ) :: ucl,ucr
     real(dp), dimension(nvar), intent(out) :: f
     real(dp)                               :: dl,ul,vl,vvl,pl,dr,ur,vr,vvr,pr,vc,g1
     real(dp)                               :: sqrtdl,sqrtdr,vxroe,vyroe,vzroe,hroe
     real(dp)                               :: el,er,sqrtd
     real(dp), dimension(nvar     )         :: ulocl,ulocr,lambda,a,fl,fr,ff
     real(dp), dimension(nvar,nvar)         :: lem,rem
     integer                                :: n,h
     
     do n = 1,nvar
        ulocl(n) = ucl(ivar(n,ic))
        ulocr(n) = ucr(ivar(n,ic))
     enddo
     g1 = gamma - one
   
   ! Convert to primitive quantities
     dl  = ulocl(1)
     ul  = ulocl(2)!/ulocl(1)
     vl  = ulocl(3)!/ulocl(1)
     vvl = ulocl(4)!/ulocl(1)
     pl  = ulocl(5)!g1*(ulocl(5) - half*(ulocl(2)*ulocl(2)+ulocl(3)*ulocl(3)+ulocl(4)*ulocl(4))/ulocl(1))
     el  = pl*gamma7 + half*dl*(ul*ul + vl*vl + vvl*vvl)
     
     dr  = ulocr(1)
     ur  = ulocr(2)!/ulocr(1)
     vr  = ulocr(3)!/ulocr(1)
     vvr = ulocr(4)!/ulocr(1)
     pr  = ulocr(5)!g1*(ulocr(5) - half*(ulocr(2)*ulocr(2)+ulocr(3)*ulocr(3)+ulocr(4)*ulocr(4))/ulocr(1))
     er  = pr*gamma7 + half*dr*(ur*ur + vr*vr + vvr*vvr)
     
   ! Step 1 : Compute Roe-averaged data from left and right states
   !   These averages will be the input variables to the eigen problem
     sqrtdl = sqrt(dl)
     sqrtdr = sqrt(dr)
     sqrtd  = one/(sqrtdl + sqrtdr)
     vxroe  = (sqrtdl*ul  + sqrtdr*ur )*sqrtd
     vyroe  = (sqrtdl*vl  + sqrtdr*vr )*sqrtd
     vzroe  = (sqrtdl*vvl + sqrtdr*vvr)*sqrtd
     hroe   = ((el+pl)/sqrtdl+(er+pr)/sqrtdr)*sqrtd
   
     ! Step 2 : Compute eigenvalues and eigenmatrices from Roe-averaged values
     call eigen_cons(nvar,vxroe,vyroe,vzroe,hroe,lambda,rem,lem)
   
     ! Step 3 : Create intermediate states from eigenmatrices
     a(1:nvar) = (ulocr(1)-ulocl(1))*lem(1,1:nvar)
     do n = 2,nvar
        a(1:nvar) = a(1:nvar) + (ulocr(n)-ulocl(n))*lem(n,1:nvar)
     enddo
   
     ! Step 4 : Compute L/R fluxes
     !  These are computed from the left and right input primitive variables
   
     fl(1) = ulocl(2)
     fr(1) = ulocr(2)
   
     fl(2) = ulocl(2)*ul + pl
     fr(2) = ulocr(2)*ur + pr
   
     fl(3) = ulocl(2)*vl
     fr(3) = ulocr(2)*vr
   
     fl(4) = ulocl(2)*vvl
     fr(4) = ulocr(2)*vvr
   
     fl(5) = (ulocl(5)+pl)*ul
     fr(5) = (ulocr(5)+pr)*ur
   
     ! Compute Roe intermediate states
     do n = 1,nvar
        f(ivar(n,ic)) = half*(fl(n)+fr(n))
     enddo
     ! now add in eignevalue decomposition...
     do n = 1,nvar
        do h = 1,nvar
           f(ivar(h,ic)) = f(ivar(h,ic)) - half*abs(lambda(n))*a(n)*rem(n,h)
        enddo
     enddo
   
   end subroutine roe
   
   !> @brief solve
   subroutine hlle(nvar, ic, ucl,ucr,f)
    integer                  , intent(in ) :: nvar,ic
    real(dp), dimension(nvar), intent(in ) :: ucl,ucr
    real(dp), dimension(nvar), intent(out) :: f
    real(dp)                               :: dl,ul,vl,vvl,pl,dr,ur,vr,vvr,pr
    real(dp)                               :: vc,bp,bm,cfl,cfr,g1
    real(dp)                               :: sqrtdl,sqrtdr,vxroe,vyroe,vzroe,hroe
    real(dp), dimension(nvar     )         :: ulocl,ulocr,lambda,a,fl,fr,ff
    real(dp), dimension(nvar,nvar)         :: lem,rem
    integer                                :: n
   
    do n = 1,nvar
       ulocl(n) = ucl(ivar(n,ic))
       ulocr(n) = ucr(ivar(n,ic))
    enddo
    g1 = gamma - one
    
    ! Convert to primitive quantities
    dl  = ulocl(1)
    ul  = ulocl(2)/ulocl(1)
    vl  = ulocl(3)/ulocl(1)
    vvl = ulocl(4)/ulocl(1)
    pl  = g1*(ulocl(5) - half*(ulocl(2)*ulocl(2)+ulocl(3)*ulocl(3)+ulocl(4)*ulocl(4))/ulocl(1))
    
    dr  = ulocr(1)
    ur  = ulocr(2)/ulocr(1)
    vr  = ulocr(3)/ulocr(1)
    vvr = ulocr(4)/ulocr(1)
    pr  = g1*(ulocr(5) - half*(ulocr(2)*ulocr(2)+ulocr(3)*ulocr(3)+ulocr(4)*ulocr(4))/ulocr(1))
   
    ! Step 1 : Compute Roe-averaged data from left and right states
    !   These averages will be the input variables to the eigen problem
    sqrtdl = sqrt(dl)
    sqrtdr = sqrt(dr)
    vxroe  = (sqrtdl*ul  + sqrtdr*ur )/(sqrtdl+sqrtdr)
    vyroe  = (sqrtdl*vl  + sqrtdr*vr )/(sqrtdl+sqrtdr)
    vzroe  = (sqrtdl*vvl + sqrtdr*vvr)/(sqrtdl+sqrtdr)
    hroe   = ((ulocl(5)+pl)/sqrtdl+(ulocr(5)+pr)/sqrtdr)/(sqrtdl+sqrtdr)
   
    ! Step 2 : Compute eigenvalues and eigenmatrices from Roe-averaged values
    call eigen_cons(nvar,vxroe,vyroe,vzroe,hroe,lambda,rem,lem)
   
    ! Step 3 : Create intermediate states from eigenmatrices
    a(:) = (ulocr(1)-ulocl(1))*lem(1,:)
    do n = 2,nvar
       a(:) = a(:) + (ulocr(n)-ulocl(n))*lem(n,:)
    enddo
   
    ! Step 4 : Compute L/R fluxes
    !  These are computed from the left and right input primitive variables
   
    fl(1) = ulocl(2)
    fr(1) = ulocr(2)
   
    fl(2) = ulocl(2)*ul + pl
    fr(2) = ulocr(2)*ur + pr
   
    fl(3) = ulocl(2)*vl
    fr(3) = ulocr(2)*vr
   
    fl(4) = ulocl(2)*vvl
    fr(4) = ulocr(2)*vvr
   
    fl(5) = (ulocl(5)+pl)*ul
    fr(5) = (ulocr(5)+pr)*ur
   
    ! Compute HLLE wave speeds
    cfl = sqrt(gamma*pl/dl)
    cfr = sqrt(gamma*pr/dr)
    bp  = max( max(lambda(5), (ur+cfr)), zero )
    bm  = min( min(lambda(1), (ul-cfl)), zero )
   
    do n = 1,nvar
       f(ivar(n,ic)) = ((bp*fl(n)-bm*fr(n))+(bp*bm)*(ulocr(n)-ulocl(n)))/(bp-bm)
    enddo
   
   end subroutine hlle
   
   
   !> @brief solve the eigenvalue problem
   subroutine eigen_cons(nvar,vx,vy,vz,hr,lambda,rem,lem)
     integer, intent(in)            :: nvar
     real(dp)                       :: vx,vy,vz,hr,cs,vsq,norm,g1
     real(dp), dimension(nvar     ) :: lambda
     real(dp), dimension(nvar,nvar) :: rem,lem
   
     vsq  = vx*vx + vy*vy + vz*vz
     cs   = sqrt(g1*max((hr-half*vsq),epsilon(vsq)))
     norm = half/(cs*cs)
     g1 = gamma - one
   
     ! eigenvalues
     lambda(1) = vx - cs
     lambda(2) = vx
     lambda(3) = vx
     lambda(4) = vx
     lambda(5) = vx + cs
   
     ! right eigenmatrix  
     rem(1,1) = one
     rem(1,2) = vx - cs
     rem(1,3) = vy
     rem(1,4) = vz
     rem(1,5) = hr - vx*cs
   
     rem(2,1) = zero
     rem(2,2) = zero
     rem(2,3) = one
     rem(2,4) = zero
     rem(2,5) = vy
   
     rem(3,1) = zero
     rem(3,2) = zero
     rem(3,3) = zero
     rem(3,4) = one
     rem(3,5) = vz
   
     rem(4,1) = one
     rem(4,2) = vx
     rem(4,3) = vy
     rem(4,4) = vz
     rem(4,5) = half * vsq
   
     rem(5,1) = one
     rem(5,2) = vx + cs
     rem(5,3) = vy
     rem(5,4) = vz
     rem(5,5) = hr + vx*cs
   
     ! left eigenmatrix
     lem(1,1) =  norm*(half*g1*vsq + vx*cs)
     lem(2,1) = -norm*(g1*vx + cs)
     lem(3,1) = -norm*g1*vy
     lem(4,1) = -norm*g1*vz
     lem(5,1) =  norm*g1
   
     lem(1,2) = -vy
     lem(2,2) =  zero
     lem(3,2) =  one
     lem(4,2) =  zero
     lem(5,2) =  zero
   
     lem(1,3) = -vz
     lem(2,3) =  zero
     lem(3,3) =  zero
     lem(4,3) =  one
     lem(5,3) =  zero
   
     lem(1,4) =  one-norm*g1*vsq
     lem(2,4) =  g1*vx/(cs*cs)
     lem(3,4) =  g1*vy/(cs*cs)
     lem(4,4) =  g1*vz/(cs*cs)
     lem(5,4) = -g1/(cs*cs)
   
     lem(1,5) =  norm*(half*g1*vsq - vx*cs)
     lem(2,5) = -norm*(g1*vx - cs)
     lem(3,5) =  lem(3,1)
     lem(4,5) =  lem(4,1)
     lem(5,5) =  lem(5,1)
   end subroutine eigen_cons
end module solvers
