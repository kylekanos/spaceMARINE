module solvers
   use physics_declarations
   implicit none
   integer, dimension(:,:), allocatable, save :: ivar
   
   private
   
   public :: Energy, slope_lim, press, solver, swapy
contains
   !> @brief minmod limiter
   elemental real(dp_t) function minmod(a,b)
     real(dp_t), intent(in) :: a,b
     minmod = (sign(half,a)+sign(half,b))*min(abs(a), abs(b))
   end function minmod
   
   !> @brief impose slope limiter
   elemental function slope_lim(y1, y2, y3) result(sl)
     real(dp_t), intent(in) :: y1, y2, y3
     real(dp_t) :: dqm, dqp, dqc, s, sl
     real(dp_t) :: adqm, adqp, adqc
     

     dqm = two*(y2 - y1)			! (mid - left)
     dqp = two*(y3 - y2)			! (right - mid)
     dqc = half*(dqm + dqp)			! (left - right)
     s = dqm*dqp
     adqm = abs(dqm)
     adqp = abs(dqp)
     adqc = abs(dqc)
     if(s <= zero) then
        sl = 0d0
     else
        if(adqm < adqp .and. adqm < adqc) then
           sl = dqm
        elseif(adqp < adqc) then
           sl = dqp
        else
           sl = dqc
        endif
     endif
   end function slope_lim

   !> @brief HLL solver (3-wave state, MHD characteristiscs)
   function solver(NP,pl,pr,eta) result(flux)
     integer, intent(in) :: NP
     real(dp_t), intent(in) :: pl(1:NP), pr(1:NP)
     real(dp_t), optional :: eta
     real(dp_t) :: flux(1:NP)
     real(dp_t), allocatable, dimension(:) :: ur, ul, fl, fr
     real(dp_t) :: csl, csr, sl, sr, den
     integer :: n

     allocate(ur(1:NP),ul(1:NP),fl(1:NP),fr(1:NP))
     
   ! speed of sounds
     csl = fast_wavespeed(pl(1:NP))
     csr = fast_wavespeed(pr(1:NP))
     flux = zero
     
     
   ! find max characteristic speed
     sr = max(pl(2) + csl, pr(2) + csr)
     sl = min(pl(2) - csl, pr(2) - csr)
     if(present(eta)) then
        sl = sign(abs(sl)+eta, sl)
        sr = sign(abs(sr)+eta, sr)
     endif

     fl(1:NP) = prim_flux(pl(1:NP))
     fr(1:NP) = prim_flux(pr(1:NP))
!~      ul(1:NP) = conservative(pl(1:NP))
!~      ur(1:NP) = conservative(pr(1:NP))
!~      print '(a,5(f12.4,1x))', 'ul=',ul
!~      print '(a,5(f12.4,1x))', 'ur=',ur
!~      print '(2(a,1x,f16.8))','sl=',sl,', sr=',sr
     
   ! check states
     if(sl > zero) then
        flux(1:NP) = fl(1:NP)
     elseif(sr < zero) then
        flux(1:NP) = fr(1:NP)
     else
        ul(1:NP) = conservative(pl(1:NP))
        ur(1:NP) = conservative(pr(1:NP))
        den = one / (sl - sr + 1e-30_dp_t)
        do n=1,NP
           flux(n) = ( sl*sr*(ur(n) - ul(n)) + sr*fl(n) - sl*fr(n) ) * den
        enddo !- n
     endif
   end function solver

   !> @brief find the fast hydro wave speed from primitives
!DEC$ ATTRIBUTES INLINE::fast_wavespeed
   pure function fast_wavespeed(w) result(fast)
     real(dp_t), intent(in) :: w(:)
     real(dp_t) :: fast
     fast = sqrt(gamma*w(5)/w(1))
   end function fast_wavespeed

    !> @brief swap x & y terms
   subroutine swapy(v)
     real(dp_t), intent(inout) :: v(:)
     real(dp_t) :: aux
     
   ! velocity 1=y, 2=x
     aux  = v(2)
     v(2) = v(3)
     v(3) = aux
   end subroutine swapy

   !> @brief swap x, y & z terms
   subroutine swapxyz(v)
     real(dp_t), intent(inout) :: v(:)
     real(dp_t) :: aux(2)
     
   ! velocity 1=z, 2=x, 3=y
     aux(:) = v(2:3)
     v(2)   = v(4)
     v(3)   = aux(2)
     v(4)   = aux(1)
   end subroutine swapxyz

   !> @brief swap x, y & z terms
   subroutine swapzxy(v)
     real(dp_t), intent(inout) :: v(:)
     real(dp_t) :: aux(3)
     
   ! velocity 1=x, 2=y, 3=z
     aux(1:3) = v(2:4)
     v(2)   = aux(2)
     v(3)   = aux(3)
     v(4)   = aux(1)
   end subroutine swapzxy

    !> @brief compute flux based on primitive input 
   pure function prim_flux(w) result(fx)
     real(dp_t), intent(in) :: w(:)
     real(dp_t) :: fx(1:size(w))
     real(dp_t) :: p, E
     
     p = w(5)
     E = energy(w)
     
     fx(1)  =  w(1)*w(2)			! rho*vx
     fx(2)  = fx(1)*w(2) + p		! rho*vx*vx + p
     fx(3)  = fx(1)*w(3)  			! rho*vx*vy
     fx(4)  = fx(1)*w(4)  			! rho*vx*vz
     fx(5)  = (E + p)*w(2)			! (E + p)*vx
   end function prim_flux

   !> @brief compute fluid energy from primitives
   pure function energy(p) result(E)
     real(dp_t), intent(in) :: p(:)
     real(dp_t) :: E
     E = p(5)*gamma7 + half*sum(p(2:4)**2)/p(1)
   end function energy

   !> @brief compute fluid pressure from conservatives
   pure function press(q) result (p)
     real(dp_t), intent(in) :: q(:)
     real(dp_t) :: p
     p = gamma1*(q(5) - half*sum(q(2:4)**2)/q(1))
   end function press
   
end module solvers
