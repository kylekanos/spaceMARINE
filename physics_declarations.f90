module physics_declarations
   use bl_types, only: dp_t
   
   public
   
   ! adjust as necessary, be sure it is PHYSICALLY CORRECT!
   real(dp_t), parameter ::  TIMESCALE =   3155692600.00000     , &
                             LSCALE    =  3.085680249999995E+018, &
                             MSCALE    =  6.077143661962447E+031, &
                             RSCALE    =  2.068458160421109E-024, &
                             VELSCALE  =   977813951.206778     , &
                             PSCALE    =  1.977694471123365E-006, &
                             NSCALE    =   1.00000000000000     , &
                             BSCALE    =  4.985222330659412E-003, &
                             TEMPSCALE =   14324369256.4538     , &
                             SCALEGRAV =  1.374560173618349E-012, &
                             SCALEC    =   30.6594580318688     , &
                             SCALECOOL =  1.595642120700027E+015
   real(dp_t), parameter :: gamma  = 1.6666666666666666666667_dp_t,&
                            gamma1 = gamma - 1.0_dp_t,&
                            gamma7 = 1.0_dp_t/gamma1, &
                            half = 0.5_dp_t, &
                            one = 1.0_dp_t, &
                            zero = 0.0_dp_t
 contains
 
   !> @brief convert conservative to primitive
   pure function primitive(u) result(q)
     real(dp_t), intent(in) :: u(:)
     real(dp_t) :: q(1:size(u))
     q = u
     q(2:4) = u(2:4)/u(1)
     q(5) = gamma7*(u(5)-half*SUM(u(2:4)**2)/u(1))
   end function primitive
   
   !> @brief convert primitive to conservative
   pure function conservative(q) result(u)
     real(dp_t), intent(in) :: q(:)
     real(dp_t) :: u(1:size(q))
     u = q
     u(2:4) = q(2:4)*q(1)
     u(5) = q(5)/gamma1 + half*sum(q(2:4)**2)/q(1)
   end function conservative
end module physics_declarations
