! 12 Jan 2025 - FDS

subroutine oneway_coupling() 
    use variables 
    implicit none


!----------------------------------------------------------------------!
!      update the bulk solid velocity and position                     !
!----------------------------------------------------------------------!

	!---------------For spinning cylinder(s)-----------
	! blade_omega = blade_alpha*U_inf/blade_r			! Cylinder spin speed
	
	! sol_speed = U_inf + ABS(blade_alpha)
	!--------------------------------------------------
	

	!---------------For Magnus effect VAWT (1 blade)-----------
	! blade_omega = blade_alpha*U_inf/blade_r			! Cylinder spin speed
	! rotor_omega = rotor_tsr*U_inf/(rotor_r+blade_r)
	! AOA = AOA1 + rotor_omega*( time - StartDynamic_time )*180.d0/PI

	! T_gen = -totalTorq
	! P_gen = T_gen*rotor_omega		
	
	! ! Rotate the 1st blade center with AOA (needed to define the blade rotation)
	! ! 1st blade center initial position for AOA=0 is at Z+ 

		! z0_t= COS((AOA)*PI/180.d0)*(rotor_r)
		! y0_t= SIN((AOA)*PI/180.d0)*(rotor_r)
		
		! y0_t = y0_t + y0
		! z0_t = z0_t + z0

	  ! sol_speed = ABS(blade_alpha)	
	
	!----------------------------------------------------------



	!---------------For Magnus effect VAWT (2 blades)-----------
	! blade_omega = blade_alpha*U_inf/blade_r			! Cylinder spin speed
	! rotor_omega = rotor_tsr*U_inf/(rotor_r+blade_r)
	! AOA = AOA1 + rotor_omega*( time - StartDynamic_time )*180.d0/PI

	! T_gen = -totalTorq
	! P_gen = T_gen*rotor_omega	
	
	! ! Rotate the 1st blade center with AOA (needed to define the blade rotation)
	! ! 1st blade center initial position for AOA=0 is at Z+ 

		! z0_t1= COS((AOA)*PI/180.d0)*(rotor_r)
		! y0_t1= SIN((AOA)*PI/180.d0)*(rotor_r)
		
		! y0_t1 = y0_t1 + y0
		! z0_t1 = z0_t1 + z0
		
	! ! Rotate the 2nd blade center with AOA (needed to define the blade rotation)

		! z0_t2= COS((AOA+120.d0)*PI/180.d0)*(rotor_r)
		! y0_t2= SIN((AOA+120.d0)*PI/180.d0)*(rotor_r)
		
		! y0_t2 = y0_t2 + y0
		! z0_t2 = z0_t2 + z0

	  ! sol_speed = ABS(blade_alpha)	
	
	!----------------------------------------------------------


	!---------------For Magnus effect VAWT (3 blades)-----------
	! blade_omega = blade_alpha*U_inf/blade_r			! Cylinder spin speed
	! rotor_omega = rotor_tsr*U_inf/(rotor_r+blade_r)
	! AOA = AOA1 + rotor_omega*( time - StartDynamic_time )*180.d0/PI

	! T_gen = -totalTorq
	! P_gen = T_gen*rotor_omega	
	
	! ! Rotate the 1st blade center with AOA (needed to define the blade rotation)
	! ! 1st blade center initial position for AOA=0 is at Z+ 

		! z0_t1= COS((AOA)*PI/180.d0)*(rotor_r)
		! y0_t1= SIN((AOA)*PI/180.d0)*(rotor_r)
		
		! y0_t1 = y0_t1 + y0
		! z0_t1 = z0_t1 + z0
		
	! ! Rotate the 2nd blade center with AOA (needed to define the blade rotation)

		! z0_t2= COS((AOA+120.d0)*PI/180.d0)*(rotor_r)
		! y0_t2= SIN((AOA+120.d0)*PI/180.d0)*(rotor_r)
		
		! y0_t2 = y0_t2 + y0
		! z0_t2 = z0_t2 + z0

	! ! ! Rotate the 3nd blade center with AOA (needed to define the blade rotation)

		! z0_t3= COS((AOA+240.d0)*PI/180.d0)*(rotor_r)
		! y0_t3= SIN((AOA+240.d0)*PI/180.d0)*(rotor_r)
		
		! y0_t3 = y0_t3 + y0
		! z0_t3 = z0_t3 + z0

	  ! sol_speed = ABS(blade_alpha)	
	
	!----------------------------------------------------------



!----------------------------------------------------------------------!
!          Wall Damping Function (approximated values)                 !
!----------------------------------------------------------------------!	  

	Re_t = MAX(sol_speed*L_ch/nu,3.d0)
    Cf =  (2.d0*LOG10(Re_t)-0.65d0)**(-2.3)

	
end subroutine oneway_coupling





subroutine twoway_coupling() 
    use variables 
    implicit none

real*8		:: dy_dt,dv_dt
real*8		:: k1_u_solid,k2_u_solid,k3_u_solid,k4_u_solid
real*8		:: k1_x_solid,k2_x_solid,k3_x_solid,k4_x_solid
real*8		:: k1_v_solid,k2_v_solid,k3_v_solid,k4_v_solid
real*8		:: k1_y_solid,k2_y_solid,k3_y_solid,k4_y_solid
real*8		:: k1_w_solid,k2_w_solid,k3_w_solid,k4_w_solid
real*8		:: k1_z_solid,k2_z_solid,k3_z_solid,k4_z_solid

real*8		:: dAOA_dt,domega_dt
real*8		:: k1_om,k2_om,k3_om,k4_om
real*8		:: k1_AOA,k2_AOA,k3_AOA,k4_AOA


!******************************************************
! NOTES: Integral of virtual forces are:              !
!                                                     !
! 	totalFX_   , totalFY_   , totalFZ_                !
! 	totalTorqx , totalTorqy , totalTorqz , totalTorq  !
!                                                     !
!	*according to virtualForceIntegrator.f90          !
!******************************************************


!----------------------------------------------------------------------!
!      update the bulk solid velocity and position caused by 2way      !
!----------------------------------------------------------------------!

	!---------- For VIV of a cylinder (Euler)---------
	! v_solid = v_solid + dt*(((-1.d0 * totalFY_) / (den_flu * (den_sol / den_flu) * solid_volume)) - k_spring * (y0_t - y0) - c_damper * v_solid)

	! y0_t = y0_t + v_solid*dt
	
	! sol_speed = sqrt(u_solid**2 + v_solid**2 + (w_solid - U_inf)**2)	
	!-------------------------------------------------



	!---------- For VIV of a cylinder (RK4)-----------

	!! --> dv/dt = f(y0_t,v_solid) = ((-1.d0 * totalFY_) / (den_flu * (den_sol / den_flu) * solid_volume)) - k_spring * (y0_t - y0) - c_damper * v_solid
	!! --> dy/dt = f(v_solid)) = v_solid

		! k1_v_solid = dv_dt(y0_t , v_solid)
		! k1_y_solid = dy_dt(v_solid)
	
		! k2_v_solid = dv_dt(y0_t+0.5d0*k1_y_solid*dt , v_solid+0.5d0*k1_v_solid*dt)
		! k2_y_solid = dy_dt(v_solid+0.5d0*k1_v_solid*dt)

		! k3_v_solid = dv_dt(y0_t+0.5d0*k2_y_solid*dt , v_solid+0.5d0*k2_v_solid*dt)
		! k3_y_solid = dy_dt(v_solid+0.5d0*k2_v_solid*dt)

		! k4_v_solid = dv_dt(y0_t+k3_y_solid*dt , v_solid+k3_v_solid*dt)
		! k4_y_solid = dy_dt(v_solid+k3_v_solid*dt)
	
	! v_solid = v_solid + (1.d0 / 6.d0) * (k1_v_solid + 2.d0 * k2_v_solid + 2.d0 * k3_v_solid + k4_v_solid) * dt	
	! y0_t = y0_t + (1.d0 / 6.d0) * (k1_y_solid + 2.d0 * k2_y_solid + 2.d0 * k3_y_solid + k4_y_solid) * dt
	! !-------------------------------------------------
	
	
	! sol_speed = sqrt(u_solid**2 + v_solid**2 + (w_solid - U_inf)**2)	

	!-------------------------------------------------




	!---------------For free falling sphere (Euler)-----------
		! du_solid = dt*(-totalFX_/den_sol/solid_volume)
	! u_solid = u_solid + du_solid

		! dv_solid = dt*(-totalFY_/den_sol/solid_volume)		
	! v_solid = v_solid + dv_solid
	
		! dw_solid = dt*(-g+(-totalFZ_/den_sol/solid_volume)+(den_flu*g/den_sol))
	! w_solid = w_solid + dw_solid
   
	! rotate_sx = rotate_sx + dt*(-totalTorqx*2.5d0/den_sol/solid_volume/r**2)
	! rotate_sy = rotate_sy + dt*(-totalTorqy*2.5d0/den_sol/solid_volume/r**2)
	! rotate_sz = rotate_sz + dt*(-totalTorqz*2.5d0/den_sol/solid_volume/r**2)


	! x0_t = x0_t + u_solid*dt
	! y0_t = y0_t + v_solid*dt
	! z0_t = z0_t + w_solid*dt
      
	! sol_speed = sqrt(u_solid**2 + v_solid**2 + w_solid**2)
		
	!-------------------------------------------------



	!---------------For free falling sphere (RK4)-----------
		du_solid = dt*(-totalFX_/den_sol/solid_volume)
	u_solid = u_solid + du_solid

	!! --> du/dt = f(Constant) = -totalFX_/den_sol/solid_volume
	!! --> dx/dt = f(u_solid)) = u_solid

		! --> k1_u_solid = f(Constant)
		! --> k1_x_solid = fuvw(u_solid)
	
		! --> k2_u_solid = f(Constant)
		! --> k2_x_solid = fuvw(u_solid+0.5d0*k1_u_solid*dt)

		! --> k3_u_solid = f(Constant)
		! --> k3_x_solid = fuvw(u_solid+0.5d0*k2_u_solid*dt)

		! --> k4_u_solid = f(Constant)
		! --> k4_x_solid = fuvw(u_solid+k3_u_solid*dt)
		
	x0_t = x0_t + (1.d0 / 6.d0) * (6.d0*u_solid + 3.d0*du_solid) * dt	! RK4 final equation


		dv_solid = dt*(-totalFY_/den_sol/solid_volume)		
	v_solid = v_solid + dv_solid	
	y0_t = y0_t + (1.d0 / 6.d0) * (6.d0*v_solid + 3.d0*dv_solid) * dt	! RK4 final equation



		dw_solid = dt*(-g+(-totalFZ_/den_sol/solid_volume)+(den_flu*g/den_sol))
	w_solid = w_solid + dw_solid	
	z0_t = z0_t + (1.d0 / 6.d0) * (6.d0*w_solid + 3.d0*dw_solid) * dt	! RK4 final equation


   
	rotate_sx = rotate_sx + dt*(-totalTorqx*2.5d0/den_sol/solid_volume/r**2)
	rotate_sy = rotate_sy + dt*(-totalTorqy*2.5d0/den_sol/solid_volume/r**2)
	rotate_sz = rotate_sz + dt*(-totalTorqz*2.5d0/den_sol/solid_volume/r**2)
		

	sol_speed = sqrt(u_solid**2 + v_solid**2 + w_solid**2)	
	!-------------------------------------------------



	!---------------For Translation-----------
	! u_solid = u_solid
	! v_solid = v_solid - totalFY_*dt/den_sol/solid_volume
	! w_solid = w_solid

	! x0_t = x0_t + u_solid*dt
	! y0_t = y0_t + v_solid*dt
	! z0_t = z0_t + w_solid*dt
      
	! sol_speed = sqrt(u_solid**2 + v_solid**2 + (w_solid - U_inf)**2)

	!-------------------------------------------------
   



	!------------------------For Magnus effect VAWT (2 blades)--------------------------
	! blade_omega = blade_alpha*U_inf/blade_r			! Cylinder spin speed


	! !For constant T_gen (Euler)---------------------- 
	! ! rotor_omega = rotor_omega + (- totalTorq - T_gen)*dt/i_sol
	! ! AOA = AOA + rotor_omega*dt*180.d0/PI
	! ! P_gen = T_gen*rotor_omega ! for constant T_gen
	! !-------------------------------------------------


	! !For constant P_gen (RK4)--------------------- 
	! ! ! --> d_omega/dt = f(rotor_omega) = (- totalTorq - P_gen/rotor_omega)/i_sol
	! ! ! --> d_AOA/dt = f(rotor_omega)) = rotor_omega*180.d0/PI

		! ! k1_om  = domega_dt(rotor_omega)
		! ! k1_AOA = dAOA_dt(rotor_omega)
	
		! ! k2_om  = domega_dt(rotor_omega+0.5d0*k1_om*dt)
		! ! k2_AOA = dAOA_dt(rotor_omega+0.5d0*k1_om*dt)

		! ! k3_om  = domega_dt(rotor_omega+0.5d0*k2_om*dt)
		! ! k3_AOA = dAOA_dt(rotor_omega+0.5d0*k2_om*dt)

		! ! k4_om  = domega_dt(rotor_omega+k3_om*dt)
		! ! k4_AOA = dAOA_dt(rotor_omega+k3_om*dt)

	! ! rotor_omega = rotor_omega + (1.d0 / 6.d0) * (k1_om + 2.d0 * k2_om + 2.d0 * k3_om + k4_om) * dt	
	! ! AOA = AOA + (1.d0 / 6.d0) * (k1_AOA + 2.d0 * k2_AOA + 2.d0 * k3_AOA + k4_AOA) * dt

	! ! T_gen = P_gen/rotor_omega
	! !-------------------------------------------------



	! ! Rotate the 1st blade center with AOA (needed to define the blade rotation)
	! ! 1st blade center initial position for AOA=0 is at Z+ 

		! z0_t1= COS((AOA)*PI/180.d0)*(rotor_r)
		! y0_t1= SIN((AOA)*PI/180.d0)*(rotor_r)
		
		! y0_t1 = y0_t1 + y0
		! z0_t1 = z0_t1 + z0
		
	! ! Rotate the 2nd blade center with AOA (needed to define the blade rotation)

		! z0_t2= COS((AOA+120.d0)*PI/180.d0)*(rotor_r)
		! y0_t2= SIN((AOA+120.d0)*PI/180.d0)*(rotor_r)
		
		! y0_t2 = y0_t2 + y0
		! z0_t2 = z0_t2 + z0
		
	  ! sol_speed = ABS(blade_alpha)
	!-----------------------------------------------------------------------------------





	!---------------------For Magnus effect VAWT (3 blades)-----------------------------
	! blade_omega = blade_alpha*U_inf/blade_r			! Cylinder spin speed


	! !For constant T_gen (Euler)---------------------- 
	! ! rotor_omega = rotor_omega + (- totalTorq - T_gen)*dt/i_sol
	! ! AOA = AOA + rotor_omega*dt*180.d0/PI
	! ! P_gen = T_gen*rotor_omega ! for constant T_gen
	! !-------------------------------------------------


	! !For constant P_gen (RK4)--------------------- 
	! ! --> d_omega/dt = f(rotor_omega) = (- totalTorq - P_gen/rotor_omega)/i_sol
	! ! --> d_AOA/dt = f(rotor_omega)) = rotor_omega*180.d0/PI

		! k1_om  = domega_dt(rotor_omega)
		! k1_AOA = dAOA_dt(rotor_omega)
	
		! k2_om  = domega_dt(rotor_omega+0.5d0*k1_om*dt)
		! k2_AOA = dAOA_dt(rotor_omega+0.5d0*k1_om*dt)

		! k3_om  = domega_dt(rotor_omega+0.5d0*k2_om*dt)
		! k3_AOA = dAOA_dt(rotor_omega+0.5d0*k2_om*dt)

		! k4_om  = domega_dt(rotor_omega+k3_om*dt)
		! k4_AOA = dAOA_dt(rotor_omega+k3_om*dt)

	! rotor_omega = rotor_omega + (1.d0 / 6.d0) * (k1_om + 2.d0 * k2_om + 2.d0 * k3_om + k4_om) * dt	
	! AOA = AOA + (1.d0 / 6.d0) * (k1_AOA + 2.d0 * k2_AOA + 2.d0 * k3_AOA + k4_AOA) * dt
	
	! T_gen = P_gen/rotor_omega
	! !-------------------------------------------------




	! ! Rotate the 1st blade center with AOA (needed to define the blade rotation)
	! ! 1st blade center initial position for AOA=0 is at Z+ 

		! z0_t1= COS((AOA)*PI/180.d0)*(rotor_r)
		! y0_t1= SIN((AOA)*PI/180.d0)*(rotor_r)
		
		! y0_t1 = y0_t1 + y0
		! z0_t1 = z0_t1 + z0
		
	! ! Rotate the 2nd blade center with AOA (needed to define the blade rotation)

		! z0_t2= COS((AOA+120.d0)*PI/180.d0)*(rotor_r)
		! y0_t2= SIN((AOA+120.d0)*PI/180.d0)*(rotor_r)
		
		! y0_t2 = y0_t2 + y0
		! z0_t2 = z0_t2 + z0

	! ! Rotate the 3nd blade center with AOA (needed to define the blade rotation)

		! z0_t3= COS((AOA+240.d0)*PI/180.d0)*(rotor_r)
		! y0_t3= SIN((AOA+240.d0)*PI/180.d0)*(rotor_r)
		
		! y0_t3 = y0_t3 + y0
		! z0_t3 = z0_t3 + z0
		
	  ! sol_speed = ABS(blade_alpha)
		
	!-----------------------------------------------------------------------------------



	! !-------------------For Darrius VAWT ------------------
	! rotor_omega = rotor_omega - totalTorq*dt/i_sol

	! AOA = AOA + rotor_omega*dt*180.d0/PI

	  ! sol_speed = rotor_omega*rotor_r
	!--------------------------------------------------------


	  
	  
	  
!----------------------------------------------------------------------!
!          Wall Damping Function (approximated values)                 !
!----------------------------------------------------------------------!	  

	Re_t = MAX(sol_speed*L_ch/nu,3.d0)
    Cf =  (2.d0*LOG10(Re_t)-0.65d0)**(-2.3)
	

	
end subroutine twoway_coupling


!--------------------Runge-Kutta functions for VIV of a cylinder (RK4)----------------------------------------
! FUNCTION dv_dt(y0_t_,v_solid_) result(dv_dt_)
! USE variables
! IMPLICIT NONE
! real*8, intent(in) :: y0_t_,v_solid_
! real*8		:: dv_dt_

! dv_dt_ =  ((-1.d0 * totalFY_) / (den_flu * (den_sol / den_flu) * solid_volume)) - k_spring * (y0_t_ - y0) - c_damper * v_solid_

! END FUNCTION


! FUNCTION dy_dt(v_solid_) result(dy_dt_)
! USE variables
! IMPLICIT NONE
! real*8, intent(in) :: v_solid_
! real*8		:: dy_dt_

! dy_dt_ = v_solid_

! END FUNCTION

!-------------------------------------------------------------------------------------------------------------


!--------------------Runge-Kutta functions for Magnus effect VAWT (2 or 3 blades)----------------------------------------
! FUNCTION domega_dt(rotor_omega_) result(domega_dt_)
! USE variables
! IMPLICIT NONE
! real*8, intent(in) :: rotor_omega_
! real*8		:: domega_dt_

! domega_dt_ =  (- totalTorq - P_gen/rotor_omega_)/i_sol

! END FUNCTION


! FUNCTION dAOA_dt(rotor_omega_) result(dAOA_dt_)
! USE variables
! IMPLICIT NONE
! real*8, intent(in) :: rotor_omega_
! real*8		:: dAOA_dt_

! dAOA_dt_ = rotor_omega_*180.d0/PI

! END FUNCTION

!-------------------------------------------------------------------------------------------------------------






! subroutine dynamic_AOA() 
    ! use variables
    ! implicit none

    ! !AOA = AOA1 + AOA_amp * SIN( 2.d0 * reduce_frequency_k * ( time - dt*StartDynamic ) )
    ! AOA = AOA1 + AOA_amp * ( SIN( 2.d0 * reduce_frequency_k * ( time - dt*StartDynamic ) + 1.5d0*PI ) + 1.d0 )      !FROM 0

    ! angular_vel = 2.d0*reduce_frequency_k*AOA_amp*COS( 2.d0 * reduce_frequency_k * ( time - dt*StartDynamic ) + 1.5d0*PI ) / 180.d0 * PI


! end subroutine dynamic_AOA