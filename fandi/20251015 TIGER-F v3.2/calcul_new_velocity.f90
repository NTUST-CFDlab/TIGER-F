! 12 Jan 2025 - FDS

subroutine calcul_new_velocity()
   use variables
   implicit none
	real*8 :: ETA_sub


   !-------------------------------------------------!
   !    Calculation of velocity field at t = dt*n+1  !
   !												 !
   !	Select the correct block according to motion !
   !-------------------------------------------------!

! Sending data to accelerator
!$acc data present(Xs,Ys,Zs,Dxs,Dys,Dzs,ETA,p(:,:,istart-2:iend+2), &
!$acc	u0(:,:,istart-2:iend+2),v0(:,:,istart-2:iend+2),w0(:,:,istart-2:iend+2), &
!$acc	u2(:,:,istart-2:iend+2),v2(:,:,istart-2:iend+2),w2(:,:,istart-2:iend+2), &
!$acc	FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2)) &
!$acc create(u1(:,:,istart-2:iend+2),v1(:,:,istart-2:iend+2),w1(:,:,istart-2:iend+2))

   !!! In x direction !!!

  !$OMP PARALLEL
   !$OMP DO PRIVATE(i,j,ETA_sub,u_solid1) collapse(nclps)
  !$acc parallel
   !$acc loop independent private(ETA_sub,u_solid1) collapse(3) gang vector
   do k=istart,iend 
   do j=1,ny;do i=1,nx

      u1(i,j,k) = u0(i,j,k) - dt*(p(i+1,j,k)-p(i,j,k))/den_flu / Dxs(i)
		!u1(i,j,k) = u0(i,j,k) - dt*(p(i+1,j,k)-p(i,j,k))/den_flu / Dxs(i) - du_solid			! moving X reference frame

		ETA_sub = 0.5d0*(ETA(i,j,k)+ETA(i+1,j,k))

	!free-falling sphere-----------------------------------------------------------------
      !u_solid1 = u_solid + (1.d0*rotate_sy)*(Zs(k)-z0_t) + (-1.d0*rotate_sz)*(Ys(j)-y0_t)
      !u_solid1 = u_solid + (1.d0*rotate_sy)*(Zs(k)-z0) + (-1.d0*rotate_sz)*(Ys(j)-y0_t)	! moving Z reference frame
	  !u_solid1 = (1.d0*rotate_sy)*(Zs(k)-z0) + (-1.d0*rotate_sz)*(Ys(j)-y0_t)				! moving XZ reference frame
	!------------------------------------------------------------------------------------
	
	  u_solid1 = u_solid

      u2(i,j,k) = ETA_sub * u_solid1 + (1.d0- ETA_sub) * u1(i,j,k)
      
      FX(i,j,k) = (u2(i,j,k) - u1(i,j,k)) * inv_dt

   end do; end do; end do
   !$OMP END DO

   !!! In y direction !!!

   !$OMP DO PRIVATE(i,j,ETA_sub,v_solid1) collapse(nclps)
   !$acc loop independent private(ETA_sub,v_solid1) collapse(3) gang vector
   do k=istart,iend
   do j=1,ny;do i=1,nx

      v1(i,j,k) = v0(i,j,k) - dt*(p(i,j+1,k)-p(i,j,k))/den_flu / Dys(j)

		ETA_sub = 0.5d0*(ETA(i,j,k)+ETA(i,j+1,k))

	!Set the total imposed and 2way velocity (rotor + 1blade rotation)---------------------------------------------------

	   if( SQRT((Zs(k)-z0_t)**2 + (Ys(j)-y0_t)**2) .LT. (blade_r+dzSml)) then
			 v_solid1= v_solid + rotor_omega * (Zs(k)-z0) + blade_omega * (Zs(k)-z0_t) 	! blade 1, Positive omega is CCW
	   else
			 v_solid1= v_solid + rotor_omega * (Zs(k)-z0)
	   endif
	!--------------------------------------------------------------------------------------------------------------------


	!Set the total imposed and 2way velocity (rotor + 2blade rotation)---------------------------------------------------

	   ! if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (blade_r+dzSml)) then
			 ! v_solid1= v_solid + rotor_omega * (Zs(k)-z0) + blade_omega * (Zs(k)-z0_t1) 	! blade 1, Positive omega is CCW
	   ! else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (blade_r+dzSml)) then
			 ! v_solid1= v_solid + rotor_omega * (Zs(k)-z0) + blade_omega * (Zs(k)-z0_t2) 	! blade 2, Positive omega is CCW 
	   ! else
			 ! v_solid1= v_solid + rotor_omega * (Zs(k)-z0)
	   ! endif
	!--------------------------------------------------------------------------------------------------------------------


	!Set the total imposed and 2way velocity (rotor + 3blade rotation)---------------------------------------------------

	   ! if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (blade_r+dzSml)) then
			 ! v_solid1= v_solid + rotor_omega * (Zs(k)-z0) + blade_omega * (Zs(k)-z0_t1) 	! blade 1, Positive omega is CCW
	   ! else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (blade_r+dzSml)) then
			 ! v_solid1= v_solid + rotor_omega * (Zs(k)-z0) + blade_omega * (Zs(k)-z0_t2) 	! blade 2, Positive omega is CCW 
	   ! else if( SQRT((Zs(k)-z0_t3)**2 + (Ys(j)-y0_t3)**2) .LT. (blade_r+dzSml)) then
			 ! v_solid1= v_solid + rotor_omega * (Zs(k)-z0) + blade_omega * (Zs(k)-z0_t3) 	! blade 3, Positive omega is CCW 
	   ! else
			 ! v_solid1= v_solid + rotor_omega * (Zs(k)-z0)
	   ! endif
	!--------------------------------------------------------------------------------------------------------------------

	!Set the 2way velocity for darrius 3 blade---------------------------------------------------

			 ! v_solid1= v_solid + rotor_omega * (Zs(k)-z0)  	! Positive omega is CCW

	!--------------------------------------------------------------------------------------------
	
	!free-falling sphere-----------------------------------------------------------------
	  !v_solid1= v_solid + (-1.d0*rotate_sx)*(Zs(k)-z0_t) + ( 1.d0*rotate_sz)*(Xs(i)-x0_t)
	  !v_solid1= v_solid + (-1.d0*rotate_sx)*(Zs(k)-z0) + ( 1.d0*rotate_sz)*(Xs(i)-x0_t)	! moving Z reference frame
	  !v_solid1= v_solid + (-1.d0*rotate_sx)*(Zs(k)-z0) + ( 1.d0*rotate_sz)*(Xs(i)-x0)		! moving XZ reference frame
	!------------------------------------------------------------------------------------

      v2(i,j,k) = ETA_sub * v_solid1 + (1.d0- ETA_sub) * v1(i,j,k)

      FY(i,j,k) = (v2(i,j,k) - v1(i,j,k)) * inv_dt

   end do; end do; end do
   !$OMP END DO   


   !!! In w direction !!!

   !$OMP DO PRIVATE(i,j,ETA_sub,w_solid1) collapse(nclps)
   !$acc loop independent private(ETA_sub,w_solid1) collapse(3) gang vector
   do k=istart,iend 
   do j=1,ny;do i=1,nx

      w1(i,j,k) = w0(i,j,k) - dt*(p(i,j,k+1)-p(i,j,k))/den_flu  / Dzs(k)
		!w1(i,j,k) = w0(i,j,k) - dt*(p(i,j,k+1)-p(i,j,k))/den_flu  / Dzs(k) - dw_solid			! moving Z reference frame

		ETA_sub = 0.5d0*(ETA(i,j,k)+ETA(i,j,k+1))

	!Set the total imposed and 2way velocity (rotor + 1blade rotation)---------------------------------------------------

	   if( SQRT((Zs(k)-z0_t)**2 + (Ys(j)-y0_t)**2) .LT. (blade_r+dySml)) then
			 w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0) + (-1.d0 * blade_omega) * (Ys(j)-y0_t)	! blade 1, Positive omega is CCW
	   else
			 w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0)
	   endif
	!--------------------------------------------------------------------------------------------------------------------

	
	!Set the total imposed and 2way velocity (rotor + 2blade rotation)---------------------------------------------------

	   ! if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (blade_r+dySml)) then
			 ! w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0) + (-1.d0 * blade_omega) * (Ys(j)-y0_t1)	! blade 1, Positive omega is CCW
	   ! else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (blade_r+dySml)) then
			 ! w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0) + (-1.d0 * blade_omega) * (Ys(j)-y0_t2)	! blade 2, Positive omega is CCW
	   ! else
			 ! w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0)
	   ! endif
	!--------------------------------------------------------------------------------------------------------------------


	!Set the total imposed and 2way velocity (rotor + 3blade rotation)---------------------------------------------------

	   ! if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (blade_r+dySml)) then
			 ! w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0) + (-1.d0 * blade_omega) * (Ys(j)-y0_t1)	! blade 1, Positive omega is CCW
	   ! else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (blade_r+dySml)) then
			 ! w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0) + (-1.d0 * blade_omega) * (Ys(j)-y0_t2)	! blade 2, Positive omega is CCW
	   ! else if( SQRT((Zs(k)-z0_t3)**2 + (Ys(j)-y0_t3)**2) .LT. (blade_r+dySml)) then
			 ! w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0) + (-1.d0 * blade_omega) * (Ys(j)-y0_t3)	! blade 3, Positive omega is CCW
	   ! else
			 ! w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0)
	   ! endif
	!--------------------------------------------------------------------------------------------------------------------


	!Set the 2way velocity for darrius 3 blade---------------------------------------------------

			 ! w_solid1= w_solid + (-1.d0 * rotor_omega) * (Ys(j)-y0) 	! Positive omega is CCW

	!--------------------------------------------------------------------------------------------

	!free-falling sphere-----------------------------------------------------------------
	!w_solid1= w_solid + (-1.d0*rotate_sy)*(Xs(i)-x0_t) + ( 1.d0*rotate_sx)*(Ys(j)-y0_t)
	!w_solid1= (-1.d0*rotate_sy)*(Xs(i)-x0_t) + ( 1.d0*rotate_sx)*(Ys(j)-y0_t)				! moving Z reference frame
	!w_solid1= (-1.d0*rotate_sy)*(Xs(i)-x0) + ( 1.d0*rotate_sx)*(Ys(j)-y0_t)				! moving XZ reference frame
	!------------------------------------------------------------------------------------
      
      w2(i,j,k) = ETA_sub * w_solid1 + (1.d0- ETA_sub) * w1(i,j,k)
      
      FZ(i,j,k) = (w2(i,j,k) - w1(i,j,k)) * inv_dt

   end do;end do
   end do
   !$acc end parallel 
   !$OMP END DO
  !$OMP END PARALLEL
!$acc end data

end subroutine calcul_new_velocity


subroutine Updating_velocity()
   use variables
   implicit none
   !---------------------------------------------------------!
   !       loops to update the fluid velocity fields         !
   !---------------------------------------------------------!

!$acc data present(u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2), &
!$acc	u2(:,:,istart-2:iend+2),v2(:,:,istart-2:iend+2),w2(:,:,istart-2:iend+2))

   !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps) 
   !$acc parallel loop independent collapse(3) gang vector
   do k=istart,iend  
   do j=1,ny;do i=1,nx

      u(i,j,k) = u2(i,j,k)    		!!! In x direction !!!
      v(i,j,k) = v2(i,j,k) 			!!! In y direction !!!
      w(i,j,k) = w2(i,j,k) 			!!! In z direction !!!

   end do; end do; end do
   !$acc end parallel
   !$OMP END PARALLEL DO
!$acc end data   
   
end subroutine Updating_velocity