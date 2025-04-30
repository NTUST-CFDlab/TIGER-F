! 22 Aug 2023 - FDS

 !    only cpu 0 answer correct
! subroutine virtualForceIntegrator()
	! use variables
	! implicit none
   ! integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   ! !---------------------------------------------------!
   ! !    LOCAL VARIABLES                                !
   ! !---------------------------------------------------!

	! iBgnVOS = 1
	! iEndVOS = nx
	! jBgnVOS = 1
	! jEndVOS = ny
	! kBgnVOS = 1
	! kEndVOS = nz




	! !Integrating components in z-direction using Simpson's 1/3 rule
	! do j=jBgnVOS,jEndVOS ; do i= iBgnVOS,iEndVOS

		! FXz(i,j) = 0.0
		! FYz(i,j) = 0.0
		
		! do k= kBgnVOS,kEndVOS

			
				! FXz(i,j) = FXz(i,j) + ( FZ(i,j,k-1)*Dzs(k-1) + 4.0*FZ(i,j,k)*Dzs(k) + FZ(i,j,k+1)*Dzs(k+1) ) / 6.0
				! FYz(i,j) = FYz(i,j) + ( FY(i,j,k-1)*Dzs(k-1) + 4.0*FY(i,j,k)*Dzs(k) + FY(i,j,k+1)*Dzs(k+1) ) / 6.0
			
			
		! end do

	! end do; end do



	! !Integrating components in y-direction using Simpson's 1/3 rule
	! do i= iBgnVOS,iEndVOS

		! FXy(i) = 0.0
		! FYy(i) = 0.0
		
		! do j=jBgnVOS,jEndVOS

			! if (FXz(i,j) /= 0.0 ) then
				! FXy(i) = FXy(i) + (FXz(i,j-1)*Dys(j-1) + 4.0*FXz(i,j)*Dys(j) + FXz(i,j+1)*Dys(j+1)) / 6.0
			! end if
			
			! if (FYz(i,j) /= 0.0 ) then
				! FYy(i) = FYy(i) + (FYz(i,j-1)*Dys(j-1) + 4.0*FYz(i,j)*Dys(j) + FYz(i,j+1)*Dys(j+1)) / 6.0
			! end if
			
		! end do
		
	! end do



	! !Integrating components in x-direction using Simpson's 1/3 rule
	! totalFX = 0.0
	! totalFY = 0.0
	! do i= iBgnVOS,iEndVOS

		! if (FXy(i) /= 0.0 ) then
			! totalFX = totalFX + (FXy(i-1)*Dxs(i-1) + 4.0*FXy(i)*Dxs(i) + FXy(i+1)*Dxs(i+1) ) / 6.0
		! end if
		
		! if (FYy(i) /= 0.0 ) then
			! totalFY = totalFy + (FYy(i-1)*Dxs(i-1) + 4.0*FYy(i)*Dxs(i) + FYy(i+1)*Dxs(i+1) ) / 6.0
		! end if

	! end do




	! cDrag = (-2.d0) * totalFX / (ref_area*U_inf*U_inf*1.d0) 
	! cLift = (-2.d0) * totalFY / (ref_area*U_inf*U_inf*1.d0) 



! end subroutine virtualForceIntegrator





subroutine virtualForceIntegrator_nima()
	use variables
	implicit none

	totalFX = 0.0
	totalFY = 0.0

  !$OMP PARALLEL
   !$OMP DO PRIVATE(j,i)  reduction( +: totalFX)
   do k=1,nz;do j=1,ny;do i=1,nx

	totalFX = totalFX + ( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )

   end do;end do;end do
   !$OMP END DO


   !$OMP DO PRIVATE(j,i) reduction( +: totalFY)
   do k=1,nz;do j=1,ny;do i=1,nx

	totalFY = totalFY + ( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )

   end do;end do;end do
   !$OMP END DO
  !$OMP END PARALLEL
    
	cDrag = (-2.d0) * totalFX / (ref_area*U_inf*U_inf*1.d0)  ! because for dimensional: (-2.d0)*totalFX*den_flu / (den_flu*ref_area*U_inf*U_inf) 
	cLift = (-2.d0) * totalFY / (ref_area*U_inf*U_inf*1.d0)

	!cDrag = (-2.0) * totalFX                 
	!cLift = (-2.0) * totalFY


end subroutine virtualForceIntegrator_nima


subroutine virtualForceTorqueIntegrator()
	use variables
	implicit none

	totalFX = 0.d0
	totalFY = 0.d0
    totalTX = 0.d0
	totalTY = 0.d0


   !$OMP PARALLEL	
	!$OMP DO PRIVATE(j,i)  reduction( +: totalFX)
	do k=1,nz;do j=1,ny;do i=1,nx

		totalFX = totalFX + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )

	end do;end do;end do
	!$OMP END DO

    

    !$OMP DO PRIVATE(j,i)  reduction( +: totalTX)
    do k=1,nz;do j=1,ny;do i=1,nx
        

		totalTX = totalTX + den_flu*(-FZ(i,j,k)*(ys(j)-yc )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    


    end do;end do;end do
    !$OMP END DO



	!$OMP DO PRIVATE(j,i) reduction( +: totalFY)
	do k=1,nz;do j=1,ny;do i=1,nx

		totalFY = totalFY + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )

	end do;end do;end do
	!$OMP END DO

            
     
	!$OMP DO PRIVATE(j,i)  reduction( +: totalTY)
    do k=1,nz;do j=1,ny;do i=1,nx
	 

		totalTY = totalTY + den_flu*(FY(i,j,k)*(zs(k)-zc )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0
     

    end do;end do;end do
    !$OMP END DO
   !$OMP END PARALLEL
	
	
	
	totalTorq =   (totalTX + totalTY)
	
	cDrag = -2.d0 * totalFX / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
	cLift = -2.d0 * totalFY / (ref_area*U_inf*U_inf*den_flu*1.d0)
    cTorq = -2.d0 * totalTorq /(2.d0*(rotor_r+blade_r)*U_inf*U_inf*den_flu*rotor_r*1.d0)     ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																				 ! Positive torque = CCW
    
	 
	cPower   =   cTorq*dabs(rotor_tsr)
	

end subroutine virtualForceTorqueIntegrator
