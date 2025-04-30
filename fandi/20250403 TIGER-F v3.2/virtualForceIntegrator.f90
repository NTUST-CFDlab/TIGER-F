! 12 Jan 2025 - FDS

subroutine virtualForceIntegrator_nima()
	use variables
	implicit none

	totalFZ = 0.d0
	totalFY = 0.d0

!$acc data present(iDx,iDy,iDz,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))

   !$OMP PARALLEL DO PRIVATE(j,i)  reduction(+: totalFZ, totalFY) collapse(nclps)
   !$acc parallel
   !$acc loop independent reduction(+: totalFZ, totalFY) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	totalFZ = totalFZ + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
	totalFY = totalFY + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
	
	end do;end do;end do
   !$acc end parallel 
   !$OMP END PARALLEL DO
!$acc end data 

    call MPI_ALLREDUCE( totalFZ, totalFZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY, totalFY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )



    !if(myid==master) then  
		cDrag = (-2.d0) * totalFZ_ / (den_flu*ref_area*U_inf*U_inf*1.d0)  ! because for dimensional: (-2.d0)*totalFX*den_flu / (den_flu*ref_area*U_inf*U_inf) 
		cLift = (-2.d0) * totalFY_ / (den_flu*ref_area*U_inf*U_inf*1.d0)
	!endif

end subroutine virtualForceIntegrator_nima



subroutine virtualForceTorqueIntegrator()
	use variables
	implicit none

	totalFZ = 0.d0
	totalFY = 0.d0
    totalTZ = 0.d0
	totalTY = 0.d0

!$acc data present(iDx,iDy,iDz,Xs,Ys,Zs,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))

	!$OMP PARALLEL	
	!$OMP DO PRIVATE(j,i)  reduction(+: totalFZ, totalFY) collapse(nclps)
	!$acc parallel
	!$acc loop independent reduction(+: totalFZ, totalFY) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

		totalFZ = totalFZ + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY = totalFY + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )

	end do;end do;end do
	!$OMP END DO

    
	! Calculate torque wrt (z0,y0)

    !$OMP DO PRIVATE(j,i)  reduction(+: totalTZ, totalTY) collapse(nclps)
	!$acc loop independent reduction(+: totalTZ, totalTY) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

		totalTZ = totalTZ + den_flu*(-FZ(i,j,k)*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY = totalTY + den_flu*( FY(i,j,k)*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0

    end do;end do;end do
    !$OMP END DO
	!$OMP END PARALLEL
   !$acc end parallel 
!$acc end data 
	

    call MPI_ALLREDUCE( totalFZ, totalFZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY, totalFY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ, totalTZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY, totalTY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )	

    !if(myid==master) then 	
		totalTorq =   (totalTZ_ + totalTY_)
	
		cDrag = -2.d0 * totalFZ_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift = -2.d0 * totalFY_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq = -2.d0 * totalTorq/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower   =   cTorq*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)
	!endif
	

end subroutine virtualForceTorqueIntegrator


subroutine virtualForceTorqueIntegrator3D()
	use variables
	implicit none

	totalFX = 0.d0
	totalFZ = 0.d0
	totalFY = 0.d0
	
	totalT_XZ = 0.d0
	totalT_XY = 0.d0
	totalT_YZ = 0.d0
	totalT_YX = 0.d0
	totalT_ZY = 0.d0
	totalT_ZX = 0.d0

!$acc data present(iDx,iDy,iDz,Xs,Ys,Zs,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))


	!$OMP PARALLEL	
	!$OMP DO PRIVATE(j,i)  reduction(+: totalFX, totalFZ, totalFY) collapse(nclps)
	!$acc parallel
	!$acc loop independent reduction(+: totalFX, totalFZ, totalFY) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

		totalFX = totalFX + den_flu*( FX(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFZ = totalFZ + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY = totalFY + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )

	end do;end do;end do
	!$OMP END DO

    
	! Calculate torque wrt (z0,y0)

    !$OMP DO PRIVATE(j,i)  reduction(+: totalT_XZ, totalT_XY, totalT_YZ, totalT_YX, totalT_ZY, totalT_ZX) collapse(nclps)
	!$acc loop independent reduction(+: totalT_XZ, totalT_XY, totalT_YZ, totalT_YX, totalT_ZY, totalT_ZX) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx


		totalT_XZ = totalT_XZ + den_flu*( FZ(i,j,k)*(ys(j)-y0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0   
		totalT_XY = totalT_XY + den_flu*(-FY(i,j,k)*(zs(k)-z0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0		
			!totalT_XY = totalT_XY + den_flu*(-FY(i,j,k)*(zs(k)-z0)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0		moving Z reference
		
		
		totalT_YZ = totalT_YZ + den_flu*(-FZ(i,j,k)*(xs(i)-x0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
			!totalT_YZ = totalT_YZ + den_flu*(-FZ(i,j,k)*(xs(i)-x0)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0   	moving X reference
		totalT_YX = totalT_YX + den_flu*( FX(i,j,k)*(zs(k)-z0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0	
			!totalT_YX = totalT_YX + den_flu*( FX(i,j,k)*(zs(k)-z0)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0		moving Z reference
		
		
		totalT_ZY = totalT_ZY + den_flu*( FY(i,j,k)*(xs(i)-x0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
			!totalT_ZY = totalT_ZY + den_flu*( FY(i,j,k)*(xs(i)-x0)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    	moving X reference
		totalT_ZX = totalT_ZX + den_flu*(-FX(i,j,k)*(ys(j)-y0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0


    end do;end do;end do
    !$OMP END DO
	!$OMP END PARALLEL
   !$acc end parallel 
!$acc end data 
	

    call MPI_ALLREDUCE( totalFX, totalFX_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFZ, totalFZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY, totalFY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
	
	call MPI_ALLREDUCE( totalT_XZ, totalT_XZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalT_XY, totalT_XY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

	
	call MPI_ALLREDUCE( totalT_YZ, totalT_YZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalT_YX, totalT_YX_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
	
	call MPI_ALLREDUCE( totalT_ZY, totalT_ZY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalT_ZX, totalT_ZX_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )


    !if(myid==master) then 	
	    totalTorqx =   (totalT_XZ_ + totalT_XY_)
		totalTorqy =   (totalT_YZ_ + totalT_YX_)
		totalTorqz =   (totalT_ZY_ + totalT_ZX_)
	
		cDrag = -2.d0 * totalFZ_ / (ref_area*w_solid*w_solid*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift = -2.d0 * totalFY_ / (ref_area*w_solid*w_solid*den_flu*1.d0)
		cTorq = -2.d0 * totalTorq/(2.d0*(rotor_r+blade_r)*w_solid*w_solid*den_flu*rotor_r*1.d0) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower   =   cTorq*ABS(rotor_tsr)
	!endif
	

end subroutine virtualForceTorqueIntegrator3D


subroutine virtualForceTorque_frozen()
	use variables
	implicit none

	totalFZ = 0.d0
	totalFY = 0.d0
    totalTZ = 0.d0
	totalTY = 0.d0


	!$OMP PARALLEL	
	!$OMP DO PRIVATE(j,i)  reduction(+: totalFZ, totalFY) collapse(nclps)
	do k=istart,iend 
	do j=1,ny;do i=1,nx

		totalFZ = totalFZ + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY = totalFY + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )

	end do;end do;end do
	!$OMP END DO

    
	! Calculate torque wrt (z0,y0)

    !$OMP DO PRIVATE(j,i)  reduction(+: totalTZ, totalTY) collapse(nclps)
	do k=istart,iend 
	do j=1,ny;do i=1,nx

		!totalTZ = totalTZ + den_flu*(-FZ(i,j,k)*(ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		!totalTY = totalTY + den_flu*( FY(i,j,k)*(zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0

		!------------------------- Frozen rotor analysis --------------------------------
		totalTZ = totalTZ + den_flu*(-FZ(i,j,k)*(Ys(j)-y0+rotor_r*dsin((-AOA)*PI/180.d0))*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY = totalTY + den_flu*( FY(i,j,k)*(Zs(k)-z0-rotor_r*dcos((-AOA)*PI/180.d0))*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0
		!------------------------- Frozen rotor analysis --------------------------------

    end do;end do;end do
    !$OMP END DO
	!$OMP END PARALLEL
	

    call MPI_ALLREDUCE( totalFZ, totalFZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY, totalFY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ, totalTZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY, totalTY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )	

    !if(myid==master) then 	
		totalTorq =   (totalTZ_ + totalTY_)
	
		cDrag = -2.d0 * totalFZ_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift = -2.d0 * totalFY_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq = -2.d0 * totalTorq/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower   =   cTorq*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)
	!endif
	

end subroutine virtualForceTorque_frozen






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

