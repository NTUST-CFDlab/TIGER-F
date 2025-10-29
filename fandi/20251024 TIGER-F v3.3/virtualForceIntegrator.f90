! 30 Oct 2025 - FDS


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


subroutine virtualForce_tandemVIV()
	use variables
	implicit none

	totalFZ1 = 0.d0 ; totalFZ2 = 0.d0
	totalFY1 = 0.d0 ; totalFY2 = 0.d0

!$acc data present(iDx,iDy,iDz,Zs,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))

   !$OMP PARALLEL DO PRIVATE(j,i)  reduction(+: totalFZ1, totalFY1, totalFZ2, totalFY2) collapse(nclps)
   !$acc parallel
   !$acc loop independent reduction(+: totalFZ1, totalFY1, totalFZ2, totalFY2) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

		 if(Zs(k) .LT. Zsplit_pos) then
			totalFZ1 = totalFZ1 + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
			totalFY1 = totalFY1 + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		 else
			totalFZ2 = totalFZ2 + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
			totalFY2 = totalFY2 + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		 endif
	
	end do;end do;end do
   !$acc end parallel 
   !$OMP END PARALLEL DO
!$acc end data 

    call MPI_ALLREDUCE( totalFZ1, totalFZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY1, totalFY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

    call MPI_ALLREDUCE( totalFZ2, totalFZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY2, totalFY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

  
		cDrag1 = -2.d0 * totalFZ1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift1 = -2.d0 * totalFY1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cPower1   =   cLift1 * v_solid_1 / U_inf

		cDrag2 = -2.d0 * totalFZ2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift2 = -2.d0 * totalFY2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cPower2   =   cLift2 * v_solid_2 / U_inf		
		

end subroutine virtualForce_tandemVIV



subroutine virtualForce_Magnus2blade()
	use variables
	implicit none
	
	totalFZ1 = 0.d0
	totalFY1 = 0.d0
    totalTZ1 = 0.d0
	totalTY1 = 0.d0
	totalTZ1_in = 0.d0
	totalTY1_in = 0.d0 
	
	totalFZ2 = 0.d0
	totalFY2 = 0.d0
    totalTZ2 = 0.d0
	totalTY2 = 0.d0
	totalTZ2_in = 0.d0
	totalTY2_in = 0.d0
	

!$acc data present(iDx,iDy,iDz,Xs,Ys,Zs,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))

	!$OMP PARALLEL	
	!$OMP DO PRIVATE(j,i)  reduction(+: totalFZ1, totalFY1, totalFZ2, totalFY2) collapse(nclps)
	!$acc parallel
	!$acc loop independent reduction(+: totalFZ1, totalFY1, totalFZ2, totalFY2) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	   if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (3.d0*blade_r)) then
		totalFZ1 = totalFZ1 + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY1 = totalFY1 + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )	
	   else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (3.d0*blade_r)) then
		totalFZ2 = totalFZ2 + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY2 = totalFY2 + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
	   endif

	end do;end do;end do
	!$OMP END DO

    
	! Calculate torque wrt (z0,y0)

    !$OMP DO PRIVATE(j,i)  reduction(+: totalTZ1, totalTY1, totalTZ2, totalTY2) collapse(nclps)
	!$acc loop independent reduction(+: totalTZ1, totalTY1, totalTZ2, totalTY2) collapse(3) gang vector
	do k=istart,iend
	do j=1,ny;do i=1,nx

	   if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (3.d0*blade_r)) then
		totalTZ1 = totalTZ1 + den_flu*(-FZ(i,j,k)*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY1 = totalTY1 + den_flu*( FY(i,j,k)*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0		

	   else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (3.d0*blade_r)) then
		totalTZ2 = totalTZ2 + den_flu*(-FZ(i,j,k)*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY2 = totalTY2 + den_flu*( FY(i,j,k)*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0

	   endif
	  
    end do;end do;end do
    !$OMP END DO
	!$OMP END PARALLEL



	! Calculate input torque wrt each cylinder center

    !$OMP DO PRIVATE(j,i)  reduction(+: totalTZ1_in, totalTY1_in, totalTZ2_in, totalTY2_in) collapse(nclps)
	!$acc loop independent reduction(+: totalTZ1_in, totalTY1_in, totalTZ2_in, totalTY2_in) collapse(3) gang vector
	do k=istart,iend
	do j=1,ny;do i=1,nx

	   if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (blade_r+dzSml)) then
		totalTZ1_in = totalTZ1_in + den_flu*(-FZ(i,j,k)*(Ys(j)-y0_t1 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY1_in = totalTY1_in + den_flu*( FY(i,j,k)*(Zs(k)-z0_t1 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0		

	   else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (blade_r+dzSml)) then
		totalTZ2_in = totalTZ2_in + den_flu*(-FZ(i,j,k)*(Ys(j)-y0_t2 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY2_in = totalTY2_in + den_flu*( FY(i,j,k)*(Zs(k)-z0_t2 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0

	   endif
	  
    end do;end do;end do
    !$OMP END DO
	!$OMP END PARALLEL
   !$acc end parallel 
   
   
   
!$acc end data 
	

    call MPI_ALLREDUCE( totalFZ1, totalFZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY1, totalFY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ1, totalTZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY1, totalTY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ1_in, totalTZ1_in_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY1_in, totalTY1_in_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )	

    call MPI_ALLREDUCE( totalFZ2, totalFZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY2, totalFY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ2, totalTZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY2, totalTY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ2_in, totalTZ2_in_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY2_in, totalTY2_in_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

		totalFZ_ = totalFZ1_ + totalFZ2_ 
		totalFY_ = totalFY1_ + totalFY2_


		totalTorq1 =   (totalTZ1_ + totalTY1_)
		totalTorq2 =   (totalTZ2_ + totalTY2_)

		totalTorq  = totalTorq1 + totalTorq2


		totalTorq1_in =   (totalTZ1_in_ + totalTY1_in_)
		totalTorq2_in =   (totalTZ2_in_ + totalTY2_in_)
		
	
		cDrag = -2.d0 * totalFZ_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift = -2.d0 * totalFY_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq = -2.d0 * totalTorq/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower   =   cTorq*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cPower_in   =   -2.d0 * (totalTorq1_in * blade_omega + totalTorq2_in * blade_omega) / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
		
		cPower_nett = cPower - cPower_in



		cDrag1 = -2.d0 * totalFZ1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift1 = -2.d0 * totalFY1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq1 = -2.d0 * totalTorq1/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower1   =   cTorq1*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cPower1_in   =   -2.d0 * (totalTorq1_in * blade_omega) / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
		
		cDrag2 = -2.d0 * totalFZ2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift2 = -2.d0 * totalFY2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq2 = -2.d0 * totalTorq2/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower2   =   cTorq2*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cPower2_in   =   -2.d0 * (totalTorq2_in * blade_omega) / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
	

end subroutine virtualForce_Magnus2blade


subroutine virtualForce_Magnus3blade()
	use variables
	implicit none
	
	totalFZ1 = 0.d0
	totalFY1 = 0.d0
    totalTZ1 = 0.d0
	totalTY1 = 0.d0
	
	totalFZ2 = 0.d0
	totalFY2 = 0.d0
    totalTZ2 = 0.d0
	totalTY2 = 0.d0

	totalFZ3 = 0.d0
	totalFY3 = 0.d0
    totalTZ3 = 0.d0
	totalTY3 = 0.d0

!$acc data present(iDx,iDy,iDz,Xs,Ys,Zs,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))

	!$OMP PARALLEL	
	!$OMP DO PRIVATE(j,i)  reduction(+: totalFZ1, totalFY1, totalFZ2, totalFY2, totalFZ3, totalFY3) collapse(nclps)
	!$acc parallel
	!$acc loop independent reduction(+: totalFZ1, totalFY1, totalFZ2, totalFY2, totalFZ3, totalFY3) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	   if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (3.d0*blade_r)) then
		totalFZ1 = totalFZ1 + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY1 = totalFY1 + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )	
	   else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (3.d0*blade_r)) then
		totalFZ2 = totalFZ2 + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY2 = totalFY2 + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
	   else if( SQRT((Zs(k)-z0_t3)**2 + (Ys(j)-y0_t3)**2) .LT. (3.d0*blade_r)) then
		totalFZ3 = totalFZ3 + den_flu*( FZ(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY3 = totalFY3 + den_flu*( FY(i,j,k) * ( iDx(i)*iDy(j)*iDz(k) ) )
	   endif

	end do;end do;end do
	!$OMP END DO

    
	! Calculate torque wrt (z0,y0)

    !$OMP DO PRIVATE(j,i)  reduction(+: totalTZ1, totalTY1, totalTZ2, totalTY2, totalTZ3, totalTY3) collapse(nclps)
	!$acc loop independent reduction(+: totalTZ1, totalTY1, totalTZ2, totalTY2, totalTZ3, totalTY3) collapse(3) gang vector
	do k=istart,iend
	do j=1,ny;do i=1,nx

	   if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (3.d0*blade_r)) then
		totalTZ1 = totalTZ1 + den_flu*(-FZ(i,j,k)*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY1 = totalTY1 + den_flu*( FY(i,j,k)*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0
	   else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (3.d0*blade_r)) then
		totalTZ2 = totalTZ2 + den_flu*(-FZ(i,j,k)*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY2 = totalTY2 + den_flu*( FY(i,j,k)*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0
	   else if( SQRT((Zs(k)-z0_t3)**2 + (Ys(j)-y0_t3)**2) .LT. (3.d0*blade_r)) then
		totalTZ3 = totalTZ3 + den_flu*(-FZ(i,j,k)*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY3 = totalTY3 + den_flu*( FY(i,j,k)*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0
	   endif

    end do;end do;end do
    !$OMP END DO
	!$OMP END PARALLEL
   !$acc end parallel 
!$acc end data 
	

    call MPI_ALLREDUCE( totalFZ1, totalFZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY1, totalFY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ1, totalTZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY1, totalTY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )	

    call MPI_ALLREDUCE( totalFZ2, totalFZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY2, totalFY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ2, totalTZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY2, totalTY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

    call MPI_ALLREDUCE( totalFZ3, totalFZ3_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY3, totalFY3_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ3, totalTZ3_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY3, totalTY3_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

		totalFZ_ = totalFZ1_ + totalFZ2_ + totalFZ3_ 
		totalFY_ = totalFY1_ + totalFY2_ + totalFY3_


		totalTorq1 =   (totalTZ1_ + totalTY1_)
		totalTorq2 =   (totalTZ2_ + totalTY2_)
		totalTorq3 =   (totalTZ3_ + totalTY3_)
		
		totalTorq  = totalTorq1 + totalTorq2 + totalTorq3
		
	
		cDrag = -2.d0 * totalFZ_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift = -2.d0 * totalFY_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq = -2.d0 * totalTorq/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower   =   cTorq*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cDrag1 = -2.d0 * totalFZ1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift1 = -2.d0 * totalFY1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq1 = -2.d0 * totalTorq1/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower1   =   cTorq1*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)
		
		cDrag2 = -2.d0 * totalFZ2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift2 = -2.d0 * totalFY2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq2 = -2.d0 * totalTorq2/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower2   =   cTorq2*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cDrag3 = -2.d0 * totalFZ3_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift3 = -2.d0 * totalFY3_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq3 = -2.d0 * totalTorq3/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower3   =   cTorq3*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)
	

end subroutine virtualForce_Magnus3blade


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

