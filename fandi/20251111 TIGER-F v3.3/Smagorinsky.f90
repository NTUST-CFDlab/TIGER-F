! 11 Nov 2025 - FDS



subroutine CalculateSmagorinskyViscosity()
    use variables 
    implicit none
	real*8		:: du_dx,dv_dy,dw_dz,du_dy,dv_dx,du_dz,dw_dx,dv_dz,dw_dy,Sxy,Sxz,Syz
	real*8		:: denom_Dxs,denom_Dys,denom_Dzs

!$acc data present(iDx,iDy,iDz,Dxs,Dys,Dzs,ETA,nut(:,:,istart:iend+1), &
!$acc	u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2))

    !$OMP PARALLEL DO               PRIVATE(du_dx,dv_dy,dw_dz,du_dy,dv_dx,du_dz,dw_dx,dv_dz,dw_dy,Sxy,Sxz,Syz, &
	!$OMP							denom_Dxs,denom_Dys,denom_Dzs,strainRateTensor,delta,tau,ufric,Yplus,fwall) COLLAPSE(3)
    !$acc parallel loop independent private(du_dx,dv_dy,dw_dz,du_dy,dv_dx,du_dz,dw_dx,dv_dz,dw_dy,Sxy,Sxz,Syz, &
	!$acc							denom_Dxs,denom_Dys,denom_Dzs,strainRateTensor,delta,tau,ufric,Yplus,fwall) collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx

	!Normal components (Central Diffrence Scheme at cell center)
	du_dx = (u(i,j,k) - u(i-1,j,k)) / iDx(i)
	dv_dy = (v(i,j,k) - v(i,j-1,k)) / iDy(j)
	dw_dz = (w(i,j,k) - w(i,j,k-1)) / iDz(k)
	
	!Shear components (Central Diffrence Scheme at cell center)
	
	denom_Dxs = 1.d0/(Dxs(i-1)+Dxs(i))
	denom_Dys = 1.d0/(Dys(j-1)+Dys(j))
	denom_Dzs = 1.d0/(Dzs(k-1)+Dzs(k))
	
	du_dy = (u(i,j+1,k) - u(i,j-1,k) + u(i-1,j+1,k) - u(i-1,j-1,k)) * denom_Dys
	dv_dx = (v(i+1,j,k) - v(i-1,j,k) + v(i+1,j-1,k) - v(i-1,j-1,k)) * denom_Dxs

	du_dz = (u(i,j,k+1) - u(i,j,k-1) + u(i-1,j,k+1) - u(i-1,j,k-1)) * denom_Dzs
	dw_dx = (w(i+1,j,k) - w(i-1,j,k) + w(i+1,j,k-1) - w(i-1,j,k-1)) * denom_Dxs

	dv_dz = (v(i,j,k+1) - v(i,j,k-1) + v(i,j-1,k+1) - v(i,j-1,k-1)) * denom_Dzs
	dw_dy = (w(i,j+1,k) - w(i,j-1,k) + w(i,j+1,k-1) - w(i,j-1,k-1)) * denom_Dys

	Sxy = 0.25d0*(du_dy + dv_dx)	! 0.5d0*(0.5*du_dy + 0.5*dv_dx)
	Sxz = 0.25d0*(du_dz + dw_dx)
	Syz = 0.25d0*(dv_dz + dw_dy)
	

	!Calculation of Strain rate-----------------------------------
    strainRateTensor =	du_dx*du_dx + dv_dy*dv_dy + dw_dz*dw_dz &
						+ 2.d0*(Sxy*Sxy + Sxz*Sxz + Syz*Syz)

    delta = ( iDx(i)*iDy(j)*iDz(k) ) ** (1.d0/3.d0)       

! Applying damping function to solid boundary cells	

		if (0.d0 < ETA(i,j,k) .AND. ETA(i,j,k) < 1.d0) then
			tau = 0.5d0*Cf*den_flu * 0.25d0*((u(i,j,k)+u(i-1,j,k))**2 + (v(i,j,k)+v(i,j-1,k))**2 + (w(i,j,k)+w(i,j,k-1))**2)
			ufric = SQRT(tau/den_flu)
			Yplus =  dySml*ufric/nu
			fwall = (1.d0 - EXP(-Yplus/25.d0))

			nut(i,j,k) = (Cs*delta*fwall)**2*SQRT(2.d0*strainRateTensor) 
		else
			nut(i,j,k) = (Cs*delta)**2*SQRT(2.d0*strainRateTensor)   
		end if 


    end do; end do; end do
    !$acc end parallel
    !$OMP END PARALLEL DO
!$acc end data


end subroutine CalculateSmagorinskyViscosity
