! 11 Nov 2025 - FDS

subroutine CalculateWALEViscosity()
    use variables 
    implicit none
	real*8		:: du_dx,dv_dy,dw_dz,du_dy,dv_dx,du_dz,dw_dx,dv_dz,dw_dy,Sxy,Sxz,Syz
	real*8		:: inv_dx,inv_dy,inv_dz
	real*8		:: u_jp,w_jp,u_jn,w_jn,v_ip,w_ip,v_in,w_in,u_kp,v_kp,u_kn,v_kn,v_jp,v_jn,u_ip,u_in,w_kp,w_kn
	real*8		:: guu,guv,guw,gvu,gvv,gvw,gwu,gwv,gww
	real*8		:: SDxx,SDyy,SDzz,SDxy,SDxz,SDyz,traceG,tensorInvariant

!$acc data present(iDx,iDy,iDz,Dxs,Dys,Dzs,ETA,nut(:,:,istart:iend+1), &
!$acc	u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2))

    !$OMP PARALLEL DO               PRIVATE(du_dx,dv_dy,dw_dz,du_dy,dv_dx,du_dz,dw_dx,dv_dz,dw_dy,Sxy,Sxz,Syz, &
	!$OMP							inv_dx,inv_dy,inv_dz,strainRateTensor,tensorInvariant,delta, &
	!$OMP							u_jp,w_jp,u_jn,w_jn,v_ip,w_ip,v_in,w_in,u_kp,v_kp,u_kn,v_kn,v_jp,v_jn,u_ip,u_in,w_kp,w_kn, &	
	!$OMP							guu,guv,guw,gvu,gvv,gvw,gwu,gwv,gww,SDxy,SDxz,SDyz,traceG,SDxx,SDyy,SDzz) COLLAPSE(3)
    !$acc parallel loop independent private(du_dx,dv_dy,dw_dz,du_dy,dv_dx,du_dz,dw_dx,dv_dz,dw_dy,Sxy,Sxz,Syz, &
	!$acc							inv_dx,inv_dy,inv_dz,strainRateTensor,tensorInvariant,delta, &
	!$acc							u_jp,w_jp,u_jn,w_jn,v_ip,w_ip,v_in,w_in,u_kp,v_kp,u_kn,v_kn,v_jp,v_jn,u_ip,u_in,w_kp,w_kn, &
	!$acc							guu,guv,guw,gvu,gvv,gvw,gwu,gwv,gww,SDxy,SDxz,SDyz,traceG,SDxx,SDyy,SDzz) collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx

	inv_dx = 1.d0/(Dxs(i-1) + Dxs(i))
	inv_dy = 1.d0/(Dys(j-1) + Dys(j))
	inv_dz = 1.d0/(Dzs(k-1) + Dzs(k))	

  if (ETA(i,j,k) > 0.d0) then		! Ommitting solid & interface cells
  
	nut(i,j,k) = 0.d0
	cycle
	
  else if(ANY(ETA(i-1:i+1, j-1:j+1, k-1:k+1) > 0.d0)) then		! Suppressing spurious velocity gradients due to step-wise interface.
	if(ETA(i,j+1,k) > 0.d0) then
		u_jp = u(i-1,j+2,k) + u(i,j+2,k)
		w_jp = w(i,j+2,k-1) + w(i,j+2,k)
		v_jp = v(i,j+1,k)
	else
		u_jp = u(i-1,j+1,k) + u(i,j+1,k)
		w_jp = w(i,j+1,k-1) + w(i,j+1,k)
		v_jp = v(i,j,k)
	endif

	if(ETA(i,j-1,k) > 0.d0) then
		u_jn = u(i-1,j-2,k) + u(i,j-2,k)
		w_jn = w(i,j-2,k-1) + w(i,j-2,k)
		v_jn = v(i,j-2,k)
	else
		u_jn = u(i-1,j-1,k) + u(i,j-1,k)
		w_jn = w(i,j-1,k-1) + w(i,j-1,k)
		v_jn = v(i,j-1,k)
	endif	

	if(ETA(i+1,j,k) > 0.d0) then
		v_ip = v(i+2,j-1,k) + v(i+2,j,k)
		w_ip = w(i+2,j,k-1) + w(i+2,j,k)
		u_ip = u(i+1,j,k)
	else
		v_ip = v(i+1,j-1,k) + v(i+1,j,k)
		w_ip = w(i+1,j,k-1) + w(i+1,j,k)
		u_ip = u(i,j,k)
	endif	

	if(ETA(i-1,j,k) > 0.d0) then
		v_in = v(i-2,j-1,k) + v(i-2,j,k)
		w_in = w(i-2,j,k-1) + w(i-2,j,k)
		u_in = u(i-2,j,k)
	else
		v_in = v(i-1,j-1,k) + v(i-1,j,k)
		w_in = w(i-1,j,k-1) + w(i-1,j,k)
		u_in = u(i-1,j,k)
	endif

	if(ETA(i,j,k+1) > 0.d0) then
		u_kp = u(i-1,j,k+2) + u(i,j,k+2)
		v_kp = v(i,j-1,k+2) + v(i,j,k+2)
		w_kp = w(i,j,k+1)
	else
		u_kp = u(i-1,j,k+1) + u(i,j,k+1)
		v_kp = v(i,j-1,k+1) + v(i,j,k+1)
		w_kp = w(i,j,k)
	endif

	if(ETA(i,j,k-1) > 0.d0) then
		u_kn = u(i-1,j,k-2) + u(i,j,k-2)
		v_kn = v(i,j-1,k-2) + v(i,j,k-2)
		w_kn = w(i,j,k-2)
	else
		u_kn = u(i-1,j,k-1) + u(i,j,k-1)
		v_kn = v(i,j-1,k-1) + v(i,j,k-1)
		w_kn = w(i,j,k-1)
	endif
	

	!Normal components (central difference)
	du_dx = (u_ip - u_in) / iDx(i)
	dv_dy = (v_jp - v_jn) / iDy(j)
	dw_dz = (w_kp - w_kn) / iDz(k)
	
	!Shear components (central difference)
	du_dy = 0.5d0*(u_jp - u_jn) * inv_dy
	dv_dx = 0.5d0*(v_ip - v_in) * inv_dx

	du_dz = 0.5d0*(u_kp - u_kn) * inv_dz
	dw_dx = 0.5d0*(w_ip - w_in) * inv_dx

	dv_dz = 0.5d0*(v_kp - v_kn) * inv_dz
	dw_dy = 0.5d0*(w_jp - w_jn) * inv_dy

  else

	!Normal components (central difference)
	du_dx = 0.5d0*(u(i+1,j,k) - u(i-1,j,k)) / Dxs(i)
	dv_dy = 0.5d0*(v(i,j+1,k) - v(i,j-1,k)) / Dys(j)
	dw_dz = 0.5d0*(w(i,j,k+1) - w(i,j,k-1)) / Dzs(k)
	
	!Shear components (central difference)
	du_dy = (u(i,j+1,k) - u(i,j-1,k)) * inv_dy
	dv_dx = (v(i+1,j,k) - v(i-1,j,k)) * inv_dx

	du_dz = (u(i,j,k+1) - u(i,j,k-1)) * inv_dz
	dw_dx = (w(i+1,j,k) - w(i-1,j,k)) * inv_dx

	dv_dz = (v(i,j,k+1) - v(i,j,k-1)) * inv_dz
	dw_dy = (w(i,j+1,k) - w(i,j-1,k)) * inv_dy

  endif

	Sxy = 0.5d0*(du_dy + dv_dx)	! 0.5d0*(du_dy + dv_dx)
	Sxz = 0.5d0*(du_dz + dw_dx)
	Syz = 0.5d0*(dv_dz + dw_dy)
	

	!Calculation of Strain rate-----------------------------------
    strainRateTensor =	du_dx*du_dx + dv_dy*dv_dy + dw_dz*dw_dz &
						+ 2.d0*(Sxy*Sxy + Sxz*Sxz + Syz*Syz)
     

	! Calculation of WALE tensor

		! Square of Gradient Tensor
		guu = du_dx*du_dx + du_dy*dv_dx + du_dz*dw_dx
		guv = du_dx*du_dy + du_dy*dv_dy + du_dz*dw_dy
		guw = du_dx*du_dz + du_dy*dv_dz + du_dz*dw_dz
	
		gvu = dv_dx*du_dx + dv_dy*dv_dx + dv_dz*dw_dx
		gvv = dv_dx*du_dy + dv_dy*dv_dy + dv_dz*dw_dy
		gvw = dv_dx*du_dz + dv_dy*dv_dz + dv_dz*dw_dz
	
		gwu = dw_dx*du_dx + dw_dy*dv_dx + dw_dz*dw_dx
		gwv = dw_dx*du_dy + dw_dy*dv_dy + dw_dz*dw_dy
		gww = dw_dx*du_dz + dw_dy*dv_dz + dw_dz*dw_dz

		SDxy = 0.5d0*(guv+gvu)
		SDxz = 0.5d0*(guw+gwu)
		SDyz = 0.5d0*(gvw+gwv)
	
		traceG = (guu + gvv + gww)/3.d0
	
		SDxx = guu - traceG
		SDyy = gvv - traceG
		SDzz = gww - traceG

	! Calculation of second tensor invariant
	tensorInvariant = SDxx*SDxx + SDyy*SDyy + SDzz*SDzz &
					  + 2.d0*(SDxy*SDxy + SDxz*SDxz + SDyz*SDyz)

    delta = ( iDx(i)*iDy(j)*iDz(k) ) ** (1.d0/3.d0)  


	nut(i,j,k) = (Cw*delta)**2 * (tensorInvariant**1.5) / &
				 (strainRateTensor**2.5 + tensorInvariant**1.25 + 1.d-16)

    end do; end do; end do
    !$acc end parallel
    !$OMP END PARALLEL DO
!$acc end data
end subroutine CalculateWALEViscosity



subroutine CalculateSmagorinskyViscosity()
    use variables 
    implicit none
	real*8		:: du_dx,dv_dy,dw_dz,du_dy,dv_dx,du_dz,dw_dx,dv_dz,dw_dy,Sxy,Sxz,Syz
	real*8		:: inv_dx,inv_dy,inv_dz

!$acc data present(iDx,iDy,iDz,Dxs,Dys,Dzs,ETA,nut(:,:,istart:iend+1), &
!$acc	u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2))

    !$OMP PARALLEL DO               PRIVATE(du_dx,dv_dy,dw_dz,du_dy,dv_dx,du_dz,dw_dx,dv_dz,dw_dy,Sxy,Sxz,Syz, &
	!$OMP							inv_dx,inv_dy,inv_dz,strainRateTensor,delta,tau_d,ufric_d,Yplus_d,fwall_d) COLLAPSE(3)
    !$acc parallel loop independent private(du_dx,dv_dy,dw_dz,du_dy,dv_dx,du_dz,dw_dx,dv_dz,dw_dy,Sxy,Sxz,Syz, &
	!$acc							inv_dx,inv_dy,inv_dz,strainRateTensor,delta,tau_d,ufric_d,Yplus_d,fwall_d) collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx

	inv_dx = 1.d0/(Dxs(i-1) + Dxs(i))
	inv_dy = 1.d0/(Dys(j-1) + Dys(j))
	inv_dz = 1.d0/(Dzs(k-1) + Dzs(k))

	!Normal components (central difference)
	du_dx = 0.5d0*(u(i+1,j,k) - u(i-1,j,k)) / Dxs(i)
	dv_dy = 0.5d0*(v(i,j+1,k) - v(i,j-1,k)) / Dys(j)
	dw_dz = 0.5d0*(w(i,j,k+1) - w(i,j,k-1)) / Dzs(k)
	
	!Shear components (central difference)
	du_dy = (u(i,j+1,k) - u(i,j-1,k)) * inv_dy
	dv_dx = (v(i+1,j,k) - v(i-1,j,k)) * inv_dx

	du_dz = (u(i,j,k+1) - u(i,j,k-1)) * inv_dz
	dw_dx = (w(i+1,j,k) - w(i-1,j,k)) * inv_dx

	dv_dz = (v(i,j,k+1) - v(i,j,k-1)) * inv_dz
	dw_dy = (w(i,j+1,k) - w(i,j-1,k)) * inv_dy

	Sxy = 0.5d0*(du_dy + dv_dx)	! 0.5d0*(du_dy + dv_dx)
	Sxz = 0.5d0*(du_dz + dw_dx)
	Syz = 0.5d0*(dv_dz + dw_dy)
	

	!Calculation of Strain rate-----------------------------------
    strainRateTensor =	du_dx*du_dx + dv_dy*dv_dy + dw_dz*dw_dz &
						+ 2.d0*(Sxy*Sxy + Sxz*Sxz + Syz*Syz)

    delta = ( iDx(i)*iDy(j)*iDz(k) ) ** (1.d0/3.d0)       

! Applying damping function to solid boundary cells	

		if (0.d0 < ETA(i,j,k) .AND. ETA(i,j,k) < 1.d0) then
			tau_d = 0.5d0*Cf_d*den_flu * 0.25d0*((u(i,j,k)+u(i-1,j,k))**2 + (v(i,j,k)+v(i,j-1,k))**2 + (w(i,j,k)+w(i,j,k-1))**2)
			ufric_d = SQRT(tau_d/den_flu)
			Yplus_d =  dySml*ufric_d/nu
			fwall_d = (1.d0 - EXP(-Yplus_d/25.d0))

			nut(i,j,k) = (Cs*delta*fwall_d)**2*SQRT(2.d0*strainRateTensor) 
		else
			nut(i,j,k) = (Cs*delta)**2*SQRT(2.d0*strainRateTensor)   
		end if 


    end do; end do; end do
    !$acc end parallel
    !$OMP END PARALLEL DO
!$acc end data

end subroutine CalculateSmagorinskyViscosity
