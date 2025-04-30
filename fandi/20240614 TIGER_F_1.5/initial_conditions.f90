! 22 Aug 2023 - FDS
subroutine initial_conditions()
use variables
implicit none

!---------------------------------------------------!
!         Initial conditions calculation            !
!---------------------------------------------------!

	!----------- Solid dynamics --------------------------------------!
	xc_t = xc ; yc_t = yc ;	zc_t = zc
    u_solid=0.d0;v_solid=0.d0;w_solid=0.d0

	rotor_omega=0.d0
	AOA = AOA1
	
	totalFX_=0.d0;totalFY_=0.d0;totalFZ_=0.d0
	totalTorq=0.d0

	!------------ Turbulent viscosity for LES = OFF -------------------!
	Viseff = 0.d0
	
	!-------------Damping Wall Function (stationary)-------------------!

   if (Re .GT. 10.0) then
    Cf=  (2.d0*dlog10(Re)-0.65d0)**(-2.3)
    !Cf = 0.026/ (Re)**(1/7)
    tau = 0.5d0*Cf*den_flu *U_inf**2
    ufric = dsqrt(tau/den_flu)
    Yplus =  (MINVAL(iDy)*ufric)/nu
    fwall = (1.d0 - dexp(-Yplus/25.d0))
   endif


!$OMP PARALLEL DO PRIVATE(j,i)  
do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2

   !u(i,j,k) = 0.0
   !v(i,j,k) = 0.0
   w(i,j,k) = U_inf
   !u1(i,j,k) = u(i,j,k)
   !v1(i,j,k) = v(i,j,k)
   !w1(i,j,k) = w(i,j,k)
   !u2(i,j,k) = u(i,j,k)
   !v2(i,j,k) = v(i,j,k)
   !w2(i,j,k) = w(i,j,k)
   !p(i,j,k) = 0.0
   !u_star(i,j,k) = 0.0
   !v_star(i,j,k) = 0.0
   !w_star(i,j,k) = 0.0
   !ETA(i,j,k) = 0.0
   !F_tavex(i,j,k) = 0.0
   !F_tavey(i,j,k) = 0.0

end do; end do; end do
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(i,j)
do k=1,nz 
do j=1,ny; do i=1,nx

   Den(i,j,k) = ( iDy(j) * iDz(k) / Dxs(i) + iDy(j) * iDz(k) / Dxs(i-1) &
                + iDx(i) * iDz(k) / Dys(j) + iDx(i) * iDz(k) / Dys(j-1) &
                + iDx(i) * iDy(j) / Dzs(k) + iDx(i) * iDy(j) / Dzs(k-1) )

   Den_inv(i,j,k) = 1.d0/Den(i,j,k)

   P_Den(i,j,k,1) = iDy(j) * iDz(k) / Dxs(i) 
   P_Den(i,j,k,2) = iDy(j) * iDz(k) / Dxs(i-1)
   P_Den(i,j,k,3) = iDx(i) * iDz(k) / Dys(j) 
   P_Den(i,j,k,4) = iDx(i) * iDz(k) / Dys(j-1)
   P_Den(i,j,k,5) = iDx(i) * iDy(j) / Dzs(k) 
   P_Den(i,j,k,6) = iDx(i) * iDy(j) / Dzs(k-1)

enddo; enddo; enddo
!$OMP END PARALLEL DO




end subroutine initial_conditions



