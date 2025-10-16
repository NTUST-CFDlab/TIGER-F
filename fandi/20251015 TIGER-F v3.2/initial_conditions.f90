! 12 Jan 2025 - FDS

subroutine initial_conditions()
use variables
implicit none

!---------------------------------------------------!
!         Initial conditions calculation            !
!---------------------------------------------------!

	!----------- Solid dynamics --------------------------------------!
	x0_t = x0 ; y0_t = y0 ;	z0_t = z0
    u_solid=0.d0;v_solid=0.d0;w_solid=0.d0
	du_solid=0.d0;dw_solid=0.d0		! for falling sphere

	rotor_omega=0.d0
	AOA = AOA1

	rotate_sx=0.d0; rotate_sy=0.d0; rotate_sz=0.d0

	! ------------ 2D KE calculation -------------------------------------!
    ! KE_timeavg(1, :, :) = 0.0d0
	! KE_timestd(1, :, :) = 0.0d0	

	!------------ Turbulent viscosity for LES = OFF -------------------!
	nut = 0.d0
	
	!-------------Initialize Wall Damping Function-------------------!

	Re_t = MAX((U_inf*L_ch)/nu,3.d0)
	Cf = (2.d0*LOG10(Re_t)-0.65d0)**(-2.3)
	

!--------- pinned memory------------------
allocate(   p(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(   u(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(   v(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(   w(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate( ETA(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(  u0(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(  v0(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(  w0(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(  FX(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(  FY(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(  FZ(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate( nut( 1:nx  , 1:ny  , 1:nz+1) )

if (wall_model==1) then
allocate(yplus_print(1:nx,1:ny,1:nz) )
allocate(uplus_print(1:nx,1:ny,1:nz) )
allocate( dist_print(1:nx,1:ny,1:nz) )
end if 

if (filer_running_avg==1) then
allocate(p_TimeAvg(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(u_TimeAvg(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(v_TimeAvg(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate(w_TimeAvg(-1:nx+2,-1:ny+2,-1:nz+2) )
allocate( 	   TKE(-1:nx+2,-1:ny+2,-1:nz+2) )
end if 

!$OMP PARALLEL DO PRIVATE(j,i)  
!$acc parallel device_type(host)
!$acc loop independent collapse(3)
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
!$acc end parallel
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(i,j)
!$acc parallel device_type(host)
!$acc loop independent collapse(3)
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
!$acc end parallel
!$OMP END PARALLEL DO




end subroutine initial_conditions



