! 20 Nov 2023 - FDS



subroutine twoway_dynamic() 
    use variables 
    implicit none


   !---------------------------------------------------------!
   !      update the bulk solid velocity caused by 2way      !
   !---------------------------------------------------------!

	!---------------For free falling sphere-----------
	u_solid = u_solid + dt*(-totalFX_/den_sol/solid_volume)
	v_solid = v_solid + dt*(-totalFY_/den_sol/solid_volume)
	w_solid = w_solid + dt*(-g+(-totalFZ_/den_sol/solid_volume)+(den_flu*g/den_sol))
   
	rotate_sx = rotate_sx + dt*(-totalTorqx*2.5d0/den_sol/solid_volume/r**2)
	rotate_sy = rotate_sy + dt*(-totalTorqy*2.5d0/den_sol/solid_volume/r**2)
	rotate_sz = rotate_sz + dt*(-totalTorqz*2.5d0/den_sol/solid_volume/r**2)
		
	!---------------For free falling sphere-----------

	!---------------For Translation-----------
	! u_solid = u_solid
	! v_solid = v_solid - totalFY_*dt/den_sol/solid_volume
	! w_solid = w_solid
	!---------------For Translation-----------
   

	!---------------For Magnus effect VAWT-----------
	!rotor_omega = rotor_omega - totalTorq*dt/i_sol
	!---------------For Magnus effect VAWT-----------

	! !---------------For Darrius VAWT-----------
	! rotor_omega = rotor_omega - totalTorq*dt/i_sol
	! !---------------For Darrius VAWT-----------


   !---------------------------------------------------------!
   !    update the bulk solid position caused by 2way        !
   !---------------------------------------------------------!

	!---------------For free falling sphere-----------
	xc_t = xc_t + u_solid*dt
	yc_t = yc_t + v_solid*dt
	zc_t = zc_t + w_solid*dt
      
	sol_speed = sqrt(u_solid**2 + v_solid**2 + (w_solid - U_inf)**2)
	!---------------For free falling sphere-----------

	! !---------------------For Translation-----------------------------
	! xc_t = xc_t + u_solid*dt
	! yc_t = yc_t + v_solid*dt
	! zc_t = zc_t + w_solid*dt
      
	! sol_speed = sqrt(u_solid**2 + v_solid**2 + (w_solid - U_inf)**2)
	! !---------------------For Translation-----------------------------	
	

	! !---------------------For Magnus effect VAWT (2 blades)-----------------------------
	! AOA = AOA + rotor_omega*dt*180.d0/PI

	! ! Rotate the 1st blade center with AOA (needed to define the blade rotation)
	! ! Blades' center is shifted from (0,0) to (-rotor_r, 0)

		! zc_t1= dCOS((AOA)*PI/180.d0)*(-rotor_r) - dSIN((AOA)*PI/180.d0)*0.d0
		! yc_t1= dSIN((AOA)*PI/180.d0)*(-rotor_r) + dCOS((AOA)*PI/180.d0)*0.d0
		
		! yc_t1 = yc_t1 + yc
		! zc_t1 = zc_t1 + zc
		
	! ! Rotate the 2nd blade center with AOA (needed to define the blade rotation)
	! ! Blades' center is shifted from (0,0) to (-rotor_r, 0)

		! zc_t2= dCOS((AOA+180.d0)*PI/180.d0)*(-rotor_r) - dSIN((AOA+180.d0)*PI/180.d0)*0.d0
		! yc_t2= dSIN((AOA+180.d0)*PI/180.d0)*(-rotor_r) + dCOS((AOA+180.d0)*PI/180.d0)*0.d0
		
		! yc_t2 = yc_t2 + yc
		! zc_t2 = zc_t2 + zc
		
	  ! sol_speed = U_inf
		
	! !---------------------For Magnus effect VAWT (2 blades)-----------------------------


	! !---------------------For Darrius VAWT (3 blades)-----------------------------
	! AOA = AOA + rotor_omega*dt*180.d0/PI

	  ! sol_speed = U_inf
		
	! !---------------------For Darrius VAWT (3 blades)-----------------------------

	  
	  Re_t= (sol_speed*L_ch)/nu
	  
	  !-------------Damping Wall Function (two-way)-------------------!

	if (Re_t .GT. 10.0) then
      Cf_d=  (2.d0*dlog10(Re_t*(1.d0+DABS(blade_alpha)))-0.65d0)**(-2.3)
      tau_d = 0.5d0*Cf_d*den_flu *(sol_speed*(1.d0+DABS(blade_alpha)))**2
      ufric_d = dsqrt(tau_d/den_flu)
      Yplus_d =  (MINVAL(iDy)*ufric_d)/nu
      fwall_d = (1.d0 - dexp(-Yplus_d/25.d0))
	endif
	
end subroutine twoway_dynamic


subroutine rotor_dynamic() 
    use variables 
    implicit none

	rotor_omega = rotor_tsr*U_inf/(rotor_r+blade_r)
	AOA = AOA1 + rotor_omega*( time - StartDynamic_time )*180.d0/PI

	!---------------------For Magnus effect VAWT (1 blade)-----------------------------
	! Rotate the blade center with AOA (needed to define the blade rotation)
	! Blades' center is shifted from (0,0) to (-rotor_r, 0)

		zc_t= dCOS((AOA)*PI/180.d0)*(-rotor_r) - dSIN((AOA)*PI/180.d0)*0.d0
		yc_t= dSIN((AOA)*PI/180.d0)*(-rotor_r) + dCOS((AOA)*PI/180.d0)*0.d0
		
		yc_t = yc_t + yc
		zc_t = zc_t + zc
		
	!---------------------For Magnus effect VAWT (1 blade)-----------------------------


	! !---------------------For Magnus effect VAWT (2 blades)-----------------------------
	! Rotate the 1st blade center with AOA (needed to define the blade rotation)
	! Blades' center is shifted from (0,0) to (-rotor_r, 0)

		! zc_t1= dCOS((AOA)*PI/180.d0)*(-rotor_r) - dSIN((AOA)*PI/180.d0)*0.d0
		! yc_t1= dSIN((AOA)*PI/180.d0)*(-rotor_r) + dCOS((AOA)*PI/180.d0)*0.d0
		
		! yc_t1 = yc_t1 + yc
		! zc_t1 = zc_t1 + zc
		
	! Rotate the 2nd blade center with AOA (needed to define the blade rotation)
	! Blades' center is shifted from (0,0) to (-rotor_r, 0)

		! zc_t2= dCOS((AOA+180.d0)*PI/180.d0)*(-rotor_r) - dSIN((AOA+180.d0)*PI/180.d0)*0.d0
		! yc_t2= dSIN((AOA+180.d0)*PI/180.d0)*(-rotor_r) + dCOS((AOA+180.d0)*PI/180.d0)*0.d0
		
		! yc_t2 = yc_t2 + yc
		! zc_t2 = zc_t2 + zc
	! !---------------------For Magnus effect VAWT (2 blades)-----------------------------
	
end subroutine rotor_dynamic


! subroutine blade_dynamic() 
    ! use variables
    ! implicit none

    ! blade_omega = blade_alpha*U_inf/blade_r

! end subroutine blade_dynamic



subroutine func_darrius_3blade()

use variables
implicit none
real*8 :: c_p, t_p, s_p
real*8 :: z_trans1, y_trans1, z_trans2, y_trans2, z_trans3, y_trans3
real*8 :: z_aoa1, y_aoa1, z_aoa2, y_aoa2, z_aoa3, y_aoa3
INTEGER:: az_min,az_max,ay_min,ay_max
integer :: l ,m ,n
real*8 :: xi
real*8 :: dxg, dyg, dzg
real*8 ,dimension(1:nSubGrids_f+1) :: SX
real*8 ,dimension(1:nSubGrids_f+1) :: SY
real*8 ,dimension(1:nSubGrids_f+1) :: SZ

   integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

c_p = L_ch  ! blade's chord
t_p = 0.22d0 ! blade's thickness (NACA 00XX --> t_p = XX/100
s_p = 0.022  ! CM offset from leading edge


! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,ETA) create(SY,SZ)

!$OMP PARALLEL DO PRIVATE(j,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel loop independent private(z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) &
!$acc 						  reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


do k=0,nz+1
	do j=0,ny+1

	! rotation wrt. AOA for blade 1
	z_aoa1 =  (Zs(k) - zc)*dcos((-AOA)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA)*PI/180.d0)
	y_aoa1 =  (Zs(k) - zc)*dsin((-AOA)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA)*PI/180.d0)

	! z and y transformation
	z_trans1 = z_aoa1 + s_p
	y_trans1 = y_aoa1 + rotor_r


	! rotation wrt. AOA for blade 2
	z_aoa2 =  (Zs(k) - zc)*dcos((-AOA+120.d0)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA+120.d0)*PI/180.d0)
	y_aoa2 =  (Zs(k) - zc)*dsin((-AOA+120.d0)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA+120.d0)*PI/180.d0)

	! z and y transformation
	z_trans2 = z_aoa2 + s_p
	y_trans2 = y_aoa2 + rotor_r


	! rotation wrt. AOA for blade 3
	z_aoa3 =  (Zs(k) - zc)*dcos((-AOA+240.d0)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA+240.d0)*PI/180.d0)
	y_aoa3 =  (Zs(k) - zc)*dsin((-AOA+240.d0)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA+240.d0)*PI/180.d0)

	! z and y transformation
	z_trans3 = z_aoa3 + s_p
	y_trans3 = y_aoa3 + rotor_r

	
	if( ((y_trans1 - 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .LE. 0.) .AND. &
	     (y_trans1 + 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .GE. 0.)) .OR. &
		((y_trans2 - 5.*t_p*c_p*(0.2969*sqrt(z_trans2/c_p) - 0.126*z_trans2/c_p - 0.3516*(z_trans2/c_p)**2 + 0.2843*(z_trans2/c_p)**3 - 0.1015*(z_trans2/c_p)**4) .LE. 0.) .AND. &
		 (y_trans2 + 5.*t_p*c_p*(0.2969*sqrt(z_trans2/c_p) - 0.126*z_trans2/c_p - 0.3516*(z_trans2/c_p)**2 + 0.2843*(z_trans2/c_p)**3 - 0.1015*(z_trans2/c_p)**4) .GE. 0.)) .OR. &
		((y_trans3 - 5.*t_p*c_p*(0.2969*sqrt(z_trans3/c_p) - 0.126*z_trans3/c_p - 0.3516*(z_trans3/c_p)**2 + 0.2843*(z_trans3/c_p)**3 - 0.1015*(z_trans3/c_p)**4) .LE. 0.) .AND. &
		 (y_trans3 + 5.*t_p*c_p*(0.2969*sqrt(z_trans3/c_p) - 0.126*z_trans3/c_p - 0.3516*(z_trans3/c_p)**2 + 0.2843*(z_trans3/c_p)**3 - 0.1015*(z_trans3/c_p)**4) .GE. 0.)) ) then
		ETA(i,j,k)=1.D0
		az_min = MIN0(k,az_min)
		az_max = MAX0(k,az_max)
		ay_min = MIN0(j,ay_min)
		ay_max = MAX0(j,ay_max)
	else
		ETA(i,j,k)=0.D0		
	endif
	
	end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

! Setting the boundary for subgrids
	iBgnVOS = 0
	iEndVOS = nx+1
	jBgnVOS = ay_min - 5
	jEndVOS = ay_max + 5
	kBgnVOS = az_min - 5
	kEndVOS = az_max + 5

! Creating subgrids
i=iBgnVOS
!$OMP PARALLEL DO PRIVATE(j,m,n,dzg,dyg,SY,SZ,xi,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) collapse(2)
!$acc parallel loop independent private(m,n,dzg,dyg,SY,SZ,xi,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) collapse(2) gang vector

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	! Finding potential grids for subgrids
	if( (abs(ETA(i,j,k)-ETA(i,j-1,k))+abs(ETA(i,j,k)-ETA(i,j+1,k))+ &
		abs(ETA(i,j,k)-ETA(i,j,k-1))+abs(ETA(i,j,k)-ETA(i,j,k+1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k-1))+abs(ETA(i,j,k)-ETA(i,j+1,k-1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k+1))+abs(ETA(i,j,k)-ETA(i,j+1,k+1))) .GT. 0.5d0) then
			dzg = iDz(k) / nSubGrids_f
			dyg = iDy(j) / nSubGrids_f
			do m=1,nSubGrids_f+1
				SY(m) = Y(j) + (m-1) * dyg		
			end do
			do n=1,nSubGrids_f+1
				SZ(n) = Z(k) + (n-1) * dzg	
			end do

		    xi = 0.0
			do m=1,nSubGrids_f; do n=1,nSubGrids_f

				! rotation wrt. AOA for blade 1			
				z_aoa1 =  (SZ(n) - zc)*dcos((-AOA)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA)*PI/180.d0)
				y_aoa1 =  (SZ(n) - zc)*dsin((-AOA)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA)*PI/180.d0)

				! z and y transformation
				z_trans1 = z_aoa1 + s_p
				y_trans1 = y_aoa1 + rotor_r
				
				! rotation wrt. AOA for blade 2			
				z_aoa2 =  (SZ(n) - zc)*dcos((-AOA+120)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA+120)*PI/180.d0)
				y_aoa2 =  (SZ(n) - zc)*dsin((-AOA+120)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA+120)*PI/180.d0)

				! z and y transformation
				z_trans2 = z_aoa2 + s_p
				y_trans2 = y_aoa2 + rotor_r
					
				! rotation wrt. AOA for blade 3			
				z_aoa3 =  (SZ(n) - zc)*dcos((-AOA+240)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA+240)*PI/180.d0)
				y_aoa3 =  (SZ(n) - zc)*dsin((-AOA+240)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA+240)*PI/180.d0)

				! z and y transformation
				z_trans3 = z_aoa3 + s_p
				y_trans3 = y_aoa3 + rotor_r					
					
		if( ((y_trans1 - 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .LE. 0.) .AND. &
			 (y_trans1 + 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .GE. 0.)) .OR. &
			((y_trans2 - 5.*t_p*c_p*(0.2969*sqrt(z_trans2/c_p) - 0.126*z_trans2/c_p - 0.3516*(z_trans2/c_p)**2 + 0.2843*(z_trans2/c_p)**3 - 0.1015*(z_trans2/c_p)**4) .LE. 0.) .AND. &
			 (y_trans2 + 5.*t_p*c_p*(0.2969*sqrt(z_trans2/c_p) - 0.126*z_trans2/c_p - 0.3516*(z_trans2/c_p)**2 + 0.2843*(z_trans2/c_p)**3 - 0.1015*(z_trans2/c_p)**4) .GE. 0.)) .OR. &
			((y_trans3 - 5.*t_p*c_p*(0.2969*sqrt(z_trans3/c_p) - 0.126*z_trans3/c_p - 0.3516*(z_trans3/c_p)**2 + 0.2843*(z_trans3/c_p)**3 - 0.1015*(z_trans3/c_p)**4) .LE. 0.) .AND. &
			 (y_trans3 + 5.*t_p*c_p*(0.2969*sqrt(z_trans3/c_p) - 0.126*z_trans3/c_p - 0.3516*(z_trans3/c_p)**2 + 0.2843*(z_trans3/c_p)**3 - 0.1015*(z_trans3/c_p)**4) .GE. 0.)) ) then							
				xi = xi + 1
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!$acc parallel loop independent collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)

end do
end do
end do    
!$acc end parallel 

!$acc end data
        
    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, -1*AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine func_darrius_3blade




subroutine func_Cylplate_2blade()

use variables
implicit none
real*8 :: L_p, t_p, S_p, g_p
real*8 :: z_trans1, y_trans1, z_trans2, y_trans2
real*8 :: z_aoa1, y_aoa1, z_aoa2, y_aoa2
INTEGER:: az_min,az_max,ay_min,ay_max
integer :: l ,m ,n
real*8 :: xi
real*8 :: dxg, dyg, dzg
real*8 ,dimension(1:nSubGrids_f+1) :: SX
real*8 ,dimension(1:nSubGrids_f+1) :: SY
real*8 ,dimension(1:nSubGrids_f+1) :: SZ

   integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

L_p = 1.d0 ! Plate's length
t_p = 0.05d0 ! Plate's thickness
S_p = (0.d0/180.d0)*PI ! Plate's shift angle
g_p = 0.1d0 ! Plate's gap


! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,ETA) create(SY,SZ)

!$OMP PARALLEL DO PRIVATE(j,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel loop independent private(z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) &
!$acc 				 reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector

do k=0,nz+1
	do j=0,ny+1

	! rotation wrt. AOA for blade 1
	z_aoa1 =  (Zs(k) - zc)*dcos((-AOA)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA)*PI/180.d0)
	y_aoa1 =  (Zs(k) - zc)*dsin((-AOA)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans1 = (z_aoa1 + rotor_r)           + r*dsin(S_p)
	y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - dcos(S_p)) + 0.5d0*t_p*dsin(S_p)


	! rotation wrt. AOA for blade 2
	z_aoa2 =  (Zs(k) - zc)*dcos((-AOA+180.d0)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA+180.d0)*PI/180.d0)
	y_aoa2 =  (Zs(k) - zc)*dsin((-AOA+180.d0)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA+180.d0)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans2 = (z_aoa2 + rotor_r)           + r*dsin(S_p)
	y_trans2 = (y_aoa2 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - dcos(S_p)) + 0.5d0*t_p*dsin(S_p)
	

	if( (((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (dabs(z_trans1/t_p + y_trans1/L_p) + dabs(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
	    (((z_aoa2+rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (dabs(z_trans2/t_p + y_trans2/L_p) + dabs(z_trans2/t_p - y_trans2/L_p)) .LE. 1)) then
		ETA(i,j,k)=1.D0
		az_min = MIN0(k,az_min)
		az_max = MAX0(k,az_max)
		ay_min = MIN0(j,ay_min)
		ay_max = MAX0(j,ay_max)
	else
		ETA(i,j,k)=0.D0		
	endif
	
	end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

! Setting the boundary for subgrids
	iBgnVOS = 0
	iEndVOS = nx+1
	jBgnVOS = ay_min - 5
	jEndVOS = ay_max + 5
	kBgnVOS = az_min - 5
	kEndVOS = az_max + 5


! Creating subgrids
i=iBgnVOS
!$OMP PARALLEL DO PRIVATE(j,m,n,dzg,dyg,SY,SZ,xi,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) collapse(2)
!$acc parallel loop independent private(m,n,dzg,dyg,SY,SZ,xi,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) collapse(2) gang vector

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	! Finding potential grids for subgrids
	if( (abs(ETA(i,j,k)-ETA(i,j-1,k))+abs(ETA(i,j,k)-ETA(i,j+1,k))+ &
		abs(ETA(i,j,k)-ETA(i,j,k-1))+abs(ETA(i,j,k)-ETA(i,j,k+1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k-1))+abs(ETA(i,j,k)-ETA(i,j+1,k-1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k+1))+abs(ETA(i,j,k)-ETA(i,j+1,k+1))) .GT. 0.5d0) then
			dzg = iDz(k) / nSubGrids_f
			dyg = iDy(j) / nSubGrids_f
			do m=1,nSubGrids_f+1
				SY(m) = Y(j) + (m-1) * dyg		
			end do
			do n=1,nSubGrids_f+1
				SZ(n) = Z(k) + (n-1) * dzg	
			end do

		    xi = 0.0
			do m=1,nSubGrids_f; do n=1,nSubGrids_f

				! rotation wrt. AOA for blade 1			
				z_aoa1 =  (SZ(n) - zc)*dcos((-AOA)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA)*PI/180.d0)
				y_aoa1 =  (SZ(n) - zc)*dsin((-AOA)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans1 = (z_aoa1 + rotor_r)           + r*dsin(S_p)
				y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - dcos(S_p)) + 0.5d0*t_p*dsin(S_p)
				
				! rotation wrt. AOA for blade 2			
				z_aoa2 =  (SZ(n) - zc)*dcos((-AOA+180)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA+180)*PI/180.d0)
				y_aoa2 =  (SZ(n) - zc)*dsin((-AOA+180)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA+180)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans2 = (z_aoa2 + rotor_r)           + r*dsin(S_p)
				y_trans2 = (y_aoa2 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - dcos(S_p)) + 0.5d0*t_p*dsin(S_p)
					
				if( (((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (dabs(z_trans1/t_p + y_trans1/L_p) + dabs(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
					(((z_aoa2+rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (dabs(z_trans2/t_p + y_trans2/L_p) + dabs(z_trans2/t_p - y_trans2/L_p)) .LE. 1)) then								
						xi = xi + 1
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	endif
end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

!$acc parallel loop independent collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)

end do
end do
end do   
!$acc end parallel 

!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, -1*AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine func_Cylplate_2blade


subroutine func_Cylplate_1blade()

use variables
implicit none
real*8 :: L_p, t_p, S_p, g_p
real*8 :: z_trans1, y_trans1
real*8 :: z_aoa1, y_aoa1
INTEGER:: az_min,az_max,ay_min,ay_max
integer :: l ,m ,n
real*8 :: xi
real*8 :: dxg, dyg, dzg
real*8 ,dimension(1:nSubGrids_f+1) :: SX
real*8 ,dimension(1:nSubGrids_f+1) :: SY
real*8 ,dimension(1:nSubGrids_f+1) :: SZ

   integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

L_p = 0.5313d0 ! Plate's length
t_p = 0.05d0 ! Plate's thickness
S_p = (0.d0/180.d0)*PI ! Plate's Shift angle
g_p = 0.1d0 ! Plate's gap


! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,ETA) create(SY,SZ)

!$OMP PARALLEL DO PRIVATE(j,z_aoa1,y_aoa1,z_trans1,y_trans1) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel loop independent private(z_aoa1,y_aoa1,z_trans1,y_trans1) &
!$acc reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


do k=0,nz+1
	do j=0,ny+1

	! rotation wrt. AOA for blade 1
	z_aoa1 =  (Zs(k) - zc)*dcos((-AOA)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA)*PI/180.d0)
	y_aoa1 =  (Zs(k) - zc)*dsin((-AOA)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans1 = (z_aoa1 + rotor_r)           + r*dsin(S_p)
	y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - dcos(S_p)) + 0.5d0*t_p*dsin(S_p)
	

	if( ((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. &
		(dabs(z_trans1/t_p + y_trans1/L_p) + dabs(z_trans1/t_p - y_trans1/L_p)) .LE. 1) then
		ETA(i,j,k)=1.D0
		az_min = MIN0(k,az_min)
		az_max = MAX0(k,az_max)
		ay_min = MIN0(j,ay_min)
		ay_max = MAX0(j,ay_max)
	else
		ETA(i,j,k)=0.D0		
	endif
	
	end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

! Setting the boundary for subgrids
	iBgnVOS = 0
	iEndVOS = nx+1
	jBgnVOS = ay_min - 5
	jEndVOS = ay_max + 5
	kBgnVOS = az_min - 5
	kEndVOS = az_max + 5

! Creating subgrids
i=iBgnVOS
!$OMP PARALLEL DO PRIVATE(j,m,n,dzg,dyg,SY,SZ,xi,z_aoa1,y_aoa1,z_trans1,y_trans1) collapse(2)
!$acc parallel loop independent private(m,n,dzg,dyg,SY,SZ,xi,z_aoa1,y_aoa1,z_trans1,y_trans1) collapse(2) gang vector

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	! Finding potential grids for subgrids
	if( (abs(ETA(i,j,k)-ETA(i,j-1,k))+abs(ETA(i,j,k)-ETA(i,j+1,k))+ &
		abs(ETA(i,j,k)-ETA(i,j,k-1))+abs(ETA(i,j,k)-ETA(i,j,k+1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k-1))+abs(ETA(i,j,k)-ETA(i,j+1,k-1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k+1))+abs(ETA(i,j,k)-ETA(i,j+1,k+1))) .GT. 0.5d0) then
			dzg = iDz(k) / nSubGrids_f
			dyg = iDy(j) / nSubGrids_f
			do m=1,nSubGrids_f+1
				SY(m) = Y(j) + (m-1) * dyg		
			end do
			do n=1,nSubGrids_f+1
				SZ(n) = Z(k) + (n-1) * dzg	
			end do

		    xi = 0.0
			do m=1,nSubGrids_f; do n=1,nSubGrids_f

				! rotation wrt. AOA				
				z_aoa1 =  (SZ(n) - zc)*dcos((-AOA)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA)*PI/180.d0)
				y_aoa1 =  (SZ(n) - zc)*dsin((-AOA)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans1 = (z_aoa1 + rotor_r)          + r*dsin(S_p)
				y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - dcos(S_p)) + 0.5d0*t_p*dsin(S_p)

				if( ((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. &
					(dabs(z_trans1/t_p + y_trans1/L_p) + dabs(z_trans1/t_p - y_trans1/L_p)) .LE. 1) then
						xi = xi + 1
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!$acc parallel loop independent collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)

end do
end do
end do 
!$acc end parallel    

!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, -1*AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine func_Cylplate_1blade


subroutine func_Cylplate_noblade()

use variables
implicit none
real*8 :: z_aoa1, y_aoa1
INTEGER:: az_min,az_max,ay_min,ay_max
integer :: l ,m ,n
real*8 :: xi
real*8 :: dxg, dyg, dzg
real*8 ,dimension(1:nSubGrids_f+1) :: SX
real*8 ,dimension(1:nSubGrids_f+1) :: SY
real*8 ,dimension(1:nSubGrids_f+1) :: SZ

   integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!



! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,ETA) create(SY,SZ)

!$OMP PARALLEL DO PRIVATE(j,z_aoa1,y_aoa1) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel loop independent private(z_aoa1,y_aoa1) &
!$acc reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


do k=0,nz+1
	do j=0,ny+1

	! rotation wrt. AOA for blade 1
	z_aoa1 =  (Zs(k) - zc)*dcos((-AOA)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA)*PI/180.d0)
	y_aoa1 =  (Zs(k) - zc)*dsin((-AOA)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA)*PI/180.d0)


	if( ((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2) then
		ETA(i,j,k)=1.D0
		az_min = MIN0(k,az_min)
		az_max = MAX0(k,az_max)
		ay_min = MIN0(j,ay_min)
		ay_max = MAX0(j,ay_max)
	else
		ETA(i,j,k)=0.D0		
	endif
	
	end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

! Setting the boundary for subgrids
	iBgnVOS = 0
	iEndVOS = nx+1
	jBgnVOS = ay_min - 5
	jEndVOS = ay_max + 5
	kBgnVOS = az_min - 5
	kEndVOS = az_max + 5

! Creating subgrids
i=iBgnVOS
!$OMP PARALLEL DO PRIVATE(j,m,n,dzg,dyg,SY,SZ,xi,z_aoa1,y_aoa1) collapse(2)
!$acc parallel loop independent private(m,n,dzg,dyg,SY,SZ,xi,z_aoa1,y_aoa1) collapse(2) gang vector

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	! Finding potential grids for subgrids
	if( (abs(ETA(i,j,k)-ETA(i,j-1,k))+abs(ETA(i,j,k)-ETA(i,j+1,k))+ &
		abs(ETA(i,j,k)-ETA(i,j,k-1))+abs(ETA(i,j,k)-ETA(i,j,k+1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k-1))+abs(ETA(i,j,k)-ETA(i,j+1,k-1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k+1))+abs(ETA(i,j,k)-ETA(i,j+1,k+1))) .GT. 0.5d0) then
			dzg = iDz(k) / nSubGrids_f
			dyg = iDy(j) / nSubGrids_f
			do m=1,nSubGrids_f+1
				SY(m) = Y(j) + (m-1) * dyg		
			end do
			do n=1,nSubGrids_f+1
				SZ(n) = Z(k) + (n-1) * dzg	
			end do

		    xi = 0.0
			do m=1,nSubGrids_f; do n=1,nSubGrids_f

				! rotation wrt. AOA				
				z_aoa1 =  (SZ(n) - zc)*dcos((-AOA)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA)*PI/180.d0)
				y_aoa1 =  (SZ(n) - zc)*dsin((-AOA)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA)*PI/180.d0)


				if( ((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2) then
						xi = xi + 1
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!$acc parallel loop independent collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)

end do
end do
end do 
!$acc end parallel    

!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, -1*AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine func_Cylplate_noblade


subroutine func_Cylinder()

use variables
implicit none
real*8 :: L_p, t_p, B_p, g_p
real*8 :: z_trans, y_trans
INTEGER:: az_min,az_max,ay_min,ay_max
real*8 :: DISTANCE,SDIST
real*8 :: DIAGONAL
integer :: l ,m ,n
real*8 :: xi
real*8 :: dxg, dyg, dzg
real*8 ,dimension(1:nSubGrids_f+1) :: SX
real*8 ,dimension(1:nSubGrids_f+1) :: SY
real*8 ,dimension(1:nSubGrids_f+1) :: SZ

   integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ)

!$OMP PARALLEL DO PRIVATE(j) REDUCTION(max : az_max,ay_max) &
!$OMP			  REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel
!$acc loop independent reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector

do k=0,nz+1
	do j=0,ny+1


	if( ((Zs(k) - zc_t)**2 + (Ys(j) - yc_t)**2) .LE. r**2) then
		ETA(i,j,k)=1.D0
		az_min = MIN0(k,az_min)
		az_max = MAX0(k,az_max)
		ay_min = MIN0(j,ay_min)
		ay_max = MAX0(j,ay_max)
	else
		ETA(i,j,k)=0.D0		
	endif
	
	end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

! Setting the boundary for subgrids
	iBgnVOS = 0
	iEndVOS = nx+1
	jBgnVOS = ay_min - 5
	jEndVOS = ay_max + 5
	kBgnVOS = az_min - 5
	kEndVOS = az_max + 5

! Creating subgrids
i=iBgnVOS
!$OMP PARALLEL DO PRIVATE(j,m,n,dzg,dyg,SY,SZ,xi) collapse(2)
!$acc parallel loop independent private(m,n,dzg,dyg,SY,SZ,xi) collapse(2) gang vector

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	! Finding potential grids for subgrids
	if( (abs(ETA(i,j,k)-ETA(i,j-1,k))+abs(ETA(i,j,k)-ETA(i,j+1,k))+ &
		abs(ETA(i,j,k)-ETA(i,j,k-1))+abs(ETA(i,j,k)-ETA(i,j,k+1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k-1))+abs(ETA(i,j,k)-ETA(i,j+1,k-1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k+1))+abs(ETA(i,j,k)-ETA(i,j+1,k+1))) .GT. 0.5d0) then
			!$acc loop seq
			do m=1,nSubGrids_f+1
			!$acc loop seq
			do n=1,nSubGrids_f+1
				dzg = iDz(k) / nSubGrids_f
				dyg = iDy(j) / nSubGrids_f
				SZ(n) = Z(k) + (n-1) * dzg
				SY(m) = Y(j) + (m-1) * dyg		
			end do; end do

		    xi = 0.0
			!$acc loop seq
			do m=1,nSubGrids_f
			!$acc loop seq
			do n=1,nSubGrids_f
				if( ((SZ(n) - zc_t)**2 + (SY(m) - yc_t)**2) .LE. r**2) then
						xi = xi + 1
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!$acc parallel loop independent collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)

end do
end do
end do    
!$acc end parallel    

!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,x,y,z,u_solid,v_solid,w_solid'
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,6(3X,F12.7))') time, xc_t, yc_t, zc_t, u_solid, v_solid, w_solid
	endif
    close(61)

end subroutine func_Cylinder


subroutine func_Sphere()
use variables
implicit none
real*4 :: DISTANCE,SDIST
real*4 :: DIAGONAL
integer :: l ,m ,n
real*4 :: xi
real*4 :: dxg, dyg, dzg
real*4 ,dimension(1:nSubGrids_f+1) :: SX
real*4 ,dimension(1:nSubGrids_f+1) :: SY
real*4 ,dimension(1:nSubGrids_f+1) :: SZ

   integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

	iBgnVOS = 1
	iEndVOS = nx
	jBgnVOS = 1
	jEndVOS = ny
	kBgnVOS = 1
	kEndVOS = nz

!$acc data present(X,Y,Z,iDx,iDy,iDz,ETA) create(SX,SY,SZ)

!$OMP PARALLEL DO PRIVATE(i,j,DIAGONAL,DISTANCE,SDIST,l,m,n,dxg,dyg,dzg,SX,SY,SZ,xi) collapse(3)
!$acc parallel loop independent private(DIAGONAL,DISTANCE,SDIST,l,m,n,dxg,dyg,dzg,SX,SY,SZ,xi) collapse(3) gang vector
do k=kBgnVOS,kEndVOS 
  do j=jBgnVOS,jEndVOS
    do i=iBgnVOS,iEndVOS

         DIAGONAL = sqrt( iDx(i)*iDx(i) + iDy(j)*iDy(j) + iDz(k)*iDz(k) ) / 2.D0
         DISTANCE = sqrt((((X(i)+X(i+1))/2.D0)-xc_t)**2.D0+ &
                         (((Y(j)+Y(j+1))/2.D0)-yc_t)**2.D0+ &
                         (((Z(k)+Z(k+1))/2.D0)-zc_t)**2.D0  )

           if( abs(DISTANCE - r) < DIAGONAL ) then
              !$acc loop seq
              do l=1,nSubGrids_f+1,1
			   !$acc loop seq
               do m=1,nSubGrids_f+1,1
			    !$acc loop seq
                do n=1,nSubGrids_f+1,1
					dxg = iDx(i) / nSubGrids_f
					dyg = iDy(j) / nSubGrids_f
					dzg = iDz(k) / nSubGrids_f
					SX(n)=X(i)+(n-1)*dxg
					SY(m)=Y(j)+(m-1)*dyg
					SZ(l)=Z(k)+(l-1)*dzg
                end do
               end do
              end do

              xi = 0.0
              !$acc loop seq
			  do l=1,nSubGrids_f
			    !$acc loop seq
                do m=1,nSubGrids_f
				  !$acc loop seq
                  do n=1,nSubGrids_f
                    
                    SDIST=SQRT((((SX(n)+SX(n+1))/2.D0)-xc_t)**2.D0+ &
                               (((SY(m)+SY(m+1))/2.D0)-yc_t)**2.D0+ &
                               (((SZ(l)+SZ(l+1))/2.D0)-zc_t)**2.D0  )
                    if(SDIST <= r) then
                      xi = xi + 1
                    end if

                  end do
                end do
              end do   
              ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f*nSubGrids_f)
              
              !print*,'ETA(i,j,k) = ',ETA(i,j,k),'xi = ',xi,i,j,k
           
           else if (DISTANCE <= r) then
             ETA(i,j,k)=1.D0
           else
             ETA(i,j,k)=0.D0
           end if
                
    end do
  end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,x,y,z,u_solid,v_solid,w_solid'
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,6(3X,F12.7))') time, xc_t, yc_t, zc_t, u_solid, v_solid, w_solid
	endif
    close(61)

end subroutine func_Sphere



! subroutine func_Cylplate_2blade()

! use variables
! implicit none
! real*8 :: L_p, t_p, B_p, g_p
! real*8 :: z_trans_, y_trans_, z_trans1, y_trans1, z_trans2, y_trans2
! real*8 :: z_aoa1, y_aoa1, z_aoa2, y_aoa2
! INTEGER, DIMENSION(0:1) :: az, ay
! integer :: l ,m ,n
! real*8 :: xi
! real*8 :: dxg, dyg, dzg
! real*8 ,dimension(1:nSubGrids_f+1) :: SX
! real*8 ,dimension(1:nSubGrids_f+1) :: SY
! real*8 ,dimension(1:nSubGrids_f+1) :: SZ

   ! integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   ! !---------------------------------------------------!
   ! !    LOCAL VARIABLES                                !
   ! !---------------------------------------------------!

! L_p = 1.d0 ! Plate's length
! t_p = 0.2d0 ! Plate's thickness
! B_p = (0.d0/180.d0)*PI ! Plate's Beta angle
! g_p = 0.1d0 ! Plate's gap


! ! Create coarse ETA (only 0 or 1)
! i=0
! az(0) = nz ; az(1) = 0
! ay(0) = ny ; ay(1) = 0
! do k=0,nz+1
	! do j=0,ny+1

	! ! rotation wrt. AOA for blade 1
	! z_aoa1 =  (Zs(k) - zc)*dcos((-AOA)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA)*PI/180.d0)
	! y_aoa1 =  (Zs(k) - zc)*dsin((-AOA)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA)*PI/180.d0)

	! ! z and y transformation for flat plate
	! z_trans_ = z_aoa1 + rotor_r
	! y_trans_ = y_aoa1 + (g_p+r)+0.5d0*L_p
	
	! z_trans1 = z_trans_*dcos(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*dsin(B_p)
	! y_trans1 = z_trans_*dsin(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*dcos(B_p) + (g_p+r)+0.5d0*L_p


	! ! rotation wrt. AOA for blade 2
	! z_aoa2 =  (Zs(k) - zc)*dcos((-AOA+180.d0)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA+180.d0)*PI/180.d0)
	! y_aoa2 =  (Zs(k) - zc)*dsin((-AOA+180.d0)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA+180.d0)*PI/180.d0)

	! ! z and y transformation for flat plate
	! z_trans_ = z_aoa2 + rotor_r
	! y_trans_ = y_aoa2 + (g_p+r)+0.5d0*L_p
	
	! z_trans2 = z_trans_*dcos(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*dsin(B_p)
	! y_trans2 = z_trans_*dsin(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*dcos(B_p) + (g_p+r)+0.5d0*L_p
	

	! if( (((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (dabs(z_trans1/t_p + y_trans1/L_p) + dabs(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
	    ! (((z_aoa2+rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (dabs(z_trans2/t_p + y_trans2/L_p) + dabs(z_trans2/t_p - y_trans2/L_p)) .LE. 1)) then
		! ETA(i,j,k)=1.D0
		! if (k < az(0)) then
			! az(0) = k
		! else if (k > az(1)) then
			! az(1) = k
		! endif
		! if (j < ay(0)) then
			! ay(0) = j
		! else if (j > ay(1)) then
			! ay(1) = j
		! endif

	! else
		! ETA(i,j,k)=0.D0		
	! endif
	
	! end do
! end do

! ! Setting the boundary for subgrids
	! iBgnVOS = 0
	! iEndVOS = nx+1
	! jBgnVOS = ay(0) - 5
	! jEndVOS = ay(1) + 5
	! kBgnVOS = az(0) - 5
	! kEndVOS = az(1) + 5

! ! Creating subgrids
! i=iBgnVOS
! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
	! ! Finding potential grids for subgrids
	! if( (abs(ETA(i,j,k)-ETA(i,j-1,k))+abs(ETA(i,j,k)-ETA(i,j+1,k))+ &
		! abs(ETA(i,j,k)-ETA(i,j,k-1))+abs(ETA(i,j,k)-ETA(i,j,k+1))+ &
		! abs(ETA(i,j,k)-ETA(i,j-1,k-1))+abs(ETA(i,j,k)-ETA(i,j+1,k-1))+ &
		! abs(ETA(i,j,k)-ETA(i,j-1,k+1))+abs(ETA(i,j,k)-ETA(i,j+1,k+1))) .GT. 0.5d0) then
			! dzg = iDz(k) / nSubGrids_f
			! dyg = iDy(j) / nSubGrids_f
			! do m=1,nSubGrids_f+1
				! SY(m) = Y(j) + (m-1) * dyg		
			! end do
			! do n=1,nSubGrids_f+1
				! SZ(n) = Z(k) + (n-1) * dzg	
			! end do

		    ! xi = 0.0
			! do m=1,nSubGrids_f; do n=1,nSubGrids_f

				! ! rotation wrt. AOA for blade 1			
				! z_aoa1 =  (SZ(n) - zc)*dcos((-AOA)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA)*PI/180.d0)
				! y_aoa1 =  (SZ(n) - zc)*dsin((-AOA)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA)*PI/180.d0)

				! ! z and y transformation for flat plate
				! z_trans_ = z_aoa1 + rotor_r
				! y_trans_ = y_aoa1 + (g_p+r)+0.5d0*L_p
	
				! z_trans1 = z_trans_*dcos(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*dsin(B_p)
				! y_trans1 = z_trans_*dsin(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*dcos(B_p) + (g_p+r)+0.5d0*L_p
				
				! ! rotation wrt. AOA for blade 2			
				! z_aoa2 =  (SZ(n) - zc)*dcos((-AOA+180)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA+180)*PI/180.d0)
				! y_aoa2 =  (SZ(n) - zc)*dsin((-AOA+180)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA+180)*PI/180.d0)

				! ! z and y transformation for flat plate
				! z_trans_ = z_aoa2 + rotor_r
				! y_trans_ = y_aoa2 + (g_p+r)+0.5d0*L_p
	
				! z_trans2 = z_trans_*dcos(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*dsin(B_p)
				! y_trans2 = z_trans_*dsin(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*dcos(B_p) + (g_p+r)+0.5d0*L_p
					
				! if( (((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (dabs(z_trans1/t_p + y_trans1/L_p) + dabs(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
					! (((z_aoa2+rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (dabs(z_trans2/t_p + y_trans2/L_p) + dabs(z_trans2/t_p - y_trans2/L_p)) .LE. 1)) then								
						! xi = xi + 1
				! endif
			! end do; end do
		   
		  
		  ! ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	! endif
! end do
! end do		


! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
! do i=iBgnVOS+1,iEndVOS
	! ETA(i,j,k) = ETA(iBgnVOS,j,k)

! end do
! end do
! end do    

    ! if(myid==master .AND. istep == 0)then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,*)'                 '
        ! write(61,*)'SOLID MOTION'
		! write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        ! write(61,*)'                 '
    ! else if (myid==master) then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, -1*AOA, rotor_omega  ! for rotation
	! endif
    ! close(61)

! end subroutine func_Cylplate_2blade


! subroutine func_Cylplate_1blade()

! use variables
! implicit none
! real*8 :: L_p, t_p, B_p, g_p
! real*8 :: z_trans_, y_trans_, z_trans1, y_trans1
! real*8 :: z_aoa1, y_aoa1
! INTEGER, DIMENSION(0:1) :: az, ay
! integer :: l ,m ,n
! real*8 :: xi
! real*8 :: dxg, dyg, dzg
! real*8 ,dimension(1:nSubGrids_f+1) :: SX
! real*8 ,dimension(1:nSubGrids_f+1) :: SY
! real*8 ,dimension(1:nSubGrids_f+1) :: SZ

   ! integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   ! !---------------------------------------------------!
   ! !    LOCAL VARIABLES                                !
   ! !---------------------------------------------------!

! L_p = 1.d0 ! Plate's length
! t_p = 0.05d0 ! Plate's thickness
! B_p = (-30.d0/180.d0)*PI ! Plate's Beta angle
! g_p = 0.1d0 ! Plate's gap


! ! Create coarse ETA (only 0 or 1)
! i=0
! az(0) = nz ; az(1) = 0
! ay(0) = ny ; ay(1) = 0
! do k=0,nz+1
	! do j=0,ny+1

	! ! rotation wrt. AOA for blade 1
	! z_aoa1 =  (Zs(k) - zc)*dcos((-AOA)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA)*PI/180.d0)
	! y_aoa1 =  (Zs(k) - zc)*dsin((-AOA)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA)*PI/180.d0)

	! ! z and y transformation for flat plate
	! z_trans_ = z_aoa1 + rotor_r
	! y_trans_ = y_aoa1 + (g_p+r)+0.5d0*L_p
	
	! z_trans1 = z_trans_*dcos(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*dsin(B_p)
	! y_trans1 = z_trans_*dsin(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*dcos(B_p) + (g_p+r)+0.5d0*L_p
	

	! if( ((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. &
		! (dabs(z_trans1/t_p + y_trans1/L_p) + dabs(z_trans1/t_p - y_trans1/L_p)) .LE. 1) then
		! ETA(i,j,k)=1.D0
		! if (k < az(0)) then
			! az(0) = k
		! else if (k > az(1)) then
			! az(1) = k
		! endif
		! if (j < ay(0)) then
			! ay(0) = j
		! else if (j > ay(1)) then
			! ay(1) = j
		! endif

	! else
		! ETA(i,j,k)=0.D0		
	! endif
	
	! end do
! end do

! ! Setting the boundary for subgrids
	! iBgnVOS = 0
	! iEndVOS = nx+1
	! jBgnVOS = ay(0) - 5
	! jEndVOS = ay(1) + 5
	! kBgnVOS = az(0) - 5
	! kEndVOS = az(1) + 5

! ! Creating subgrids
! i=iBgnVOS
! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
	! ! Finding potential grids for subgrids
	! if( (abs(ETA(i,j,k)-ETA(i,j-1,k))+abs(ETA(i,j,k)-ETA(i,j+1,k))+ &
		! abs(ETA(i,j,k)-ETA(i,j,k-1))+abs(ETA(i,j,k)-ETA(i,j,k+1))+ &
		! abs(ETA(i,j,k)-ETA(i,j-1,k-1))+abs(ETA(i,j,k)-ETA(i,j+1,k-1))+ &
		! abs(ETA(i,j,k)-ETA(i,j-1,k+1))+abs(ETA(i,j,k)-ETA(i,j+1,k+1))) .GT. 0.5d0) then
			! dzg = iDz(k) / nSubGrids_f
			! dyg = iDy(j) / nSubGrids_f
			! do m=1,nSubGrids_f+1
				! SY(m) = Y(j) + (m-1) * dyg		
			! end do
			! do n=1,nSubGrids_f+1
				! SZ(n) = Z(k) + (n-1) * dzg	
			! end do

		    ! xi = 0.0
			! do m=1,nSubGrids_f; do n=1,nSubGrids_f

				! ! rotation wrt. AOA				
				! z_aoa1 =  (SZ(n) - zc)*dcos((-AOA)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA)*PI/180.d0)
				! y_aoa1 =  (SZ(n) - zc)*dsin((-AOA)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA)*PI/180.d0)

				! ! z and y transformation for flat plate
				! z_trans_ = z_aoa1 + rotor_r
				! y_trans_ = y_aoa1 + (g_p+r)+0.5d0*L_p
	
				! z_trans1 = z_trans_*dcos(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*dsin(B_p)
				! y_trans1 = z_trans_*dsin(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*dcos(B_p) + (g_p+r)+0.5d0*L_p

				! if( ((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. &
					! (dabs(z_trans1/t_p + y_trans1/L_p) + dabs(z_trans1/t_p - y_trans1/L_p)) .LE. 1) then
						! xi = xi + 1
				! endif
			! end do; end do
		   
		  
		  ! ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	! endif
! end do
! end do		


! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
! do i=iBgnVOS+1,iEndVOS
	! ETA(i,j,k) = ETA(iBgnVOS,j,k)

! end do
! end do
! end do    

    ! if(myid==master .AND. istep == 0)then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,*)'                 '
        ! write(61,*)'SOLID MOTION'
		! write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        ! write(61,*)'                 '
    ! else if (myid==master) then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, -1*AOA, rotor_omega  ! for rotation
	! endif
    ! close(61)

! end subroutine func_Cylplate_1blade



! subroutine func_Cylinder()

! use variables
! implicit none
! real*4 :: DISTANCE,SDIST
! real*4 :: DIAGONAL
! integer :: l ,m ,n
! real*4 :: xi
! real*4 :: dxg, dyg, dzg
! real*4 :: solid_volume
! real*4 ,dimension(1:nSubGrids_f+1) :: SX
! real*4 ,dimension(1:nSubGrids_f+1) :: SY
! real*4 ,dimension(1:nSubGrids_f+1) :: SZ

   ! integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   ! !---------------------------------------------------!
   ! !    LOCAL VARIABLES                                !
   ! !---------------------------------------------------!

	! iBgnVOS = 0
	! iEndVOS = nx+1
	! jBgnVOS = 0
	! jEndVOS = ny+1
	! kBgnVOS = 0
	! kEndVOS = nz+1


! i=iBgnVOS
! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS

	 ! DIAGONAL = sqrt(  iDx(i)*iDx(i) + iDy(j)*iDy(j) + iDz(k)*iDz(k) ) / 2.D0

	 ! DISTANCE = sqrt( (( (Z(k)+Z(k+1)) / 2.D0 ) -   zc )  ** 2.D0 &
                 ! +  (( (Y(j)+Y(j+1)) / 2.D0 ) - yc ) ** 2.D0 )  

	   ! if( abs(DISTANCE - r) < DIAGONAL ) then
		  
		  ! xi = 0.0
		  
		 
			! do m=1,nSubGrids_f+1; do n=1,nSubGrids_f+1
				! dzg = iDz(k) / nSubGrids_f
				! dyg = iDy(j) / nSubGrids_f
				! SZ(n) = Z(k) + (n-1) * dzg
				! SY(m) = Y(j) + (m-1) * dyg		
			! end do; end do
		  
		  

		  
			! do m=1,nSubGrids_f; do n=1,nSubGrids_f
				  ! SDIST = sqrt( (( (SZ(n)+SZ(n+1)) / 2.D0 ) - zc )  ** 2.D0 &
                    ! +   (( (SY(m)+SY(m+1)) / 2.D0 ) - yc ) ** 2.D0 )
				  ! if(SDIST <= r) then
					! xi = xi + 1
				  ! end if
			! end do; end do
		   
		  
		  ! ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
		  
		  ! !print*,'ETA(i,j,k) = ',ETA(i,j,k),'xi = ',xi,i,j,k
	   
	   ! else if (DISTANCE > r) then
		 ! ETA(i,j,k)=0.D0
	   ! else
		 ! ETA(i,j,k)=1.D0
	   ! end if
			
! end do
! end do


! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
! do i=iBgnVOS+1,iEndVOS
	! ETA(i,j,k) = ETA(iBgnVOS,j,k)

! end do
! end do
! end do    


! end subroutine func_Cylinder






