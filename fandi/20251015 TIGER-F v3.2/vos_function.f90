! 3 Nov 2024 - FDS


!********************************** Extruded solids ***********************************



subroutine func_darrius_3blade()

use variables
implicit none

integer:: az_min,az_max,ay_min,ay_max
integer:: m, n
real*8 :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

real*8 :: c_p, t_p, s_p
real*8 :: z_trans1, y_trans1, z_trans2, y_trans2, z_trans3, y_trans3
real*8 :: z_aoa1, y_aoa1, z_aoa2, y_aoa2, z_aoa3, y_aoa3

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

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) &
!$acc 	   reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


do k=1,nz
	do j=1,ny

		Zf = Zs(k) - z0
		Yf = Ys(j) - y0

	! rotation wrt. AOA for blade 1
	z_aoa1 =  Zf*COS((AOA-90.d0)*PI/180.d0)+Yf*SIN((AOA-90.d0)*PI/180.d0)
	y_aoa1 = -Zf*SIN((AOA-90.d0)*PI/180.d0)+Yf*COS((AOA-90.d0)*PI/180.d0)

	! z and y transformation
	z_trans1 = z_aoa1 + s_p
	y_trans1 = y_aoa1 - rotor_r


	! rotation wrt. AOA for blade 2
	z_aoa2 =  Zf*COS((AOA-90.d0+120.d0)*PI/180.d0)+Yf*SIN((AOA-90.d0+120.d0)*PI/180.d0)
	y_aoa2 = -Zf*SIN((AOA-90.d0+120.d0)*PI/180.d0)+Yf*COS((AOA-90.d0+120.d0)*PI/180.d0)

	! z and y transformation
	z_trans2 = z_aoa2 + s_p
	y_trans2 = y_aoa2 - rotor_r


	! rotation wrt. AOA for blade 3
	z_aoa3 =  Zf*COS((AOA-90.d0+240.d0)*PI/180.d0)+Yf*SIN((AOA-90.d0+240.d0)*PI/180.d0)
	y_aoa3 = -Zf*SIN((AOA-90.d0+240.d0)*PI/180.d0)+Yf*COS((AOA-90.d0+240.d0)*PI/180.d0)

	! z and y transformation
	z_trans3 = z_aoa3 + s_p
	y_trans3 = y_aoa3 - rotor_r

	
	if( ((y_trans1 - 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .LE. 0.) .AND. &
	     (y_trans1 + 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .GE. 0.)) .OR. &
		((y_trans2 - 5.*t_p*c_p*(0.2969*sqrt(z_trans2/c_p) - 0.126*z_trans2/c_p - 0.3516*(z_trans2/c_p)**2 + 0.2843*(z_trans2/c_p)**3 - 0.1015*(z_trans2/c_p)**4) .LE. 0.) .AND. &
		 (y_trans2 + 5.*t_p*c_p*(0.2969*sqrt(z_trans2/c_p) - 0.126*z_trans2/c_p - 0.3516*(z_trans2/c_p)**2 + 0.2843*(z_trans2/c_p)**3 - 0.1015*(z_trans2/c_p)**4) .GE. 0.)) .OR. &
		((y_trans3 - 5.*t_p*c_p*(0.2969*sqrt(z_trans3/c_p) - 0.126*z_trans3/c_p - 0.3516*(z_trans3/c_p)**2 + 0.2843*(z_trans3/c_p)**3 - 0.1015*(z_trans3/c_p)**4) .LE. 0.) .AND. &
		 (y_trans3 + 5.*t_p*c_p*(0.2969*sqrt(z_trans3/c_p) - 0.126*z_trans3/c_p - 0.3516*(z_trans3/c_p)**2 + 0.2843*(z_trans3/c_p)**3 - 0.1015*(z_trans3/c_p)**4) .GE. 0.)) ) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1


! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) collapse(2) gang

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then
            
			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    ETA(i,j,k) = 0.d0
            !$acc loop private(n,Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

					Zf = SZ(n) - z0
					Yf = SY(m) - y0

				! rotation wrt. AOA for blade 1
				z_aoa1 =  Zf*COS((AOA-90.d0)*PI/180.d0)+Yf*SIN((AOA-90.d0)*PI/180.d0)
				y_aoa1 = -Zf*SIN((AOA-90.d0)*PI/180.d0)+Yf*COS((AOA-90.d0)*PI/180.d0)

				! z and y transformation
				z_trans1 = z_aoa1 + s_p
				y_trans1 = y_aoa1 - rotor_r


				! rotation wrt. AOA for blade 2
				z_aoa2 =  Zf*COS((AOA-90.d0+120.d0)*PI/180.d0)+Yf*SIN((AOA-90.d0+120.d0)*PI/180.d0)
				y_aoa2 = -Zf*SIN((AOA-90.d0+120.d0)*PI/180.d0)+Yf*COS((AOA-90.d0+120.d0)*PI/180.d0)

				! z and y transformation
				z_trans2 = z_aoa2 + s_p
				y_trans2 = y_aoa2 - rotor_r


				! rotation wrt. AOA for blade 3
				z_aoa3 =  Zf*COS((AOA-90.d0+240.d0)*PI/180.d0)+Yf*SIN((AOA-90.d0+240.d0)*PI/180.d0)
				y_aoa3 = -Zf*SIN((AOA-90.d0+240.d0)*PI/180.d0)+Yf*COS((AOA-90.d0+240.d0)*PI/180.d0)

				! z and y transformation
				z_trans3 = z_aoa3 + s_p
				y_trans3 = y_aoa3 - rotor_r				
					
		if( ((y_trans1 - 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .LE. 0.) .AND. &
			 (y_trans1 + 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .GE. 0.)) .OR. &
			((y_trans2 - 5.*t_p*c_p*(0.2969*sqrt(z_trans2/c_p) - 0.126*z_trans2/c_p - 0.3516*(z_trans2/c_p)**2 + 0.2843*(z_trans2/c_p)**3 - 0.1015*(z_trans2/c_p)**4) .LE. 0.) .AND. &
			 (y_trans2 + 5.*t_p*c_p*(0.2969*sqrt(z_trans2/c_p) - 0.126*z_trans2/c_p - 0.3516*(z_trans2/c_p)**2 + 0.2843*(z_trans2/c_p)**3 - 0.1015*(z_trans2/c_p)**4) .GE. 0.)) .OR. &
			((y_trans3 - 5.*t_p*c_p*(0.2969*sqrt(z_trans3/c_p) - 0.126*z_trans3/c_p - 0.3516*(z_trans3/c_p)**2 + 0.2843*(z_trans3/c_p)**3 - 0.1015*(z_trans3/c_p)**4) .LE. 0.) .AND. &
			 (y_trans3 + 5.*t_p*c_p*(0.2969*sqrt(z_trans3/c_p) - 0.126*z_trans3/c_p - 0.3516*(z_trans3/c_p)**2 + 0.2843*(z_trans3/c_p)**3 - 0.1015*(z_trans3/c_p)**4) .GE. 0.)) ) then							
					!$acc atomic update
					ETA(i,j,k) = ETA(i,j,k) + 1.d0
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!	If there are traces of ETA, it means that the initial condition was not set correctly, or
!	the solid diplacement between timesteps is too large
!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS-2,kEndVOS+2
do j=jBgnVOS-2,jEndVOS+2
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
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
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine func_darrius_3blade



subroutine func_Cylplate_3blade()

use variables
implicit none

integer:: az_min,az_max,ay_min,ay_max
integer:: m, n
real*8 :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

real*8 :: L_p, t_p, S_p, g_p
real*8 :: z_trans1, y_trans1, z_trans2, y_trans2, z_trans3, y_trans3
real*8 :: z_aoa1, y_aoa1, z_aoa2, y_aoa2, z_aoa3, y_aoa3

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

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) &
!$acc 			   reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector

do k=1,nz
	do j=1,ny
	
		Zf = Zs(k) - z0
		Yf = Ys(j) - y0
	
	! 1st blade initial position for AOA=0 is at Z+
	! rotation wrt. AOA for blade 1
	z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
	y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans1 = (z_aoa1 - rotor_r)           + r*SIN(S_p)
	y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)


	! rotation wrt. AOA for blade 2
	z_aoa2 =  Zf*COS((AOA+120.d0)*PI/180.d0)+Yf*SIN((AOA+120.d0)*PI/180.d0)
	y_aoa2 = -Zf*SIN((AOA+120.d0)*PI/180.d0)+Yf*COS((AOA+120.d0)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans2 = (z_aoa2 - rotor_r)           + r*SIN(S_p)
	y_trans2 = (y_aoa2 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)
	

	! rotation wrt. AOA for blade 3
	z_aoa3 =  Zf*COS((AOA+240.d0)*PI/180.d0)+Yf*SIN((AOA+240.d0)*PI/180.d0)
	y_aoa3 = -Zf*SIN((AOA+240.d0)*PI/180.d0)+Yf*COS((AOA+240.d0)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans3 = (z_aoa3 - rotor_r)           + r*SIN(S_p)
	y_trans3 = (y_aoa3 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)


	if( (((z_aoa1-rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
	    (((z_aoa2-rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (ABS(z_trans2/t_p + y_trans2/L_p) + ABS(z_trans2/t_p - y_trans2/L_p)) .LE. 1) .OR. &
	    (((z_aoa3-rotor_r)**2 + y_aoa3**2) .LE. r**2 .OR. (ABS(z_trans3/t_p + y_trans3/L_p) + ABS(z_trans3/t_p - y_trans3/L_p)) .LE. 1)) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1


! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) collapse(2) gang

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then

			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    ETA(i,j,k) = 0.d0
            !$acc loop private(Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2,z_aoa3,y_aoa3,z_trans3,y_trans3) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

					Zf = SZ(n) - z0
					Yf = SY(m) - y0

				! rotation wrt. AOA for blade 1			
				z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
				y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans1 = (z_aoa1 - rotor_r)           + r*SIN(S_p)
				y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)


				! rotation wrt. AOA for blade 2
				z_aoa2 =  Zf*COS((AOA+120.d0)*PI/180.d0)+Yf*SIN((AOA+120.d0)*PI/180.d0)
				y_aoa2 = -Zf*SIN((AOA+120.d0)*PI/180.d0)+Yf*COS((AOA+120.d0)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans2 = (z_aoa2 - rotor_r)           + r*SIN(S_p)
				y_trans2 = (y_aoa2 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)
	

				! rotation wrt. AOA for blade 3
				z_aoa3 =  Zf*COS((AOA+240.d0)*PI/180.d0)+Yf*SIN((AOA+240.d0)*PI/180.d0)
				y_aoa3 = -Zf*SIN((AOA+240.d0)*PI/180.d0)+Yf*COS((AOA+240.d0)*PI/180.d0)
	
				! z and y transformation for flat plate
				z_trans3 = (z_aoa3 - rotor_r)           + r*SIN(S_p)
				y_trans3 = (y_aoa3 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)
					
			if( (((z_aoa1-rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
				(((z_aoa2-rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (ABS(z_trans2/t_p + y_trans2/L_p) + ABS(z_trans2/t_p - y_trans2/L_p)) .LE. 1) .OR. &
				(((z_aoa3-rotor_r)**2 + y_aoa3**2) .LE. r**2 .OR. (ABS(z_trans3/t_p + y_trans3/L_p) + ABS(z_trans3/t_p - y_trans3/L_p)) .LE. 1)) then							
					!$acc atomic update
					ETA(i,j,k) = ETA(i,j,k) + 1.d0
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

!	If there are traces of ETA, it means that the initial condition was not set correctly, or
!	the solid diplacement between timesteps is too large
!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS-2,kEndVOS+2
do j=jBgnVOS-2,jEndVOS+2
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
end do
end do
end do   
!$acc end parallel
!$OMP END PARALLEL DO


	!-------------- Calculate moment of inertia -------------------

		if( istep == 0 )then
		
			I_moment = 0.d0
		
			!$OMP PARALLEL DO PRIVATE(k) REDUCTION(+: I_moment) collapse(3)
			!$acc parallel loop reduction(+: I_moment) collapse(3) gang vector
			do k = kBgnVOS-2,kEndVOS+2
				do j = jBgnVOS-2,jEndVOS+2
					do i=1,nx
						I_moment = I_moment+den_sol*ETA(i,j,k)*((ys(j) - y0)**2 +(zs(k)-z0)**2)*iDx(i)*iDy(j)*iDz(k)
					end do
				end do
			end do
			!$acc end parallel
			!$OMP END PARALLEL DO 
		
		end if 

		!write(*,*) 'Moment of inertia = ', I_moment



!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine func_Cylplate_3blade


subroutine func_Cylplate_2blade()

use variables
implicit none

integer:: az_min,az_max,ay_min,ay_max
integer:: m, n
real*8 :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

real*8 :: L_p, t_p, S_p, g_p
real*8 :: z_trans1, y_trans1, z_trans2, y_trans2
real*8 :: z_aoa1, y_aoa1, z_aoa2, y_aoa2

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

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) &
!$acc 	   reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector

do k=1,nz
	do j=1,ny

		Zf = Zs(k) - z0
		Yf = Ys(j) - y0
	
	! 1st blade initial position for AOA=0 is at Z+
	! rotation wrt. AOA for blade 1
	z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
	y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans1 = (z_aoa1 - rotor_r)           + r*SIN(S_p)
	y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)


	! rotation wrt. AOA for blade 2
	z_aoa2 =  Zf*COS((AOA+180.d0)*PI/180.d0)+Yf*SIN((AOA+180.d0)*PI/180.d0)
	y_aoa2 = -Zf*SIN((AOA+180.d0)*PI/180.d0)+Yf*COS((AOA+180.d0)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans2 = (z_aoa2 - rotor_r)           + r*SIN(S_p)
	y_trans2 = (y_aoa2 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)
	

	if( (((z_aoa1-rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
	    (((z_aoa2-rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (ABS(z_trans2/t_p + y_trans2/L_p) + ABS(z_trans2/t_p - y_trans2/L_p)) .LE. 1)) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1


! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) collapse(2) gang

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then

			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    ETA(i,j,k) = 0.d0
            !$acc loop private(n,Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

					Zf = SZ(n) - z0
					Yf = SY(m) - y0

				! rotation wrt. AOA for blade 1			
				z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
				y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans1 = (z_aoa1 - rotor_r)           + r*SIN(S_p)
				y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)


				! rotation wrt. AOA for blade 2
				z_aoa2 =  Zf*COS((AOA+180.d0)*PI/180.d0)+Yf*SIN((AOA+180.d0)*PI/180.d0)
				y_aoa2 = -Zf*SIN((AOA+180.d0)*PI/180.d0)+Yf*COS((AOA+180.d0)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans2 = (z_aoa2 - rotor_r)           + r*SIN(S_p)
				y_trans2 = (y_aoa2 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)
					
				if( (((z_aoa1-rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
					(((z_aoa2-rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (ABS(z_trans2/t_p + y_trans2/L_p) + ABS(z_trans2/t_p - y_trans2/L_p)) .LE. 1)) then								
						!$acc atomic update
						ETA(i,j,k) = ETA(i,j,k) + 1.d0
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

!	If there are traces of ETA, it means that the initial condition was not set correctly, or
!	the solid diplacement between timesteps is too large
!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS-2,kEndVOS+2
do j=jBgnVOS-2,jEndVOS+2
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
end do
end do
end do   
!$acc end parallel
!$OMP END PARALLEL DO

	!-------------- Calculate moment of inertia -------------------

		if( istep == 0 )then
		
			I_moment = 0.d0
		
			!$OMP PARALLEL DO PRIVATE(k) REDUCTION(+: I_moment) collapse(3)
			!$acc parallel loop reduction(+: I_moment) collapse(3) gang vector
			do k = kBgnVOS-2,kEndVOS+2
				do j = jBgnVOS-2,jEndVOS+2
					do i=1,nx
						I_moment = I_moment+den_sol*ETA(i,j,k)*((ys(j) - y0)**2 +(zs(k)-z0)**2)*iDx(i)*iDy(j)*iDz(k)
					end do
				end do
			end do
			!$acc end parallel
			!$OMP END PARALLEL DO 
		
		end if 

		!write(*,*) 'Moment of inertia = ', I_moment

!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine func_Cylplate_2blade


subroutine func_Cylplate_1blade()

use variables
implicit none

integer:: az_min,az_max,ay_min,ay_max
integer:: m, n
real*8 :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

real*8 :: L_p, t_p, S_p, g_p
real*8 :: z_trans1, y_trans1
real*8 :: z_aoa1, y_aoa1

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

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1) &
!$acc      reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


do k=1,nz
	do j=1,ny

		Zf = Zs(k) - z0
		Yf = Ys(j) - y0
	
	! 1st blade initial position for AOA=0 is at Z+
	! rotation wrt. AOA for blade 1
	z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
	y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans1 = (z_aoa1 - rotor_r)           + r*SIN(S_p)
	y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)
	

	if( ((z_aoa1-rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. &
		(ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1

! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1) collapse(2) gang

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then
			
			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    ETA(i,j,k)=0.d0
            !$acc loop private(n,Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

					Zf = SZ(n) - z0
					Yf = SY(m) - y0

				! rotation wrt. AOA for blade 1			
				z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
				y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans1 = (z_aoa1 - rotor_r)           + r*SIN(S_p)
				y_trans1 = (y_aoa1 + (g_p+r)+0.5d0*L_p) - r*(1.d0 - COS(S_p)) + 0.5d0*t_p*SIN(S_p)

				if( ((z_aoa1-rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. &
					(ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) then
						!$acc atomic update
						ETA(i,j,k) = ETA(i,j,k) + 1.d0
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!	If there are traces of ETA, it means that the initial condition was not set correctly, or
!	the solid diplacement between timesteps is too large
!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS-2,kEndVOS+2
do j=jBgnVOS-2,jEndVOS+2
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
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
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine func_Cylplate_1blade


subroutine func_Cylplate_noblade()

use variables
implicit none

integer:: az_min,az_max,ay_min,ay_max
integer:: m, n
real*8 :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

real*8 :: z_aoa1, y_aoa1

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!



! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_aoa1,y_aoa1) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(Yf,Zf,z_aoa1,y_aoa1) &
!$acc      reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


do k=1,nz
	do j=1,ny

		Zf = Zs(k) - z0
		Yf = Ys(j) - y0
	
	! 1st blade initial position for AOA=0 is at Z+
	! rotation wrt. AOA for blade 1
	z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
	y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)


	if( ((z_aoa1-rotor_r)**2 + y_aoa1**2) .LE. r**2) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1

! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1) collapse(2) gang

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then
			
			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo


		    ETA(i,j,k) = 0.d0
            !$acc loop private(n,z_aoa1,y_aoa1) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

					Zf = SZ(n) - z0
					Yf = SY(m) - y0

				! rotation wrt. AOA for blade 1			
				z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
				y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)


				if( ((z_aoa1-rotor_r)**2 + y_aoa1**2) .LE. r**2) then
					!$acc atomic update
					ETA(i,j,k) = ETA(i,j,k) + 1.d0
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!	If there are traces of ETA, it means that the initial condition was not set correctly, or
!	the solid diplacement between timesteps is too large
!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS-2,kEndVOS+2
do j=jBgnVOS-2,jEndVOS+2
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
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
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine func_Cylplate_noblade


subroutine func_Cylinder()

use variables
implicit none

integer:: az_min,az_max,ay_min,ay_max
integer:: m, n
real*8 :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j,Yf,Zf) REDUCTION(max : az_max,ay_max) &
!$OMP			  REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(Yf,Zf) reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector

do k=1,nz
	do j=1,ny

		Zf = Zs(k) - z0_t
		Yf = Ys(j) - y0_t

	if( (Zf**2 + Yf**2) .LE. r**2) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1

! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,SY,SZ,Yf,Zf) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,SY,SZ,Yf,Zf) collapse(2) gang
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then
            
			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    ETA(i,j,k) = 0.d0
            !$acc loop private(n,Yf,Zf) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

					Zf = SZ(n) - z0_t
					Yf = SY(m) - y0_t
					
				if( (Zf**2 + Yf**2) .LE. r**2) then
					!$acc atomic update
					ETA(i,j,k) = ETA(i,j,k) + 1.d0
				endif
			end do; end do
		   
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!	If there are traces of ETA, it means that the initial condition was not set correctly, or
!	the solid diplacement between timesteps is too large
!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS-2,kEndVOS+2
do j=jBgnVOS-2,jEndVOS+2
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
end do
end do
end do   
!$acc end parallel
!$OMP END PARALLEL DO  


	!-------------- Calculate moment of inertia -------------------

		if( istep == 0 )then
		
			I_moment = 0.d0
		
			!$OMP PARALLEL DO PRIVATE(k) REDUCTION(+: I_moment) collapse(3)
			!$acc parallel loop reduction(+: I_moment) collapse(3) gang vector
			do k = kBgnVOS-2,kEndVOS+2
				do j = jBgnVOS-2,jEndVOS+2
					do i=1,nx
						I_moment = I_moment+den_sol*ETA(i,j,k)*((ys(j) - y0_t)**2 +(zs(k)-z0_t)**2)*iDx(i)*iDy(j)*iDz(k)
					end do
				end do
			end do
			!$acc end parallel
			!$OMP END PARALLEL DO 
		
		end if 

		!write(*,*) 'Moment of inertia = ', I_moment

!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,x,y,z,u_solid,v_solid,w_solid'
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,6(3X,F12.7))') time, x0_t, y0_t, z0_t, u_solid, v_solid, w_solid
	endif
    close(61)

end subroutine func_Cylinder



subroutine airfoil_00xx()

use variables
implicit none

integer:: az_min,az_max,ay_min,ay_max
integer:: m, n
real*8 :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

real*8 :: c_p, t_p, s_p
real*8 :: z_trans1, y_trans1

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

c_p = L_ch  ! blade's chord
t_p = 0.12d0 ! blade's thickness (NACA 00XX --> t_p = XX/100
s_p = 0.d0  ! CM offset from leading edge


! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_trans1) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(Yf,Zf,z_trans1) &
!$acc 	   reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


do k=1,nz
	do j=1,ny

		Zf = Zs(k) - z0
		Yf = Ys(j) - y0

	! rotation wrt. AOA for blade 1
	z_trans1 =  Zf*COS(AOA*PI/180.d0)+Yf*SIN(AOA*PI/180.d0) + s_p
	y_trans1 = -Zf*SIN(AOA*PI/180.d0)+Yf*COS(AOA*PI/180.d0)

	
	if( ((y_trans1 - 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .LE. 0.) .AND. &
	     (y_trans1 + 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .GE. 0.))) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1

! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_trans1,y_trans1) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,Yf,Zf,SY,SZ,z_trans1,y_trans1) collapse(2) gang

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then
			
			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    ETA(i,j,k) = 0.d0
            !$acc loop private(n,Yf,Zf,z_trans1) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

					Zf = SZ(n) - z0
					Yf = SY(m) - y0

				! rotation wrt. AOA for blade 1
				z_trans1 =  Zf*COS(AOA*PI/180.d0)+Yf*SIN(AOA*PI/180.d0) + s_p
				y_trans1 = -Zf*SIN(AOA*PI/180.d0)+Yf*COS(AOA*PI/180.d0)
								
					
		if( ((y_trans1 - 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .LE. 0.) .AND. &
			 (y_trans1 + 5.*t_p*c_p*(0.2969*sqrt(z_trans1/c_p) - 0.126*z_trans1/c_p - 0.3516*(z_trans1/c_p)**2 + 0.2843*(z_trans1/c_p)**3 - 0.1015*(z_trans1/c_p)**4) .GE. 0.))) then							
					!$acc atomic update
					ETA(i,j,k) = ETA(i,j,k) + 1.d0
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!	If there are traces of ETA, it means that the initial condition was not set correctly, or
!	the solid diplacement between timesteps is too large
!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS-2,kEndVOS+2
do j=jBgnVOS-2,jEndVOS+2
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
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
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine airfoil_00xx


subroutine NLR7301()

use variables
implicit none

integer:: az_min,az_max,ay_min,ay_max
integer:: m, n
real*8 :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

real*8 :: z_aoa1, y_aoa1

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!


! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_aoa1,y_aoa1) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(Yf,Zf,z_aoa1,y_aoa1) &
!$acc 	   reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


do k=1,nz
	do j=1,ny

		Zf = Zs(k) - z0_t
		Yf = Ys(j) - y0_t

	z_aoa1 =  Zf * COS((AOA)*PI/180.d0) + Yf * SIN((AOA)*PI/180.d0) + offset_z * (1 - COS((AOA)*PI/180.d0))
	y_aoa1 = -Zf * SIN((AOA)*PI/180.d0) + Yf * COS((AOA)*PI/180.d0) - offset_z * SIN((AOA)*PI/180.d0)

	if(((-16040082377189.02539062d0*z_aoa1**8 + 2076613785850.65429688d0*z_aoa1**7 - 111151138382.01135254d0*z_aoa1**6 + 3181591509.27247143d0*z_aoa1**5 - 52688585.24160857d0*z_aoa1**4 + 512560.45074645d0*z_aoa1**3 - 2914.70429386d0*z_aoa1**2 + 11.29370989d0*z_aoa1 + 0.00145690d0 .GE. y_aoa1) .AND. &
		(14317365091714.06640625d0*z_aoa1**8 - 1874506187987.77319336d0*z_aoa1**7 + 101472830998.53988647d0*z_aoa1**6 - 2935929084.82835007d0*z_aoa1**5 + 49049520.99299136d0*z_aoa1**4 - 478830.97776009d0*z_aoa1**3 + 2686.68293892d0*z_aoa1**2 - 9.59170081d0*z_aoa1 - 0.00327512d0 .LE. y_aoa1) .AND. &
		(z_aoa1 .GE. 0.d0) .AND. (z_aoa1 .LE. 0.0305844d0)) .OR. &
	   ((-13.91507664d0*z_aoa1**8 + 58.76283187d0*z_aoa1**7 - 102.84730271d0*z_aoa1**6 + 97.27949815d0*z_aoa1**5 - 54.52243098d0*z_aoa1**4 + 18.65274402d0*z_aoa1**3 - 4.03367889d0*z_aoa1**2 + 0.59917692d0*z_aoa1 + 0.03937969d0 .GE. y_aoa1) .AND. &
		(11.85968221d0*z_aoa1**8 - 57.91030005d0*z_aoa1**7 + 112.51994610d0*z_aoa1**6 - 114.85681566d0*z_aoa1**5 + 67.58720217d0*z_aoa1**4 - 23.39537083d0*z_aoa1**3 + 4.91637212d0*z_aoa1**2 - 0.68204380d0*z_aoa1 - 0.02573834d0 .LE. y_aoa1) .AND. &
		(z_aoa1 .GE. 0.0305844d0) .AND. (z_aoa1 .LE. 1.d0))	.OR. &
	   ((1122225947.02511907d0*z_aoa1**8 - 2927239004.64574718d0*z_aoa1**7 + 1018097536.27265966d0*z_aoa1**6 + 2139895321.05184412d0*z_aoa1**5 - 7054343.26825002d0*z_aoa1**4 - 1922291256.07341337d0*z_aoa1**3 - 801320538.51426792d0*z_aoa1**2 + 2102797818.81258583d0*z_aoa1 - 725108770.23475766d0 .GE. y_aoa1) .AND. &
		(-578107206.20548809d0*z_aoa1**8 + 1508464378.49922228d0*z_aoa1**7 - 524344172.27626121d0*z_aoa1**6 - 1105147540.37096477d0*z_aoa1**5 + 6155172.39750056d0*z_aoa1**4 + 989134008.15751648d0*z_aoa1**3 + 416413656.77780855d0*z_aoa1**2 - 1087571331.23178387d0*z_aoa1 + 375001688.30047375d0 .LE. y_aoa1) .AND. &
		(z_aoa1 .GE. 0.943491d0) .AND. (z_aoa1 .LE. 0.950344d0)) .OR. &
	   ((-21872.36792795d0*z_aoa1**8 + 195265.86529625d0*z_aoa1**7 - 761872.49234961d0*z_aoa1**6 + 1696856.18512821d0*z_aoa1**5 - 2359563.53647953d0*z_aoa1**4 + 2097698.47387663d0*z_aoa1**3 - 1164347.44946958d0*z_aoa1**2 + 368925.52822860d0*z_aoa1 - 51090.21897113d0 .GE. y_aoa1) .AND. &
		(6197.08005899d0*z_aoa1**8 - 54151.40746550d0*z_aoa1**7 + 206734.78917549d0*z_aoa1**6 - 450381.06427308d0*z_aoa1**5 + 612385.90625673d0*z_aoa1**4 - 532166.55309417d0*z_aoa1**3 + 288634.23821150d0*z_aoa1**2 - 89332.77797694d0*z_aoa1 + 12079.74827825d0 .LE. y_aoa1) .AND. &
		(z_aoa1 .GE. 0.950344d0) .AND. (z_aoa1 .LE. 1.27339d0))	&
		 ) then	
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1

! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1) collapse(2) gang

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then
			
			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    ETA(i,j,k) = 0.d0
            !$acc loop private(n,Yf,Zf,z_aoa1,y_aoa1) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

					Zf = SZ(n) - z0_t
					Yf = SY(m) - y0_t

				z_aoa1 =  Zf * COS((AOA)*PI/180.d0) + Yf * SIN((AOA)*PI/180.d0) + offset_z * (1 - COS((AOA)*PI/180.d0))
				y_aoa1 = -Zf * SIN((AOA)*PI/180.d0) + Yf * COS((AOA)*PI/180.d0) - offset_z * SIN((AOA)*PI/180.d0)			
					
	if(((-16040082377189.02539062d0*z_aoa1**8 + 2076613785850.65429688d0*z_aoa1**7 - 111151138382.01135254d0*z_aoa1**6 + 3181591509.27247143d0*z_aoa1**5 - 52688585.24160857d0*z_aoa1**4 + 512560.45074645d0*z_aoa1**3 - 2914.70429386d0*z_aoa1**2 + 11.29370989d0*z_aoa1 + 0.00145690d0 .GE. y_aoa1) .AND. &
		(14317365091714.06640625d0*z_aoa1**8 - 1874506187987.77319336d0*z_aoa1**7 + 101472830998.53988647d0*z_aoa1**6 - 2935929084.82835007d0*z_aoa1**5 + 49049520.99299136d0*z_aoa1**4 - 478830.97776009d0*z_aoa1**3 + 2686.68293892d0*z_aoa1**2 - 9.59170081d0*z_aoa1 - 0.00327512d0 .LE. y_aoa1) .AND. &
		(z_aoa1 .GE. 0.d0) .AND. (z_aoa1 .LE. 0.0305844d0)) .OR. &
	   ((-13.91507664d0*z_aoa1**8 + 58.76283187d0*z_aoa1**7 - 102.84730271d0*z_aoa1**6 + 97.27949815d0*z_aoa1**5 - 54.52243098d0*z_aoa1**4 + 18.65274402d0*z_aoa1**3 - 4.03367889d0*z_aoa1**2 + 0.59917692d0*z_aoa1 + 0.03937969d0 .GE. y_aoa1) .AND. &
		(11.85968221d0*z_aoa1**8 - 57.91030005d0*z_aoa1**7 + 112.51994610d0*z_aoa1**6 - 114.85681566d0*z_aoa1**5 + 67.58720217d0*z_aoa1**4 - 23.39537083d0*z_aoa1**3 + 4.91637212d0*z_aoa1**2 - 0.68204380d0*z_aoa1 - 0.02573834d0 .LE. y_aoa1) .AND. &
		(z_aoa1 .GE. 0.0305844d0) .AND. (z_aoa1 .LE. 1.d0))	.OR. &
	   ((1122225947.02511907d0*z_aoa1**8 - 2927239004.64574718d0*z_aoa1**7 + 1018097536.27265966d0*z_aoa1**6 + 2139895321.05184412d0*z_aoa1**5 - 7054343.26825002d0*z_aoa1**4 - 1922291256.07341337d0*z_aoa1**3 - 801320538.51426792d0*z_aoa1**2 + 2102797818.81258583d0*z_aoa1 - 725108770.23475766d0 .GE. y_aoa1) .AND. &
		(-578107206.20548809d0*z_aoa1**8 + 1508464378.49922228d0*z_aoa1**7 - 524344172.27626121d0*z_aoa1**6 - 1105147540.37096477d0*z_aoa1**5 + 6155172.39750056d0*z_aoa1**4 + 989134008.15751648d0*z_aoa1**3 + 416413656.77780855d0*z_aoa1**2 - 1087571331.23178387d0*z_aoa1 + 375001688.30047375d0 .LE. y_aoa1) .AND. &
		(z_aoa1 .GE. 0.943491d0) .AND. (z_aoa1 .LE. 0.950344d0)) .OR. &
	   ((-21872.36792795d0*z_aoa1**8 + 195265.86529625d0*z_aoa1**7 - 761872.49234961d0*z_aoa1**6 + 1696856.18512821d0*z_aoa1**5 - 2359563.53647953d0*z_aoa1**4 + 2097698.47387663d0*z_aoa1**3 - 1164347.44946958d0*z_aoa1**2 + 368925.52822860d0*z_aoa1 - 51090.21897113d0 .GE. y_aoa1) .AND. &
		(6197.08005899d0*z_aoa1**8 - 54151.40746550d0*z_aoa1**7 + 206734.78917549d0*z_aoa1**6 - 450381.06427308d0*z_aoa1**5 + 612385.90625673d0*z_aoa1**4 - 532166.55309417d0*z_aoa1**3 + 288634.23821150d0*z_aoa1**2 - 89332.77797694d0*z_aoa1 + 12079.74827825d0 .LE. y_aoa1) .AND. &
		(z_aoa1 .GE. 0.950344d0) .AND. (z_aoa1 .LE. 1.27339d0))	&
		 ) then					
				!$acc atomic update
				ETA(i,j,k) = ETA(i,j,k) + 1.d0
	endif
			end do; end do
		   
		  
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!	If there are traces of ETA, it means that the initial condition was not set correctly, or
!	the solid diplacement between timesteps is too large
!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS-2,kEndVOS+2
do j=jBgnVOS-2,jEndVOS+2
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
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
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine NLR7301


! subroutine func_Cylplate_2blade()

! use variables
! implicit none

! integer:: az_min,az_max,ay_min,ay_max
! integer:: m, n
! real*8 :: yf, zf
! real*8 ,dimension(1:nSubGrids_f) :: SY
! real*8 ,dimension(1:nSubGrids_f) :: SZ
! real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
! real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

! real*8 :: L_p, t_p, B_p, g_p
! real*8 :: z_trans_, y_trans_, z_trans1, y_trans1, z_trans2, y_trans2
! real*8 :: z_aoa1, y_aoa1, z_aoa2, y_aoa2

   ! !---------------------------------------------------!
   ! !    LOCAL VARIABLES                                !
   ! !---------------------------------------------------!

! L_p = 1.d0 ! Plate's length
! t_p = 0.2d0 ! Plate's thickness
! B_p = (0.d0/180.d0)*PI ! Plate's Beta angle
! g_p = 0.1d0 ! Plate's gap


! ! Create coarse ETA (only 0 or 1)
! i=0
! az_min = nz ; az_max = 0
! ay_min = ny ; ay_max = 0

! !$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

! !$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) &
! !$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
! !$acc parallel vector_length(32)
! !$acc loop private(Yf,Zf,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) &
! !$acc 	   reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


! do k=1,nz
	! do j=1,ny

		! Zf = Zs(k) - z0
		! Yf = Ys(j) - y0
	
	! ! 1st blade initial position for AOA=0 is at Z+
	! ! rotation wrt. AOA for blade 1
	! z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
	! y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

	! ! z and y transformation for flat plate
	! z_trans_ = z_aoa1 + rotor_r
	! y_trans_ = y_aoa1 + (g_p+r)+0.5d0*L_p
	
	! z_trans1 = z_trans_*COS(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*SIN(B_p)
	! y_trans1 = z_trans_*SIN(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*COS(B_p) + (g_p+r)+0.5d0*L_p


	! ! rotation wrt. AOA for blade 2
	! z_aoa2 =  Zf*COS((AOA+180.d0)*PI/180.d0)+Yf*SIN((AOA+180.d0)*PI/180.d0)
	! y_aoa2 = -Zf*SIN((AOA+180.d0)*PI/180.d0)+Yf*COS((AOA+180.d0)*PI/180.d0)

	! ! z and y transformation for flat plate
	! z_trans_ = z_aoa2 + rotor_r
	! y_trans_ = y_aoa2 + (g_p+r)+0.5d0*L_p
	
	! z_trans2 = z_trans_*COS(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*SIN(B_p)
	! y_trans2 = z_trans_*SIN(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*COS(B_p) + (g_p+r)+0.5d0*L_p
	

	! if( (((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
	    ! (((z_aoa2+rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (ABS(z_trans2/t_p + y_trans2/L_p) + ABS(z_trans2/t_p - y_trans2/L_p)) .LE. 1)) then
		! ETA(i,j,k)=1.D0
		! az_min = MIN(k,az_min)
		! az_max = MAX(k,az_max)
		! ay_min = MIN(j,ay_min)
		! ay_max = MAX(j,ay_max)
	! else
		! ETA(i,j,k)=0.D0		
	! endif
	
	! end do
! end do
! !$acc end parallel
! !$OMP END PARALLEL DO

! ! Setting the boundary for subgrids
	! iBgnVOS = 0
	! iEndVOS = nx+1
	! jBgnVOS = ay_min - 1
	! jEndVOS = ay_max + 1
	! kBgnVOS = az_min - 1
	! kEndVOS = az_max + 1

! ! Finding potential grids for subgrids

! i=iBgnVOS
! !$OMP PARALLEL DO COLLAPSE(2)
! !$acc parallel loop collapse(2) gang vector
! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
	! if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		! ETA_2D(k,j) = 1
	! else
		! ETA_2D(k,j) = 0
	! endif
! enddo;enddo
! !$acc end parallel
! !$OMP END PARALLEL DO


! ! Creating subgrids

! !$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) collapse(2)
! !$acc parallel vector_length(32)
! !$acc loop private(m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) collapse(2) gang

! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS

	! if(ETA_2D(k,j) == 1) then
			
			! !$acc loop vector
            ! do n=1,nSubGrids_f
				! SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				! SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			! enddo

		    ! ETA(i,j,k) = 0.d0
            ! !$acc loop private(n,Yf,Zf,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1,z_aoa2,y_aoa2,z_trans2,y_trans2) collapse(2) vector
            ! do m=1,nSubGrids_f
            ! do n=1,nSubGrids_f

					! Zf = SZ(n) - z0
					! Yf = SY(m) - y0

				! ! 1st blade initial position for AOA=0 is at Z+
				! ! rotation wrt. AOA for blade 1
				! z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
				! y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

				! ! z and y transformation for flat plate
				! z_trans_ = z_aoa1 + rotor_r
				! y_trans_ = y_aoa1 + (g_p+r)+0.5d0*L_p
	
				! z_trans1 = z_trans_*COS(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*SIN(B_p)
				! y_trans1 = z_trans_*SIN(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*COS(B_p) + (g_p+r)+0.5d0*L_p


				! ! rotation wrt. AOA for blade 2
				! z_aoa2 =  Zf*COS((AOA+180.d0)*PI/180.d0)+Yf*SIN((AOA+180.d0)*PI/180.d0)
				! y_aoa2 = -Zf*SIN((AOA+180.d0)*PI/180.d0)+Yf*COS((AOA+180.d0)*PI/180.d0)

				! ! z and y transformation for flat plate
				! z_trans_ = z_aoa2 + rotor_r
				! y_trans_ = y_aoa2 + (g_p+r)+0.5d0*L_p
	
				! z_trans2 = z_trans_*COS(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*SIN(B_p)
				! y_trans2 = z_trans_*SIN(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*COS(B_p) + (g_p+r)+0.5d0*L_p
					
				! if( (((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2 .OR. (ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) .OR. &
					! (((z_aoa2+rotor_r)**2 + y_aoa2**2) .LE. r**2 .OR. (ABS(z_trans2/t_p + y_trans2/L_p) + ABS(z_trans2/t_p - y_trans2/L_p)) .LE. 1)) then								
						! !$acc atomic update
						! ETA(i,j,k) = ETA(i,j,k) + 1.d0
				! endif
			! end do; end do
		   
		  
		  ! ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	! endif
! end do
! end do		
! !$acc end parallel
! !$OMP END PARALLEL DO

! !	If there are traces of ETA, it means that the initial condition was not set correctly, or
! !	the solid diplacement between timesteps is too large
! !$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
! !$acc parallel loop collapse(3) gang vector
! do k=kBgnVOS-2,kEndVOS+2
! do j=jBgnVOS-2,jEndVOS+2
! do i=iBgnVOS+1,iEndVOS
	! ETA(i,j,k) = ETA(iBgnVOS,j,k)
! end do
! end do
! end do   
! !$acc end parallel
! !$OMP END PARALLEL DO

! !$acc end data   

    ! if(myid==master .AND. istep == 0)then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,*)'                 '
        ! write(61,*)'SOLID MOTION'
		! write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        ! write(61,*)'                 '
    ! else if (myid==master) then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	! endif
    ! close(61)

! end subroutine func_Cylplate_2blade


! subroutine func_Cylplate_1blade()

! use variables
! implicit none

! integer:: az_min,az_max,ay_min,ay_max
! integer:: m, n
! real*8 :: yf, zf
! real*8 ,dimension(1:nSubGrids_f) :: SY
! real*8 ,dimension(1:nSubGrids_f) :: SZ
! real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
! real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

! real*8 :: L_p, t_p, B_p, g_p
! real*8 :: z_trans_, y_trans_, z_trans1, y_trans1
! real*8 :: z_aoa1, y_aoa1

   ! !---------------------------------------------------!
   ! !    LOCAL VARIABLES                                !
   ! !---------------------------------------------------!

! L_p = 1.d0 ! Plate's length
! t_p = 0.2d0 ! Plate's thickness
! B_p = (0.d0/180.d0)*PI ! Plate's Beta angle
! g_p = 0.1d0 ! Plate's gap


! ! Create coarse ETA (only 0 or 1)
! i=0
! az_min = nz ; az_max = 0
! ay_min = ny ; ay_max = 0

! !$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

! !$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1) &
! !$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
! !$acc parallel vector_length(32)
! !$acc loop private(Yf,Zf,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1) &
! !$acc 	   reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


! do k=1,nz
	! do j=1,ny

		! Zf = Zs(k) - z0
		! Yf = Ys(j) - y0
	
	! ! 1st blade initial position for AOA=0 is at Z+
	! ! rotation wrt. AOA for blade 1
	! z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
	! y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

	! ! z and y transformation for flat plate
	! z_trans_ = z_aoa1 + rotor_r
	! y_trans_ = y_aoa1 + (g_p+r)+0.5d0*L_p
	
	! z_trans1 = z_trans_*COS(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*SIN(B_p)
	! y_trans1 = z_trans_*SIN(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*COS(B_p) + (g_p+r)+0.5d0*L_p


	! if( (((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2) .OR. ((ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1)) then
		! ETA(i,j,k)=1.D0
		! az_min = MIN(k,az_min)
		! az_max = MAX(k,az_max)
		! ay_min = MIN(j,ay_min)
		! ay_max = MAX(j,ay_max)
	! else
		! ETA(i,j,k)=0.D0		
	! endif
	
	! end do
! end do
! !$acc end parallel
! !$OMP END PARALLEL DO

! ! Setting the boundary for subgrids
	! iBgnVOS = 0
	! iEndVOS = nx+1
	! jBgnVOS = ay_min - 1
	! jEndVOS = ay_max + 1
	! kBgnVOS = az_min - 1
	! kEndVOS = az_max + 1

! ! Finding potential grids for subgrids

! i=iBgnVOS
! !$OMP PARALLEL DO COLLAPSE(2)
! !$acc parallel loop collapse(2) gang vector
! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
	! if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		! ETA_2D(k,j) = 1
	! else
		! ETA_2D(k,j) = 0
	! endif
! enddo;enddo
! !$acc end parallel
! !$OMP END PARALLEL DO


! ! Creating subgrids

! !$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1) collapse(2)
! !$acc parallel vector_length(32)
! !$acc loop private(m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1) collapse(2) gang

! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS

	! if(ETA_2D(k,j) == 1) then
			
			! !$acc loop vector
            ! do n=1,nSubGrids_f
				! SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				! SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			! enddo

		    ! ETA(i,j,k) = 0.d0
            ! !$acc loop private(n,Yf,Zf,z_aoa1,y_aoa1,z_trans_,y_trans_,z_trans1,y_trans1) collapse(2) vector
            ! do m=1,nSubGrids_f
            ! do n=1,nSubGrids_f

					! Zf = SZ(n) - z0
					! Yf = SY(m) - y0

				! ! 1st blade initial position for AOA=0 is at Z+
				! ! rotation wrt. AOA for blade 1
				! z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
				! y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

				! ! z and y transformation for flat plate
				! z_trans_ = z_aoa1 + rotor_r
				! y_trans_ = y_aoa1 + (g_p+r)+0.5d0*L_p
	
				! z_trans1 = z_trans_*COS(B_p) - (y_trans_-(g_p+r)-0.5d0*L_p)*SIN(B_p)
				! y_trans1 = z_trans_*SIN(B_p) + (y_trans_-(g_p+r)-0.5d0*L_p)*COS(B_p) + (g_p+r)+0.5d0*L_p

					
				! if( (((z_aoa1+rotor_r)**2 + y_aoa1**2) .LE. r**2) .OR. ((ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1)) then							
						! !$acc atomic update
						! ETA(i,j,k) = ETA(i,j,k) + 1.d0
				! endif
			! end do; end do
		   
		  
		  ! ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	! endif
! end do
! end do		
! !$acc end parallel
! !$OMP END PARALLEL DO

! !	If there are traces of ETA, it means that the initial condition was not set correctly, or
! !	the solid diplacement between timesteps is too large
! !$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
! !$acc parallel loop collapse(3) gang vector
! do k=kBgnVOS-2,kEndVOS+2
! do j=jBgnVOS-2,jEndVOS+2
! do i=iBgnVOS+1,iEndVOS
	! ETA(i,j,k) = ETA(iBgnVOS,j,k)
! end do
! end do
! end do   
! !$acc end parallel
! !$OMP END PARALLEL DO

! !$acc end data   

    ! if(myid==master .AND. istep == 0)then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,*)'                 '
        ! write(61,*)'SOLID MOTION'
		! write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        ! write(61,*)'                 '
    ! else if (myid==master) then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
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

	 ! DISTANCE = sqrt( (( (Z(k)+Z(k+1)) / 2.D0 ) -   z0 )  ** 2.D0 &
                 ! +  (( (Y(j)+Y(j+1)) / 2.D0 ) - y0 ) ** 2.D0 )  

	   ! if( abs(DISTANCE - r) < DIAGONAL ) then
		  
		  ! xi = 0.0
		  
		 
			! do m=1,nSubGrids_f+1; do n=1,nSubGrids_f+1
				! dzg = iDz(k) / nSubGrids_f
				! dyg = iDy(j) / nSubGrids_f
				! SZ(n) = Z(k) + (n-1) * dzg
				! SY(m) = Y(j) + (m-1) * dyg		
			! end do; end do
		  
		  

		  
			! do m=1,nSubGrids_f; do n=1,nSubGrids_f
				  ! SDIST = sqrt( (( (SZ(n)+SZ(n+1)) / 2.D0 ) - z0 )  ** 2.D0 &
                    ! +   (( (SY(m)+SY(m+1)) / 2.D0 ) - y0 ) ** 2.D0 )
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


subroutine square_cyl()

use variables
implicit none

integer:: az_min,az_max,ay_min,ay_max
integer:: m, n
real*8 :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

real*8 :: L_p, t_p, S_p, g_p
real*8 :: z_trans1, y_trans1
real*8 :: z_aoa1, y_aoa1

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

L_p = 1.d0 ! Plate's length
t_p = 1.d0 ! Plate's thickness


! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j,Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1) &
!$acc      reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector


do k=1,nz
	do j=1,ny

		Zf = Zs(k) - z0
		Yf = Ys(j) - y0

	! rotation wrt. AOA for blade 1
	z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
	y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans1 = z_aoa1
	y_trans1 = y_aoa1
	

	if( (ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1

! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,Yf,Zf,SY,SZ,z_aoa1,y_aoa1,z_trans1,y_trans1) collapse(2) gang

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then
			
			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    ETA(i,j,k) = 0.d0
            !$acc loop private(n,Yf,Zf,z_aoa1,y_aoa1,z_trans1,y_trans1) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

					Zf = SZ(n) - z0
					Yf = SY(m) - y0

				! rotation wrt. AOA for blade 1
				z_aoa1 =  Zf*COS((AOA)*PI/180.d0)+Yf*SIN((AOA)*PI/180.d0)
				y_aoa1 = -Zf*SIN((AOA)*PI/180.d0)+Yf*COS((AOA)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans1 = z_aoa1
				y_trans1 = y_aoa1

				if( (ABS(z_trans1/t_p + y_trans1/L_p) + ABS(z_trans1/t_p - y_trans1/L_p)) .LE. 1) then
					!$acc atomic update
					ETA(i,j,k) = ETA(i,j,k) + 1.d0
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS-2,kEndVOS+2
do j=jBgnVOS-2,jEndVOS+2
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
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
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine square_cyl


subroutine flat_plate()

use variables
implicit none

integer	:: az_min,az_max,ay_min,ay_max
integer	:: m ,n
real*8  :: yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid2 = 1.d0/ (nSubGrids_f**2*1.d0)

real*8 	:: z_trans1, y_trans1
real*8 	:: z_aoa1, y_aoa1

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!


! Create coarse ETA (only 0 or 1)
i=0
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0

!$acc data present(Y,Z,Ys,Zs,iDy,iDz,ETA) create(SY,SZ,ETA_2D)

!$OMP PARALLEL DO PRIVATE(j) &
!$OMP			  REDUCTION(max : az_max,ay_max) REDUCTION(min : az_min,ay_min) collapse(2)
!$acc parallel vector_length(32)
!$acc loop reduction(max : az_max,ay_max) reduction(min : az_min,ay_min) collapse(2) gang vector

do k=1,nz
	do j=1,ny
	
	if( ((Ys(j) - y0) .LE. 0.d0) .OR. ((Ys(j) .LE. (y0+0.005)) .AND. (Zs(k) .GE. 0.5d0) .AND. (Zs(k) .LE. 0.508d0)) ) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
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
	jBgnVOS = 1 !ay_min - 5
	jEndVOS = ay_max + 10
	kBgnVOS = 1 !az_min - 5
	kEndVOS = nz !az_max + 5

! Finding potential grids for subgrids

i=iBgnVOS
!$OMP PARALLEL DO COLLAPSE(2)
!$acc parallel loop collapse(2) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i, j-1:j+1, k-1:k+1))) then
		ETA_2D(k,j) = 1
	else
		ETA_2D(k,j) = 0
	endif
enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,m,n,SY,SZ) collapse(2)
!$acc parallel vector_length(32)
!$acc loop private(m,n,SY,SZ) collapse(2) gang

do k=kBgnVOS,kEndVOS
	do j=jBgnVOS,jEndVOS

	if(ETA_2D(k,j) == 1) then
            
			!$acc loop vector
            do n=1,nSubGrids_f
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    ETA(i,j,k) = 0.d0
            !$acc loop private(n) collapse(2) vector
            do m=1,nSubGrids_f
            do n=1,nSubGrids_f

				if( ((SY(m) - y0) .LE. 0.d0) .OR. ((SY(m) .LE. (y0+0.005)) .AND. (SZ(n) .GE. 0.5d0) .AND. (SZ(n) .LE. 0.508d0))) then
					!$acc atomic update
					ETA(i,j,k) = ETA(i,j,k) + 1.d0
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = ETA(i,j,k)*inv_subgrid2
	endif
end do
end do		
!$acc end parallel
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(j,i) collapse(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)
end do
end do
end do   
!$acc end parallel
!$OMP END PARALLEL DO    


do j=jBgnVOS,jEndVOS
	write(*,*) Y(j), ETA(iBgnVOS,j,159)
end do


!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,AOA,rotor_omega'  ! for rotation
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,3X,F14.6,3X,F12.7)') time, AOA, rotor_omega  ! for rotation
	endif
    close(61)

end subroutine flat_plate


!********************************** Pure 3D solids ***********************************


! subroutine func_Sphere()
! use variables
! implicit none
! real*4 :: DISTANCE,SDIST
! real*4 :: DIAGONAL
! integer :: l ,m ,n
! real*4 :: xi
! real*4 :: dxg, dyg, dzg
! real*4 ,dimension(1:nSubGrids_f+1) :: SX
! real*4 ,dimension(1:nSubGrids_f+1) :: SY
! real*4 ,dimension(1:nSubGrids_f+1) :: SZ


   ! !---------------------------------------------------!
   ! !    LOCAL VARIABLES                                !
   ! !---------------------------------------------------!

	! iBgnVOS = 1
	! iEndVOS = nx
	! jBgnVOS = 1
	! jEndVOS = ny
	! kBgnVOS = 1
	! kEndVOS = nz

! !$acc data present(X,Y,Z,iDx,iDy,iDz,ETA) create(SX,SY,SZ)

! !$OMP PARALLEL DO PRIVATE(i,j,DIAGONAL,DISTANCE,SDIST,l,m,n,dxg,dyg,dzg,SX,SY,SZ,xi) collapse(3)
! !$acc parallel loop independent private(DIAGONAL,DISTANCE,SDIST,l,m,n,dxg,dyg,dzg,SX,SY,SZ,xi) collapse(3) gang vector
! do k=kBgnVOS,kEndVOS 
  ! do j=jBgnVOS,jEndVOS
    ! do i=iBgnVOS,iEndVOS

         ! DIAGONAL = sqrt( iDx(i)*iDx(i) + iDy(j)*iDy(j) + iDz(k)*iDz(k) ) / 2.D0
         ! DISTANCE = sqrt((((X(i)+X(i+1))/2.D0)-x0_t)**2.D0+ &
                         ! (((Y(j)+Y(j+1))/2.D0)-y0_t)**2.D0+ &
                         ! (((Z(k)+Z(k+1))/2.D0)-z0_t)**2.D0  )

           ! if( abs(DISTANCE - r) < DIAGONAL ) then
              ! !$acc loop seq
              ! do l=1,nSubGrids_f+1,1
			   ! !$acc loop seq
               ! do m=1,nSubGrids_f+1,1
			    ! !$acc loop seq
                ! do n=1,nSubGrids_f+1,1
					! dxg = iDx(i) / nSubGrids_f
					! dyg = iDy(j) / nSubGrids_f
					! dzg = iDz(k) / nSubGrids_f
					! SX(n)=X(i)+(n-1)*dxg
					! SY(m)=Y(j)+(m-1)*dyg
					! SZ(l)=Z(k)+(l-1)*dzg
                ! end do
               ! end do
              ! end do

              ! xi = 0.0
              ! !$acc loop seq
			  ! do l=1,nSubGrids_f
			    ! !$acc loop seq
                ! do m=1,nSubGrids_f
				  ! !$acc loop seq
                  ! do n=1,nSubGrids_f
                    
                    ! SDIST=SQRT((((SX(n)+SX(n+1))/2.D0)-x0_t)**2.D0+ &
                               ! (((SY(m)+SY(m+1))/2.D0)-y0_t)**2.D0+ &
                               ! (((SZ(l)+SZ(l+1))/2.D0)-z0_t)**2.D0  )
                    ! if(SDIST <= r) then
                      ! xi = xi + 1
                    ! end if

                  ! end do
                ! end do
              ! end do   
              ! ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f*nSubGrids_f)
              
              ! !print*,'ETA(i,j,k) = ',ETA(i,j,k),'xi = ',xi,i,j,k
           
           ! else if (DISTANCE <= r) then
             ! ETA(i,j,k)=1.D0
           ! else
             ! ETA(i,j,k)=0.D0
           ! end if
                
    ! end do
  ! end do
! end do
! !$acc end parallel
! !$OMP END PARALLEL DO

! !$acc end data

    ! if(myid==master .AND. istep == 0)then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,*)'                 '
        ! write(61,*)'SOLID MOTION'
		! write(61,*) ' VARIABLES = t*,x,y,z,u_solid,v_solid,w_solid'
        ! write(61,*)'                 '
    ! else if (myid==master) then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,'(F12.7,6(3X,F12.7))') time, x0_t, y0_t, z0_t, u_solid, v_solid, w_solid
	! endif
    ! close(61)

! end subroutine func_Sphere


subroutine func_Sphere()

use variables
implicit none

integer	:: az_min,az_max,ay_min,ay_max,ax_min,ax_max
integer :: l ,m ,n, xi
real*8 	:: xf, yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SX
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid3 = 1.d0/ (nSubGrids_f*nSubGrids_f**2*1.d0)


   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

! Create coarse ETA (only 0 or 1)
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0
ax_min = nx ; ax_max = 0

!$acc data present(X,Y,Z,Xs,Ys,Zs,iDx,iDy,iDz,ETA) create(SX,SY,SZ,ETA_3D)

!$OMP PARALLEL DO PRIVATE(j,i,Xf,Yf,Zf) REDUCTION(max : az_max,ay_max,ax_max) &
!$OMP			  REDUCTION(min : az_min,ay_min,ax_min) collapse(3)
!$acc parallel vector_length(32)
!$acc loop independent private(Xf,Yf,Zf) reduction(max : az_max,ay_max,ax_max) reduction(min : az_min,ay_min,ax_min) collapse(3) gang vector

do k=1,nz
	do j=1,ny
		do i=1,nx

		Xf = Xs(i) - x0_t
			!Xf = Xs(i) - x0		! moving X reference frame		
		Yf = Ys(j) - y0_t
		Zf = Zs(k) - z0_t
			!Zf = Zs(k) - z0		! moving Z reference frame

	if( Zf**2 + Yf**2 + Xf**2 .LE. r**2) then				! Geometry function
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
		ax_min = MIN(i,ax_min)
		ax_max = MAX(i,ax_max)
	else
		ETA(i,j,k)=0.D0		
	endif
		end do
	end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

! Setting the boundary for subgrids
	iBgnVOS = ax_min - 1
	iEndVOS = ax_max + 1
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1


! Finding potential grids for subgrids
!$OMP PARALLEL DO COLLAPSE(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS,iEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i-1:i+1, j-1:j+1, k-1:k+1))) then
		ETA_3D(k,j,i) = 1
	else
		ETA_3D(k,j,i) = 0
	endif
enddo;enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,i,l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3)
!$acc parallel vector_length(32)
!$acc loop independent private(l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3) gang 
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS,iEndVOS

	if(ETA_3D(k,j,i) == 1) then
            
            !$acc loop
            do n=1,nSubGrids_f,1
				SX(n) = X(i) + (n - 0.5d0)*iDx(i)*inv_subgrid1
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    xi = 0		
            !$acc loop private(m,n,Xf,Yf,Zf) reduction(+:xi)
            do l=1,nSubGrids_f,1	
            !$acc loop private(n,Xf,Yf,Zf) reduction(+:xi)
            do m=1,nSubGrids_f,1
            !$acc loop private(Xf,Yf,Zf) reduction(+:xi)
            do n=1,nSubGrids_f,1
			
				Xf = SX(l) - x0_t
					!Xf = SX(l) - x0		! moving X reference frame		
				Yf = SY(m) - y0_t
				Zf = SZ(n) - z0_t
					!Zf = SZ(n) - z0		! moving Z reference frame
				
				if( Zf**2 + Yf**2 + Xf**2 .LE. r**2) then					! Geometry function
						xi = xi + 1
				endif

			end do; end do; end do

		   
		  ETA(i,j,k) = xi*inv_subgrid3
	endif
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
        write(61,'(F12.7,6(3X,F12.7))') time, x0_t, y0_t, z0_t, u_solid, v_solid, w_solid
	endif
    close(61)

    if(myid==master .AND. istep == 0)then
		open (62,file='force_rotate.dat',position='append')
        write(62,*)'                 '
        write(62,*)'FORCE ROTATE'
		write(62,*) ' VARIABLES = t*,totalFX_,totalFY_,totalFZ_,rotate_sx,rotate_sy,rotate_sz'
		!write(62,*) ' VARIABLES = t*,totalFX_,totalFY_,totalFZ_'
        write(62,*)'                 '
    else if (myid==master) then
		open (62,file='force_rotate.dat',position='append')
        write(62,'(F12.7,6(3X,F14.7))') time, totalFX_, totalFY_, totalFZ_, rotate_sx, rotate_sy, rotate_sz
		!write(62,'(F12.7,3(3X,F14.7))') time, totalFX_, totalFY_, totalFZ_
	endif
    close(62)

end subroutine func_Sphere


! subroutine func_Sphere() !with MPI

! use variables
! implicit none

! integer	:: az_min,az_max,ay_min,ay_max,ax_min,ax_max
! integer 	:: l ,m ,n, xi
! real*8 	:: xf, yf, zf
! real*8 ,dimension(1:nSubGrids_f) :: SX
! real*8 ,dimension(1:nSubGrids_f) :: SY
! real*8 ,dimension(1:nSubGrids_f) :: SZ
! real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
! real*8	:: inv_subgrid3 = 1.d0/ (nSubGrids_f*nSubGrids_f**2*1.d0)

   ! !---------------------------------------------------!
   ! !    LOCAL VARIABLES                                !
   ! !---------------------------------------------------!

! ! Create coarse ETA (only 0 or 1)
! az_min = nz ; az_max = 0
! ay_min = ny ; ay_max = 0
! ax_min = nx ; ax_max = 0

! !$acc data present(X,Y,Z,Xs,Ys,Zs,iDx,iDy,iDz,ETA) create(SX,SY,SZ,ETA_3D)

! !$OMP PARALLEL DO PRIVATE(j,i,Xf,Yf,Zf) REDUCTION(max : az_max,ay_max,ax_max) &
! !$OMP			  REDUCTION(min : az_min,ay_min,ax_min) collapse(3)
! !$acc parallel vector_length(32)
! !$acc loop private(Xf,Yf,Zf) reduction(max : az_max,ay_max,ax_max) reduction(min : az_min,ay_min,ax_min) collapse(3) gang vector

! do k=1,nz
	! do j=1,ny
		! do i=1,nx

		! Xf = Xs(i) - x0_t
		! Yf = Ys(j) - y0_t
		! Zf = Zs(k) - z0_t

	! if( Zf**2 + Yf**2 + Xf**2 .LE. r**2) then				! Geometry function
		! ETA(i,j,k)=1.D0
		! az_min = MIN(k,az_min)
		! az_max = MAX(k,az_max)
		! ay_min = MIN(j,ay_min)
		! ay_max = MAX(j,ay_max)
		! ax_min = MIN(i,ax_min)
		! ax_max = MAX(i,ax_max)
	! else
		! ETA(i,j,k)=0.D0		
	! endif
		! end do
	! end do
! end do
! !$acc end parallel
! !$OMP END PARALLEL DO

! ! Setting the boundary for subgrids
	! iBgnVOS = ax_min - 1
	! iEndVOS = ax_max + 1
	! jBgnVOS = ay_min - 1
	! jEndVOS = ay_max + 1
	! kBgnVOS = az_min - 1
	! kEndVOS = az_max + 1



    ! !-----------------------MPI DIVISION-------------------------!
    ! Zdv_vos = (kEndVOS-kBgnVOS+1) / nproc
    ! Zr_vos  = (kEndVOS-kBgnVOS+1) - Zdv_vos * nproc 
    ! !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! !i = myid
    ! do i=0,(nproc-1)

        ! if(i < Zr_vos) then
            ! gstart_vos(i) = kBgnVOS + i * (Zdv_vos+1)
            ! gend0_vos(i) = gstart_vos(i) + Zdv_vos
        ! else
            ! gstart_vos(i) = kBgnVOS + i * Zdv_vos + Zr_vos
            ! gend0_vos(i) = gstart_vos(i) + Zdv_vos - 1
        ! end if
        
        ! gcount_vos(i) = gend0_vos(i) - gstart_vos(i) + 1
        ! gend_vos(i) = gcount_vos(i) + 2

    ! end do

    ! !----------for nz vos----------!
    ! istart_vos = gstart_vos(myid)  !
    ! iend_vos = gend0_vos(myid)     !
    ! igcount_vos = gcount_vos(myid) !
    ! !----------for nz vos----------!

    ! !-----------------------MPI DIVISION-------------------------!


! ! Finding potential grids for subgrids
! !$OMP PARALLEL DO COLLAPSE(3)
! !$acc parallel loop collapse(3) gang vector
! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
! do i=iBgnVOS,iEndVOS
	! if(ANY(ETA(i,j,k) /= ETA(i-1:i+1, j-1:j+1, k-1:k+1))) then
		! ETA_3D(k,j,i) = 1
	! else
		! ETA_3D(k,j,i) = 0
	! endif
! enddo;enddo;enddo
! !$acc end parallel
! !$OMP END PARALLEL DO


! ! Creating subgrids

! !$OMP PARALLEL DO PRIVATE(j,i,l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3)
! !$acc parallel vector_length(32)
! !$acc loop private(l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3) gang
! do k=istart_vos,iend_vos
! do j=jBgnVOS,jEndVOS
! do i=iBgnVOS,iEndVOS

	! if(ETA_3D(k,j,i) == 1) then

            ! !$acc loop
            ! do n=1,nSubGrids_f,1
				! SX(n) = X(i) + (n - 0.5d0)*iDx(i)*inv_subgrid1
				! SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				! SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			! enddo

		    ! xi = 0		
            ! !$acc loop private(m,n,Xf,Yf,Zf) reduction(+:xi)
            ! do l=1,nSubGrids_f,1	
            ! !$acc loop private(n,Xf,Yf,Zf) reduction(+:xi)
            ! do m=1,nSubGrids_f,1
            ! !$acc loop private(Xf,Yf,Zf) reduction(+:xi)
            ! do n=1,nSubGrids_f,1

					! Xf = SX(l) - x0_t
					! Yf = SY(m) - y0_t
					! Zf = SZ(n) - z0_t
				
				! if( Zf**2 + Yf**2 + Xf**2 .LE. r**2) then					! Geometry function
						! xi = xi + 1
				! endif

			! end do; end do; end do


		  ! ETA(i,j,k) = xi*inv_subgrid3
	! endif
! end do
! end do	
! end do	
! !$acc end parallel
! !$OMP END PARALLEL DO


		   
      ! !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
	  ! !$acc update self(ETA(:,:,istart_vos:iend_vos)) if(nproc>1)
      ! icount = igcount_vos*(nx+4)*(ny+4)
      ! !Send my results back to the master
      ! if(myid>master)then
         ! itag = 401
         ! call MPI_SEND( ETA(-1,-1,istart_vos), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      ! end if
! !      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      ! !Wait to receive results from each task
      ! if(myid==master)then
         ! do i = 1, (nproc-1)
            ! icount = gcount_vos(i)*(nx+4)*(ny+4)
            ! itag = 401
            ! call MPI_RECV( ETA(-1,-1,gstart_vos(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )    
         ! end do
      ! end if
	  ! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! !>>>>>>>>> data transformation from master to all nodes <<<<<<<<<<<<<<<<
	! icount= (nz+4)*(nx+4)*(ny+4)
	! call MPI_BCAST ( ETA, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
	! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	! !$acc update device(ETA) if(nproc>1)
	! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


! !$acc end data

    ! if(myid==master .AND. istep == 0)then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,*)'                 '
        ! write(61,*)'SOLID MOTION'
		! write(61,*) ' VARIABLES = t*,x,y,z,u_solid,v_solid,w_solid'
        ! write(61,*)'                 '
    ! else if (myid==master) then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,'(F12.7,6(3X,F12.7))') time, x0_t, y0_t, z0_t, u_solid, v_solid, w_solid
	! endif
    ! close(61)

    ! if(myid==master .AND. istep == 0)then
		! open (62,file='force_rotate.dat',position='append')
        ! write(62,*)'                 '
        ! write(62,*)'FORCE ROTATE'
		! write(62,*) ' VARIABLES = t*,totalFX_,totalFY_,totalFZ_,rotate_sx,rotate_sy,rotate_sz'
		! !write(62,*) ' VARIABLES = t*,totalFX_,totalFY_,totalFZ_'
        ! write(62,*)'                 '
    ! else if (myid==master) then
		! open (62,file='force_rotate.dat',position='append')
        ! write(62,'(F12.7,6(3X,F14.7))') time, totalFX_, totalFY_, totalFZ_, rotate_sx, rotate_sy, rotate_sz
		! !write(62,'(F12.7,3(3X,F14.7))') time, totalFX_, totalFY_, totalFZ_
	! endif
    ! close(62)

! end subroutine func_Sphere




subroutine func_torus3d() !with MPI

use variables
implicit none

integer	:: az_min,az_max,ay_min,ay_max,ax_min,ax_max
integer :: l ,m ,n, xi
real*8 	:: xf, yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SX
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid3 = 1.d0/ (nSubGrids_f*nSubGrids_f**2*1.d0)

real*8 	:: rin,rout

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

rin= 1.0
rout= 2.0


! Create coarse ETA (only 0 or 1)
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0
ax_min = nx ; ax_max = 0

!$acc data present(X,Y,Z,Xs,Ys,Zs,iDx,iDy,iDz,ETA) create(SX,SY,SZ,ETA_3D)

!$OMP PARALLEL DO PRIVATE(j,i,Xf,Yf,Zf) REDUCTION(max : az_max,ay_max,ax_max) &
!$OMP			  REDUCTION(min : az_min,ay_min,ax_min) collapse(3)
!$acc parallel vector_length(32)
!$acc loop private(Xf,Yf,Zf) reduction(max : az_max,ay_max,ax_max) reduction(min : az_min,ay_min,ax_min) collapse(3) gang vector

do k=1,nz
	do j=1,ny
		do i=1,nx

		Xf = Xs(i) - x0_t
		Yf = Ys(j) - y0_t
		Zf = Zs(k) - z0_t

		if( Zf**2 + (sqrt(Xf**2+Yf**2) - 0.5*(rin+rout))**2 .LE. 0.25*(rin-rout)**2) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
		ax_min = MIN(i,ax_min)
		ax_max = MAX(i,ax_max)
	else
		ETA(i,j,k)=0.D0		
	endif
		end do
	end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

! Setting the boundary for subgrids
	iBgnVOS = ax_min - 1
	iEndVOS = ax_max + 1
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1



    !-----------------------MPI DIVISION-------------------------!
    Zdv_vos = (kEndVOS-kBgnVOS+1) / nproc
    Zr_vos  = (kEndVOS-kBgnVOS+1) - Zdv_vos * nproc 
    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !i = myid
    do i=0,(nproc-1)

        if(i < Zr_vos) then
            gstart_vos(i) = kBgnVOS + i * (Zdv_vos+1)
            gend0_vos(i) = gstart_vos(i) + Zdv_vos
        else
            gstart_vos(i) = kBgnVOS + i * Zdv_vos + Zr_vos
            gend0_vos(i) = gstart_vos(i) + Zdv_vos - 1
        end if
        
        gcount_vos(i) = gend0_vos(i) - gstart_vos(i) + 1
        gend_vos(i) = gcount_vos(i) + 2

    end do

    !----------for nz vos----------!
    istart_vos = gstart_vos(myid)  !
    iend_vos = gend0_vos(myid)     !
    igcount_vos = gcount_vos(myid) !
    !----------for nz vos----------!

    !-----------------------MPI DIVISION-------------------------!


! Finding potential grids for subgrids
!$OMP PARALLEL DO COLLAPSE(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS,iEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i-1:i+1, j-1:j+1, k-1:k+1))) then
		ETA_3D(k,j,i) = 1
	else
		ETA_3D(k,j,i) = 0
	endif
enddo;enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,i,l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3)
!$acc parallel vector_length(32)
!$acc loop private(l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3) gang
do k=istart_vos,iend_vos
do j=jBgnVOS,jEndVOS
do i=iBgnVOS,iEndVOS

	if(ETA_3D(k,j,i) == 1) then

            !$acc loop
            do n=1,nSubGrids_f,1
				SX(n) = X(i) + (n - 0.5d0)*iDx(i)*inv_subgrid1
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    xi = 0		
            !$acc loop private(m,n,Xf,Yf,Zf) reduction(+:xi)
            do l=1,nSubGrids_f,1	
            !$acc loop private(n,Xf,Yf,Zf) reduction(+:xi)
            do m=1,nSubGrids_f,1
            !$acc loop private(Xf,Yf,Zf) reduction(+:xi)
            do n=1,nSubGrids_f,1

					Xf = SX(l) - x0_t
					Yf = SY(m) - y0_t
					Zf = SZ(n) - z0_t
				
				if( Zf**2 + (sqrt(Xf**2+Yf**2) - 0.5*(rin+rout))**2 .LE. 0.25*(rin-rout)**2) then
						xi = xi + 1
				endif

			end do; end do; end do


		  ETA(i,j,k) = xi*inv_subgrid3
	endif
end do
end do	
end do	
!$acc end parallel
!$OMP END PARALLEL DO


		   
      !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
	  !$acc update self(ETA(:,:,istart_vos:iend_vos)) if(nproc>1)
      icount = igcount_vos*(nx+4)*(ny+4)
      !Send my results back to the master
      if(myid>master)then
         itag = 401
         call MPI_SEND( ETA(-1,-1,istart_vos), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      !Wait to receive results from each task
      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount_vos(i)*(nx+4)*(ny+4)
            itag = 401
            call MPI_RECV( ETA(-1,-1,gstart_vos(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )    
         end do
      end if
	  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	!>>>>>>>>> data transformation from master to all nodes <<<<<<<<<<<<<<<<
	icount= (nz+4)*(nx+4)*(ny+4)
	call MPI_BCAST ( ETA, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	!$acc update device(ETA) if(nproc>1)
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,x,y,z,u_solid,v_solid,w_solid'
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,6(3X,F12.7))') time, x0_t, y0_t, z0_t, u_solid, v_solid, w_solid
	endif
    close(61)

end subroutine func_torus3d



! subroutine func_rbc()

! use variables
! implicit none

! integer	:: az_min,az_max,ay_min,ay_max,ax_min,ax_max
! integer 	:: l ,m ,n, xi
! real*8 	:: xf, yf, zf
! real*8 ,dimension(1:nSubGrids_f) :: SX
! real*8 ,dimension(1:nSubGrids_f) :: SY
! real*8 ,dimension(1:nSubGrids_f) :: SZ
! real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
! real*8	:: inv_subgrid3 = 1.d0/ (nSubGrids_f*nSubGrids_f**2*1.d0)

! real*8 	:: od,bo,ho,pp,qq,rr


   ! !---------------------------------------------------!
   ! !    LOCAL VARIABLES                                !
   ! !---------------------------------------------------!

! od = 8.d0
! bo = 1.d0
! ho = 2.12d0

! pp = -0.5d0*od**2 + 0.5d0*ho**2 * ((od/bo)**2 - 1.d0) - 0.5d0*ho**2 * ((od/bo)**2 - 1.d0)*DSQRT(1.d0 - (bo/ho)**2)
! qq = pp*od**2/bo**2 + 0.25d0*bo**2 * ((od/bo)**4 - 1.d0)
! rr = -1.d0*pp*0.25d0*od**2 - 0.0625d0*od**4 


! ! Create coarse ETA (only 0 or 1)
! az_min = nz ; az_max = 0
! ay_min = ny ; ay_max = 0
! ax_min = nx ; ax_max = 0

! !$acc data present(X,Y,Z,Xs,Ys,Zs,iDx,iDy,iDz,ETA) create(SX,SY,SZ,ETA_3D)

! !$OMP PARALLEL DO PRIVATE(j,i,Xf,Yf,Zf) REDUCTION(max : az_max,ay_max,ax_max) &
! !$OMP			  REDUCTION(min : az_min,ay_min,ax_min) collapse(3)
! !$acc parallel vector_length(32)
! !$acc loop private(Xf,Yf,Zf) reduction(max : az_max,ay_max,ax_max) reduction(min : az_min,ay_min,ax_min) collapse(3) gang vector

! do k=1,nz
	! do j=1,ny
		! do i=1,nx

		! Xf = Xs(i) - x0_t
		! Yf = Ys(j) - y0_t
		! Zf = Zs(k) - z0_t

	! if( (xf**2 + yf**2 + zf**2)**2 + pp*(yf**2 + zf**2) + qq*xf**2 + rr .LE.0.d0 ) then
		! ETA(i,j,k)=1.D0
		! az_min = MIN(k,az_min)
		! az_max = MAX(k,az_max)
		! ay_min = MIN(j,ay_min)
		! ay_max = MAX(j,ay_max)
		! ax_min = MIN(i,ax_min)
		! ax_max = MAX(i,ax_max)
	! else
		! ETA(i,j,k)=0.D0		
	! endif
		! end do
	! end do
! end do
! !$acc end parallel
! !$OMP END PARALLEL DO

! ! Setting the boundary for subgrids
	! iBgnVOS = ax_min - 1
	! iEndVOS = ax_max + 1
	! jBgnVOS = ay_min - 1
	! jEndVOS = ay_max + 1
	! kBgnVOS = az_min - 1
	! kEndVOS = az_max + 1


! ! Finding potential grids for subgrids
! !$OMP PARALLEL DO COLLAPSE(3)
! !$acc parallel loop collapse(3) gang vector
! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
! do i=iBgnVOS,iEndVOS
	! if(ANY(ETA(i,j,k) /= ETA(i-1:i+1, j-1:j+1, k-1:k+1))) then
		! ETA_3D(k,j,i) = 1
	! else
		! ETA_3D(k,j,i) = 0
	! endif
! enddo;enddo;enddo
! !$acc end parallel
! !$OMP END PARALLEL DO


! ! Creating subgrids

! !$OMP PARALLEL DO PRIVATE(j,i,l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3)
! !$acc parallel vector_length(32)
! !$acc loop private(l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3) gang
! do k=kBgnVOS,kEndVOS
! do j=jBgnVOS,jEndVOS
! do i=iBgnVOS,iEndVOS

	! if(ETA_3D(k,j,i) == 1) then

            ! !$acc loop
            ! do n=1,nSubGrids_f,1
				! SX(n) = X(i) + (n - 0.5d0)*iDx(i)*inv_subgrid1
				! SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				! SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			! enddo

		    ! xi = 0		
            ! !$acc loop private(m,n,Xf,Yf,Zf) reduction(+:xi)
            ! do l=1,nSubGrids_f,1	
            ! !$acc loop private(n,Xf,Yf,Zf) reduction(+:xi)
            ! do m=1,nSubGrids_f,1
            ! !$acc loop private(Xf,Yf,Zf) reduction(+:xi)
            ! do n=1,nSubGrids_f,1

					! Xf = SX(l) - x0_t
					! Yf = SY(m) - y0_t
					! Zf = SZ(n) - z0_t
				
				! if( (xf**2 + yf**2 + zf**2)**2 + pp*(yf**2 + zf**2) + qq*xf**2 + rr .LE.0.d0 ) then
						! xi = xi + 1
				! endif

			! end do; end do; end do
		   
		  ! ETA(i,j,k) = xi*inv_subgrid3
	! endif
! end do
! end do	
! end do	
! !$acc end parallel
! !$OMP END PARALLEL DO
   

! !$acc end data

    ! if(myid==master .AND. istep == 0)then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,*)'                 '
        ! write(61,*)'SOLID MOTION'
		! write(61,*) ' VARIABLES = t*,x,y,z,u_solid,v_solid,w_solid'
        ! write(61,*)'                 '
    ! else if (myid==master) then
		! open (61,file='solid_motion.dat',position='append')
        ! write(61,'(F12.7,6(3X,F12.7))') time, x0_t, y0_t, z0_t, u_solid, v_solid, w_solid
	! endif
    ! close(61)

! end subroutine func_rbc


subroutine func_rbc() !with MPI

use variables
implicit none

integer	:: az_min,az_max,ay_min,ay_max,ax_min,ax_max
integer :: l ,m ,n, xi
real*8 	:: xf, yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SX
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid3 = 1.d0/ (nSubGrids_f*nSubGrids_f**2*1.d0)

real*8 	:: od,bo,ho,pp,qq,rr

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

od = 8.d0
bo = 1.d0
ho = 2.12d0

pp = -0.5d0*od**2 + 0.5d0*ho**2 * ((od/bo)**2 - 1.d0) - 0.5d0*ho**2 * ((od/bo)**2 - 1.d0)*DSQRT(1.d0 - (bo/ho)**2)
qq = pp*od**2/bo**2 + 0.25d0*bo**2 * ((od/bo)**4 - 1.d0)
rr = -1.d0*pp*0.25d0*od**2 - 0.0625d0*od**4 


! Create coarse ETA (only 0 or 1)
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0
ax_min = nx ; ax_max = 0

!$acc data present(X,Y,Z,Xs,Ys,Zs,iDx,iDy,iDz,ETA) create(SX,SY,SZ,ETA_3D)

!$OMP PARALLEL DO PRIVATE(j,i,Xf,Yf,Zf) REDUCTION(max : az_max,ay_max,ax_max) &
!$OMP			  REDUCTION(min : az_min,ay_min,ax_min) collapse(3)
!$acc parallel vector_length(32)
!$acc loop private(Xf,Yf,Zf) reduction(max : az_max,ay_max,ax_max) reduction(min : az_min,ay_min,ax_min) collapse(3) gang vector

do k=1,nz
	do j=1,ny
		do i=1,nx

		Xf = Xs(i) - x0_t
		Yf = Ys(j) - y0_t
		Zf = Zs(k) - z0_t

	if( (xf**2 + yf**2 + zf**2)**2 + pp*(yf**2 + zf**2) + qq*xf**2 + rr .LE.0.d0 ) then
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
		ax_min = MIN(i,ax_min)
		ax_max = MAX(i,ax_max)
	else
		ETA(i,j,k)=0.D0		
	endif
		end do
	end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

! Setting the boundary for subgrids
	iBgnVOS = ax_min - 1
	iEndVOS = ax_max + 1
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1



    !-----------------------MPI DIVISION-------------------------!
    Zdv_vos = (kEndVOS-kBgnVOS+1) / nproc
    Zr_vos  = (kEndVOS-kBgnVOS+1) - Zdv_vos * nproc 
    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !i = myid
    do i=0,(nproc-1)

        if(i < Zr_vos) then
            gstart_vos(i) = kBgnVOS + i * (Zdv_vos+1)
            gend0_vos(i) = gstart_vos(i) + Zdv_vos
        else
            gstart_vos(i) = kBgnVOS + i * Zdv_vos + Zr_vos
            gend0_vos(i) = gstart_vos(i) + Zdv_vos - 1
        end if
        
        gcount_vos(i) = gend0_vos(i) - gstart_vos(i) + 1
        gend_vos(i) = gcount_vos(i) + 2

    end do

    !----------for nz vos----------!
    istart_vos = gstart_vos(myid)  !
    iend_vos = gend0_vos(myid)     !
    igcount_vos = gcount_vos(myid) !
    !----------for nz vos----------!

    !-----------------------MPI DIVISION-------------------------!


! Finding potential grids for subgrids
!$OMP PARALLEL DO COLLAPSE(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS,iEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i-1:i+1, j-1:j+1, k-1:k+1))) then
		ETA_3D(k,j,i) = 1
	else
		ETA_3D(k,j,i) = 0
	endif
enddo;enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids 

!$OMP PARALLEL DO PRIVATE(j,i,l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3)
!$acc parallel vector_length(32)
!$acc loop private(l,m,n,SX,SY,SZ,Xf,Yf,Zf,xi) collapse(3) gang
do k=istart_vos,iend_vos
do j=jBgnVOS,jEndVOS
do i=iBgnVOS,iEndVOS

	if(ETA_3D(k,j,i) == 1) then

            !$acc loop
            do n=1,nSubGrids_f,1
				SX(n) = X(i) + (n - 0.5d0)*iDx(i)*inv_subgrid1
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    xi = 0		
            !$acc loop private(m,n,Xf,Yf,Zf) reduction(+:xi)
            do l=1,nSubGrids_f,1	
            !$acc loop private(n,Xf,Yf,Zf) reduction(+:xi)
            do m=1,nSubGrids_f,1
            !$acc loop private(Xf,Yf,Zf) reduction(+:xi)
            do n=1,nSubGrids_f,1

					Xf = SX(l) - x0_t
					Yf = SY(m) - y0_t
					Zf = SZ(n) - z0_t
				
				if( (xf**2 + yf**2 + zf**2)**2 + pp*(yf**2 + zf**2) + qq*xf**2 + rr .LE.0.d0 ) then
						xi = xi + 1
				endif

			end do; end do; end do


		  ETA(i,j,k) = xi*inv_subgrid3
	endif
end do
end do	
end do	
!$acc end parallel
!$OMP END PARALLEL DO


	  !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
	  !$acc update self(ETA(:,:,istart_vos:iend_vos)) if(nproc>1)
      icount = igcount_vos*(nx+4)*(ny+4)
      !Send my results back to the master
      if(myid>master)then
         itag = 401
         call MPI_SEND( ETA(-1,-1,istart_vos), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      !Wait to receive results from each task
      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount_vos(i)*(nx+4)*(ny+4)
            itag = 401
            call MPI_RECV( ETA(-1,-1,gstart_vos(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )    
         end do
      end if
	  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	  

	!>>>>>>>>> data transformation from master to all nodes <<<<<<<<<<<<<<<<
	icount= (nz+4)*(nx+4)*(ny+4)
	call MPI_BCAST ( ETA, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	!$acc update device(ETA) if(nproc>1)
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,x,y,z,u_solid,v_solid,w_solid'
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,6(3X,F12.7))') time, x0_t, y0_t, z0_t, u_solid, v_solid, w_solid
	endif
    close(61)

end subroutine func_rbc



subroutine func_64Spheres() !with MPI

use variables
implicit none

integer	:: az_min,az_max,ay_min,ay_max,ax_min,ax_max
integer :: l ,m ,n, xi
real*8 	:: xf, yf, zf
real*8 ,dimension(1:nSubGrids_f) :: SX
real*8 ,dimension(1:nSubGrids_f) :: SY
real*8 ,dimension(1:nSubGrids_f) :: SZ
real*8	:: inv_subgrid1 = 1.d0/ (nSubGrids_f*1.d0)
real*8	:: inv_subgrid3 = 1.d0/ (nSubGrids_f*nSubGrids_f**2*1.d0)


integer :: nn,ii
real*8 ,dimension(1:64) :: x0_,y0_,z0_,sphere


   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

		

! Create coarse ETA (only 0 or 1)
az_min = nz ; az_max = 0
ay_min = ny ; ay_max = 0
ax_min = nx ; ax_max = 0

!$acc data present(X,Y,Z,Xs,Ys,Zs,iDx,iDy,iDz,ETA) create(SX,SY,SZ,x0_,y0_,z0_,sphere,ETA_3D)


	ii = 1

	do i=1,4
		do j=1,4
			do k=1,4
				x0_(ii) = (i-1)*1.25d0+x0
				y0_(ii) = (j-1)*1.25d0+y0
				z0_(ii) = (k-1)*1.25d0+z0
				ii = ii+1				
	enddo;enddo;enddo


!$acc update device(x0_,y0_,z0_)


!$OMP PARALLEL DO PRIVATE(j,i,nn,sphere) REDUCTION(max : az_max,ay_max,ax_max) &
!$OMP			  REDUCTION(min : az_min,ay_min,ax_min) collapse(3)
!$acc parallel vector_length(32)
!$acc loop private(nn,sphere) reduction(max : az_max,ay_max,ax_max) reduction(min : az_min,ay_min,ax_min) collapse(3) gang vector

do k=1,nz
	do j=1,ny
		do i=1,nx

		!$acc loop
		do nn=1,64
			sphere(nn) = (Xs(i) - x0_(nn))**2 + (Ys(j) - y0_(nn))**2 + (Zs(k) - z0_(nn))**2
		enddo


	if( ANY(sphere(1:64) .LE. r**2) )then				! Geometry function
		ETA(i,j,k)=1.D0
		az_min = MIN(k,az_min)
		az_max = MAX(k,az_max)
		ay_min = MIN(j,ay_min)
		ay_max = MAX(j,ay_max)
		ax_min = MIN(i,ax_min)
		ax_max = MAX(i,ax_max)
	else
		ETA(i,j,k)=0.D0		
	endif
		end do
	end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

! Setting the boundary for subgrids
	iBgnVOS = ax_min - 1
	iEndVOS = ax_max + 1
	jBgnVOS = ay_min - 1
	jEndVOS = ay_max + 1
	kBgnVOS = az_min - 1
	kEndVOS = az_max + 1



    !-----------------------MPI DIVISION-------------------------!
    Zdv_vos = (kEndVOS-kBgnVOS+1) / nproc
    Zr_vos  = (kEndVOS-kBgnVOS+1) - Zdv_vos * nproc 
    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !i = myid
    do i=0,(nproc-1)

        if(i < Zr_vos) then
            gstart_vos(i) = kBgnVOS + i * (Zdv_vos+1)
            gend0_vos(i) = gstart_vos(i) + Zdv_vos
        else
            gstart_vos(i) = kBgnVOS + i * Zdv_vos + Zr_vos
            gend0_vos(i) = gstart_vos(i) + Zdv_vos - 1
        end if
        
        gcount_vos(i) = gend0_vos(i) - gstart_vos(i) + 1
        gend_vos(i) = gcount_vos(i) + 2

    end do

    !----------for nz vos----------!
    istart_vos = gstart_vos(myid)  !
    iend_vos = gend0_vos(myid)     !
    igcount_vos = gcount_vos(myid) !
    !----------for nz vos----------!

    !-----------------------MPI DIVISION-------------------------!


! Finding potential grids for subgrids
!$OMP PARALLEL DO COLLAPSE(3)
!$acc parallel loop collapse(3) gang vector
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS,iEndVOS
	if(ANY(ETA(i,j,k) /= ETA(i-1:i+1, j-1:j+1, k-1:k+1))) then
		ETA_3D(k,j,i) = 1
	else
		ETA_3D(k,j,i) = 0
	endif
enddo;enddo;enddo
!$acc end parallel
!$OMP END PARALLEL DO


! Creating subgrids

!$OMP PARALLEL DO PRIVATE(j,i,l,m,n,nn,SX,SY,SZ,xi,sphere) collapse(3)
!$acc parallel vector_length(32)
!$acc loop private(l,m,n,nn,SX,SY,SZ,xi,sphere) collapse(3) gang
do k=istart_vos,iend_vos
do j=jBgnVOS,jEndVOS
do i=iBgnVOS,iEndVOS

	if(ETA_3D(k,j,i) == 1) then

            !$acc loop
            do n=1,nSubGrids_f,1
				SX(n) = X(i) + (n - 0.5d0)*iDx(i)*inv_subgrid1
				SY(n) = Y(j) + (n - 0.5d0)*iDy(j)*inv_subgrid1
				SZ(n) = Z(k) + (n - 0.5d0)*iDz(k)*inv_subgrid1
			enddo

		    xi = 0		
            !$acc loop private(m,n) reduction(+:xi)
            do l=1,nSubGrids_f,1	
            !$acc loop private(n) reduction(+:xi)
            do m=1,nSubGrids_f,1
            !$acc loop reduction(+:xi)
            do n=1,nSubGrids_f,1

            !$acc loop
			do nn=1,64
				sphere(nn) = (SX(l) - x0_(nn))**2 + (SY(m) - y0_(nn))**2 + (SZ(n) - z0_(nn))**2
			enddo
				
				if( ANY(sphere(1:64) .LE. r**2) )then					! Geometry function
						xi = xi + 1
				endif

			end do; end do; end do


		  ETA(i,j,k) = xi*inv_subgrid3
	endif
end do
end do	
end do	
!$acc end parallel
!$OMP END PARALLEL DO


		   
      !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
	  !$acc update self(ETA(:,:,istart_vos:iend_vos)) if(nproc>1)
      icount = igcount_vos*(nx+4)*(ny+4)
      !Send my results back to the master
      if(myid>master)then
         itag = 401
         call MPI_SEND( ETA(-1,-1,istart_vos), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      !Wait to receive results from each task
      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount_vos(i)*(nx+4)*(ny+4)
            itag = 401
            call MPI_RECV( ETA(-1,-1,gstart_vos(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )    
         end do
      end if
	  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


	!>>>>>>>>> data transformation from master to all nodes <<<<<<<<<<<<<<<<
	icount= (nz+4)*(nx+4)*(ny+4)
	call MPI_BCAST ( ETA, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	!$acc update device(ETA) if(nproc>1)
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!$acc end data

    if(myid==master .AND. istep == 0)then
		open (61,file='solid_motion.dat',position='append')
        write(61,*)'                 '
        write(61,*)'SOLID MOTION'
		write(61,*) ' VARIABLES = t*,x,y,z,u_solid,v_solid,w_solid'
        write(61,*)'                 '
    else if (myid==master) then
		open (61,file='solid_motion.dat',position='append')
        write(61,'(F12.7,6(3X,F12.7))') time, x0_t, y0_t, z0_t, u_solid, v_solid, w_solid
	endif
    close(61)


end subroutine func_64Spheres
