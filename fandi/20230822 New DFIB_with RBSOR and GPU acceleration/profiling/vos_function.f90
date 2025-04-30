! 28 Jul 2023 - FDS


! subroutine dynamic_AOA() 
    ! use variables
    ! implicit none

    ! !AOA = AOA1 + AOA_amp * SIN( 2.d0 * reduce_frequency_k * ( time - StartDynamic_time ) )
    ! AOA = AOA1 + AOA_amp * ( SIN( 2.d0 * reduce_frequency_k * ( time - StartDynamic_time ) + 1.5d0*PI ) + 1.d0 )      !FROM 0

    ! angular_vel = 2.d0*reduce_frequency_k*AOA_amp*COS( 2.d0 * reduce_frequency_k * ( time - StartDynamic_time ) + 1.5d0*PI ) / 180.d0 * PI

! end subroutine dynamic_AOA


subroutine rotor_dynamic() 
    use variables 
    implicit none

	rotor_omega = rotor_tsr*U_inf/(rotor_r+blade_r)
	
	AOA = AOA1 + rotor_omega*( time - StartDynamic_time )*180.d0/PI
	
end subroutine rotor_dynamic


subroutine blade_dynamic() 
    use variables
    implicit none

    blade_omega = blade_alpha*U_inf/blade_r

end subroutine blade_dynamic


subroutine func_Cylplate()

use variables
implicit none
real*8 :: L_p, t_p, B_p, g_p
real*8 :: z_trans, y_trans, z_trans1, y_trans1
real*8 :: z_aoa, y_aoa
INTEGER, DIMENSION(0:1) :: az, ay
integer :: l ,m ,n
real*8 :: xi
real*8 :: dxg, dyg, dzg
real*8 :: solid_volume, ETAs
real*8 ,dimension(1:nSubGrids_f+1) :: SX
real*8 ,dimension(1:nSubGrids_f+1) :: SY
real*8 ,dimension(1:nSubGrids_f+1) :: SZ

   integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

L_p = 1.d0 ! Plate's length
t_p = 0.05d0 ! Plate's thickness
B_p = (-30.d0/180.d0)*PI ! Plate's Beta angle
g_p = 0.1d0 ! Plate's gap


	! Rotate the blade center with AOA (needed to define the blade rotation)
	! Blades' center is shifted from (0,0) to (-rotor_r, 0)

	
		yc_t= dSIN((AOA)*PI/180.d0)*(-rotor_r) + dCOS((AOA)*PI/180.d0)*0.d0
		zc_t= dCOS((AOA)*PI/180.d0)*(-rotor_r) - dSIN((AOA)*PI/180.d0)*0.d0
		
		yc_t = yc_t + yc
		zc_t = zc_t + zc

! Create coarse ETA (only 0 or 1)
i=0
az(0) = nz ; az(1) = 0
ay(0) = ny ; ay(1) = 0
do k=0,nz+1
	do j=0,ny+1

	! rotation wrt. AOA
	z_aoa =  (Zs(k) - zc)*dcos((-AOA)*PI/180.d0)-(Ys(j) - yc)*dsin((-AOA)*PI/180.d0)
	y_aoa =  (Zs(k) - zc)*dsin((-AOA)*PI/180.d0)+(Ys(j) - yc)*dcos((-AOA)*PI/180.d0)

	! z and y transformation for flat plate
	z_trans1 = z_aoa + rotor_r
	y_trans1 = y_aoa + (g_p+r)+0.5d0*L_p
	
	z_trans = z_trans1*dcos(B_p) - (y_trans1-(g_p+r)-0.5d0*L_p)*dsin(B_p)
	y_trans = z_trans1*dsin(B_p) + (y_trans1-(g_p+r)-0.5d0*L_p)*dcos(B_p) + (g_p+r)+0.5d0*L_p
	

	if( ((z_aoa+rotor_r)**2 + y_aoa**2) .LE. r**2 .OR. &
		(dabs(z_trans/t_p + y_trans/L_p) + dabs(z_trans/t_p - y_trans/L_p)) .LE. 1) then
		ETA(i,j,k)=1.D0
		if (k < az(0)) then
			az(0) = k
		else if (k > az(1)) then
			az(1) = k
		endif
		if (j < ay(0)) then
			ay(0) = j
		else if (j > ay(1)) then
			ay(1) = j
		endif

	else
		ETA(i,j,k)=0.D0		
	endif
	
	end do
end do

! Setting the boundary for subgrids
	iBgnVOS = 0
	iEndVOS = nx+1
	jBgnVOS = ay(0) - 1
	jEndVOS = ay(1) + 1
	kBgnVOS = az(0) - 1
	kEndVOS = az(1) + 1

! Creating subgrids
i=iBgnVOS
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
				z_aoa =  (SZ(n) - zc)*dcos((-AOA)*PI/180.d0)-(SY(m) - yc)*dsin((-AOA)*PI/180.d0)
				y_aoa =  (SZ(n) - zc)*dsin((-AOA)*PI/180.d0)+(SY(m) - yc)*dcos((-AOA)*PI/180.d0)

				! z and y transformation for flat plate
				z_trans1 = z_aoa + rotor_r
				y_trans1 = y_aoa + (g_p+r)+0.5d0*L_p
	
				z_trans = z_trans1*dcos(B_p) - (y_trans1-(g_p+r)-0.5d0*L_p)*dsin(B_p)
				y_trans = z_trans1*dsin(B_p) + (y_trans1-(g_p+r)-0.5d0*L_p)*dcos(B_p) + (g_p+r)+0.5d0*L_p

				if( ((z_aoa+rotor_r)**2 + y_aoa**2) .LE. r**2 .OR. &
					(dabs(z_trans/t_p + y_trans/L_p) + dabs(z_trans/t_p - y_trans/L_p)) .LE. 1) then
						xi = xi + 1
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	endif
end do
end do		


do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)

end do
end do
end do    

    !-----------------------calculating casting area----------------------!
    ! Z_castingarea=0.0
    ! do j=1,ny; do i=1,nx
    ! ETAs=0.0
    ! do k=1,nz 
        ! if ( ETA(i,j,k)>ETAs ) then
            ! ETAs=ETA(i,j,k)
        ! endif 
    ! enddo
    ! Z_castingarea=Z_castingarea+ETAs*iDx(i)*iDy(j)
    ! enddo; enddo

    ! Y_castingarea=0.0
    ! do k=1,nz; do i=1,nx
    ! ETAs=0.0
    ! do j=1,ny 
        ! if ( ETA(i,j,k)>ETAs ) then
            ! ETAs=ETA(i,j,k)
        ! endif 
    ! enddo
    ! Y_castingarea=Y_castingarea+ETAs*iDx(i)*iDz(k)
    ! enddo; enddo
	
	! X_castingarea=0.0
    ! do j=1,ny; do k=1,nz
    ! ETAs=0.0
    ! do i=1,nx 
        ! if ( ETA(i,j,k)>ETAs ) then
            ! ETAs=ETA(i,j,k)
        ! endif 
    ! enddo
    ! X_castingarea=X_castingarea+ETAs*iDy(j)*iDz(k)
    ! enddo; enddo

    ! Solid_volume=0.0
    ! do k=1,nz; do j=1,ny; do i=1,nx
        ! solid_volume=solid_volume + ETA(i,j,k)*iDx(i)*iDy(j)*iDz(k)
    ! enddo; enddo; enddo


	if (select_ref_area == 1) then						! Reference area in virtual force calculation
		ref_area =  Z_castingarea
    elseif (select_ref_area == 2) then		
		ref_area =  Y_castingarea
    elseif (select_ref_area == 3) then		
		ref_area =  X_castingarea
    elseif (select_ref_area == 4) then		
		ref_area =  user_defined
	endif


    !-----------------------calculating casting area----------------------!
        
    totalfinaltime = MPI_WTIME()
    ! if(myid==master)then
		! open (61,file='Geometry_info.dat',position='append')
        ! write(61,*)'                '
        ! write(61,*)'X_casting area = ',X_castingarea
        ! write(61,*)'Y_casting area = ',Y_castingarea
        ! write(61,*)'Z_casting area = ',Z_castingarea
        ! write(61,*)'Reference area for virtual forces = ', ref_area
		! write(61,*)'                 '
        ! write(61,*)'Casting_volume = ', solid_volume
        ! write(61,*)'                 '
        ! write(61,*)'Conversion cost = ',totalfinaltime-totalstarttime
        ! write(61,*)'                 '
        ! write(61,*)'---------------------------------------------------------'
        ! write(61,*)'                 '
        ! write(61,*)'                 '
    ! endif
    ! close(61)

end subroutine func_Cylplate


subroutine func_Cylinder()

use variables
implicit none
real*8 :: L_p, t_p, B_p, g_p
real*8 :: z_trans, y_trans
INTEGER, DIMENSION(0:1) :: az, ay
real*8 :: DISTANCE,SDIST
real*8 :: DIAGONAL
integer :: l ,m ,n
real*8 :: xi
real*8 :: dxg, dyg, dzg
real*8 :: solid_volume, ETAs
real*8 ,dimension(1:nSubGrids_f+1) :: SX
real*8 ,dimension(1:nSubGrids_f+1) :: SY
real*8 ,dimension(1:nSubGrids_f+1) :: SZ

   integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

		yc_t = yc
		zc_t = zc

! Create coarse ETA (only 0 or 1)
i=0
az(0) = nz ; az(1) = 0
ay(0) = ny ; ay(1) = 0
do k=0,nz+1
	do j=0,ny+1


	if( ((Zs(k) - zc)**2 + (Ys(j) - yc)**2) .LE. r**2) then
		ETA(i,j,k)=1.D0
		if (k < az(0)) then
			az(0) = k
		else if (k > az(1)) then
			az(1) = k
		endif
		if (j < ay(0)) then
			ay(0) = j
		else if (j > ay(1)) then
			ay(1) = j
		endif

	else
		ETA(i,j,k)=0.D0		
	endif
	
	end do
end do

! Setting the boundary for subgrids
	iBgnVOS = 0
	iEndVOS = nx+1
	jBgnVOS = ay(0) - 1
	jEndVOS = ay(1) + 1
	kBgnVOS = az(0) - 1
	kEndVOS = az(1) + 1

! Creating subgrids
i=iBgnVOS
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	! Finding potential grids for subgrids
	if( (abs(ETA(i,j,k)-ETA(i,j-1,k))+abs(ETA(i,j,k)-ETA(i,j+1,k))+ &
		abs(ETA(i,j,k)-ETA(i,j,k-1))+abs(ETA(i,j,k)-ETA(i,j,k+1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k-1))+abs(ETA(i,j,k)-ETA(i,j+1,k-1))+ &
		abs(ETA(i,j,k)-ETA(i,j-1,k+1))+abs(ETA(i,j,k)-ETA(i,j+1,k+1))) .GT. 0.5d0) then
			do m=1,nSubGrids_f+1; do n=1,nSubGrids_f+1
				dzg = iDz(k) / nSubGrids_f
				dyg = iDy(j) / nSubGrids_f
				SZ(n) = Z(k) + (n-1) * dzg
				SY(m) = Y(j) + (m-1) * dyg		
			end do; end do

		    xi = 0.0
			do m=1,nSubGrids_f; do n=1,nSubGrids_f
				if( ((SZ(n) - zc)**2 + (SY(m) - yc)**2) .LE. r**2) then
						xi = xi + 1
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	endif
end do
end do		


do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)

end do
end do
end do    

    !-----------------------calculating casting area----------------------!
    Z_castingarea=0.0
    do j=1,ny; do i=1,nx
    ETAs=0.0
    do k=1,nz 
        if ( ETA(i,j,k)>ETAs ) then
            ETAs=ETA(i,j,k)
        endif 
    enddo
    Z_castingarea=Z_castingarea+ETAs*iDx(i)*iDy(j)
    enddo; enddo

    Y_castingarea=0.0
    do k=1,nz; do i=1,nx
    ETAs=0.0
    do j=1,ny 
        if ( ETA(i,j,k)>ETAs ) then
            ETAs=ETA(i,j,k)
        endif 
    enddo
    Y_castingarea=Y_castingarea+ETAs*iDx(i)*iDz(k)
    enddo; enddo
	
	X_castingarea=0.0
    do j=1,ny; do k=1,nz
    ETAs=0.0
    do i=1,nx 
        if ( ETA(i,j,k)>ETAs ) then
            ETAs=ETA(i,j,k)
        endif 
    enddo
    X_castingarea=X_castingarea+ETAs*iDy(j)*iDz(k)
    enddo; enddo

    Solid_volume=0.0
    do k=1,nz; do j=1,ny; do i=1,nx
        solid_volume=solid_volume + ETA(i,j,k)*iDx(i)*iDy(j)*iDz(k)
    enddo; enddo; enddo


	if (select_ref_area == 1) then						! Reference area in virtual force calculation
		ref_area =  Z_castingarea
    elseif (select_ref_area == 2) then		
		ref_area =  Y_castingarea
    elseif (select_ref_area == 3) then		
		ref_area =  X_castingarea
    elseif (select_ref_area == 4) then		
		ref_area =  user_defined
	endif


    !-----------------------calculating casting area----------------------!
        
    totalfinaltime = MPI_WTIME()
    if(myid==master)then
        write(61,*)'                '
        write(61,*)'X_casting area = ',X_castingarea
        write(61,*)'Y_casting area = ',Y_castingarea
        write(61,*)'Z_casting area = ',Z_castingarea
        write(61,*)'Reference area for virtual forces = ', ref_area
		write(61,*)'                 '
        write(61,*)'Casting_volume = ', solid_volume
        write(61,*)'                 '
        write(61,*)'Conversion cost = ',totalfinaltime-totalstarttime
        write(61,*)'                 '
        write(61,*)'---------------------------------------------------------'
        write(61,*)'                 '
        write(61,*)'                 '
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
real*4 :: solid_volume, ETAs
real*4 ,dimension(1:nSubGrids_f+1) :: SX
real*4 ,dimension(1:nSubGrids_f+1) :: SY
real*4 ,dimension(1:nSubGrids_f+1) :: SZ

   integer          :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS

   !---------------------------------------------------!
   !    LOCAL VARIABLES                                !
   !---------------------------------------------------!

		yc_t = yc
		zc_t = zc

	iBgnVOS = 1
	iEndVOS = nx
	jBgnVOS = 1
	jEndVOS = ny
	kBgnVOS = 1
	kEndVOS = nz

do k=kBgnVOS,kEndVOS 
  do j=jBgnVOS,jEndVOS
    do i=iBgnVOS,iEndVOS

         DIAGONAL = sqrt( iDx(i)*iDx(i) + iDy(j)*iDy(j) + iDz(k)*iDz(k) ) / 2.D0
         DISTANCE = sqrt((((X(i)+X(i+1))/2.D0)-xc)**2.D0+ &
                         (((Y(j)+Y(j+1))/2.D0)-yc)**2.D0+ &
                         (((Z(k)+Z(k+1))/2.D0)-zc)**2.D0  )

           if( abs(DISTANCE - r) < DIAGONAL ) then
              
              do l=1,nSubGrids_f+1,1
              do m=1,nSubGrids_f+1,1
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
              do l=1,nSubGrids_f
                do m=1,nSubGrids_f
                  do n=1,nSubGrids_f
                    
                    SDIST=SQRT((((SX(n)+SX(n+1))/2.D0)-xc)**2.D0+ &
                               (((SY(m)+SY(m+1))/2.D0)-yc)**2.D0+ &
                               (((SZ(l)+SZ(l+1))/2.D0)-zc)**2.D0  )
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

    !-----------------------calculating casting area----------------------!
    Z_castingarea=0.0
    do j=1,ny; do i=1,nx
    ETAs=0.0
    do k=1,nz 
        if ( ETA(i,j,k)>ETAs ) then
            ETAs=ETA(i,j,k)
        endif 
    enddo
    Z_castingarea=Z_castingarea+ETAs*iDx(i)*iDy(j)
    enddo; enddo

    Y_castingarea=0.0
    do k=1,nz; do i=1,nx
    ETAs=0.0
    do j=1,ny 
        if ( ETA(i,j,k)>ETAs ) then
            ETAs=ETA(i,j,k)
        endif 
    enddo
    Y_castingarea=Y_castingarea+ETAs*iDx(i)*iDz(k)
    enddo; enddo
	
	X_castingarea=0.0
    do j=1,ny; do k=1,nz
    ETAs=0.0
    do i=1,nx 
        if ( ETA(i,j,k)>ETAs ) then
            ETAs=ETA(i,j,k)
        endif 
    enddo
    X_castingarea=X_castingarea+ETAs*iDy(j)*iDz(k)
    enddo; enddo

    Solid_volume=0.0
    do k=1,nz; do j=1,ny; do i=1,nx
        solid_volume=solid_volume + ETA(i,j,k)*iDx(i)*iDy(j)*iDz(k)
    enddo; enddo; enddo


	ref_area = Z_castingarea 									! Reference area in virtual force calculation


    !-----------------------calculating casting area----------------------!
        
    totalfinaltime = MPI_WTIME()
    if(myid==master)then
        write(61,*)'                '
        write(61,*)'X_casting area = ',X_castingarea
        write(61,*)'Y_casting area = ',Y_castingarea
        write(61,*)'Z_casting area = ',Z_castingarea
        write(61,*)'Reference area for virtual forces = ', ref_area
		write(61,*)'                 '
        write(61,*)'Casting_volume = ', solid_volume
        write(61,*)'                 '
        write(61,*)'Conversion cost = ',totalfinaltime-totalstarttime
        write(61,*)'                 '
        write(61,*)'---------------------------------------------------------'
        write(61,*)'                 '
        write(61,*)'                 '
    endif
    close(61)    

end subroutine func_Sphere



! subroutine func_Cylinder()

! use variables
! implicit none
! real*4 :: DISTANCE,SDIST
! real*4 :: DIAGONAL
! integer :: l ,m ,n
! real*4 :: xi
! real*4 :: dxg, dyg, dzg
! real*4 :: solid_volume, ETAs
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

    ! !-----------------------calculating casting area----------------------!
    ! Z_castingarea=0.0
    ! do j=1,ny; do i=1,nx
    ! ETAs=0.0
    ! do k=1,nz 
        ! if ( ETA(i,j,k)>ETAs ) then
            ! ETAs=ETA(i,j,k)
        ! endif 
    ! enddo
    ! Z_castingarea=Z_castingarea+ETAs*iDx(i)*iDy(j)
    ! enddo; enddo

    ! Y_castingarea=0.0
    ! do k=1,nz; do i=1,nx
    ! ETAs=0.0
    ! do j=1,ny 
        ! if ( ETA(i,j,k)>ETAs ) then
            ! ETAs=ETA(i,j,k)
        ! endif 
    ! enddo
    ! Y_castingarea=Y_castingarea+ETAs*iDx(i)*iDz(k)
    ! enddo; enddo
	
	! X_castingarea=0.0
    ! do j=1,ny; do k=1,nz
    ! ETAs=0.0
    ! do i=1,nx 
        ! if ( ETA(i,j,k)>ETAs ) then
            ! ETAs=ETA(i,j,k)
        ! endif 
    ! enddo
    ! X_castingarea=X_castingarea+ETAs*iDy(j)*iDz(k)
    ! enddo; enddo

    ! Solid_volume=0.0
    ! do k=1,nz; do j=1,ny; do i=1,nx
        ! solid_volume=solid_volume + ETA(i,j,k)*iDx(i)*iDy(j)*iDz(k)
    ! enddo; enddo; enddo


	! if (select_ref_area == 1) then						! Reference area in virtual force calculation
		! ref_area =  Z_castingarea
    ! elseif (select_ref_area == 2) then		
		! ref_area =  Y_castingarea
    ! elseif (select_ref_area == 3) then		
		! ref_area =  X_castingarea
    ! elseif (select_ref_area == 4) then		
		! ref_area =  user_defined
	! endif


    ! !-----------------------calculating casting area----------------------!
        
    ! totalfinaltime = MPI_WTIME()
    ! if(myid==master)then
        ! write(61,*)'                '
        ! write(61,*)'X_casting area = ',X_castingarea
        ! write(61,*)'Y_casting area = ',Y_castingarea
        ! write(61,*)'Z_casting area = ',Z_castingarea
        ! write(61,*)'Reference area for virtual forces = ', ref_area
		! write(61,*)'                 '
        ! write(61,*)'Casting_volume = ', solid_volume
        ! write(61,*)'                 '
        ! write(61,*)'Conversion cost = ',totalfinaltime-totalstarttime
        ! write(61,*)'                 '
        ! write(61,*)'---------------------------------------------------------'
        ! write(61,*)'                 '
        ! write(61,*)'                 '
    ! endif
    ! close(61)

! end subroutine func_Cylinder






