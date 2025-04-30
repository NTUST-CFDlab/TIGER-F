! 20 Nov 2023 - FDS



subroutine twoway_dynamic() 
    use variables 
    implicit none


   !---------------------------------------------------------!
   !      update the bulk solid velocity caused by 2way      !
   !---------------------------------------------------------!

	!---------------For Translation-----------
	! u_solid = u_solid
	! v_solid = v_solid - totalFY_*dt/den_sol/solid_volume
	! w_solid = w_solid
	!---------------For Translation-----------
   

	!---------------For Magnus effect VAWT-----------
	!rotor_omega = rotor_omega - totalTorq*dt/i_sol
	!---------------For Magnus effect VAWT-----------

	!---------------For Darrius VAWT-----------
	rotor_omega = rotor_omega - totalTorq*dt/i_sol
	!---------------For Darrius VAWT-----------

   !---------------------------------------------------------!
   !    update the bulk solid position caused by 2way        !
   !---------------------------------------------------------!

	!---------------------For Translation-----------------------------
	! xc_t = xc_t + u_solid*dt
	! yc_t = yc_t + v_solid*dt
	! zc_t = zc_t + w_solid*dt
      
	! sol_speed = sqrt(u_solid**2 + v_solid**2 + (w_solid - U_inf)**2)
	!---------------------For Translation-----------------------------	
	

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


	!---------------------For Darrius VAWT (3 blades)-----------------------------
	AOA = AOA + rotor_omega*dt*180.d0/PI

	  sol_speed = U_inf
		
	!---------------------For Darrius VAWT (3 blades)-----------------------------

	  
	  Re_t= (sol_speed*L_ch)/nu
	  

	
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






subroutine func_Cylinder()

use variables
implicit none
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

!$OMP PARALLEL DO PRIVATE(j) REDUCTION(max : az_max,ay_max) &
!$OMP			  REDUCTION(min : az_min,ay_min) collapse(2)
do k=1,nz
	do j=1,ny

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
do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
	! Finding potential grids for subgrids
	if( ETA(i,j,k) /= ETA(i,j+1,k-1) .OR. ETA(i,j,k) /= ETA(i,j+1,k  ) .OR. ETA(i,j,k) /= ETA(i,j+1,k+1) .OR. &
		ETA(i,j,k) /= ETA(i,j  ,k-1) .OR.									ETA(i,j,k) /= ETA(i,j  ,k+1) .OR. &
		ETA(i,j,k) /= ETA(i,j-1,k-1) .OR. ETA(i,j,k) /= ETA(i,j-1,k  ) .OR. ETA(i,j,k) /= ETA(i,j-1,k+1) ) then
			do m=1,nSubGrids_f+1; do n=1,nSubGrids_f+1
				dzg = iDz(k) / nSubGrids_f
				dyg = iDy(j) / nSubGrids_f
				SZ(n) = Z(k) + (n-1) * dzg
				SY(m) = Y(j) + (m-1) * dyg		
			end do; end do

		    xi = 0.0
			do m=1,nSubGrids_f; do n=1,nSubGrids_f
				if( ((SZ(n) - zc_t)**2 + (SY(m) - yc_t)**2) .LE. r**2) then
						xi = xi + 1
				endif
			end do; end do
		   
		  
		  ETA(i,j,k) = xi / (nSubGrids_f*nSubGrids_f)
	endif
end do
end do		
!$OMP END PARALLEL DO

do k=kBgnVOS,kEndVOS
do j=jBgnVOS,jEndVOS
do i=iBgnVOS+1,iEndVOS
	ETA(i,j,k) = ETA(iBgnVOS,j,k)

end do
end do
end do    

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
		! open (61,file='geometry_info.dat',position='append')
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






