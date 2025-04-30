! 28 Jul 2023 - FDS

! subroutine dynamic_AOA() 
    ! use variables
    ! implicit none

    ! !AOA = AOA1 + AOA_amp * SIN( 2.d0 * reduce_frequency_k * ( time - dt*StartDynamic ) )
    ! AOA = AOA1 + AOA_amp * ( SIN( 2.d0 * reduce_frequency_k * ( time - dt*StartDynamic ) + 1.5d0*PI ) + 1.d0 )      !FROM 0

    ! angular_vel = 2.d0*reduce_frequency_k*AOA_amp*COS( 2.d0 * reduce_frequency_k * ( time - dt*StartDynamic ) + 1.5d0*PI ) / 180.d0 * PI


! end subroutine dynamic_AOA

! subroutine rotor_dynamic() 
    ! use variables 
    ! implicit none

	! rotor_omega = rotor_tsr*U_inf/(rotor_r+blade_r)
	
	! AOA = AOA1 + rotor_omega*( time - StartDynamic_time )*180.d0/PI
	
! end subroutine rotor_dynamic


! subroutine blade_dynamic() 
    ! use variables
    ! implicit none

    ! blade_omega = blade_alpha*U_inf/blade_r

! end subroutine blade_dynamic


subroutine read_dat()
    use variables
    implicit none
	integer*4          	:: poly_max=99999999
    integer 			:: n ,code
	character(len=30)  	:: useless

	
    open(9,file=dat_file)										! calculate file length (numbers of rows)
    do n=1,poly_max
        read(9,*,IOSTAT=code) useless
        if(code<0)then
            poly=n-1
            exit 
        end if
    end do 
    close(9)

		if(myid==master)then
			open (61,file='geometry_info.dat',position='append')
			write(61,*)'Geometry_file = ',dat_file
	        write(61,*)'Number of points in polygon =',poly
			write(61,*)'                '
        endif

    allocate(iaz(poly))
    allocate(iay(poly))

    open( 24,file = dat_file, form = 'FORMATTED' )
        do n = 1, poly
            read(24,*) iaz(n), iay(n)
			iaz(n) = iaz(n)*scale
			iay(n) = iay(n)*scale
        end do    
    close(24)	

    close(61)
	! deallocate(iaz)
	! deallocate(iay)
     
end subroutine read_dat	


 
subroutine vos_ray2d() 
    use variables
    implicit none
    integer 			:: l ,m ,n
    real*8,dimension(1:poly)	:: az, ay
    real*8,dimension(nSubGrids_2d+1) :: ZZ, YY
    real*8 :: total, segment_d

    totalstarttime = MPI_WTIME()


	! Rotate the blade with AOA
    !$OMP PARALLEL DO
    do n = 1, poly
    
        !az(n) = iaz(n) - offset_z               ! iaz(i), iay(i) : point coordinates from GEOMETRY_file
		az(n) = iaz(n) - rotor_r

        ay(n)= SIN((AOA)*PI/180.d0)*az(n) + COS((AOA)*PI/180.d0)*iay(n) 
        az(n)= COS((AOA)*PI/180.d0)*az(n) - SIN((AOA)*PI/180.d0)*iay(n)     ! ROTATION

        ay(n) = ay(n) + yc
        az(n) = az(n) + zc

        !az(n) = az(n) + offset_z
		
    end do
    !$OMP END PARALLEL DO


	! Rotate the blade center with AOA (needed to define the blade rotation)
	! Blades' center is shifted from (0,0) to (-rotor_r, 0)

	
		yc_t= SIN((AOA)*PI/180.d0)*(-rotor_r) + COS((AOA)*PI/180.d0)*0.d0
		zc_t= COS((AOA)*PI/180.d0)*(-rotor_r) - SIN((AOA)*PI/180.d0)*0.d0
		
		yc_t = yc_t + yc
		zc_t = zc_t + zc

!     if(myid==master)then
!	 write(*,*) -AOA, zc_t, yc_t
!	 endif
	 
	 

    !GET THE POLYGON BORDERS BY COMPARING  MESH POINTS WITHIN DOMAIN   
    do j=1,ny
        if( Y(j) > minval(ay) )then
            A = j-2
            exit
        end if
    end do

    do j=A,ny
        if( Y(j) > maxval(ay) )then
            B = j+1
            exit
        end if
    end do

    do k=1,nz
        if( Z(k) > minval(az) )then
            C = k-2
            exit
        end if
    end do

    do k=C,nz
        if( Z(k) > maxval(az) )then
            E = k+1
            exit
        end if
    end do    

    intersection = 0

    !$OMP PARALLEL DO PRIVATE(j)  
    do k=1,nz ; do j=1,ny
        points(k,j)=0
        intersection(k,j)=0
        ETA_1(k,j)=0
    enddo; enddo
    !$OMP END PARALLEL DO




    !How many points on edge 
    !$OMP PARALLEL DO PRIVATE(j,k)  
    do m =1,poly
    !-------------------------!
        do k=C,E
            if( az(m) >= Z(k) .AND. az(m) < Z(k+1) )then
                do j=A,B
                    if( ay(m) >=Y(j) .AND. ay(m) < Y(j+1) )then
                        points(k,j) = 1
                        exit
                    end if
                end do
                exit
            end if
        end do
    !-------------------------!
    end do 
    !$OMP END PARALLEL DO
    

    !-----------------------------------------------!
    !   Define ETA which point inside the polygon   !
    !-----------------------------------------------!
    
    !--------------!
    do m=1,poly-1  !
        if( ((ay(m+1)-ay(m))**2 + (az(m+1)-az(m))**2) < min_dist**2) then
        !$OMP PARALLEL DO PRIVATE(k)  
        do j=A,B   !
    !--------------!

            !---------- AY(I) > AY(I+1)  ----------!
            if( (ay(m) >= Ys(j)) .AND. (ay(m+1) < Ys(j)))then            ! Xs, Ys, Zs : Midpoints of grid coordinates
                !-------------------------!
                    do k=C,E
                        if( Zs(k) > max( az(m),az(m+1) ) )then
                            intersection(k,j) = intersection(k,j) + 0
                        else 
                            intersection(k,j) = intersection(k,j) + 1
                        end if
                    end do
                !-------------------------!
            end if
            !---------- AY(I) > AY(I+1)  ----------!



            !---------- AY(I) < AY(I+1)  ----------!
            if( (ay(m) <= Ys(j)) .AND. (ay(m+1) > Ys(j)) )then
                !-------------------------!
                    do k=C,E
                        if( Zs(k) > max( az(m),az(m+1) ) )then
                            intersection(k,j) = intersection(k,j) + 0
                        else
                            intersection(k,j) = intersection(k,j) + 1
                        end if
                    end do
                !-------------------------!
            end if
            !---------- AY(I) < AY(I+1)  ----------!
    !---------!
        end do!
    !$OMP END PARALLEL DO
        end if
    end do    !
    !---------!

    !$OMP PARALLEL DO PRIVATE(j) 
    do k=C,E;do j=A,B

            if( mod( intersection(k,j),2 ) == 0 )then
                ETA_1(k,j) = 0
            else
                ETA_1(k,j) = 1
            end if
    
    end do;end do    
    !$OMP END PARALLEL DO


    
    !-----------------------------------------------!
    !              DEFINE THE SUBGRID               !
    !-----------------------------------------------!
    
    !------------!
    do k=C,E     !
        do j=A,B !
    !------------!


            if( points(k,j) > 0 )then
                    sub_intersection = 0
                    !$OMP PARALLEL DO
                    do l = 1,nSubGrids_2d+1
                        ZZ(l) = Z(k) + (l-1.d0) * ( Z(k+1)-Z(k) ) / DBLE(nSubGrids_2d) !SUB_GRID SIZE IN X
                        YY(l) = Y(j) + (l-1.d0) * ( Y(j+1)-Y(j) ) / DBLE(nSubGrids_2d) !SUB_GRID SIZE IN Y
                    enddo
                    !$OMP END PARALLEL DO

                    !--------------CALCULATION THE SUBGRID-----------------!
                    !$OMP PARALLEL DO PRIVATE(l,n)  
                    !-------------------!
                    do m=1,poly-1       !
                        if( ((ay(m+1)-ay(m))**2 + (az(m+1)-az(m))**2) < min_dist**2) then
                        do l=1,nSubGrids_2d!
                    !-------------------!

                        !---------- AY(I) > AY(I+1)  ----------!
                        if( (ay(m) >=  (YY(l)+YY(l+1))/2.d0) .AND. (ay(m+1) <= (YY(l)+YY(l+1))/2.d0) )then
                                !-------------------------!
                                do n=1, nSubGrids_2d
                                    if( (ZZ(n)+ZZ(n+1))/2.d0 > max(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 0
                                    else
                                        sub_intersection(n,l) = sub_intersection(n,l) + 1
                                    end if
                                end do
                                !-------------------------!
                            !end if
                        end if
                        !---------- AY(I) > AY(I+1)  ----------!


                        !---------- AY(I) < AY(I+1)  ----------!
                        if( (ay(m) <=  (YY(l)+YY(l+1))/2.d0) .AND. (ay(m+1) >= (YY(l)+YY(l+1))/2.d0) )then
                                !-------------------------!
                                do n=1, nSubGrids_2d
                                    if( (ZZ(n)+ZZ(n+1))/2.d0 > max(az(m),az(m+1)) )then
                                        sub_intersection(n,l) = sub_intersection(n,l) + 0
                                    else
                                        sub_intersection(n,l) = sub_intersection(n,l) + 1
                                    end if
                                end do
                                !-------------------------!
                            !end if
                        end if
                        !---------- AY(I) < AY(I+1)  ----------!

                    !---------!
                        end do!
                        end if
                    end do    !
                    !---------!
                    !$OMP END PARALLEL DO
                    !--------------CALCULATION THE SUBGRID-----------------!

                total = 0
                
                do m=1,nSubGrids_2d
                    do n=1,nSubGrids_2d

                        if( mod(sub_intersection(m,n),2) == 0 )then
                            sub_intersection(m,n) = 0
                        else
                            sub_intersection(m,n) = 1
                        end if
                        total = total + sub_intersection(m,n)

                    end do
                end do
                 

                ETA_1(k,j) = DBLE(total)/DBLE(nSubGrids_2d)/DBLE(nSubGrids_2d)
        
            end if


    !----------!
        end do !
    end do     !
    !----------!



    !$OMP PARALLEL DO PRIVATE(i,j)  
    do k=1,nz ;do j=1,ny; do i=0,nx+1
        ETA(i,j,k) = ETA_1(k,j)    
    end do; end do; end do
    !$OMP END PARALLEL DO



    return       
end subroutine vos_ray2d
