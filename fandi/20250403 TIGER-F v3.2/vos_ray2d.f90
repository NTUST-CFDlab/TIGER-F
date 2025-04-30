! 27 Oct 2024 - FDS


subroutine read_dat()
    use variables
    implicit none
	integer*4          	:: poly_max=99999999
    integer 			:: n ,code
	character(len=30)  	:: useless
	real*8              :: tun

	
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
        end do    
    close(24)	
    close(61)

    tun = (MAXVAL(iaz) - MINVAL(iaz))/ z_length2D

	!$OMP PARALLEL DO 
	!$acc parallel device_type(host)
	!$acc loop
    do n=1, poly 
        iaz(n) = iaz(n)/tun				! to scale the solid dimensions 
        iay(n) = iay(n)/tun	
    enddo
	!$acc end parallel
    !$OMP END PARALLEL DO

end subroutine read_dat	


subroutine vos_ray2d()
    use variables
    implicit none
    integer :: l ,m ,n , off
	real*8  :: Ybgn, Yend, Zbgn, Zend
    real*8  :: total, segment_d
	real*8  :: inv_subgrid1 = 1.d0/ (nSubGrids_2d*1.d0)
	real*8  :: inv_subgrid2 = 1.d0/ (nSubGrids_2d**2*1.d0)
    real*8,dimension(1:poly)		 :: az, ay
    real*8,dimension(1:nSubGrids_2d) :: SZ, SY
	
    totalstarttime = MPI_WTIME()

!$acc data present(Z,Y,Zs,Ys,iDy,iDz,ETA) copyin(iaz,iay) create(az,ay,points,intersection,ETA_1,sub_intersection)


	! Rotate the blade with AOA
    !$OMP PARALLEL DO
	!$acc parallel loop independent gang vector
    do n = 1, poly
    
        !az(n) = iaz(n) - offset_z               ! iaz(i), iay(i) : point coordinates from GEOMETRY_file
		az(n) = iaz(n) - rotor_r

        az(n)= COS((AOA)*PI/180.d0)*iaz(n) - SIN((AOA)*PI/180.d0)*iay(n)     ! ROTATION 
        ay(n)= SIN((AOA)*PI/180.d0)*iaz(n) + COS((AOA)*PI/180.d0)*iay(n) 

        ay(n) = ay(n) + y0_t
        az(n) = az(n) + z0_t

        !az(n) = az(n) + offset_z
		
    end do
	!$acc end parallel
    !$OMP END PARALLEL DO

!$acc update self(az,ay)
	 
! updating bounds after the solid transformation 

	!$acc kernels
		Ybgn = MINVAL(ay)
		Yend = MAXVAL(ay)
		Zbgn = MINVAL(az)
		Zend = MAXVAL(az)
	!$acc end kernels
	

!.....GET THE POLYGON BORDERS BY COMPARING  MESH POINTS WITHIN DOMAIN
!................ Bounding BOX  ......................
	off = 0
	do j = 1, ny
		if (Y(j) > Ybgn .AND. off == 0) then
			A = j - 2
			off = 1
		elseif (Y(j) > Yend .AND. off == 1) then
			B = j + 1
			exit 
		end if
	end do

	off = 0
	do k = 1, nz
		if (Z(k) > Zbgn .AND. off == 0) then
			C = k - 2
			off = 1
		elseif (Z(k) > Zend .AND. off == 1) then
			E = k + 1
			exit 
		end if
	end do  

    intersection = 0

    !$OMP PARALLEL DO PRIVATE(j)
	!$acc parallel loop independent collapse(2) gang vector
    do k=1,nz ; do j=1,ny
        points(k,j)=0
        intersection(k,j)=0
        ETA_1(k,j)=0
    enddo; enddo
	!$acc end parallel
    !$OMP END PARALLEL DO


    !How many points on edge 
    !$OMP PARALLEL DO PRIVATE(j,k) 
	!$acc parallel loop independent private(j,k) gang vector
    do l =1,poly
    !-------------------------!
        !$acc loop private(k)
		do k=C,E
            if( az(l) >= Z(k) .AND. az(l) < Z(k+1) )then
                !$acc loop
				do j=A,B
                    if( ay(l) >=Y(j) .AND. ay(l) < Y(j+1) )then
                        points(k,j) = 1
                        exit
                    end if
                end do
                exit
            end if
        end do
    !-------------------------!
    end do
	!$acc end parallel
    !$OMP END PARALLEL DO



    !-----------------------------------------------!
    !   Define ETA which point inside the polygon   !
    !-----------------------------------------------!


    !--------------!
	!$OMP PARALLEL DO PRIVATE(j,k)
	!$acc parallel num_gangs(128)
	!$acc loop independent private(j,k) gang
	do l=1,poly-1  !
        if( ((ay(l+1)-ay(l))**2 + (az(l+1)-az(l))**2) < min_dist**2) then
		!$acc loop private(k) vector
        do j=A,B   !
    !--------------!

            if( ((ay(l) >= Ys(j)) .AND. (ay(l+1) < Ys(j))) .OR. &		!	AY(I) > AY(I+1)
				((ay(l) <= Ys(j)) .AND. (ay(l+1) > Ys(j)))	)then      	!	AY(I) < AY(I+1)     ! Xs, Ys, Zs : Midpoints of grid coordinates
                !-------------------------!
                    !$acc loop
					do k=C,E
                        if( Zs(k) <= MAX( az(l),az(l+1) ) )then
							!$acc atomic update
                            intersection(k,j) = intersection(k,j) + 1
                        end if
                    end do
                !-------------------------!
            end if
			
    !---------!
        end do!
        end if
    end do    !
	!$acc end parallel
    !$OMP END PARALLEL DO
	!---------!


    !$OMP PARALLEL DO PRIVATE(j)
	!$acc parallel loop independent collapse(2) gang vector
    do k=C,E;do j=A,B

            if( mod( intersection(k,j),2 ) == 0 )then
                ETA_1(k,j) = 0
            else
                ETA_1(k,j) = 1
            end if
    
    end do;end do
	!$acc end parallel
    !$OMP END PARALLEL DO



    !-----------------------------------------------!
    !              DEFINE THE SUBGRID               !
    !-----------------------------------------------!

!$acc parallel num_gangs(128)
!$acc loop independent private(l,m,n,SY,SZ,sub_intersection,total) collapse(2) gang
    !------------!
     do k=C,E    !
        do j=A,B !
    !------------!

            if( points(k,j) > 0 )then
				!$OMP PARALLEL DO
				!$acc loop vector
                do n = 1,nSubGrids_2d
                    SY(n) = Y(j) + (n-0.5d0) * iDy(j) * inv_subgrid1 !SUB_GRID SIZE IN Y
                    SZ(n) = Z(k) + (n-0.5d0) * iDz(k) * inv_subgrid1 !SUB_GRID SIZE IN Z
                enddo
				!$OMP END PARALLEL DO

                sub_intersection = 0

                    !--------------CALCULATION THE SUBGRID-----------------!

                    !-------------------!
					!$OMP PARALLEL DO PRIVATE(m,n)
					!$acc loop private(m,n) vector
					do l=1,poly-1       !
                        if( ((ay(l+1)-ay(l))**2 + (az(l+1)-az(l))**2) < min_dist**2) then
						!$acc loop private(n)
						do m=1,nSubGrids_2d!
                    !-------------------!

                        if( ((ay(l) >=  SY(m)) .AND. (ay(l+1) <= SY(m))) .OR. &
							((ay(l) <=  SY(m)) .AND. (ay(l+1) >= SY(m)))	)then
                                !-------------------------!
								!$acc loop
								do n=1, nSubGrids_2d
                                    if( SZ(n)<= MAX(az(l),az(l+1)) )then
										!$acc atomic update
										sub_intersection(m,n) = sub_intersection(m,n) + 1
									end if
                                end do
                                !-------------------------!
                        end if


                    !---------!
                        end do!
                        end if
                    end do    !
					!$OMP END PARALLEL DO
                    !--------------CALCULATION THE SUBGRID-----------------!



                total = 0
                
				!$acc loop reduction(+: total) collapse(2) vector
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
				
                ETA_1(k,j) = 1.d0*total * inv_subgrid2

            end if


    !----------!
        end do !
    end do     !
    !----------!
!$acc end parallel



    !$OMP PARALLEL DO COLLAPSE(3)
	!$acc parallel loop independent collapse(3) gang vector
    do k=1,nz ;do j=1,ny; do i=0,nx+1
        ETA(i,j,k) = ETA_1(k,j)    
    end do; end do; end do
	!$acc end parallel
    !$OMP END PARALLEL DO

!$acc end data

    return       
end subroutine vos_ray2d




! subroutine vos_ray2d() !	only for CPU
    ! use variables
    ! implicit none
    ! integer 			:: l ,m ,n
    ! real*8,dimension(1:poly)	:: az, ay
    ! real*8,dimension(nSubGrids_2d+1) :: ZZ, YY
    ! real*8 :: total, segment_d

    ! totalstarttime = MPI_WTIME()

! !$acc data present(Z,Y,Zs,Ys,ETA) copyin(iaz,iay) create(az,ay,points,intersection,ETA_1,sub_intersection)


	! ! Rotate the blade with AOA
    ! !$OMP PARALLEL DO
	! !$acc parallel loop independent gang vector
    ! do n = 1, poly
    
        ! !az(n) = iaz(n) - offset_z               ! iaz(i), iay(i) : point coordinates from GEOMETRY_file
		! az(n) = iaz(n) - rotor_r

        ! ay(n)= SIN((AOA)*PI/180.d0)*az(n) + COS((AOA)*PI/180.d0)*iay(n) 
        ! az(n)= COS((AOA)*PI/180.d0)*az(n) - SIN((AOA)*PI/180.d0)*iay(n)     ! ROTATION

        ! ay(n) = ay(n) + y0
        ! az(n) = az(n) + z0

        ! !az(n) = az(n) + offset_z
		
    ! end do
	! !$acc end parallel
    ! !$OMP END PARALLEL DO

	 
! !$acc update self(az,ay)	 

    ! !GET THE POLYGON BORDERS BY COMPARING  MESH POINTS WITHIN DOMAIN   
    ! do j=1,ny
        ! if( Y(j) > minval(ay) )then
            ! A = j-2
            ! exit
        ! end if
    ! end do

    ! do j=A,ny
        ! if( Y(j) > maxval(ay) )then
            ! B = j+1
            ! exit
        ! end if
    ! end do

    ! do k=1,nz
        ! if( Z(k) > minval(az) )then
            ! C = k-2
            ! exit
        ! end if
    ! end do

    ! do k=C,nz
        ! if( Z(k) > maxval(az) )then
            ! E = k+1
            ! exit
        ! end if
    ! end do    

    ! intersection = 0

    ! !$OMP PARALLEL DO PRIVATE(j)
	! !$acc parallel loop independent collapse(2) gang vector
    ! do k=1,nz ; do j=1,ny
        ! points(k,j)=0
        ! intersection(k,j)=0
        ! ETA_1(k,j)=0
    ! enddo; enddo
	! !$acc end parallel
    ! !$OMP END PARALLEL DO




    ! !How many points on edge 
    ! !$OMP PARALLEL DO PRIVATE(j,k) 
	! !$acc parallel loop independent gang vector
    ! do m =1,poly
    ! !-------------------------!
        ! !$acc loop seq
		! do k=C,E
            ! if( az(m) >= Z(k) .AND. az(m) < Z(k+1) )then
                ! !$acc loop seq
				! do j=A,B
                    ! if( ay(m) >=Y(j) .AND. ay(m) < Y(j+1) )then
                        ! points(k,j) = 1
                        ! exit
                    ! end if
                ! end do
                ! exit
            ! end if
        ! end do
    ! !-------------------------!
    ! end do
	! !$acc end parallel
    ! !$OMP END PARALLEL DO
    

    ! !-----------------------------------------------!
    ! !   Define ETA which point inside the polygon   !
    ! !-----------------------------------------------!
    
    ! !--------------!
    ! do m=1,poly-1  !
        ! if( ((ay(m+1)-ay(m))**2 + (az(m+1)-az(m))**2) < min_dist**2) then
        ! !$OMP PARALLEL DO PRIVATE(k)
		! !$acc parallel loop independent private(k) gang vector		
        ! do j=A,B   !
    ! !--------------!

            ! !---------- AY(I) > AY(I+1)  ----------!
            ! if( (ay(m) >= Ys(j)) .AND. (ay(m+1) < Ys(j)))then            ! Xs, Ys, Zs : Midpoints of grid coordinates
                ! !-------------------------!
                    ! do k=C,E
                        ! if( Zs(k) > max( az(m),az(m+1) ) )then
                            ! intersection(k,j) = intersection(k,j) + 0
                        ! else 
                            ! intersection(k,j) = intersection(k,j) + 1
                        ! end if
                    ! end do
                ! !-------------------------!
            ! end if
            ! !---------- AY(I) > AY(I+1)  ----------!



            ! !---------- AY(I) < AY(I+1)  ----------!
            ! if( (ay(m) <= Ys(j)) .AND. (ay(m+1) > Ys(j)) )then
                ! !-------------------------!
                    ! do k=C,E
                        ! if( Zs(k) > max( az(m),az(m+1) ) )then
                            ! intersection(k,j) = intersection(k,j) + 0
                        ! else
                            ! intersection(k,j) = intersection(k,j) + 1
                        ! end if
                    ! end do
                ! !-------------------------!
            ! end if
            ! !---------- AY(I) < AY(I+1)  ----------!
    ! !---------!
        ! end do!
	! !$acc end parallel
    ! !$OMP END PARALLEL DO
        ! end if
    ! end do    !
    ! !---------!

    ! !$OMP PARALLEL DO PRIVATE(j)
	! !$acc parallel loop independent collapse(2) gang vector
    ! do k=C,E;do j=A,B

            ! if( mod( intersection(k,j),2 ) == 0 )then
                ! ETA_1(k,j) = 0
            ! else
                ! ETA_1(k,j) = 1
            ! end if
    
    ! end do;end do
	! !$acc end parallel
    ! !$OMP END PARALLEL DO


    
    ! !-----------------------------------------------!
    ! !              DEFINE THE SUBGRID               !
    ! !-----------------------------------------------!
    
    ! !------------!
	! !$acc parallel
	! !$acc loop independent private(l,m,n,total,ZZ,YY,sub_intersection) collapse(2) gang vector
     ! do k=C,E     !
        ! do j=A,B !
    ! !------------!


            ! if( points(k,j) > 0 )then
                    ! sub_intersection = 0
                    ! !$OMP PARALLEL DO
					! !$acc loop seq
                    ! do l = 1,nSubGrids_2d+1
                        ! ZZ(l) = Z(k) + (l-1.d0) * ( Z(k+1)-Z(k) ) / DBLE(nSubGrids_2d) !SUB_GRID SIZE IN X
                        ! YY(l) = Y(j) + (l-1.d0) * ( Y(j+1)-Y(j) ) / DBLE(nSubGrids_2d) !SUB_GRID SIZE IN Y
                    ! enddo
                    ! !$OMP END PARALLEL DO

                    ! !--------------CALCULATION THE SUBGRID-----------------!
                    ! !$OMP PARALLEL DO PRIVATE(l,n) 	
                    ! !-------------------!
                    ! !$acc loop seq
					! do m=1,poly-1       !
                        ! if( ((ay(m+1)-ay(m))**2 + (az(m+1)-az(m))**2) < min_dist**2) then
                        ! !$acc loop seq
						! do l=1,nSubGrids_2d!
                    ! !-------------------!

                        ! !---------- AY(I) > AY(I+1)  ----------!
                        ! if( (ay(m) >=  (YY(l)+YY(l+1))/2.d0) .AND. (ay(m+1) <= (YY(l)+YY(l+1))/2.d0) )then
                                ! !-------------------------!
                                ! !$acc loop seq
								! do n=1, nSubGrids_2d
                                    ! if( (ZZ(n)+ZZ(n+1))/2.d0 > max(az(m),az(m+1)) )then
                                        ! sub_intersection(n,l) = sub_intersection(n,l) + 0
                                    ! else
                                        ! sub_intersection(n,l) = sub_intersection(n,l) + 1
                                    ! end if
                                ! end do
                                ! !-------------------------!
                            ! !end if
                        ! end if
                        ! !---------- AY(I) > AY(I+1)  ----------!


                        ! !---------- AY(I) < AY(I+1)  ----------!
                        ! if( (ay(m) <=  (YY(l)+YY(l+1))/2.d0) .AND. (ay(m+1) >= (YY(l)+YY(l+1))/2.d0) )then
                                ! !-------------------------!
                                ! !$acc loop seq
								! do n=1, nSubGrids_2d
                                    ! if( (ZZ(n)+ZZ(n+1))/2.d0 > max(az(m),az(m+1)) )then
                                        ! sub_intersection(n,l) = sub_intersection(n,l) + 0
                                    ! else
                                        ! sub_intersection(n,l) = sub_intersection(n,l) + 1
                                    ! end if
                                ! end do
                                ! !-------------------------!
                            ! !end if
                        ! end if
                        ! !---------- AY(I) < AY(I+1)  ----------!

                    ! !---------!
                        ! end do!
                        ! end if
                    ! end do    !
                    ! !$OMP END PARALLEL DO
                    ! !--------------CALCULATION THE SUBGRID-----------------!

                ! total = 0
                
                ! !$acc loop seq
				! do m=1,nSubGrids_2d
                    ! !$acc loop seq
					! do n=1,nSubGrids_2d

                        ! if( mod(sub_intersection(m,n),2) == 0 )then
                            ! sub_intersection(m,n) = 0
                        ! else
                            ! sub_intersection(m,n) = 1
                        ! end if
                        ! total = total + sub_intersection(m,n)

                    ! end do
                ! end do
                 

                ! ETA_1(k,j) = DBLE(total)/DBLE(nSubGrids_2d)/DBLE(nSubGrids_2d)
        
            ! end if


    ! !----------!
        ! end do !
    ! end do     !
	! !$acc end parallel
    ! !----------!



    ! !$OMP PARALLEL DO PRIVATE(i,j)  
    ! do k=1,nz ;do j=1,ny; do i=0,nx+1
        ! ETA(i,j,k) = ETA_1(k,j)    
    ! end do; end do; end do
    ! !$OMP END PARALLEL DO

! !$acc end data

    ! return       
! end subroutine vos_ray2d