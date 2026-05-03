! 03 May 2026 by FDS


subroutine virtualForceIntegrator_nima()
	use variables
	implicit none
	real*8 :: FZs,FYs

	totalFZ = 0.d0
	totalFY = 0.d0

!$acc data present(iDx,iDy,iDz,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))

   !$OMP PARALLEL DO      PRIVATE(FZs,FYs) REDUCTION(+: totalFZ, totalFY) COLLAPSE(3)
   !$acc parallel
   !$acc loop independent private(FZs,FYs) reduction(+: totalFZ, totalFY) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx
	
	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))
	
	totalFZ = totalFZ + den_flu*( FZs * iDx(i)*iDy(j)*iDz(k) )
	totalFY = totalFY + den_flu*( FYs * iDx(i)*iDy(j)*iDz(k) )
	
	end do;end do;end do
   !$acc end parallel 
   !$OMP END PARALLEL DO
!$acc end data 

    call MPI_ALLREDUCE( totalFZ, totalFZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY, totalFY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )



    !if(myid==master) then  
		cDrag = (-2.d0) * totalFZ_ / (den_flu*ref_area*U_inf*U_inf*1.d0)  ! because for dimensional: (-2.d0)*totalFX*den_flu / (den_flu*ref_area*U_inf*U_inf) 
		cLift = (-2.d0) * totalFY_ / (den_flu*ref_area*U_inf*U_inf*1.d0)
	!endif

end subroutine virtualForceIntegrator_nima



subroutine virtualForceTorqueIntegrator()
	use variables
	implicit none
	real*8 :: FZs,FYs

	totalFZ = 0.d0
	totalFY = 0.d0
    totalTZ = 0.d0
	totalTY = 0.d0

!$acc data present(iDx,iDy,iDz,Xs,Ys,Zs,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))

	!$OMP PARALLEL	
	!$OMP DO               PRIVATE(FZs,FYs) REDUCTION(+: totalFZ, totalFY) COLLAPSE(3)
	!$acc parallel
	!$acc loop independent private(FZs,FYs) reduction(+: totalFZ, totalFY) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))

		totalFZ = totalFZ + den_flu*( FZs * iDx(i)*iDy(j)*iDz(k) )
		totalFY = totalFY + den_flu*( FYs * iDx(i)*iDy(j)*iDz(k) )

	end do;end do;end do
	!$OMP END DO

    
	! Calculate torque wrt (z0,y0)

    !$OMP DO               PRIVATE(FZs,FYs) REDUCTION(+: totalTZ, totalTY) COLLAPSE(3)
	!$acc loop independent private(FZs,FYs) reduction(+: totalTZ, totalTY) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))

		totalTZ = totalTZ + den_flu*(-FZs*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY = totalTY + den_flu*( FYs*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0

    end do;end do;end do
    !$OMP END DO
	!$OMP END PARALLEL
   !$acc end parallel 
!$acc end data 
	

    call MPI_ALLREDUCE( totalFZ, totalFZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY, totalFY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ, totalTZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY, totalTY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )	

    !if(myid==master) then 	
		totalTorq =   (totalTZ_ + totalTY_)
	
		cDrag = -2.d0 * totalFZ_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift = -2.d0 * totalFY_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq = -2.d0 * totalTorq/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower   =   cTorq*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)
	!endif
	

end subroutine virtualForceTorqueIntegrator


subroutine virtualForce_tandemVIV()
	use variables
	implicit none
	real*8 :: FZs,FYs

	totalFZ1 = 0.d0 ; totalFZ2 = 0.d0
	totalFY1 = 0.d0 ; totalFY2 = 0.d0

!$acc data present(iDx,iDy,iDz,Zs,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))

   !$OMP PARALLEL DO      PRIVATE(FZs,FYs) REDUCTION(+: totalFZ1, totalFY1, totalFZ2, totalFY2) COLLAPSE(3)
   !$acc parallel
   !$acc loop independent private(FZs,FYs) reduction(+: totalFZ1, totalFY1, totalFZ2, totalFY2) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))

		 if(Zs(k) .LT. Zsplit_pos) then
			totalFZ1 = totalFZ1 + den_flu*( FZs * iDx(i)*iDy(j)*iDz(k) )
			totalFY1 = totalFY1 + den_flu*( FYs * iDx(i)*iDy(j)*iDz(k) )
		 else
			totalFZ2 = totalFZ2 + den_flu*( FZs * iDx(i)*iDy(j)*iDz(k) )
			totalFY2 = totalFY2 + den_flu*( FYs * iDx(i)*iDy(j)*iDz(k) )
		 endif
	
	end do;end do;end do
   !$acc end parallel 
   !$OMP END PARALLEL DO
!$acc end data 

    call MPI_ALLREDUCE( totalFZ1, totalFZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY1, totalFY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

    call MPI_ALLREDUCE( totalFZ2, totalFZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY2, totalFY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

  
		cDrag1 = -2.d0 * totalFZ1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift1 = -2.d0 * totalFY1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cPower1   =   cLift1 * v_solid_1 / U_inf

		cDrag2 = -2.d0 * totalFZ2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift2 = -2.d0 * totalFY2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cPower2   =   cLift2 * v_solid_2 / U_inf		
		

end subroutine virtualForce_tandemVIV


subroutine virtualForceIntegrator_shear()			! only for cylinders
	use variables
	implicit none

	real*8	:: rad, theta_tan
	real*8  :: FZs,FYs
	real*8	:: fric_coeff_avg,moment_dist_avg, tnew, told
	integer	:: n
	integer, dimension(0:360)	:: near_Zn,near_Yn
	real*8	:: d_area,CF_sum,v_avg,w_avg,nut_avg
	real*8, dimension(0:360) 	:: Zn,Yn
	real*8, dimension(1:2,1:2) 	:: v_s,w_s,nut_s

	totalFZ = 0.d0
	totalFY = 0.d0

			!>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
			!$acc update self(w(:,:,iend),v(:,:,iend),nut(:,:,iend)) if(nproc>1)
			icount = (nx+4)*(ny+4)
			itag = 601
			call MPI_SENDRECV( w(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   w(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
			itag = 602
			call MPI_SENDRECV( v(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   v(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )      
			itag = 603
			call MPI_SENDRECV( nut(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   nut(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )

			call MPI_barrier(MPI_COMM_WORLD, ierr)
			!$acc update device(w(:,:,istart-1),v(:,:,istart-1),nut(:,:,istart-1)) if(nproc>1 .AND. myid/=0)     
			!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!$acc data present(Xs,Ys,Zs,iDx,iDy,iDz,u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2), &
!$acc			nut(:,:,istart:iend+1),FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2)) &
!$acc 	   create(Zn,Yn,near_Zn,near_Yn,moment_dist(:,:),fric_coeff(:,:))

   !$OMP PARALLEL DO	           PRIVATE(FZs,FYs) reduction(+: totalFZ, totalFY) COLLAPSE(3)
   !$acc parallel loop independent private(FZs,FYs) reduction(+: totalFZ, totalFY) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))

	totalFZ = totalFZ + den_flu*( FZs * iDx(i)*iDy(j)*iDz(k) )
	totalFY = totalFY + den_flu*( FYs * iDx(i)*iDy(j)*iDz(k) )
	
	end do;end do;end do
   !$acc end parallel 
   !$OMP END PARALLEL DO


	
	! Finding sample points
   !$OMP PARALLEL DO
   !$acc parallel loop independent gang vector
	do n = 0,360
		Zn(n) = z0 - (blade_r+2.d0*dzSml)*COS(PI*n/180.D0)
		Yn(n) = y0 - (blade_r+2.d0*dySml)*SIN(PI*n/180.D0)
	enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO

	d_area = 2.d0*blade_r*PI*dx/360.d0

	! Finding near points
   !$OMP PARALLEL DO			   COLLAPSE(2)
   !$acc parallel loop independent collapse(2) gang vector
	do n = 0,360
		do k=kBgnVOS-5,kEndVOS+5
			if(Zn(n) .GE. Zs(k) .AND. Zn(n) .LT. Zs(k+1)) near_Zn(n) = k
	enddo;enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO			   COLLAPSE(2)
   !$acc parallel loop independent collapse(2) gang vector
	do n = 0,360
		do j=jBgnVOS-5,jEndVOS+5
			if(Yn(n) .GE. Ys(j) .AND. Yn(n) .LT. Ys(j+1)) near_Yn(n) = j
	enddo;enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO

	fric_coeff = 0.d0
	moment_dist = 0.d0
	tau_moment = 0.d0

   tnew = 1.d0/(stepavg*1.d0)
   told = (stepavg-1)*1.d0

	! calculating wall parameters
   !$OMP PARALLEL DO	           PRIVATE(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) REDUCTION(+: tau_moment) COLLAPSE(2)
   !$acc parallel loop independent private(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) reduction(+: tau_moment) collapse(2) gang vector
	do n = 0,360
		do i=1,nx
			if(near_Zn(n) .GE. istart .AND. near_Zn(n) .LE. iend) then
				do k=1,2 ; do j=1,2
					v_s(j,k) = 0.5d0*( v(i,near_Yn(n)  +(j-1),near_Zn(n)  +(k-1))+ &
									   v(i,near_Yn(n)-1+(j-1),near_Zn(n)  +(k-1)) )
					w_s(j,k) = 0.5d0*( w(i,near_Yn(n)  +(j-1),near_Zn(n)  +(k-1))+ &
									   w(i,near_Yn(n)  +(j-1),near_Zn(n)-1+(k-1)) )
					nut_s(j,k) =     nut(i,near_Yn(n)  +(j-1),near_Zn(n)+  (k-1))
				enddo ; enddo
				
				v_avg = bilinear_interp(Yn(n), Zn(n), Ys(near_Yn(n)), Ys(near_Yn(n)+1), Zs(near_Zn(n)), Zs(near_Zn(n)+1), v_s)
				w_avg = bilinear_interp(Yn(n), Zn(n), Ys(near_Yn(n)), Ys(near_Yn(n)+1), Zs(near_Zn(n)), Zs(near_Zn(n)+1), w_s)
			  nut_avg = bilinear_interp(Yn(n), Zn(n), Ys(near_Yn(n)), Ys(near_Yn(n)+1), Zs(near_Zn(n)), Zs(near_Zn(n)+1), nut_s)

			!theta_tan = ATAN(-1.d0*(Zn(n) - z0)/(Yn(n) - y0))
			
			theta_tan = ATAN2( -(Zn(n)-z0), (Yn(n)-y0) )
			IF (theta_tan >  0.5d0*PI)  theta_tan = theta_tan - PI
			IF (theta_tan < -0.5d0*PI)  theta_tan = theta_tan + PI
			
			tau_moment 	= tau_moment + (nu + nut_avg) * den_flu * d_area * blade_r * SIGN(1.d0,(Yn(n) - y0)) * &
								(w_avg*COS(theta_tan) + v_avg*SIN(theta_tan) &
								+ blade_omega * blade_r * SIGN(1.d0,(Yn(n) - y0))) / (2.d0*dzSml)

			! moment_dist(n,i) 	= (nu + nut_avg) * den_flu * d_area * blade_r * SIGN(1.d0,(Yn(n) - y0)) * &
								! (w_avg*COS(theta_tan) + v_avg*SIN(theta_tan) &
								! + blade_omega * blade_r * SIGN(1.d0,(Yn(n) - y0))) / (1.5d0*dzSml)

			! fric_coeff(n,i)  	= (nu + nut_avg) * den_flu * (w_avg*COS(theta_tan) + v_avg*SIN(theta_tan) &	!  Cf
							! + blade_omega * blade_r * SIGN(1.d0,(Yn(n) - y0))) * SIGN(1.d0,(Yn(n) - y0))*2.d0/ (1.5d0*dzSml)
			! CF_TimeAvg(n,i) = (CF_TimeAvg(n,i) * told + fric_coeff(n,i)) * tnew
			endif
		enddo   	
	enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO	

	! stepavg = stepavg+1

	! if (time >= 310.d0 .AND. time < 310.d0+dt) then 
		! !$acc update self(moment_dist(:,:),CF_TimeAvg(:,:))
		! do n = 0,360
			! fric_coeff_avg=0.d0
			! moment_dist_avg=0.d0
				! do i=1,nx
					! fric_coeff_avg=fric_coeff_avg+CF_TimeAvg(n,i)
					! moment_dist_avg=moment_dist_avg+moment_dist(n,i)
				! enddo
				! fric_coeff_avg=fric_coeff_avg/(nx*1.d0)	
				! moment_dist_avg=moment_dist_avg/(nx*1.d0)
			! write(82,'(I5,3X,G18.10E3,3X,G18.10E3)') n,fric_coeff_avg,moment_dist_avg
			! flush(82)
		! enddo
	! endif


	
!$acc end data 


    call MPI_ALLREDUCE( totalFZ, totalFZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY, totalFY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( tau_moment, tau_moment_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )


    !if(myid==master) then  
		cDrag = (-2.d0) * totalFZ_ / (den_flu*ref_area*U_inf*U_inf*1.d0)  ! because for dimensional: (-2.d0)*totalFX*den_flu / (den_flu*ref_area*U_inf*U_inf) 
		cLift = (-2.d0) * totalFY_ / (den_flu*ref_area*U_inf*U_inf*1.d0)
		cTorq = -2.d0 * tau_moment_/(2.d0*blade_r*lx*U_inf*U_inf*den_flu*blade_r) 
	!endif

end subroutine virtualForceIntegrator_shear



subroutine virtualForce_Magnus2blade
	use variables
	implicit none

	real*8	:: rad, theta_tan
	real*8  :: FZs,FYs
	integer	:: n
	integer, dimension(0:360)	:: near_Zn1,near_Yn1,near_Zn2,near_Yn2
	real*8	:: d_area,v_avg,w_avg,nut_avg
	real*8, dimension(0:360) 	:: Zn1,Yn1,Zn2,Yn2
	real*8, dimension(1:2,1:2) 	:: v_s,w_s,nut_s

	totalFZ1 = 0.d0
	totalFY1 = 0.d0
    totalTZ1 = 0.d0
	totalTY1 = 0.d0

	totalFZ2 = 0.d0
	totalFY2 = 0.d0
    totalTZ2 = 0.d0
	totalTY2 = 0.d0


			!>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
			!$acc update self(w(:,:,iend),v(:,:,iend),nut(:,:,iend)) if(nproc>1)
			icount = (nx+4)*(ny+4)
			itag = 601
			call MPI_SENDRECV( w(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   w(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
			itag = 602
			call MPI_SENDRECV( v(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   v(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )      
			itag = 603
			call MPI_SENDRECV( nut(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   nut(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )

			call MPI_barrier(MPI_COMM_WORLD, ierr)
			!$acc update device(w(:,:,istart-1),v(:,:,istart-1),nut(:,:,istart-1)) if(nproc>1 .AND. myid/=0)     
			!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!$acc data present(Xs,Ys,Zs,iDx,iDy,iDz,u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2), &
!$acc			nut(:,:,istart:iend+1),FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2)) &
!$acc 	   create(Zn1,Yn1,Zn2,Yn2,near_Zn1,near_Yn1,near_Zn2,near_Yn2)

	!$OMP PARALLEL	
	!$OMP DO               PRIVATE(FZs,FYs) REDUCTION(+: totalFZ1, totalFY1, totalFZ2, totalFY2) COLLAPSE(3)
	!$acc parallel
	!$acc loop independent private(FZs,FYs) reduction(+: totalFZ1, totalFY1, totalFZ2, totalFY2) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))

	   if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (3.d0*blade_r)) then
		totalFZ1 = totalFZ1 + den_flu*( FZs * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY1 = totalFY1 + den_flu*( FYs * ( iDx(i)*iDy(j)*iDz(k) ) )	
	   else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (3.d0*blade_r)) then
		totalFZ2 = totalFZ2 + den_flu*( FZs * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY2 = totalFY2 + den_flu*( FYs * ( iDx(i)*iDy(j)*iDz(k) ) )
	   endif

	end do;end do;end do
	!$OMP END DO

    
	! Calculate torque wrt (z0,y0)

    !$OMP DO               PRIVATE(FZs,FYs) REDUCTION(+: totalTZ1, totalTY1, totalTZ2, totalTY2) COLLAPSE(3)
	!$acc loop independent private(FZs,FYs) reduction(+: totalTZ1, totalTY1, totalTZ2, totalTY2) collapse(3) gang vector
	do k=istart,iend
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))

	   if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (3.d0*blade_r)) then
		totalTZ1 = totalTZ1 + den_flu*(-FZs*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY1 = totalTY1 + den_flu*( FYs*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0		

	   else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (3.d0*blade_r)) then
		totalTZ2 = totalTZ2 + den_flu*(-FZs*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY2 = totalTY2 + den_flu*( FYs*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0

	   endif
	  
    end do;end do;end do
   !$acc end parallel 
   !$OMP END DO
   !$OMP END PARALLEL


	! Calculate input torque wrt each cylinder center
	
	! Finding sample points
   !$OMP PARALLEL DO
   !$acc parallel loop independent gang vector
	do n = 0,360
		Zn1(n) = z0_t1 - (blade_r+2.d0*dzSml)*COS(PI*n/180.D0)
		Yn1(n) = y0_t1 - (blade_r+2.d0*dySml)*SIN(PI*n/180.D0)

		Zn2(n) = z0_t2 - (blade_r+2.d0*dzSml)*COS(PI*n/180.D0)
		Yn2(n) = y0_t2 - (blade_r+2.d0*dySml)*SIN(PI*n/180.D0)
	enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO


	d_area = 2.d0*blade_r*PI*dx/360.d0

	! Finding near points
   !$OMP PARALLEL DO			   COLLAPSE(2)
   !$acc parallel loop independent collapse(2) gang vector
	do n = 0,360
		do k=kBgnVOS-5,kEndVOS+5
			if(Zn1(n) .GE. Zs(k) .AND. Zn1(n) .LT. Zs(k+1)) near_Zn1(n) = k
			if(Zn2(n) .GE. Zs(k) .AND. Zn2(n) .LT. Zs(k+1)) near_Zn2(n) = k
	enddo;enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO			   COLLAPSE(2)
   !$acc parallel loop independent collapse(2) gang vector
	do n = 0,360
		do j=jBgnVOS-5,jEndVOS+5
			if(Yn1(n) .GE. Ys(j) .AND. Yn1(n) .LT. Ys(j+1)) near_Yn1(n) = j
			if(Yn2(n) .GE. Ys(j) .AND. Yn2(n) .LT. Ys(j+1)) near_Yn2(n) = j
	enddo;enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO


	tau_moment1 = 0.d0
	tau_moment2 = 0.d0


	! calculating wall parameters of blade 1
   !$OMP PARALLEL DO	           PRIVATE(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) REDUCTION(+: tau_moment1) COLLAPSE(2)
   !$acc parallel loop independent private(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) reduction(+: tau_moment1) collapse(2) gang vector
	do n = 0,360
		do i=1,nx
			if(near_Zn1(n) .GE. istart .AND. near_Zn1(n) .LE. iend) then
				do k=1,2 ; do j=1,2
					v_s(j,k) = 0.5d0*( v(i,near_Yn1(n)  +(j-1),near_Zn1(n)  +(k-1))+ &
									   v(i,near_Yn1(n)-1+(j-1),near_Zn1(n)  +(k-1)) )
					w_s(j,k) = 0.5d0*( w(i,near_Yn1(n)  +(j-1),near_Zn1(n)  +(k-1))+ &
									   w(i,near_Yn1(n)  +(j-1),near_Zn1(n)-1+(k-1)) )
					nut_s(j,k) =     nut(i,near_Yn1(n)  +(j-1),near_Zn1(n)  +(k-1))
				enddo ; enddo
				
				v_avg = bilinear_interp(Yn1(n), Zn1(n), Ys(near_Yn1(n)), Ys(near_Yn1(n)+1), Zs(near_Zn1(n)), Zs(near_Zn1(n)+1), v_s)
				w_avg = bilinear_interp(Yn1(n), Zn1(n), Ys(near_Yn1(n)), Ys(near_Yn1(n)+1), Zs(near_Zn1(n)), Zs(near_Zn1(n)+1), w_s)
			  nut_avg = bilinear_interp(Yn1(n), Zn1(n), Ys(near_Yn1(n)), Ys(near_Yn1(n)+1), Zs(near_Zn1(n)), Zs(near_Zn1(n)+1), nut_s)
			
			theta_tan = ATAN2( -(Zn1(n)-z0_t1), (Yn1(n)-y0_t1) )
			IF (theta_tan >  0.5d0*PI)  theta_tan = theta_tan - PI
			IF (theta_tan < -0.5d0*PI)  theta_tan = theta_tan + PI

			tau_moment1 	= tau_moment1 + (nu + nut_avg) * den_flu * d_area * blade_r * SIGN(1.d0,(Yn1(n) - y0_t1)) * &
								(w_avg*COS(theta_tan) + v_avg*SIN(theta_tan) &
								+ blade_omega * blade_r * SIGN(1.d0,(Yn1(n) - y0_t1))) / (2.d0*dzSml)
			endif
		enddo   	
	enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO	


	! calculating wall parameters of blade 2
   !$OMP PARALLEL DO	           PRIVATE(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) REDUCTION(+: tau_moment2) COLLAPSE(2)
   !$acc parallel loop independent private(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) reduction(+: tau_moment2) collapse(2) gang vector
	do n = 0,360
		do i=1,nx
			if(near_Zn2(n) .GE. istart .AND. near_Zn2(n) .LE. iend) then
				do k=1,2 ; do j=1,2
					v_s(j,k) = 0.5d0*( v(i,near_Yn2(n)  +(j-1),near_Zn2(n)  +(k-1))+ &
									   v(i,near_Yn2(n)-1+(j-1),near_Zn2(n)  +(k-1)) )
					w_s(j,k) = 0.5d0*( w(i,near_Yn2(n)  +(j-1),near_Zn2(n)  +(k-1))+ &
									   w(i,near_Yn2(n)  +(j-1),near_Zn2(n)-1+(k-1)) )
					nut_s(j,k) =     nut(i,near_Yn2(n)  +(j-1),near_Zn2(n)  +(k-1))
				enddo ; enddo
				
				v_avg = bilinear_interp(Yn2(n), Zn2(n), Ys(near_Yn2(n)), Ys(near_Yn2(n)+1), Zs(near_Zn2(n)), Zs(near_Zn2(n)+1), v_s)
				w_avg = bilinear_interp(Yn2(n), Zn2(n), Ys(near_Yn2(n)), Ys(near_Yn2(n)+1), Zs(near_Zn2(n)), Zs(near_Zn2(n)+1), w_s)
			  nut_avg = bilinear_interp(Yn2(n), Zn2(n), Ys(near_Yn2(n)), Ys(near_Yn2(n)+1), Zs(near_Zn2(n)), Zs(near_Zn2(n)+1), nut_s)

			theta_tan = ATAN2( -(Zn2(n)-z0_t2), (Yn2(n)-y0_t2) )
			IF (theta_tan >  0.5d0*PI)  theta_tan = theta_tan - PI
			IF (theta_tan < -0.5d0*PI)  theta_tan = theta_tan + PI

			tau_moment2 	= tau_moment2 + (nu + nut_avg) * den_flu * d_area * blade_r * SIGN(1.d0,(Yn2(n) - y0_t2)) * &
								(w_avg*COS(theta_tan) + v_avg*SIN(theta_tan) &
								+ blade_omega * blade_r * SIGN(1.d0,(Yn2(n) - y0_t2))) / (2.d0*dzSml)
			endif
		enddo   	
	enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO	

	
!$acc end data 


    call MPI_ALLREDUCE( totalFZ1, totalFZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY1, totalFY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ1, totalTZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY1, totalTY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( tau_moment1, tau_moment1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )	

    call MPI_ALLREDUCE( totalFZ2, totalFZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY2, totalFY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ2, totalTZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY2, totalTY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( tau_moment2, tau_moment2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

		totalFZ_ = totalFZ1_ + totalFZ2_ 
		totalFY_ = totalFY1_ + totalFY2_


		totalTorq1 =   (totalTZ1_ + totalTY1_)
		totalTorq2 =   (totalTZ2_ + totalTY2_)

		totalTorq  = totalTorq1 + totalTorq2


		cDrag = -2.d0 * totalFZ_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift = -2.d0 * totalFY_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq = -2.d0 * totalTorq/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower   =   cTorq*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cPower_in   = -2.d0*blade_omega*(tau_moment1_  + tau_moment2_ ) / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
		
		cPower_nett = cPower + cPower_in


		cDrag1 = -2.d0 * totalFZ1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift1 = -2.d0 * totalFY1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq1 = -2.d0 * totalTorq1/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower1   =   cTorq1*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cPower1_in   = -2.d0*tau_moment1_*blade_omega / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
		
		cDrag2 = -2.d0 * totalFZ2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift2 = -2.d0 * totalFY2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq2 = -2.d0 * totalTorq2/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower2   =   cTorq2*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cPower2_in   = -2.d0*tau_moment2_*blade_omega / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
	

end subroutine virtualForce_Magnus2blade



subroutine virtualForce_Magnus3blade()
	use variables
	implicit none

	real*8	:: rad, theta_tan
	real*8  :: FZs,FYs
	integer	:: n
	integer, dimension(0:360)	:: near_Zn1,near_Yn1,near_Zn2,near_Yn2,near_Zn3,near_Yn3
	real*8	:: d_area,v_avg,w_avg,nut_avg
	real*8, dimension(0:360) 	:: Zn1,Yn1,Zn2,Yn2,Zn3,Yn3
	real*8, dimension(1:2,1:2) 	:: v_s,w_s,nut_s
	
	totalFZ1 = 0.d0
	totalFY1 = 0.d0
    totalTZ1 = 0.d0
	totalTY1 = 0.d0
	
	totalFZ2 = 0.d0
	totalFY2 = 0.d0
    totalTZ2 = 0.d0
	totalTY2 = 0.d0

	totalFZ3 = 0.d0
	totalFY3 = 0.d0
    totalTZ3 = 0.d0
	totalTY3 = 0.d0

			!>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
			!$acc update self(w(:,:,iend),v(:,:,iend),nut(:,:,iend)) if(nproc>1)
			icount = (nx+4)*(ny+4)
			itag = 601
			call MPI_SENDRECV( w(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   w(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
			itag = 602
			call MPI_SENDRECV( v(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   v(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )      
			itag = 603
			call MPI_SENDRECV( nut(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   nut(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )

			call MPI_barrier(MPI_COMM_WORLD, ierr)
			!$acc update device(w(:,:,istart-1),v(:,:,istart-1),nut(:,:,istart-1)) if(nproc>1 .AND. myid/=0)     
			!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!$acc data present(Xs,Ys,Zs,iDx,iDy,iDz,u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2), &
!$acc			nut(:,:,istart:iend+1),FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2)) &
!$acc 	   create(Zn1,Yn1,Zn2,Yn2,Zn3,Yn3,near_Zn1,near_Yn1,near_Zn2,near_Yn2,near_Zn3,near_Yn3)

	!$OMP PARALLEL	
	!$OMP DO               PRIVATE(FZs,FYs) REDUCTION(+: totalFZ1, totalFY1, totalFZ2, totalFY2, totalFZ3, totalFY3) COLLAPSE(3)
	!$acc parallel
	!$acc loop independent private(FZs,FYs) reduction(+: totalFZ1, totalFY1, totalFZ2, totalFY2, totalFZ3, totalFY3) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))

	   if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (3.d0*blade_r)) then
		totalFZ1 = totalFZ1 + den_flu*( FZs * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY1 = totalFY1 + den_flu*( FYs * ( iDx(i)*iDy(j)*iDz(k) ) )	
	   else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (3.d0*blade_r)) then
		totalFZ2 = totalFZ2 + den_flu*( FZs * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY2 = totalFY2 + den_flu*( FYs * ( iDx(i)*iDy(j)*iDz(k) ) )
	   else if( SQRT((Zs(k)-z0_t3)**2 + (Ys(j)-y0_t3)**2) .LT. (3.d0*blade_r)) then
		totalFZ3 = totalFZ3 + den_flu*( FZs * ( iDx(i)*iDy(j)*iDz(k) ) )
		totalFY3 = totalFY3 + den_flu*( FYs * ( iDx(i)*iDy(j)*iDz(k) ) )
	   endif

	end do;end do;end do
	!$OMP END DO

    
	! Calculate torque wrt (z0,y0)

    !$OMP DO               PRIVATE(FZs,FYs) REDUCTION(+: totalTZ1, totalTY1, totalTZ2, totalTY2, totalTZ3, totalTY3) COLLAPSE(3)
	!$acc loop independent private(FZs,FYs) reduction(+: totalTZ1, totalTY1, totalTZ2, totalTY2, totalTZ3, totalTY3) collapse(3) gang vector
	do k=istart,iend
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))

	   if( SQRT((Zs(k)-z0_t1)**2 + (Ys(j)-y0_t1)**2) .LT. (3.d0*blade_r)) then
		totalTZ1 = totalTZ1 + den_flu*(-FZs*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY1 = totalTY1 + den_flu*( FYs*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0
	   else if( SQRT((Zs(k)-z0_t2)**2 + (Ys(j)-y0_t2)**2) .LT. (3.d0*blade_r)) then
		totalTZ2 = totalTZ2 + den_flu*(-FZs*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY2 = totalTY2 + den_flu*( FYs*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0
	   else if( SQRT((Zs(k)-z0_t3)**2 + (Ys(j)-y0_t3)**2) .LT. (3.d0*blade_r)) then
		totalTZ3 = totalTZ3 + den_flu*(-FZs*(Ys(j)-y0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
		totalTY3 = totalTY3 + den_flu*( FYs*(Zs(k)-z0 )*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0
	   endif

    end do;end do;end do
    !$acc end parallel 
    !$OMP END DO
	!$OMP END PARALLEL



	! Calculate input torque wrt each cylinder center
	
	! Finding sample points
   !$OMP PARALLEL DO
   !$acc parallel loop independent gang vector
	do n = 0,360
		Zn1(n) = z0_t1 - (blade_r+2.d0*dzSml)*COS(PI*n/180.D0)
		Yn1(n) = y0_t1 - (blade_r+2.d0*dySml)*SIN(PI*n/180.D0)

		Zn2(n) = z0_t2 - (blade_r+2.d0*dzSml)*COS(PI*n/180.D0)
		Yn2(n) = y0_t2 - (blade_r+2.d0*dySml)*SIN(PI*n/180.D0)

		Zn3(n) = z0_t3 - (blade_r+2.d0*dzSml)*COS(PI*n/180.D0)
		Yn3(n) = y0_t3 - (blade_r+2.d0*dySml)*SIN(PI*n/180.D0)
	enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO


	d_area = 2.d0*blade_r*PI*dx/360.d0

	! Finding near points
   !$OMP PARALLEL DO			   COLLAPSE(2)
   !$acc parallel loop independent collapse(2) gang vector
	do n = 0,360
		do k=kBgnVOS-5,kEndVOS+5
			if(Zn1(n) .GE. Zs(k) .AND. Zn1(n) .LT. Zs(k+1)) near_Zn1(n) = k
			if(Zn2(n) .GE. Zs(k) .AND. Zn2(n) .LT. Zs(k+1)) near_Zn2(n) = k
			if(Zn3(n) .GE. Zs(k) .AND. Zn3(n) .LT. Zs(k+1)) near_Zn3(n) = k
	enddo;enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO			   COLLAPSE(2)
   !$acc parallel loop independent collapse(2) gang vector
	do n = 0,360
		do j=jBgnVOS-5,jEndVOS+5
			if(Yn1(n) .GE. Ys(j) .AND. Yn1(n) .LT. Ys(j+1)) near_Yn1(n) = j
			if(Yn2(n) .GE. Ys(j) .AND. Yn2(n) .LT. Ys(j+1)) near_Yn2(n) = j
			if(Yn3(n) .GE. Ys(j) .AND. Yn3(n) .LT. Ys(j+1)) near_Yn3(n) = j
	enddo;enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO


	tau_moment1 = 0.d0
	tau_moment2 = 0.d0
	tau_moment3 = 0.d0

	! calculating wall parameters of blade 1
   !$OMP PARALLEL DO	           PRIVATE(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) REDUCTION(+: tau_moment1) COLLAPSE(2)
   !$acc parallel loop independent private(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) reduction(+: tau_moment1) collapse(2) gang vector
	do n = 0,360
		do i=1,nx
			if(near_Zn1(n) .GE. istart .AND. near_Zn1(n) .LE. iend) then
				do k=1,2 ; do j=1,2
					v_s(j,k) = 0.5d0*( v(i,near_Yn1(n)  +(j-1),near_Zn1(n)  +(k-1))+ &
									   v(i,near_Yn1(n)-1+(j-1),near_Zn1(n)  +(k-1)) )
					w_s(j,k) = 0.5d0*( w(i,near_Yn1(n)  +(j-1),near_Zn1(n)  +(k-1))+ &
									   w(i,near_Yn1(n)  +(j-1),near_Zn1(n)-1+(k-1)) )
					nut_s(j,k) =     nut(i,near_Yn1(n)  +(j-1),near_Zn1(n)  +(k-1))
				enddo ; enddo
				
				v_avg = bilinear_interp(Yn1(n), Zn1(n), Ys(near_Yn1(n)), Ys(near_Yn1(n)+1), Zs(near_Zn1(n)), Zs(near_Zn1(n)+1), v_s)
				w_avg = bilinear_interp(Yn1(n), Zn1(n), Ys(near_Yn1(n)), Ys(near_Yn1(n)+1), Zs(near_Zn1(n)), Zs(near_Zn1(n)+1), w_s)
			  nut_avg = bilinear_interp(Yn1(n), Zn1(n), Ys(near_Yn1(n)), Ys(near_Yn1(n)+1), Zs(near_Zn1(n)), Zs(near_Zn1(n)+1), nut_s)
			
			theta_tan = ATAN2( -(Zn1(n)-z0_t1), (Yn1(n)-y0_t1) )
			IF (theta_tan >  0.5d0*PI)  theta_tan = theta_tan - PI
			IF (theta_tan < -0.5d0*PI)  theta_tan = theta_tan + PI

			tau_moment1 	= tau_moment1 + (nu + nut_avg) * den_flu * d_area * blade_r * SIGN(1.d0,(Yn1(n) - y0_t1)) * &
								(w_avg*COS(theta_tan) + v_avg*SIN(theta_tan) &
								+ blade_omega * blade_r * SIGN(1.d0,(Yn1(n) - y0_t1))) / (2.d0*dzSml)
			endif
		enddo   	
	enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO	


	! calculating wall parameters of blade 2
   !$OMP PARALLEL DO	           PRIVATE(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) REDUCTION(+: tau_moment2) COLLAPSE(2)
   !$acc parallel loop independent private(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) reduction(+: tau_moment2) collapse(2) gang vector
	do n = 0,360
		do i=1,nx
			if(near_Zn2(n) .GE. istart .AND. near_Zn2(n) .LE. iend) then
				do k=1,2 ; do j=1,2
					v_s(j,k) = 0.5d0*( v(i,near_Yn2(n)  +(j-1),near_Zn2(n)  +(k-1))+ &
									   v(i,near_Yn2(n)-1+(j-1),near_Zn2(n)  +(k-1)) )
					w_s(j,k) = 0.5d0*( w(i,near_Yn2(n)  +(j-1),near_Zn2(n)  +(k-1))+ &
									   w(i,near_Yn2(n)  +(j-1),near_Zn2(n)-1+(k-1)) )
					nut_s(j,k) =     nut(i,near_Yn2(n)  +(j-1),near_Zn2(n)  +(k-1))
				enddo ; enddo
				
				v_avg = bilinear_interp(Yn2(n), Zn2(n), Ys(near_Yn2(n)), Ys(near_Yn2(n)+1), Zs(near_Zn2(n)), Zs(near_Zn2(n)+1), v_s)
				w_avg = bilinear_interp(Yn2(n), Zn2(n), Ys(near_Yn2(n)), Ys(near_Yn2(n)+1), Zs(near_Zn2(n)), Zs(near_Zn2(n)+1), w_s)
			  nut_avg = bilinear_interp(Yn2(n), Zn2(n), Ys(near_Yn2(n)), Ys(near_Yn2(n)+1), Zs(near_Zn2(n)), Zs(near_Zn2(n)+1), nut_s)

			theta_tan = ATAN2( -(Zn2(n)-z0_t2), (Yn2(n)-y0_t2) )
			IF (theta_tan >  0.5d0*PI)  theta_tan = theta_tan - PI
			IF (theta_tan < -0.5d0*PI)  theta_tan = theta_tan + PI

			tau_moment2 	= tau_moment2 + (nu + nut_avg) * den_flu * d_area * blade_r * SIGN(1.d0,(Yn2(n) - y0_t2)) * &
								(w_avg*COS(theta_tan) + v_avg*SIN(theta_tan) &
								+ blade_omega * blade_r * SIGN(1.d0,(Yn2(n) - y0_t2))) / (2.d0*dzSml)
			endif
		enddo   	
	enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO	


	! calculating wall parameters of blade 3
   !$OMP PARALLEL DO	           PRIVATE(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) REDUCTION(+: tau_moment3) COLLAPSE(2)
   !$acc parallel loop independent private(k,j,v_s,w_s,nut_s,v_avg,w_avg,nut_avg,theta_tan) reduction(+: tau_moment3) collapse(2) gang vector
	do n = 0,360
		do i=1,nx
			if(near_Zn3(n) .GE. istart .AND. near_Zn3(n) .LE. iend) then
				do k=1,2 ; do j=1,2
					v_s(j,k) = 0.5d0*( v(i,near_Yn3(n)  +(j-1),near_Zn3(n)  +(k-1))+ &
									   v(i,near_Yn3(n)-1+(j-1),near_Zn3(n)  +(k-1)) )
					w_s(j,k) = 0.5d0*( w(i,near_Yn3(n)  +(j-1),near_Zn3(n)  +(k-1))+ &
									   w(i,near_Yn3(n)  +(j-1),near_Zn3(n)-1+(k-1)) )
					nut_s(j,k) =     nut(i,near_Yn3(n)  +(j-1),near_Zn3(n)  +(k-1))
				enddo ; enddo
				
				v_avg = bilinear_interp(Yn3(n), Zn3(n), Ys(near_Yn3(n)), Ys(near_Yn3(n)+1), Zs(near_Zn3(n)), Zs(near_Zn3(n)+1), v_s)
				w_avg = bilinear_interp(Yn3(n), Zn3(n), Ys(near_Yn3(n)), Ys(near_Yn3(n)+1), Zs(near_Zn3(n)), Zs(near_Zn3(n)+1), w_s)
			  nut_avg = bilinear_interp(Yn3(n), Zn3(n), Ys(near_Yn3(n)), Ys(near_Yn3(n)+1), Zs(near_Zn3(n)), Zs(near_Zn3(n)+1), nut_s)

			theta_tan = ATAN2( -(Zn3(n)-z0_t3), (Yn3(n)-y0_t3) )
			IF (theta_tan >  0.5d0*PI)  theta_tan = theta_tan - PI
			IF (theta_tan < -0.5d0*PI)  theta_tan = theta_tan + PI

			tau_moment3 	= tau_moment3 + (nu + nut_avg) * den_flu * d_area * blade_r * SIGN(1.d0,(Yn3(n) - y0_t3)) * &
								(w_avg*COS(theta_tan) + v_avg*SIN(theta_tan) &
								+ blade_omega * blade_r * SIGN(1.d0,(Yn3(n) - y0_t3))) / (2.d0*dzSml)
			endif
		enddo   	
	enddo
   !$acc end parallel 
   !$OMP END PARALLEL DO	
	
!$acc end data 
	

    call MPI_ALLREDUCE( totalFZ1, totalFZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY1, totalFY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ1, totalTZ1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY1, totalTY1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )	
    call MPI_ALLREDUCE( tau_moment1, tau_moment1_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

    call MPI_ALLREDUCE( totalFZ2, totalFZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY2, totalFY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ2, totalTZ2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY2, totalTY2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( tau_moment2, tau_moment2_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

    call MPI_ALLREDUCE( totalFZ3, totalFZ3_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY3, totalFY3_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTZ3, totalTZ3_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalTY3, totalTY3_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( tau_moment3, tau_moment3_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

		totalFZ_ = totalFZ1_ + totalFZ2_ + totalFZ3_ 
		totalFY_ = totalFY1_ + totalFY2_ + totalFY3_


		totalTorq1 =   (totalTZ1_ + totalTY1_)
		totalTorq2 =   (totalTZ2_ + totalTY2_)
		totalTorq3 =   (totalTZ3_ + totalTY3_)
		
		totalTorq  = totalTorq1 + totalTorq2 + totalTorq3
		
	
		cDrag = -2.d0 * totalFZ_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift = -2.d0 * totalFY_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq = -2.d0 * totalTorq/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower   =   cTorq*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cPower_in   = -2.d0*blade_omega*(tau_moment1_  + tau_moment2_  + tau_moment3_) / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
		
		cPower_nett = cPower + cPower_in


		cDrag1 = -2.d0 * totalFZ1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift1 = -2.d0 * totalFY1_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq1 = -2.d0 * totalTorq1/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower1   =   cTorq1*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cPower1_in   = -2.d0*tau_moment1_*blade_omega / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
		
		cDrag2 = -2.d0 * totalFZ2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift2 = -2.d0 * totalFY2_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq2 = -2.d0 * totalTorq2/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower2   =   cTorq2*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)

		cPower2_in   = -2.d0*tau_moment2_*blade_omega / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
	
		cDrag3 = -2.d0 * totalFZ3_ / (ref_area*U_inf*U_inf*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift3 = -2.d0 * totalFY3_ / (ref_area*U_inf*U_inf*den_flu*1.d0)
		cTorq3 = -2.d0 * totalTorq3/(2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*den_flu*(rotor_r+blade_r)) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower3   =   cTorq3*ABS(rotor_omega*(rotor_r+blade_r)/U_inf)
		
		cPower3_in   = -2.d0*tau_moment3_*blade_omega / (2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*U_inf*den_flu)
		

end subroutine virtualForce_Magnus3blade


subroutine virtualForceTorqueIntegrator3D()
	use variables
	implicit none
	real*8 :: FZs,FYs,FXs

	totalFX = 0.d0
	totalFZ = 0.d0
	totalFY = 0.d0
	
	totalT_XZ = 0.d0
	totalT_XY = 0.d0
	totalT_YZ = 0.d0
	totalT_YX = 0.d0
	totalT_ZY = 0.d0
	totalT_ZX = 0.d0

!$acc data present(iDx,iDy,iDz,Xs,Ys,Zs,FX(:,:,istart-2:iend+2),FY(:,:,istart-2:iend+2),FZ(:,:,istart-2:iend+2))


	!$OMP PARALLEL	
	!$OMP DO               PRIVATE(FZs,FYs,FXs) REDUCTION(+: totalFX, totalFZ, totalFY) COLLAPSE(3)
	!$acc parallel
	!$acc loop independent private(FZs,FYs,FXs) reduction(+: totalFX, totalFZ, totalFY) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))
	FXs = 0.5d0*(FX(i,j,k)+FX(i-1,j,k))

		totalFZ = totalFZ + den_flu*( FZs * iDx(i)*iDy(j)*iDz(k) )
		totalFY = totalFY + den_flu*( FYs * iDx(i)*iDy(j)*iDz(k) )
		totalFX = totalFX + den_flu*( FXs * iDx(i)*iDy(j)*iDz(k) )

	end do;end do;end do
	!$OMP END DO

    
	! Calculate torque wrt (z0,y0)

    !$OMP DO               PRIVATE(FZs,FYs,FXs) REDUCTION(+: totalT_XZ, totalT_XY, totalT_YZ, totalT_YX, totalT_ZY, totalT_ZX) COLLAPSE(3)
	!$acc loop independent private(FZs,FYs,FXs) reduction(+: totalT_XZ, totalT_XY, totalT_YZ, totalT_YX, totalT_ZY, totalT_ZX) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx

	! Cell center values
	FZs = 0.5d0*(FZ(i,j,k)+FZ(i,j,k-1))
	FYs = 0.5d0*(FY(i,j,k)+FY(i,j-1,k))
	FXs = 0.5d0*(FX(i,j,k)+FX(i-1,j,k))

		totalT_XZ = totalT_XZ + den_flu*( FZs*(ys(j)-y0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0   
		totalT_XY = totalT_XY + den_flu*(-FYs*(zs(k)-z0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0		
			!totalT_XY = totalT_XY + den_flu*(-FYs*(zs(k)-z0)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0		moving Z reference
		
		
		totalT_YZ = totalT_YZ + den_flu*(-FZs*(xs(i)-x0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
			!totalT_YZ = totalT_YZ + den_flu*(-FZ(i,j,k)*(xs(i)-x0)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0   	moving X reference
		totalT_YX = totalT_YX + den_flu*( FXs*(zs(k)-z0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0	
			!totalT_YX = totalT_YX + den_flu*( FXs*(zs(k)-z0)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0		moving Z reference
		
		
		totalT_ZY = totalT_ZY + den_flu*( FYs*(xs(i)-x0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    
			!totalT_ZY = totalT_ZY + den_flu*( FYs*(xs(i)-x0)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0    	moving X reference
		totalT_ZX = totalT_ZX + den_flu*(-FXs*(ys(j)-y0_t)*(iDx(i)*iDy(j)*iDz(k) ) )        !CCW > 0


    end do;end do;end do
    !$OMP END DO
	!$OMP END PARALLEL
   !$acc end parallel 
!$acc end data 
	

    call MPI_ALLREDUCE( totalFX, totalFX_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFZ, totalFZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalFY, totalFY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
	
	call MPI_ALLREDUCE( totalT_XZ, totalT_XZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalT_XY, totalT_XY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

	
	call MPI_ALLREDUCE( totalT_YZ, totalT_YZ_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalT_YX, totalT_YX_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
	
	call MPI_ALLREDUCE( totalT_ZY, totalT_ZY_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE( totalT_ZX, totalT_ZX_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )


    !if(myid==master) then 	
	    totalTorqx =   (totalT_XZ_ + totalT_XY_)
		totalTorqy =   (totalT_YZ_ + totalT_YX_)
		totalTorqz =   (totalT_ZY_ + totalT_ZX_)
	
		cDrag = -2.d0 * totalFZ_ / (ref_area*w_solid*w_solid*den_flu*1.d0)  ! because for dimensional: (-2.d0)*totalFX / (den_flu*ref_area*U_inf*U_inf) 
		cLift = -2.d0 * totalFY_ / (ref_area*w_solid*w_solid*den_flu*1.d0)
		cTorq = -2.d0 * totalTorq/(2.d0*(rotor_r+blade_r)*w_solid*w_solid*den_flu*rotor_r*1.d0) ! (-2.d0)*totalTorq / (den_flu*area*U_inf*U_inf*rotor_r)
																							! Positive torque = CCW   
		cPower   =   cTorq*ABS(rotor_tsr)
	!endif
	

end subroutine virtualForceTorqueIntegrator3D


