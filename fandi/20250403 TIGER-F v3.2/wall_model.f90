! 09 Oct 2024 by FDS


subroutine wall_mod_cyl()
use variables
implicit none

integer		:: jj,kk
real*8		:: zb_nut,yb_nut,zr_nut,yr_nut
real*8		:: wra,wrb,wrc,wrd,vra,vrb,vrc,vrd
real*8		:: wr_ceil,wr_floor,wr
real*8		:: vr_ceil,vr_floor,vr
real*8		:: uc_r,uc_b,theta_b
real*8		:: dist_r,dist_b
real*8		:: KappaRe,yplusR,yplusR_,yplusB,uPlusB


   ! open (83,file='wall_yplus.dat',position='append')
   open (84,file='wall_uplus_instant.dat',position='append')
   ! open (85,file='wall_info_instant.dat',position='append')
   ! open (86,file='wall_nut_instant.dat',position='append')
   ! open (87,file='wall_yplus_instant.dat',position='append')
   ! open (88,file='wall_wdist_instant.dat',position='append')
   ! open (89,file='wall_moddamp_instant.dat',position='append')

   ! open (83,file='face_center_u.dat',position='append')
   ! open (84,file='face_center_v.dat',position='append')


!$acc data present(ETA,Zs,Ys,v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2),nut(:,:,istart:iend+1),yplus_print(:,:,istart:iend),uplus_print(:,:,istart:iend),dist_print(:,:,istart:iend))

!-----------------------------w-component-----------------------------------



!$OMP PARALLEL DO PRIVATE(j,i,zb_nut,yb_nut,zr_nut,yr_nut,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(3)
!$acc parallel
!$acc loop independent private(zb_nut,yb_nut,zr_nut,yr_nut,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(3) gang vector
do k=istart,iend
	do j=jBgnVOS-INT(c_ref+1.),jEndVOS+INT(c_ref+1.)
	do i=1,nx

		zb_nut = Zs(k)
		yb_nut = Ys(j)
		
		! Calculate the surface parallel slope
		theta_b = ATAN(-1.d0*(zb_nut - z0_t)/(yb_nut - y0_t))

	! Find the impossing face centers, ignoring vertical surface to avoid a division by zero
	if(ETA(i,j,k) .LT. 1.d0 .AND. ((zb_nut-z0_t)**2 + (yb_nut-y0_t)**2) .LT. (r+c_ref*dySml)**2) then

		zr_nut = (r+c_ref*dySml)*(zb_nut-z0_t)/SQRT((zb_nut-z0_t)**2 + (yb_nut-y0_t)**2) + z0_t
		yr_nut = (r+c_ref*dySml)*(yb_nut-y0_t)/SQRT((zb_nut-z0_t)**2 + (yb_nut-y0_t)**2) + y0_t

		
		! Find Uc using the bilinear interpolation
		!$acc loop seq
		do kk=istart,iend
			if((zr_nut .GE. Zs(kk)) .AND. (zr_nut .LT. Zs(kk+1))) then
				!$acc loop seq
				do jj=jBgnVOS-INT(c_ref+1.),jEndVOS+INT(c_ref+1.)
					if((yr_nut .GE. Ys(jj)) .AND. (yr_nut .LT. Ys(jj+1))) then
						wra = 0.5d0*(w(i,jj,kk)     + w(i,jj,kk-1))
						wrb = 0.5d0*(w(i,jj,kk+1)   + w(i,jj,kk)) 
						wrc = 0.5d0*(w(i,jj+1,kk)   + w(i,jj+1,kk-1)) 
						wrd = 0.5d0*(w(i,jj+1,kk+1) + w(i,jj+1,kk)) 

						vra = 0.5d0*(v(i,jj,kk)     + v(i,jj-1,kk))
						vrb = 0.5d0*(v(i,jj,kk+1)   + v(i,jj-1,kk+1))
						vrc = 0.5d0*(v(i,jj+1,kk)   + v(i,jj,kk))
						vrd = 0.5d0*(v(i,jj+1,kk+1) + v(i,jj,kk+1))
						
						wr_ceil  = (zr_nut - Zs(kk))*(wrd - wrc)/(Zs(kk+1) - Zs(kk)) + wrc
						wr_floor = (zr_nut - Zs(kk))*(wrb - wra)/(Zs(kk+1) - Zs(kk)) + wra
						
						wr = (yr_nut - Ys(jj))*(wr_ceil - wr_floor)/(Ys(jj+1) - Ys(jj)) + wr_floor

						vr_ceil  = (zr_nut - Zs(kk))*(vrd - vrc)/(Zs(kk+1) - Zs(kk)) + vrc
						vr_floor = (zr_nut - Zs(kk))*(vrb - vra)/(Zs(kk+1) - Zs(kk)) + vra
							
						vr = (yr_nut - Ys(jj))*(vr_ceil - vr_floor)/(Ys(jj+1) - Ys(jj)) + vr_floor						
							
						uc_r = ABS(wr*COS(theta_b) + vr*SIN(theta_b))
						!write(83,'(7(3x,F15.10))') zb_u,yb_u,zr_u,yr_u,wr,vr,uc_r


						! Calculate uTau at point R using Newton-Raphson iteration (= uTau at point B)
						dist_r = c_ref*dySml
						dist_b = SQRT((zb_nut-z0_t)**2+(yb_nut-y0_t)**2) - r
							dist_b = MAX(0.d0,dist_b)
						
						dist_print(i,j,k) = dist_b
						
						KappaRe = Kappa*uc_r*dist_r/nu
						yplusR = yplusLam
						
						! Newton-Raphson iteration	
						NRit=0 ; NRerr=1.d0
						do while (NRtolerance < NRerr .AND. NRit < NRmaxIter)					
							yplusR_ = (KappaRe + yplusR)/(1.d0 + LOG(E_WM * yplusR))

							NRerr = ABS((yplusR - yplusR_)/yplusR)					
							NRit=NRit+1
							yplusR = yplusR_
						end do
							
						yplusR = MAX(0.d0,yplusR)
						uTau = yplusR*nu/dist_r
						
						! Calculate yPlus at point B
						yplusB = dist_b*uTau/nu
						yplus_print(i,j,k) = yplusB
						! Calculate uPlus at point B
						if(yplusB .LE. yplusLam) then
							uPlusB = uplusB
							uplus_print(i,j,k) = yplusB
							nut(i,j,k) = 0.d0
						else
							uPlusB = LOG(E_WM * yplusB)/Kappa
							uplus_print(i,j,k) = uplusB
							nut(i,j,k) = MAX(0.d0,(yplusB*nu/(uPlusB + ROOTVSMALL)) - nu)
						endif

						! Calculate the magnitude of Uc at point B
						!uc_b = uPlusB*uTau

					

!write(*,'(5(3x,F15.10))') yplusR,yplusB,uPlusB,uTau
						
! if (myid==master .AND. i==1 .AND. (istep == 10 .OR. mod(istep,1000)==0)) then
! write(84,'(I6,7(3x,F15.10))') istep,zb_nut,yb_nut,yplusR,yplusB,uPlusB,uTau,nut(i,j,k)
! endif
					endif
				enddo
			endif
		enddo						
	endif
enddo
enddo
enddo
!$acc end parallel
!$OMP END PARALLEL DO


!$acc end data 

    ! close(83)
    close(84)
    ! close(85)
    ! close(86)
    ! close(87)
    ! close(88)
	! close(89)

 
end subroutine wall_mod_cyl


! subroutine wall_mod_cyl()
! use variables
! implicit none

! integer		:: jj,kk
! real*8		:: zb_u,yb_u,zr_u,yr_u
! real*8		:: zb_v,yb_v,zr_v,yr_v
! real*8		:: wra,wrb,wrc,wrd,vra,vrb,vrc,vrd
! real*8		:: wr_ceil,wr_floor,wr
! real*8		:: vr_ceil,vr_floor,vr
! real*8		:: uc_r,uc_b,theta_b
! real*8		:: dist_r,dist_b
! real*8		:: KappaRe,yplusR,yplusR_,yplusB,uPlusB


   ! ! open (83,file='wall_yplus.dat',position='append')
   ! ! open (84,file='wall_uplus_instant.dat',position='append')
   ! ! open (85,file='wall_info_instant.dat',position='append')
   ! ! open (86,file='wall_nut_instant.dat',position='append')
   ! ! open (87,file='wall_yplus_instant.dat',position='append')
   ! ! open (88,file='wall_wdist_instant.dat',position='append')
   ! ! open (89,file='wall_moddamp_instant.dat',position='append')

   ! open (83,file='face_center_u.dat',position='append')
   ! open (84,file='face_center_v.dat',position='append')


! !$acc data present(ETA,Z,Zs,Y,Ys,v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2),v2(:,:,istart-2:iend+2),w2(:,:,istart-2:iend+2))

! !-----------------------------w-component-----------------------------------



! !$OMP PARALLEL DO PRIVATE(j,i,zb_u,yb_u,zr_u,yr_u,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(3)
! !$acc parallel
! !$acc loop independent private(j,i,zb_u,yb_u,zr_u,yr_u,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(3) gang vector
! do k=istart,iend
	! do j=jBgnVOS,jEndVOS
	! do i=1,nx

		! zb_u = Z(k+1)
		! yb_u = Ys(j)
		
		! ! Calculate the surface parallel slope
		! theta_b = ATAN(-1.d0*(zb_u - z0_t)/(yb_u - y0_t))

	! ! Find the impossing face centers, ignoring vertical surface to avoid a division by zero
	! if((0.5d0*ETA(i,j,k)+ETA(i,j,k+1)) .LT. 0.5d0 .AND. ((zb_u-z0_t)**2 + (yb_u-y0_t)**2) .LT. (r+c_ref*dySml)**2) then

		! zr_u = (r+c_ref*dySml)*(zb_u-z0_t)/SQRT((zb_u-z0_t)**2 + (yb_u-y0_t)**2) + z0_t
		! yr_u = (r+c_ref*dySml)*(yb_u-y0_t)/SQRT((zb_u-z0_t)**2 + (yb_u-y0_t)**2) + y0_t

		
		! ! Find Uc using the bilinear interpolation (w2 and v2 are used to avoid overlaping when updating w and v)
		! !$acc loop seq
		! do kk=istart,iend
			! if((zr_u .GE. Zs(kk)) .AND. (zr_u .LT. Zs(kk+1))) then
				! !$acc loop seq
				! do jj=jBgnVOS,jEndVOS
					! if((yr_u .GE. Ys(jj)) .AND. (yr_u .LT. Ys(jj+1))) then
						! wra = 0.5d0*(w2(i,jj,kk)     + w2(i,jj,kk-1))
						! wrb = 0.5d0*(w2(i,jj,kk+1)   + w2(i,jj,kk)) 
						! wrc = 0.5d0*(w2(i,jj+1,kk)   + w2(i,jj+1,kk-1)) 
						! wrd = 0.5d0*(w2(i,jj+1,kk+1) + w2(i,jj+1,kk)) 

						! vra = 0.5d0*(v2(i,jj,kk)     + v2(i,jj-1,kk))
						! vrb = 0.5d0*(v2(i,jj,kk+1)   + v2(i,jj-1,kk+1))
						! vrc = 0.5d0*(v2(i,jj+1,kk)   + v2(i,jj,kk))
						! vrd = 0.5d0*(v2(i,jj+1,kk+1) + v2(i,jj,kk+1))
						
						! wr_ceil  = (zr_u - Zs(kk))*(wrd - wrc)/(Zs(kk+1) - Zs(kk)) + wrc
						! wr_floor = (zr_u - Zs(kk))*(wrb - wra)/(Zs(kk+1) - Zs(kk)) + wra
						
						! wr = (yr_u - Ys(jj))*(wr_ceil - wr_floor)/(Ys(jj+1) - Ys(jj)) + wr_floor

						! vr_ceil  = (zr_u - Zs(kk))*(vrd - vrc)/(Zs(kk+1) - Zs(kk)) + vrc
						! vr_floor = (zr_u - Zs(kk))*(vrb - vra)/(Zs(kk+1) - Zs(kk)) + vra
							
						! vr = (yr_u - Ys(jj))*(vr_ceil - vr_floor)/(Ys(jj+1) - Ys(jj)) + vr_floor						
							
						! uc_r = ABS(wr*COS(theta_b) + vr*SIN(theta_b))
						! !write(83,'(7(3x,F15.10))') zb_u,yb_u,zr_u,yr_u,wr,vr,uc_r


						! ! Calculate uTau at point R using Newton-Raphson iteration (= uTau at point B)
						! dist_r = c_ref*dySml
						! dist_b = SQRT((zb_u-z0_t)**2+(yb_u-y0_t)**2) - r
							! dist_b = MAX(0.d0,dist_b)
						
						! KappaRe = Kappa*uc_r*dist_r/nu
						! yplusR = yplusLam
						
						! ! Newton-Raphson iteration	
						! NRit=0 ; NRerr=1.d0
						! do while (NRtolerance < NRerr .AND. NRit < NRmaxIter)					
							! yplusR_ = (KappaRe + yplusR)/(1.d0 + LOG(E_WM * yplusR))

							! NRerr = ABS((yplusR - yplusR_)/yplusR)					
							! NRit=NRit+1
							! yplusR = yplusR_
						! end do
							
						! yplusR = MAX(0.d0,yplusR)
						! uTau = yplusR*nu/dist_r
						
						! ! Calculate yPlus at point B
						! yplusB = dist_b*uTau/nu
						
						! ! Calculate uPlus at point B
						! if(yplusB .LE. yplusLam) then
							! uPlusB = yplusB
						! else
							! uPlusB = LOG(E_WM * yplusB)/Kappa
						! endif

						! ! Calculate the magnitude of Uc at point B
						! uc_b = uPlusB*uTau

! !write(*,'(5(3x,F15.10))') yplusR,yplusB,uPlusB,uTau,uc_b
				
						! ! Impose w-velocity, the last term represents fluid velocity portion on the interface cell
						! ! SIGN is used to re-imposse the same velocity direction as on R
						! w(i,j,k) = wr * (uc_b/uc_r) * (1.d0 - (0.5d0*ETA(i,j,k)+ETA(i,j,k+1)))
					
						
						! ! yplusRjk(j,k) = yplusR
						! ! yplusBjk(j,k) = yplusB
						! ! wjk(j,k) = w(1,j,k)
						
! ! if (myid==master .AND. (istep == 10 .OR. mod(istep,1000)==0)) then
! ! write(*,'(I6,5(3x,F15.10))') istep,yplusR,yplusB,zb_u,yb_u,w(1,j,k)
! ! endif
					! endif
				! enddo
			! endif
		! enddo						
	! endif
! enddo
! enddo
! enddo
! !$acc end parallel
! !$OMP END PARALLEL DO


! !-----------------------------v-component-----------------------------------

! !$OMP PARALLEL DO PRIVATE(j,i,zb_v,yb_v,zr_v,yr_v,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(3)
! !$acc parallel
! !$acc loop independent private(j,i,zb_v,yb_v,zr_v,yr_v,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(3) gang vector

! do k=istart,iend
	! do j=jBgnVOS,jEndVOS
	! do i=1,nx

		! zb_v = Zs(k)
		! yb_v = Y(j+1)

		! ! Calculate the surface parallel slope
		! theta_b = ATAN(-1.d0*(zb_v - z0_t)/(yb_v - y0_t))

	! ! Find the impossing face centers, ignoring horizontal surface to avoid a division by zero
	! if((0.5d0*ETA(i,j,k)+ETA(i,j+1,k)) .LT. 0.5d0 .AND. ((zb_v-z0_t)**2 + (yb_v-y0_t)**2) .LT. (r+c_ref*dySml)**2) then

		! zr_v = (r+c_ref*dySml)*(zb_v-z0_t)/SQRT((zb_v-z0_t)**2 + (yb_v-y0_t)**2) + z0_t
		! yr_v = (r+c_ref*dySml)*(yb_v-y0_t)/SQRT((zb_v-z0_t)**2 + (yb_v-y0_t)**2) + y0_t	


		! ! Find Uc using the bilinear interpolation (w2 and v2 are used to avoid overlaping when updating w and v)
		! !$acc loop seq

		! do kk=istart,iend
			! if((zr_v .GE. Zs(kk)) .AND. (zr_v .LT. Zs(kk+1))) then
				! !$acc loop seq
				! do jj=jBgnVOS,jEndVOS
					! if((yr_v .GE. Ys(jj)) .AND. (yr_v .LT. Ys(jj+1))) then
						! wra = 0.5d0*(w2(i,jj,kk)     + w2(i,jj,kk-1))
						! wrb = 0.5d0*(w2(i,jj,kk+1)   + w2(i,jj,kk)) 
						! wrc = 0.5d0*(w2(i,jj+1,kk)   + w2(i,jj+1,kk-1)) 
						! wrd = 0.5d0*(w2(i,jj+1,kk+1) + w2(i,jj+1,kk)) 

						! vra = 0.5d0*(v2(i,jj,kk)     + v2(i,jj-1,kk))
						! vrb = 0.5d0*(v2(i,jj,kk+1)   + v2(i,jj-1,kk+1))
						! vrc = 0.5d0*(v2(i,jj+1,kk)   + v2(i,jj,kk))
						! vrd = 0.5d0*(v2(i,jj+1,kk+1) + v2(i,jj,kk+1))

						! wr_ceil  = (zr_u - Zs(kk))*(wrd - wrc)/(Zs(kk+1) - Zs(kk)) + wrc
						! wr_floor = (zr_u - Zs(kk))*(wrb - wra)/(Zs(kk+1) - Zs(kk)) + wra
						
						! wr = (yr_u - Ys(jj))*(wr_ceil - wr_floor)/(Ys(jj+1) - Ys(jj)) + wr_floor

						! vr_ceil  = (zr_v - Zs(kk))*(vrd - vrc)/(Zs(kk+1) - Zs(kk)) + vrc
						! vr_floor = (zr_v - Zs(kk))*(vrb - vra)/(Zs(kk+1) - Zs(kk)) + vra
							
						! vr = (yr_v - Ys(jj))*(vr_ceil - vr_floor)/(Ys(jj+1) - Ys(jj)) + vr_floor
						

						! uc_r = ABS(wr*COS(theta_b) + vr*SIN(theta_b))
						! !write(84,'(7(3x,F15.10))') zb_v,yb_v,zr_v,yr_v,wr,vr,uc_r

						! ! Calculate uTau at point R using Newton-Raphson iteration (= uTau at point B)
						! dist_r = c_ref*dySml
						! dist_b = SQRT((zb_v-z0_t)**2+(yb_v-y0_t)**2) - r
							! dist_b = MAX(0.d0,dist_b)
						
						! KappaRe = Kappa*uc_r*dist_r/nu
						! yplusR = yplusLam
						
						! ! Newton-Raphson iteration
						! NRit=0 ; NRerr=1.d0
						! do while (NRtolerance < NRerr .AND. NRit < NRmaxIter)					
							! yplusR_ = (KappaRe + yplusR)/(1.d0 + LOG(E_WM * yplusR))

							! NRerr = ABS((yplusR - yplusR_)/yplusR)					
							! NRit=NRit+1
							! yplusR = yplusR_
						! end do
							
						! yplusR = MAX(0.d0,yplusR)
						! uTau = yplusR*nu/dist_r
						
						! ! Calculate yPlus at point B
						! yplusB = dist_b*uTau/nu
						
						! ! Calculate uPlus at point B
						! if(yplusB .LE. yplusLam) then
							! uPlusB = yplusB
						! else
							! uPlusB = LOG(E_WM * yplusB)/Kappa
						! endif

						! ! Calculate the magnitude of Uc at point B
						! uc_b = uPlusB*uTau

! !write(*,'(3(3x,F15.10))') zb_v,yb_v,v(1,j,k)				
						! ! Impose v-velocity, the last term represents fluid velocity portion on the interface cell
						! ! SIGN is used to re-imposse the same velocity direction as on R
						! v(i,j,k) = vr * (uc_b/uc_r) *  (1.d0 - (0.5d0*ETA(i,j,k)+ETA(i,j+1,k))) 
				
						
! !write(*,'(6(3x,F15.10))') yplusR,yplusB,uPlusB,uTau,uc_b,v(1,j,k)
! !write(*,'(6(3x,F15.10))') yplusR,yplusB,uPlusB,zb_v,yb_v,v(1,j,k)
					! endif
				! enddo
			! endif
		! enddo						
	! endif
! enddo
! enddo
! enddo
! !$acc end parallel
! !$OMP END PARALLEL DO



! !-------------------------------------------------------------------------------------------------------------------------------------------------------------------		


! !$acc end data 

    ! close(83)
    ! close(84)
    ! ! close(85)
    ! ! close(86)
    ! ! close(87)
    ! ! close(88)
	! ! close(89)

 
! end subroutine wall_mod_cyl


! Only 1 layer in spanwise, and then copied

! subroutine wall_mod_cyl()
! use variables
! implicit none

! integer		:: jj,kk
! real*8		:: zb_u,yb_u,zr_u,yr_u
! real*8		:: zb_v,yb_v,zr_v,yr_v
! real*8		:: wra,wrb,wrc,wrd,vra,vrb,vrc,vrd
! real*8		:: wr_ceil,wr_floor,wr
! real*8		:: vr_ceil,vr_floor,vr
! real*8		:: uc_r,uc_b,theta_b
! real*8		:: dist_r,dist_b
! real*8		:: KappaRe,yplusR,yplusR_,yplusB,uPlusB


   ! ! open (83,file='wall_yplus.dat',position='append')
   ! ! open (84,file='wall_uplus_instant.dat',position='append')
   ! ! open (85,file='wall_info_instant.dat',position='append')
   ! ! open (86,file='wall_nut_instant.dat',position='append')
   ! ! open (87,file='wall_yplus_instant.dat',position='append')
   ! ! open (88,file='wall_wdist_instant.dat',position='append')
   ! ! open (89,file='wall_moddamp_instant.dat',position='append')

   ! open (83,file='face_center_u.dat',position='append')
   ! open (84,file='face_center_v.dat',position='append')


! !-----------------------------u-component-----------------------------------

! !$acc data present(ETA,Z,Zs,Y,Ys,v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2))


! !$OMP PARALLEL DO PRIVATE(j,i,zb_u,yb_u,zr_u,yr_u,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(2)
! !$acc parallel
! !$acc loop independent private(j,i,zb_u,yb_u,zr_u,yr_u,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(2) gang vector
! do k=istart,iend
	! do j=jBgnVOS,jEndVOS
	
	! ! Find the impossing face centers 
	! if(((0.5d0*ETA(INT(nx/2),j,k)+ETA(INT(nx/2),j,k+1)) .LE. 0.5d0) .AND. (((Z(k+1)-z0_t)**2 + (Ys(j)-y0_t)**2) .LT. (r+c_ref*dySml)**2 )) then
		! zb_u = Z(k+1)
		! yb_u = Ys(j)
		! zr_u = (r+c_ref*dySml)*(zb_u-z0_t)/SQRT((zb_u-z0_t)**2 + (yb_u-y0_t)**2) + z0_t
		! yr_u = (r+c_ref*dySml)*(yb_u-y0_t)/SQRT((zb_u-z0_t)**2 + (yb_u-y0_t)**2) + y0_t

		! ! Calculate the surface parallel slope
		! theta_b = ATAN(-1.d0*(zb_u - z0_t)/(yb_u - y0_t))
		
		! ! Find Uc using the bilinear interpolation
		! !$acc loop seq
		! do kk=istart,iend
			! if((zr_u .GE. Zs(kk)) .AND. (zr_u .LT. Zs(kk+1))) then
				! !$acc loop seq
				! do jj=jBgnVOS,jEndVOS
					! if((yr_u .GE. Ys(jj)) .AND. (yr_u .LT. Ys(jj+1))) then
						! wra = 0.5d0*(w(INT(nx/2),jj,kk)     + w(INT(nx/2),jj,kk-1))
						! wrb = 0.5d0*(w(INT(nx/2),jj,kk+1)   + w(INT(nx/2),jj,kk)) 
						! wrc = 0.5d0*(w(INT(nx/2),jj+1,kk)   + w(INT(nx/2),jj+1,kk-1)) 
						! wrd = 0.5d0*(w(INT(nx/2),jj+1,kk+1) + w(INT(nx/2),jj+1,kk)) 
						
						! vra = 0.5d0*(v(INT(nx/2),jj,kk)     + v(INT(nx/2),jj-1,kk))
						! vrb = 0.5d0*(v(INT(nx/2),jj,kk+1)   + v(INT(nx/2),jj-1,kk+1))
						! vrc = 0.5d0*(v(INT(nx/2),jj+1,kk)   + v(INT(nx/2),jj,kk))
						! vrd = 0.5d0*(v(INT(nx/2),jj+1,kk+1) + v(INT(nx/2),jj,kk+1))
						
						! wr_ceil  = (zr_u - Zs(kk))*(wrd - wrc)/(Zs(kk+1) - Zs(kk)) + wrc
						! wr_floor = (zr_u - Zs(kk))*(wrb - wra)/(Zs(kk+1) - Zs(kk)) + wra
						
						! wr = (yr_u - Ys(jj))*(wr_ceil - wr_floor)/(Ys(jj+1) - Ys(jj)) + wr_floor
						

						! vr_ceil  = (zr_u - Zs(kk))*(vrd - vrc)/(Zs(kk+1) - Zs(kk)) + vrc
						! vr_floor = (zr_u - Zs(kk))*(vrb - vra)/(Zs(kk+1) - Zs(kk)) + vra
							
						! vr = (yr_u - Ys(jj))*(vr_ceil - vr_floor)/(Ys(jj+1) - Ys(jj)) + vr_floor
							
						! uc_r = wr*COS(theta_b) + vr*SIN(theta_b)
						! !write(83,'(7(3x,F15.10))') zb_u,yb_u,zr_u,yr_u,wr,vr,uc_r


						! ! Calculate uTau at point R using Newton-Raphson iteration (= uTau at point B)
						! dist_r = c_ref*dySml
						! dist_b = SQRT((zb_u-z0_t)**2+(yb_u-y0_t)**2) - r
							! dist_b = MAX(0.d0,dist_b)
						
						! KappaRe = ABS(Kappa*uc_r*dist_r/nu) ! Use ABS because KappaRe must be +, while uc_r may be +/-
						! yplusR = yplusLam
						
						! ! Newton-Raphson iteration	
						! NRit=0 ; NRerr=1.d0
						! do while (NRtolerance < NRerr .AND. NRit < NRmaxIter)					
							! yplusR_ = (KappaRe + yplusR)/(1.d0 + LOG(E_WM * yplusR))

							! NRerr = ABS((yplusR - yplusR_)/yplusR)					
							! NRit=NRit+1
							! yplusR = yplusR_
						! end do
							
						! yplusR = MAX(0.d0,yplusR)
						! uTau = yplusR*nu/dist_r
						
						! ! Calculate yPlus at point B
						! yplusB = dist_b*uTau/nu
						
						! ! Calculate uPlus at point B
						! if(yplusB .LE. yplusLam) then
							! uPlusB = yplusB
						! else
							! uPlusB = LOG(E_WM * yplusB)/Kappa
						! endif

						! ! Calculate Uc at point B
						! uc_b = uPlusB*uTau*SIGN(1.d0,uc_r)	! SIGN is used to re-imposse the direction of uc based on uc_r

! !write(*,'(5(3x,F15.10))') yplusR,yplusB,uPlusB,uTau,uc_b
				
						! ! Impose w-velocity
						! w(1,j,k) = uc_b * COS(theta_b)
						
						! ! yplusRjk(j,k) = yplusR
						! ! yplusBjk(j,k) = yplusB
						! ! wjk(j,k) = w(1,j,k)
						
! ! if (myid==master .AND. (istep == 10 .OR. mod(istep,1000)==0)) then
! ! write(*,'(I6,5(3x,F15.10))') istep,yplusR,yplusB,zb_u,yb_u,w(1,j,k)
! ! endif
					! endif
				! enddo
			! endif
		! enddo						
	! endif
! enddo
! enddo
! !$acc end parallel
! !$OMP END PARALLEL DO


! !-----------------------------v-component-----------------------------------

! !$OMP PARALLEL DO PRIVATE(j,i,zb_v,yb_v,zr_v,yr_v,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(2)
! !$acc parallel
! !$acc loop independent private(j,i,zb_v,yb_v,zr_v,yr_v,theta_b,kk,jj,wra,wrb,wrc,wrd,vra,vrb,vrc,vrd,wr_ceil,wr_floor,wr,vr_ceil,vr_floor,vr,uc_r,dist_r,dist_b,KappaRe,yplusR,yplusR_,NRit,NRerr,uTau,yplusB,uPlusB,uc_b) collapse(2) gang vector

! do k=istart,iend
	! do j=jBgnVOS,jEndVOS

	! ! Find the impossing face centers 	
	! if(((0.5d0*ETA(INT(nx/2),j,k)+ETA(INT(nx/2),j+1,k)) .LE. 0.5d0) .AND. (((Zs(k)-z0_t)**2 + (Y(j+1)-y0_t)**2) .LT. (r+c_ref*dySml)**2 )) then
		! zb_v = Zs(k)
		! yb_v = Y(j+1)
		! zr_v = (r+c_ref*dySml)*(zb_v-z0_t)/SQRT((zb_v-z0_t)**2 + (yb_v-y0_t)**2) + z0_t
		! yr_v = (r+c_ref*dySml)*(yb_v-y0_t)/SQRT((zb_v-z0_t)**2 + (yb_v-y0_t)**2) + y0_t	

		! ! Calculate the surface parallel slope
		! theta_b = ATAN(-1.d0*(zb_v - z0_t)/(yb_v - y0_t))

		! ! Find Uc using the bilinear interpolation
		! !$acc loop seq

		! do kk=istart,iend
			! if((zr_v .GE. Zs(kk)) .AND. (zr_v .LT. Zs(kk+1))) then
				! !$acc loop seq
				! do jj=jBgnVOS,jEndVOS
					! if((yr_v .GE. Ys(jj)) .AND. (yr_v .LT. Ys(jj+1))) then
						! wra = 0.5d0*(w(INT(nx/2),jj,kk)     + w(INT(nx/2),jj,kk-1))
						! wrb = 0.5d0*(w(INT(nx/2),jj,kk+1)   + w(INT(nx/2),jj,kk)) 
						! wrc = 0.5d0*(w(INT(nx/2),jj+1,kk)   + w(INT(nx/2),jj+1,kk-1)) 
						! wrd = 0.5d0*(w(INT(nx/2),jj+1,kk+1) + w(INT(nx/2),jj+1,kk)) 
						
						! vra = 0.5d0*(v(INT(nx/2),jj,kk)     + v(INT(nx/2),jj-1,kk))
						! vrb = 0.5d0*(v(INT(nx/2),jj,kk+1)   + v(INT(nx/2),jj-1,kk+1))
						! vrc = 0.5d0*(v(INT(nx/2),jj+1,kk)   + v(INT(nx/2),jj,kk))
						! vrd = 0.5d0*(v(INT(nx/2),jj+1,kk+1) + v(INT(nx/2),jj,kk+1))

						! vr_ceil  = (zr_v - Zs(kk))*(vrd - vrc)/(Zs(kk+1) - Zs(kk)) + vrc
						! vr_floor = (zr_v - Zs(kk))*(vrb - vra)/(Zs(kk+1) - Zs(kk)) + vra
							
						! vr = (yr_v - Ys(jj))*(vr_ceil - vr_floor)/(Ys(jj+1) - Ys(jj)) + vr_floor
						

						! wr_ceil  = (zr_v - Zs(kk))*(wrd - wrc)/(Zs(kk+1) - Zs(kk)) + wrc
						! wr_floor = (zr_v - Zs(kk))*(wrb - wra)/(Zs(kk+1) - Zs(kk)) + wra
						
						! wr = (yr_v - Ys(jj))*(wr_ceil - wr_floor)/(Ys(jj+1) - Ys(jj)) + wr_floor
							
						! uc_r = wr*COS(theta_b) + vr*SIN(theta_b)
						! !write(84,'(7(3x,F15.10))') zb_v,yb_v,zr_v,yr_v,wr,vr,uc_r

						! ! Calculate uTau at point R using Newton-Raphson iteration (= uTau at point B)
						! dist_r = c_ref*dySml
						! dist_b = SQRT((zb_v-z0_t)**2+(yb_v-y0_t)**2) - r
							! dist_b = MAX(0.d0,dist_b)
						
						! KappaRe = ABS(Kappa*uc_r*dist_r/nu) ! Use ABS because KappaRe must be +, while uc_r may be +/-
						! yplusR = yplusLam
						
						! ! Newton-Raphson iteration
						! NRit=0 ; NRerr=1.d0
						! do while (NRtolerance < NRerr .AND. NRit < NRmaxIter)					
							! yplusR_ = (KappaRe + yplusR)/(1.d0 + LOG(E_WM * yplusR))

							! NRerr = ABS((yplusR - yplusR_)/yplusR)					
							! NRit=NRit+1
							! yplusR = yplusR_
						! end do
							
						! yplusR = MAX(0.d0,yplusR)
						! uTau = yplusR*nu/dist_r
						
						! ! Calculate yPlus at point B
						! yplusB = dist_b*uTau/nu
						
						! ! Calculate uPlus at point B
						! if(yplusB .LE. yplusLam) then
							! uPlusB = yplusB
						! else
							! uPlusB = LOG(E_WM * yplusB)/Kappa
						! endif

						! ! Calculate Uc at point B
						! uc_b = uPlusB*uTau*SIGN(1.d0,uc_r) ! SIGN is used to re-imposse the direction of uc based on uc_r

! !write(*,'(3(3x,F15.10))') zb_v,yb_v,v(1,j,k)				
						! ! Impose v-velocity (positive slope)
						! v(1,j,k) = uc_b * SIN(theta_b)
						
! !write(*,'(6(3x,F15.10))') yplusR,yplusB,uPlusB,uTau,uc_b,v(1,j,k)
! !write(*,'(6(3x,F15.10))') yplusR,yplusB,uPlusB,zb_v,yb_v,v(1,j,k)
					! endif
				! enddo
			! endif
		! enddo						
	! endif
! enddo
! enddo
! !$acc end parallel
! !$OMP END PARALLEL DO

! ! Copying to spanwise cells (can not do this in the above loops because the race condition may occur)
! !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
! !$acc parallel loop independent collapse(2)
! do k=istart,iend
! do j=jBgnVOS,jEndVOS
	! if(((0.5d0*ETA(INT(nx/2),j,k)+ETA(INT(nx/2),j,k+1)) .LE. 0.5d0) .AND. (((Z(k+1)-z0_t)**2 + (Ys(j)-y0_t)**2) .LT. (r+c_ref*dySml)**2 )) then
		! do i=2,nx
			! w(i,j,k) = w(1,j,k)
		! end do
	! endif
		
! end do
! end do
! !$acc end parallel
! !$OMP END PARALLEL DO


! !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
! !$acc parallel loop independent collapse(2)
! do k=istart,iend
! do j=jBgnVOS,jEndVOS
	! if	(((0.5d0*ETA(INT(nx/2),j,k)+ETA(INT(nx/2),j+1,k)) .LE. 0.5d0) .AND. (((Zs(k)-z0_t)**2 + (Y(j+1)-y0_t)**2) .LT. (r+c_ref*dySml)**2 )) then
		! do i=2,nx
			! v(i,j,k) = v(1,j,k)
		! end do
	! endif
		
! end do
! end do
! !$acc end parallel
! !$OMP END PARALLEL DO


! !-------------------------------------------------------------------------------------------------------------------------------------------------------------------		


! !$acc end data 

    ! close(83)
    ! close(84)
    ! ! close(85)
    ! ! close(86)
    ! ! close(87)
    ! ! close(88)
	! ! close(89)

 
! end subroutine wall_mod_cyl




! This subroutine is used to calculate the shortest distance and orientation of each cell to the VOS.
! DISCLAIMER: this is only for extruded span solids!!!

subroutine wall_distance()
use variables
implicit none

real*8, dimension(:), allocatable   			 :: wp_up,wp_low,z_wp
real*8, dimension(:,:), allocatable 			 :: y_wp_up, y_wp_low
real*8											 :: dist

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------
! count how many k has wall points
! This part is enough to be executed in CPU cores
wp_k = 0

!$OMP PARALLEL DO PRIVATE(j) REDUCTION(+: wp_k,void)
!$acc parallel loop independent private(j) reduction(+: wp_k,void) device_type(host)
do k=kBgnVOS,kEndVOS
	void = 0
	do j=jBgnVOS,jEndVOS
		if ((ETA(1,j,k)>=0.5d0 .AND. ETA(1,j+1,k)<0.5d0) .OR. (ETA(1,j,k)>=0.5d0 .AND. ETA(1,j-1,k)<0.5d0))  then			
			void = void+1
		endif
	enddo
	if(void .GE. 1) then
		wp_k = wp_k + 1
	endif
enddo
!$acc end parallel
!$OMP END PARALLEL DO


allocate(wp_up(1:nz))
allocate(wp_low(1:nz))
allocate(z_wp(1:nz))
allocate(y_wp_up(1:nz,1:10))
allocate(y_wp_low(1:nz,1:10))

! count how many wall points on each k

!$acc data present(Y,Ys,Zs,Dys,ETA,wdist,wsin,wcos) &
!$acc create(wp_up,wp_low,z_wp,y_wp_up,y_wp_low)



!$OMP PARALLEL DO PRIVATE(j,void,k_wp,wp_j_up,wp_j_low,ratio)
!$acc parallel
!$acc loop independent private(j,void,k_wp,wp_j_up,wp_j_low,ratio)
do k=kBgnVOS,kEndVOS
	void = 0
	wp_j_up = 0
	wp_j_low = 0
	!$acc loop seq
	do j=jBgnVOS,jEndVOS
		if (ETA(1,j,k)>=0.5d0 .AND. ETA(1,j+1,k)<0.5d0) then			! for upper surface			
			wp_j_up = wp_j_up + 1
            if (ETA(1,j,k)>=1.d0 .AND. ETA(1,j+1,k)<=0.d0) then
				y_wp_up(k,wp_j_up) = Y(j+1)					
			else
				ratio=(ETA(1,j,k)-0.5d0)/(ETA(1,j,k)-ETA(1,j+1,k))
				y_wp_up(k,wp_j_up) = Ys(j)+ratio*Dys(j)
			endif
			void = void+1
		elseif (ETA(1,j,k)>=0.5d0 .AND. ETA(1,j-1,k)<0.5d0) then		! for lower surface
			wp_j_low = wp_j_low + 1
            if (ETA(1,j,k)>=1.d0 .AND. ETA(1,j-1,k)<=0.d0) then
				y_wp_low(k,wp_j_low) = Y(j)			
			else
				ratio=(ETA(1,j,k)-0.5d0)/(ETA(1,j,k)-ETA(1,j-1,k))	
				y_wp_low(k,wp_j_low) = Ys(j)-ratio*Dys(j)
			endif
			void = void+1
		endif
	enddo

	if(void .GE. 1) then
		z_wp(k) = Zs(k)
		wp_up(k) = wp_j_up			! number of upper wall points in each k
		wp_low(k) = wp_j_low			! number of lower wall points in each k

	endif
enddo
!$acc end parallel
!$OMP END PARALLEL DO

! Sort wp ascending sequentially from the lowest to the highest k
! This must run sequentially!
! NOTES: needs a better sorting algorithm to enable parallelization!!!!

wp_k = 1

!$acc serial
!$acc loop seq
do k=kBgnVOS,kEndVOS
	if (z_wp(k) .GT. 0.d0) then
		z_wp(wp_k) = Zs(k)
		wp_up(wp_k) = wp_up(k)
		wp_low(wp_k) = wp_low(k)
			!$acc loop seq
			do wp_j_up=1,wp_up(k)
				y_wp_up(wp_k,wp_j_up) = y_wp_up(k,wp_j_up)
			enddo
			!$acc loop seq
			do wp_j_low=1,wp_low(k)
				y_wp_low(wp_k,wp_j_low) = y_wp_low(k,wp_j_low)
			enddo
		wp_k = wp_k + 1
	endif
enddo
!$acc end serial


! Calculate the shortest distance of each nearest cell
max_dist = (dzSml**2+dySml**2)


!$OMP PARALLEL DO PRIVATE(j,wdistance,wdistance_,hyp,wsin_,wcos_,k_wp,wp_j_up,wp_j_low,j_wp,j_wpp)
!$acc parallel
!$acc loop independent private(j,wdistance,wdistance_,hyp,wsin_,wcos_,k_wp,wp_j_up,wp_j_low,j_wp,j_wpp) collapse(2) gang vector
do k=kBgnVOS,kEndVOS
	do j=jBgnVOS,jEndVOS
	wdistance = max_dist				
	wsin_ = 0.d0
	wcos_ = 0.d0
	if (ETA(1,j,k) .LT. 0.5d0) then
	 !$acc loop seq
		do k_wp = 1, wp_k-1
			!$acc loop seq
			do wp_j_up  = 1,  wp_up(k_wp)
			!$acc loop seq
			do wp_j_low = 1, wp_low(k_wp)
				if     (ABS(y_wp_up(k_wp,wp_j_up)-y_wp_low(k_wp,wp_j_low)) .LT. (0.01d0)) then
					wdistance_ = dist(Zs(k),Ys(j),z_wp(k_wp),y_wp_up(k_wp,wp_j_up),z_wp(k_wp),y_wp_low(k_wp,wp_j_low))
					if ((wdistance_ .LE. max_dist) .AND. (wdistance_ .LT. wdistance)) then
						wdistance = wdistance_
						wsin_ = 1.d0
						wcos_ = 0.d0
					endif
				endif
			enddo
			enddo			
			!$acc loop seq
			do j_wp  = 1, wp_up(k_wp)
			!$acc loop seq
			do j_wpp = 1, wp_up(k_wp-1)
				if     (ABS(y_wp_up(k_wp,j_wp)-y_wp_up(k_wp-1,j_wpp)) .LT. (0.01d0)) then
					wdistance_ = dist(Zs(k),Ys(j),z_wp(k_wp),y_wp_up(k_wp,j_wp),z_wp(k_wp-1),y_wp_up(k_wp-1,j_wpp))
					if ((wdistance_ .LE. max_dist) .AND. (wdistance_ .LT. wdistance)) then
						wdistance = wdistance_
						hyp  = SQRT((y_wp_up(k_wp,j_wp) - y_wp_up(k_wp-1,j_wpp))**2 + (z_wp(k_wp) - z_wp(k_wp-1))**2)
						wsin_ = (y_wp_up(k_wp,j_wp) - y_wp_up(k_wp-1,j_wpp))/hyp
						wcos_ = (z_wp(k_wp) - z_wp(k_wp-1))/hyp
					endif
				endif
			enddo
			enddo
		
			!$acc loop seq
			do j_wp  = 1, wp_low(k_wp)
			!$acc loop seq
			do j_wpp = 1, wp_low(k_wp-1)
				if     (ABS(y_wp_low(k_wp,j_wp)-y_wp_low(k_wp-1,j_wpp)) .LT. (0.01d0)) then
					wdistance_ = dist(Zs(k),Ys(j),z_wp(k_wp),y_wp_low(k_wp,j_wp),z_wp(k_wp-1),y_wp_low(k_wp-1,j_wpp))
					if ((wdistance_ .LE. max_dist) .AND. (wdistance_ .LT. wdistance)) then
						wdistance = wdistance_
						hyp  = SQRT((y_wp_low(k_wp,j_wp) - y_wp_low(k_wp-1,j_wpp))**2 + (z_wp(k_wp) - z_wp(k_wp-1))**2)
						wsin_ = (y_wp_low(k_wp,j_wp) - y_wp_low(k_wp-1,j_wpp))/hyp
						wcos_ = (z_wp(k_wp) - z_wp(k_wp-1))/hyp
					endif
				endif
			enddo
			enddo
		enddo
	endif	
			wdist(j,k) = SQRT(wdistance)	
			wsin(j,k) = wsin_
			wcos(j,k) = wcos_
	enddo
enddo
!$acc end parallel
!$OMP END PARALLEL DO

!$acc end data 


end subroutine wall_distance



! This subroutine is used to apply a wall function and/or a damping function. (BETA version)
subroutine wall_mod()
use variables
implicit none

   open (83,file='wall_yplus.dat',position='append')
   open (84,file='wall_uplus_instant.dat',position='append')
   open (85,file='wall_info_instant.dat',position='append')
   open (86,file='wall_nut_instant.dat',position='append')
   open (87,file='wall_yplus_instant.dat',position='append')
   open (88,file='wall_wdist_instant.dat',position='append')
   open (89,file='wall_moddamp_instant.dat',position='append')

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------

!$acc data present(ETA,wdist,wsin,wcos,fdamp(:,istart:iend),ypluss(:,istart:iend),upluss(:,istart:iend), &
!$acc              v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2),nut(:,:,istart:iend))


! Calculate friction velocity uTau using the Newton-Raphson iteration method, calculate nut, modify nut

countYplus = 0 ; yplusss = 0.d0 ; yplusMAX = 0.d0 ; yplusMIN = 300.d0

!$OMP PARALLEL DO PRIVATE(j,Ucp,Uwp,Urel,uTau,uTau_,NRit,NRerr,kUu,fkUu,NR_f,NR_df) reduction(+: yplusss,countYplus) reduction(MAX: yplusMAX) reduction(MIN: yplusMIN)
!$acc parallel loop independent private(j,Ucp,Uwp,Urel,uTau,uTau_,NRit,NRerr,kUu,fkUu,NR_f,NR_df) reduction(+: yplusss,countYplus) reduction(MAX: yplusMAX) reduction(MIN: yplusMIN) collapse(2) gang vector
do k=istart,iend
	do j=jBgnVOS,jEndVOS
		if(ETA(1,j,k) .GE. 0.5d0) then
			nut(1,j,k) = 0.d0
		elseif(wdist(j,k) .LE. SQRT(max_dist)) then 					! Damping
			Ucp = 0.5d0*( w(INT(nx/2),j,k)+w(INT(nx/2),j,k-1) )*wcos(j,k)+0.5d0*( v(INT(nx/2),j,k)+v(INT(nx/2),j-1,k) )*wsin(j,k)	! Approximanted using the values at span center (need a further check)
			Uwp = 0.d0 																				!(for stationary wall)
			Urel = ABS(Ucp - Uwp)
	
			uTau = SQRT((nut(1,j,k) + nu) * Urel/wdist(j,k))
			if(uTau .GT. ROOTVSMALL) then			! to avoid a division by zero, and wall function is not needed when Urel = 0
				NRit=0 ; NRerr=1.d0
				do while (NRtolerance < NRerr .AND. NRit < NRmaxIter)					
					kUu = MIN(Kappa*Urel/uTau,50.d0)
					fkUu = EXP(kUu) - 1.d0 - kUu*(1.d0+0.5d0*kUu)
					
					NR_f = - uTau*wdist(j,k)/nu + Urel/uTau + invE*(fkUu-kUu*kUu*kUu/6.d0)
					NR_df= wdist(j,k)/nu + Urel/(uTau*uTau) + invE*kUu*fkUu/uTau
				
					uTau_ = uTau + NR_f/NR_df
					
					NRerr = ABS((uTau - uTau_)/uTau)					
					NRit=NRit+1
					uTau = uTau_
			
				end do
				uTau = MAX(0.d0,uTau)
				ypluss(j,k) = uTau*wdist(j,k)/nu
				upluss(j,k) = Urel/uTau
				
				if(wall_model == 1 .AND. wdist(j,k) .LE. SQRT(dzSml**2+dySml**2)) then			! Applied only on the nearest cells
					nut(1,j,k) = MAX(0.d0,(uTau*uTau/(Urel/wdist(j,k) + ROOTVSMALL)) - nu)
				elseif(Damping_F == 1 .AND. ypluss(j,k) .LT. 300.d0) then
					fdamp(j,k) = 1.d0 - dexp(-ypluss(j,k)/25.d0)
					nut(1,j,k) = fdamp(j,k)**2*nut(1,j,k)				! Damping
				endif

				if(wdist(j,k) .LE. SQRT(dzSml**2+dySml**2)) then
					yplusss = yplusss + ypluss(j,k)
					yplusMAX = MAX(yplusMAX,ypluss(j,k))
					yplusMIN = MIN(yplusMIN,ypluss(j,k))
					countYplus = countYplus+1
				endif

			endif	

			!if(istep == INT(1.d0/dt)) then
			!        write(84,'(10(3x,F15.10),3x,I4)') Ys(j),Zs(k),wdist(j,k),ypluss(j,k),upluss(j,k),wcos(j,k),Urel,uTau,kUu,NRerr,NRit
			
			!endif
			
		endif
    end do
end do
!$acc end parallel
!$OMP END PARALLEL DO
 
call MPI_ALLREDUCE( yplusss, yplusss_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
call MPI_ALLREDUCE( yplusMAX, yplusMAX_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr )
call MPI_ALLREDUCE( yplusMIN, yplusMIN_, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr )
call MPI_ALLREDUCE( countYplus, countYplus_, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
 
 
! Copying to spanwise cells
!$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
!$acc parallel loop independent collapse(3)
do k=istart,iend
do j=jBgnVOS,jEndVOS
do i=2,nx
	nut(i,j,k) = nut(1,j,k)
end do
end do
end do
!$acc end parallel
!$OMP END PARALLEL DO


 if(myid==master) then
 write(83,'(4(3x,F15.7))') istep*dt,yplusMIN_,yplusMAX_,yplusss_/countYplus_*1.d0
 endif
 

! Transfering data back to master for write out
if (istep == INT(100.d0/dt)) then     
	!$acc update self(wdist,fdamp(:,istart:iend),ypluss(:,istart:iend),upluss(:,istart:iend),nut(:,:,istart:iend))
          
      !>>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<<
      if(myid>master)then
         itag = 501
         call MPI_SEND( ypluss(-1,istart), igcount*(ny+4), MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 502
         call MPI_SEND( upluss(-1,istart), igcount*(ny+4), MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 503
         call MPI_SEND( fdamp(-1,istart), igcount*(ny+4), MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 504
         call MPI_SEND( nut(-1,-1,istart), igcount*(nx+4)*(ny+4), MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if

      if(myid==master)then
         do i = 1, (nproc-1)
            itag = 501
            call MPI_RECV( ypluss(-1,gstart(i)), igcount*(ny+4), MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 502
            call MPI_RECV( upluss(-1,gstart(i)), igcount*(ny+4), MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 503
            call MPI_RECV( fdamp(-1,gstart(i)), igcount*(ny+4), MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 504
            call MPI_RECV( nut(-1,-1,gstart(i)), igcount*(nx+4)*(ny+4), MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

           
      if(myid==master)then
		do k=kBgnVOS,kEndVOS
			do j=jBgnVOS,jEndVOS
				write(84,'(2(3x,F15.10))') ypluss(j,k),upluss(j,k)
				write(85,'(5(3x,F15.10))') Ys(j),Zs(k),wdist(j,k),ypluss(j,k),nut(1,j,k)
				write(86,'(3(3x,F15.10))') Ys(j),Zs(k),nut(1,j,k)
				write(87,'(3(3x,F15.10))') Ys(j),Zs(k),ypluss(j,k)
				write(88,'(3(3x,F15.10))') Ys(j),Zs(k),wdist(j,k)
				write(89,'(3(3x,F15.10))') Ys(j),Zs(k),fdamp(j,k)
			end do
		end do
	  endif
endif


!-------------------------------------------------------------------------------------------------------------------------------------------------------------------		

    close(83)
    close(84)
    close(85)
    close(86)
    close(87)
    close(88)
	close(89)

!$acc end data 
 
end subroutine wall_mod




function dist(x_w, y_w, x1_w, y1_w, x2_w, y2_w)
    implicit none
!$acc routine

    real*8, intent(in) :: x_w, y_w, x1_w, y1_w, x2_w, y2_w
    real*8 :: a_w, b_w, c_w, d_w, lenSq, param, dot, xx_w, yy_w, dx_w, dy_w, dist

    a_w = x_w - x1_w
    b_w = y_w - y1_w
    c_w = x2_w - x1_w
    d_w = y2_w - y1_w

    lenSq = c_w * c_w + d_w * d_w
    if (lenSq /= 0.0d0) then
        dot = a_w * c_w + b_w * d_w
        param = dot / lenSq
    else
        param = -1.0d0
    end if

    if (param < 0.0d0) then
        xx_w = x1_w
        yy_w = y1_w
    else if (param > 1.0d0) then
        xx_w = x2_w
        yy_w = y2_w
    else
        xx_w = x1_w + param * c_w
        yy_w = y1_w + param * d_w
    end if

    dx_w = x_w - xx_w
    dy_w = y_w - yy_w

    dist = dx_w**2 + dy_w**2
	
	RETURN
	
end function dist