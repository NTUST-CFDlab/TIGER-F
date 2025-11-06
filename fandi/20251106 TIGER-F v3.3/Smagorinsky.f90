! 02 Nov 2025 - FDS



subroutine CalculateSmagorinskyViscosity()
    use variables 
    implicit none

!$acc data present(iDx,iDy,iDz,Dxs,Dys,Dzs,ETA,nut(:,:,istart:iend+1), &
!$acc	u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2))

    !$OMP PARALLEL DO               PRIVATE(strainRateTensor,delta,tau,ufric,Yplus,fwall) COLLAPSE(3)
    !$acc parallel loop independent private(strainRateTensor,delta,tau,ufric,Yplus,fwall) collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx

	!Calculation of Strain rate (Uses Central Diffrence Scheme for gradient)--------------------------
    strainRateTensor =										&
			  ( 0.5d0*(u(i-1,j,k)-u(i+1,j,k))/Dxs(i) )**2 	&
			+ ( 0.5d0*(v(i,j-1,k)-v(i,j+1,k))/Dys(j) )**2 	&
			+ ( 0.5d0*(w(i,j,k-1)-w(i,j,k+1))/Dzs(k) )**2 	&
			+ 0.5d0*( 0.5d0*(u(i,j-1,k)-u(i,j+1,k))/iDy(j) + 0.5d0*(v(i-1,j,k)-v(i+1,j,k))/iDx(i) )**2				&
            + 0.5d0*( 0.5d0*(u(i,j,k-1)-u(i,j,k+1))/iDz(k) + 0.5d0*(w(i-1,j,k)-w(i+1,j,k))/iDx(i) )**2				&
            + 0.5d0*( 0.5d0*(v(i,j,k-1)-v(i,j,k+1))/iDz(k) + 0.5d0*(w(i,j-1,k)-w(i,j+1,k))/iDy(j) )**2

    delta = ( iDx(i)*iDy(j)*iDz(k) ) ** (1.d0/3.d0)       

! Applying damping function to solid boundary cells	

		if (0.d0 < ETA(i,j,k) .AND. ETA(i,j,k) < 1.d0) then
			tau = 0.5d0*Cf*den_flu * 0.25d0*((u(i,j,k)+u(i-1,j,k))**2 + (v(i,j,k)+v(i,j-1,k))**2 + (w(i,j,k)+w(i,j,k-1))**2)
			ufric = SQRT(tau/den_flu)
			Yplus =  dySml*ufric/nu
			fwall = (1.d0 - EXP(-Yplus/25.d0))

			nut(i,j,k) = (Cs*delta*fwall)**2*SQRT(2.d0*strainRateTensor) 
		else
			nut(i,j,k) = (Cs*delta)**2*SQRT(2.d0*strainRateTensor)   
		end if 


    end do; end do; end do
    !$acc end parallel
    !$OMP END PARALLEL DO
!$acc end data


end subroutine CalculateSmagorinskyViscosity



! subroutine CalculateSmagorinskyViscosity()
    ! use variables 
    ! implicit none
	! real*8		:: dPhidX

! !$acc data present(iDx,iDy,iDz,Dxs,Dys,Dzs,ETA,nut(:,:,istart:iend+1), &
! !$acc	u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2))

    ! !$OMP PARALLEL DO               PRIVATE(strainRateTensor,delta,tau,ufric,Yplus,fwall) COLLAPSE(3)
    ! !$acc parallel loop independent private(strainRateTensor,delta,tau,ufric,Yplus,fwall) collapse(3) gang vector
    ! do k=istart,iend
    ! do j=1,ny; do i=1,nx

    ! strainRateTensor =												&
			  ! dPhidX(u(i-1,j,k),u(i+1,j,k),Dxs(i))**2 	&
			! + dPhidX(v(i,j-1,k),v(i,j+1,k),Dys(j))**2 	&
			! + dPhidX(w(i,j,k-1),w(i,j,k+1),Dzs(k))**2 	&
			! + 0.5d0*(dPhidX(u(i,j-1,k),u(i,j+1,k),iDy(j))+dPhidX(v(i-1,j,k),v(i+1,j,k),iDx(i)) )**2				&
            ! + 0.5d0*(dPhidX(u(i,j,k-1),u(i,j,k+1),iDz(k))+dPhidX(w(i-1,j,k),w(i+1,j,k),iDx(i)) )**2				&
            ! + 0.5d0*(dPhidX(v(i,j,k-1),v(i,j,k+1),iDz(k))+dPhidX(w(i,j-1,k),w(i,j+1,k),iDy(j)) )**2

    ! delta = ( iDx(i)*iDy(j)*iDz(k) ) ** (1.d0/3.d0)       

! ! Applying damping function to solid boundary cells	

		! if (0.d0 < ETA(i,j,k) .AND. ETA(i,j,k) < 1.d0) then
			! tau = 0.5d0*Cf*den_flu * 0.25d0*((u(i,j,k)+u(i-1,j,k))**2 + (v(i,j,k)+v(i,j-1,k))**2 + (w(i,j,k)+w(i,j,k-1))**2)
			! ufric = SQRT(tau/den_flu)
			! Yplus =  dySml*ufric/nu
			! fwall = (1.d0 - EXP(-Yplus/25.d0))

			! nut(i,j,k) = (Cs*delta*fwall)**2*SQRT(2.d0*strainRateTensor) 
		! else
			! nut(i,j,k) = (Cs*delta)**2*SQRT(2.d0*strainRateTensor)   
		! end if 


    ! end do; end do; end do
    ! !$acc end parallel
    ! !$OMP END PARALLEL DO
! !$acc end data


! end subroutine CalculateSmagorinskyViscosity



! function dPhidX(Phi_1,Phi_2,delta_x)
    ! implicit none
! !$acc routine

    ! real*8, intent(in) :: Phi_1,Phi_2,delta_x
    ! real*8 :: dPhidX

	! !----------- Central Diffrence Scheme--------------------------------
	! dPhidX  = 0.5d0*(Phi_1 - Phi_2) / delta_x
	
	! RETURN
	
! end function dPhidX