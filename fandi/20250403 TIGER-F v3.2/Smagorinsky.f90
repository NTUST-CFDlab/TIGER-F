! 12 Jan 2025 - FDS

subroutine CalculateSmagorinskyViscosity()
    use variables 
    implicit none

!$acc data present(iDx,iDy,iDz,Dxs,Dys,Dzs,ETA,nut(:,:,istart:iend+1), &
!$acc	u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2)) &
!$acc create(dudx(:,:,istart:iend,:),dvdx(:,:,istart:iend,:),dwdx(:,:,istart:iend,:))

   !$OMP PARALLEL
    call GradientPhiGauss(u,dudx,Dxs,iDy,iDz)
    call GradientPhiGauss(v,dvdx,iDx,Dys,iDz)
    call GradientPhiGauss(w,dwdx,iDx,iDy,Dzs)

    !$OMP DO PRIVATE(i,j,strainRateTensor,delta,tau,ufric,Yplus,fwall) collapse(nclps)
    !$acc parallel loop independent private(strainRateTensor,delta,tau,ufric,Yplus,fwall) collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx

    strainRateTensor =												&
			dudx(i,j,k,1)**2 + dvdx(i,j,k,2)**2 + dwdx(i,j,k,3)**2	&
			+ 0.5d0*(dudx(i,j,k,2)+dvdx(i,j,k,1) )**2				&
            + 0.5d0*(dudx(i,j,k,3)+dwdx(i,j,k,1) )**2				&
            + 0.5d0*(dvdx(i,j,k,3)+dwdx(i,j,k,2) )**2

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
    !$OMP END DO    
    !$OMP END PARALLEL 
!$acc end data


end subroutine CalculateSmagorinskyViscosity




subroutine GradientPhiGauss(Phi,dPhidX,Dxm,Dym,Dzm)
    use variables 
    implicit none

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: Phi
    real*8, dimension(nx,ny,nz,3) :: dPhidX
    real*8 ,dimension (-1:nx+2) :: Dxm
    real*8 ,dimension (-1:ny+2) :: Dym
    real*8 ,dimension (-1:nz+2) :: Dzm

    !--------Local variables--------!
    !real*8 :: fact
    !real*8 :: Phiface_w, Phiface_e, Phiface_n, Phiface_s, Phiface_b, Phiface_f
    !--------Local variables--------!
!$acc data create(Phi(:,:,istart-2:iend+2),dPhidX(:,:,istart:iend,:),Dxm,Dym,Dzm)

    !$OMP DO PRIVATE(i,j) collapse(nclps)
    !$acc parallel loop independent collapse(3) gang vector
    do k=istart,iend 
    do j=1,ny; do i=1,nx

       ! Phiface_e = 0.5 * ( Phi(i+1,j,k) + Phi(i,j,k) )
       ! Phiface_w = 0.5 * ( Phi(i-1,j,k) + Phi(i,j,k) )

       ! Phiface_n = 0.5 * ( Phi(i,j+1,k) + Phi(i,j,k) )
       ! Phiface_s = 0.5 * ( Phi(i,j-1,k) + Phi(i,j,k) )

       ! Phiface_f = 0.5 * ( Phi(i,j,k+1) + Phi(i,j,k) )
       ! Phiface_b = 0.5 * ( Phi(i,j,k-1) + Phi(i,j,k) )
		
       ! dPhidX(i,j,k,1) = Phiface_w*(-1)*(Dym(j)*Dzm(k)) + Phiface_e*(1)*(Dym(j)*Dzm(k)) 

       ! dPhidX(i,j,k,2) = Phiface_s*(-1)*(Dxm(i)*Dzm(k)) + Phiface_n*(1)*(Dxm(i)*Dzm(k)) 

       ! dPhidX(i,j,k,3) = Phiface_b*(-1)*(Dxm(i)*Dym(j)) + Phiface_f*(1)*(Dxm(i)*Dym(j))

       ! fact = 1.0/(Dxm(i)*Dym(j)*Dzm(k))

       ! dPhidX(i,j,k,1) =  fact * dPhidX(i,j,k,1) 
       ! dPhidX(i,j,k,2) =  fact * dPhidX(i,j,k,2) 
       ! dPhidX(i,j,k,3) =  fact * dPhidX(i,j,k,3) 

        !----------- Central Diffrence Scheme--------------------------------

          dPhidX(i,j,k,1) = 0.5d0/Dxm(i)*(Phi(i-1,j,k) - Phi(i+1,j,k))                   

          dPhidX(i,j,k,2) = 0.5d0/Dym(j)*(Phi(i,j-1,k) - Phi(i,j+1,k)) 

          dPhidX(i,j,k,3) = 0.5d0/Dzm(k)*(Phi(i,j,k-1) - Phi(i,j,k+1)) 

    end do;end do
    end do
	!$acc end parallel
    !$OMP END DO
!$acc end data

end subroutine GradientPhiGauss