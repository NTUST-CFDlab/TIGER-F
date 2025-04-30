! 22 Aug 2023 - FDS
subroutine CalculateSmagorinskyViscosity()
    use variables 
    implicit none

   !$OMP PARALLEL
    call GradientPhiGauss(u,dudx,Dxs,iDy,iDz)
    call GradientPhiGauss(v,dvdx,iDx,Dys,iDz)
    call GradientPhiGauss(w,dwdx,iDx,iDy,Dzs)


    !$OMP DO PRIVATE(i,j,mutsgs,delta) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx

    mutsgs = abs(&
                 2*dudx(i,j,k,1)**2 + 2*dvdx(i,j,k,2)**2 + 2*dwdx(i,j,k,3)**2 + &
                 ( dudx(i,j,k,2)+dvdx(i,j,k,1) )**2 + &
                 ( dudx(i,j,k,3)+dwdx(i,j,k,1) )**2 + &
                 ( dvdx(i,j,k,3)+dwdx(i,j,k,2) )**2         )

    delta = ( iDx(i)*iDy(j)*iDz(k) ) ** (1./3.)       

! Applying damping function to solid boundary cells	

	  if ((0.5d0*(ETA(i,j,k)+ETA(i,j+1,k))) > 0.0d0 .AND. (0.5d0*(ETA(i,j,k)+ETA(i,j+1,k))) < 1.d0) then

      mutsgs = (Cs*delta*fwall)**2*sqrt(mutsgs) 

      else

      mutsgs = (Cs*delta)**2*sqrt(mutsgs)    
    
	  end if 

	Viseff(i,j,k) = mutsgs


    end do; end do; end do
    !$OMP END DO    
   !$OMP END PARALLEL 


end subroutine CalculateSmagorinskyViscosity




subroutine GradientPhiGauss(Phi,dPhidX,Dxm,Dym,Dzm)
    use variables 
    implicit none

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: Phi
    real*8, dimension(nx,ny,nz,3) :: dPhidX
    real*8 ,dimension (-1:nx+2) :: Dxm, Dym, Dzm 

    !--------Local variables--------!
    !real*8 :: fact
    !real*8 :: Phiface_w, Phiface_e, Phiface_n, Phiface_s, Phiface_b, Phiface_f
    !--------Local variables--------!


    !$OMP DO PRIVATE(i,j) collapse(nclps)
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
    !$OMP END DO


end subroutine GradientPhiGauss