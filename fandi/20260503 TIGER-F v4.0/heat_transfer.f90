! 03 May 2026 by FDS


subroutine discretisation_QUICK_HT()
  use variables
  implicit none


!$acc data present(iDx,iDy,iDz,Dxs,Dys,Dzs,nut(:,:,istart:iend+1), &
!$acc	     u(:,:,istart-2:iend+2), v(:,:,istart-2:iend+2), w(:,:,istart-2:iend+2), &
!$acc	theta(:,:,istart-2:iend+2),thetastar(:,:,istart-2:iend+2))


   !$OMP PARALLEL DO      PRIVATE(tte,ttw,ttn,tts,ttf,ttb,alpha) COLLAPSE(3)
   !$acc parallel
   !$acc loop independent private(tte,ttw,ttn,tts,ttf,ttb,alpha) collapse(3) gang vector
   do k=istart,iend
   do j=1,ny; do i=1,nx

        if (u(i,j,k) > 0) then 

        tte = 0.5d0*(theta(i,j,k)     + theta(i+1,j,k)  ) -0.125d0*Dxs(i)*Dxs(i)/iDx(i)* &
				 (  (theta(i+1,j,k)   - theta(i,j,k)    ) / Dxs(i) &
				  - (theta(i,j,k)     - theta(i-1,j,k)  ) / Dxs(i)   )
                
        else

        tte = 0.5d0*(theta(i,j,k)     + theta(i+1,j,k)  ) -0.125d0*Dxs(i)*Dxs(i)/iDx(i+1)* &
				 (  (theta(i+2,j,k)   - theta(i+1,j,k)  ) / Dxs(i+1) &
				  - (theta(i+1,j,k)   - theta(i,j,k)    ) / Dxs(i) )

        end if

        if (u(i-1,j,k) > 0) then 

        ttw = 0.5d0*(theta(i-1,j,k)   + theta(i,j,k)    ) -0.125d0*Dxs(i-1)*Dxs(i-1)/iDx(i-1)* &
				 (  (theta(i,j,k)     - theta(i-1,j,k)  ) / Dxs(i-1)   &
				  - (theta(i-1,j,k)   - theta(i-2,j,k)  ) / Dxs(i-2) ) 
                
        else

        ttw = 0.5d0*(theta(i-1,j,k)   + theta(i,j,k)    ) -0.125d0*Dxs(i-1)*Dxs(i-1)/iDx(i)* &
				 (  (theta(i+1,j,k)   - theta(i,j,k)    ) / Dxs(i) &
				  - (theta(i,j,k)     - theta(i-1,j,k)  ) / Dxs(i-1)   ) 

        end if

        if (v(i,j,k) > 0) then 

        ttn = 0.5d0*(theta(i,j,k)     + theta(i,j+1,k)  ) -0.125d0*Dys(j)*Dys(j)/iDy(j)* &
				 (  (theta(i,j+1,k)   - theta(i,j,k)    ) / Dys(j)  &
				  - (theta(i,j,k)     - theta(i,j-1,k)  ) / Dys(j-1))  
                
        else

        ttn = 0.5d0*(theta(i,j,k)     + theta(i,j+1,k) )-0.125d0*Dys(j)*Dys(j)/iDy(j+1)* &
				 (  (theta(i,j+2,k)   - theta(i,j+1,k)  ) / Dys(j+1)&
				  - (theta(i,j+1,k)   - theta(i,j,k)    ) / Dys(j)  )

        end if

        if (v(i,j-1,k) > 0) then 

        tts = 0.5d0*(theta(i,j-1,k)   + theta(i,j,k)    ) -0.125d0*Dys(j-1)*Dys(j-1)/iDy(j-1)* &
				 (  (theta(i,j,k)     - theta(i,j-1,k)  ) / Dys(j-1)&
				  - (theta(i,j-1,k)   - theta(i,j-2,k)  ) / Dys(j-2)) 
                
        else

        tts = 0.5d0*(theta(i,j-1,k)   + theta(i,j,k)    ) -0.125d0*Dys(j-1)*Dys(j-1)/iDy(j)* &
				 (  (theta(i,j+1,k)   - theta(i,j,k)    ) / Dys(j)  &
				  - (theta(i,j,k)     - theta(i,j-1,k)  ) / Dys(j-1))  

        end if
 

       if (w(i,j,k) > 0) then 

        ttf = 0.5d0*(theta(i,j,k)     + theta(i,j,k+1)  ) -0.125d0*Dzs(k)*Dzs(k)/iDz(k)* &
				 (  (theta(i,j,k+1)   - theta(i,j,k)    ) / Dzs(k)  &
				  - (theta(i,j,k)     - theta(i,j,k-1)  ) / Dzs(k-1))  
                
        else

        ttf = 0.5d0*(theta(i,j,k)     + theta(i,j,k+1) )-0.125d0*Dzs(k)*Dzs(k)/iDz(k+1)* &
				 (  (theta(i,j,k+2)   - theta(i,j,k+1)  ) / Dzs(k+1)&
				  - (theta(i,j,k+1)   - theta(i,j,k)    ) / Dzs(k)  )

        end if

        if (w(i,j,k-1) > 0) then 

        ttb = 0.5d0*(theta(i,j,k-1)   + theta(i,j,k)    ) -0.125d0*Dzs(k-1)*Dzs(k-1)/iDz(k-1)* &
				 (  (theta(i,j,k)     - theta(i,j,k-1)  ) / Dzs(k-1)&
				  - (theta(i,j,k-1)   - theta(i,j,k-2)  ) / Dzs(k-2)) 
                
        else

        ttb = 0.5d0*(theta(i,j,k-1)   + theta(i,j,k)    ) -0.125d0*Dzs(k-1)*Dzs(k-1)/iDz(k)* &
				 (  (theta(i,j,k+1)   - theta(i,j,k)    ) / Dzs(k)  &
				  - (theta(i,j,k)     - theta(i,j,k-1)  ) / Dzs(k-1))  

        end if

		alpha = (1.d0-ETA(i,j,k)) * (alpha_f + nut(i,j,k)/Prt) + ETA(i,j,k)*alpha_r*alpha_f

        thetastar(i,j,k) = - (u(i,j,k)*tte-u(i-1,j,k)*ttw) / iDx(i) 	&
						   - (v(i,j,k)*ttn-v(i,j-1,k)*tts) / iDy(j) 	&
						   - (w(i,j,k)*ttf-w(i,j,k-1)*ttb) / iDz(k) 	&

                       + alpha*( (theta(i+1,j,k)-theta(i,j,k)) / Dxs(i) - (theta(i,j,k)-theta(i-1,j,k)) / Dxs(i-1) ) / iDx(i) &
                       + alpha*( (theta(i,j+1,k)-theta(i,j,k)) / Dys(j) - (theta(i,j,k)-theta(i,j-1,k)) / Dys(j-1) ) / iDy(j) &
                       + alpha*( (theta(i,j,k+1)-theta(i,j,k)) / Dzs(k) - (theta(i,j,k)-theta(i,j,k-1)) / Dzs(k-1) ) / iDz(k)
                        
  	end do;end do;end do
	!$acc end parallel   
    !$OMP END PARALLEL DO


!$acc end data

end subroutine discretisation_QUICK_HT


subroutine AdamsBashforth_HT()
    use variables
    implicit none
    real*8      :: r_dt1, r_dt2, c2_a, c2_b


!$acc data present(theta(:,:,istart-2:iend+2),theta0(:,:,istart-2:iend+2),theta1(:,:,istart-2:iend+2), &
!$acc			   thetastar(:,:,istart-2:iend+2),thetastar1(:,:,istart-2:iend+2))

	r_dt1 = dt_1/dt			!2nd order

    if( istep >= 2) then

		c2_a  = 1.d0 + 0.5d0/r_dt1		!2nd order
		c2_b  = -0.5d0/r_dt1			!2nd order

		!$OMP PARALLEL DO COLLAPSE(3) 
		!$acc parallel loop independent collapse(3) gang vector 
		do k=istart,iend 
		do j=1,ny; do i=1,nx
		
            theta0(i,j,k) = theta(i,j,k) + ( c2_a * thetastar(i,j,k) + c2_b * thetastar1(i,j,k) ) * dt
			theta1(i,j,k) = theta0(i,j,k)

			thetastar1(i,j,k) = thetastar(i,j,k)		!2nd order

		enddo; enddo; enddo
		!$acc end parallel
		!$OMP END PARALLEL DO


    else if( istep == 1) then

		!$OMP PARALLEL DO COLLAPSE(3) 
		!$acc parallel loop independent collapse(3) gang vector 
		do k=istart,iend 
		do j=1,ny; do i=1,nx

            theta0(i,j,k) = theta(i,j,k) + thetastar(i,j,k) * dt
			theta1(i,j,k) = theta0(i,j,k)

			thetastar1(i,j,k) = thetastar(i,j,k)		!2nd order

		enddo; enddo; enddo
		!$acc end parallel
		!$OMP END PARALLEL DO

    end if


!$acc end data	

	dt_1  = dt			!2nd order
	
end subroutine AdamsBashforth_HT



subroutine calcul_new_theta()
   use variables
   implicit none
   
! Sending data to accelerator
!$acc data present(ETA,theta0(:,:,istart-2:iend+2),theta1(:,:,istart-2:iend+2),Ftheta(:,:,istart-2:iend+2))

  !$OMP PARALLEL DO   PRIVATE(tsolid1) COLLAPSE(3)
  !$acc parallel loop private(tsolid1) collapse(3) gang vector
   do k=istart,iend 
   do j=1,ny;do i=1,nx
      
	  tsolid1 = tsolid

      theta1(i,j,k) = ETA(i,j,k) * tsolid1 + (1.d0- ETA(i,j,k)) * theta0(i,j,k)
      
      Ftheta(i,j,k) = (theta1(i,j,k) - theta0(i,j,k)) * inv_dt
	  
   end do; end do; end do
   !$acc end parallel 
   !$OMP END PARALLEL DO

!$acc end data

end subroutine calcul_new_theta



subroutine Updating_theta()
   use variables
   implicit none

!$acc data present(theta(:,:,istart-2:iend+2),theta1(:,:,istart-2:iend+2),dtheta_dt(:,:,istart-2:iend+2))


   !$OMP PARALLEL DO COLLAPSE(3) 
   !$acc parallel loop collapse(3) gang vector
   do k=istart,iend  
   do j=1,ny;do i=1,nx

	  dtheta_dt(i,j,k) = (theta1(i,j,k) - theta(i,j,k))  * inv_dt

      theta(i,j,k) = theta1(i,j,k)    		

   end do; end do; end do
   !$acc end parallel
   !$OMP END PARALLEL DO


!$acc end data   
   
end subroutine Updating_theta



subroutine thermalForceIntegrator()
	use variables
	implicit none

	totalFtheta = 0.d0

!$acc data present(iDx,iDy,iDz,Ftheta(:,:,istart-2:iend+2))

   !$OMP PARALLEL DO 			   REDUCTION(+: totalFtheta) COLLAPSE(3)
   !$acc parallel loop independent reduction(+: totalFtheta) collapse(3) gang vector
	do k=istart,iend 
	do j=1,ny;do i=1,nx
	
	  if(ETA(i,j,k) .GT. 0.d0) then
		totalFtheta = totalFtheta + ( Ftheta(i,j,k) * iDx(i)*iDy(j)*iDz(k) )
	  endif
	end do;end do;end do
   !$acc end parallel 
   !$OMP END PARALLEL DO
!$acc end data 

    call MPI_ALLREDUCE( totalFtheta, totalFtheta_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

	Qsolid = - den_flu*cp_f*totalFtheta_		! Solid Cooling/Heating Load

end subroutine thermalForceIntegrator