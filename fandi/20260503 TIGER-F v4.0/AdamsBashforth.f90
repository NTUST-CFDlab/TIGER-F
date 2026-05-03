! 03 May 2026 by FDS

subroutine AdamsBashforth()
    use variables
    implicit none
    real*8      :: r_dt1, r_dt2, c2_a, c2_b, c3_a, c3_b, c3_c


!$acc data present(u(:,:,istart-2:iend+2),      v(:,:,istart-2:iend+2),      w(:,:,istart-2:iend+2), &
!$acc	          u0(:,:,istart-2:iend+2),     v0(:,:,istart-2:iend+2),     w0(:,:,istart-2:iend+2), &
!$acc	      u_star(:,:,istart-2:iend+2), v_star(:,:,istart-2:iend+2), w_star(:,:,istart-2:iend+2), &
!$acc	     u_star2(:,:,istart-2:iend+2),v_star2(:,:,istart-2:iend+2),w_star2(:,:,istart-2:iend+2), &
!$acc	     u_star1(:,:,istart-2:iend+2),v_star1(:,:,istart-2:iend+2),w_star1(:,:,istart-2:iend+2))

	r_dt1 = dt_1/dt			!2nd order
	r_dt2 = dt_2/dt			!3rd order


    if( istep >= 3) then	

		c3_a  = (2.d0 + 6.d0*r_dt1 + 3.d0*r_dt2)/(6.d0*r_dt1*(r_dt1+r_dt2)) + 1.d0	!3rd order
		c3_b  = - (2.d0 + 3.d0*r_dt1 + 3.d0*r_dt2)/(6.d0*r_dt1*r_dt2)				!3rd order
		c3_c  = (2.d0 + 3.d0*r_dt1)/(6.d0*r_dt2*(r_dt1+r_dt2))						!3rd order

		!$OMP PARALLEL DO COLLAPSE(3) 
		!$acc parallel loop independent collapse(3) gang vector 
		do k=istart,iend 
		do j=1,ny; do i=1,nx

			u0(i,j,k) = u(i,j,k) + ( c3_a * u_star(i,j,k) + c3_b * u_star1(i,j,k) + c3_c * u_star2(i,j,k) ) * dt
            v0(i,j,k) = v(i,j,k) + ( c3_a * v_star(i,j,k) + c3_b * v_star1(i,j,k) + c3_c * v_star2(i,j,k) ) * dt 
            w0(i,j,k) = w(i,j,k) + ( c3_a * w_star(i,j,k) + c3_b * w_star1(i,j,k) + c3_c * w_star2(i,j,k) ) * dt

			u_star2(i,j,k) = u_star1(i,j,k)		!3rd order
			u_star1(i,j,k) = u_star(i,j,k)		!2nd order

			v_star2(i,j,k) = v_star1(i,j,k)		!3rd order
			v_star1(i,j,k) = v_star(i,j,k)		!2nd order

			w_star2(i,j,k) = w_star1(i,j,k)		!3rd order
			w_star1(i,j,k) = w_star(i,j,k)		!2nd order

		enddo; enddo; enddo
		!$acc end parallel
		!$OMP END PARALLEL DO


    else if( istep == 2) then		! else if( istep >= 2) then

		c2_a  = 1.d0 + 0.5d0/r_dt1		!2nd order
		c2_b  = -0.5d0/r_dt1			!2nd order

		!$OMP PARALLEL DO COLLAPSE(3) 
		!$acc parallel loop independent collapse(3) gang vector 
		do k=istart,iend 
		do j=1,ny; do i=1,nx
		
            u0(i,j,k) = u(i,j,k) + ( c2_a * u_star(i,j,k) + c2_b * u_star1(i,j,k) ) * dt
            v0(i,j,k) = v(i,j,k) + ( c2_a * v_star(i,j,k) + c2_b * v_star1(i,j,k) ) * dt
            w0(i,j,k) = w(i,j,k) + ( c2_a * w_star(i,j,k) + c2_b * w_star1(i,j,k) ) * dt

			u_star2(i,j,k) = u_star1(i,j,k)		!3rd order
			u_star1(i,j,k) = u_star(i,j,k)		!2nd order

			v_star2(i,j,k) = v_star1(i,j,k)		!3rd order
			v_star1(i,j,k) = v_star(i,j,k)		!2nd order

			w_star2(i,j,k) = w_star1(i,j,k)		!3rd order
			w_star1(i,j,k) = w_star(i,j,k)		!2nd order

		enddo; enddo; enddo
		!$acc end parallel
		!$OMP END PARALLEL DO


    else if( istep == 1) then

		!$OMP PARALLEL DO COLLAPSE(3) 
		!$acc parallel loop independent collapse(3) gang vector 
		do k=istart,iend 
		do j=1,ny; do i=1,nx

            u0(i,j,k) = u(i,j,k) + u_star(i,j,k) * dt
            v0(i,j,k) = v(i,j,k) + v_star(i,j,k) * dt
            w0(i,j,k) = w(i,j,k) + w_star(i,j,k) * dt

			u_star1(i,j,k) = u_star(i,j,k)		!2nd order

			v_star1(i,j,k) = v_star(i,j,k)		!2nd order

			w_star1(i,j,k) = w_star(i,j,k)		!2nd order

		enddo; enddo; enddo
		!$acc end parallel
		!$OMP END PARALLEL DO
			

    !else if( istep == 4) then	! DO NOT USE FOR CONSTANT CFL

		! !$OMP PARALLEL DO COLLAPSE(3) 
		! !$acc parallel loop independent collapse(3) gang vector 
		! do k=istart,iend 
		! do j=1,ny; do i=1,nx

        !    u0(i,j,k) = u(i,j,k) + ( 55.0 * u_star(i,j,k) - 59.0 * u_star1(i,j,k) + 37.0 * u_star2(i,j,k) - 9.0 * u_star3(i,j,k) ) * dt / 24.0
        !    v0(i,j,k) = v(i,j,k) + ( 55.0 * v_star(i,j,k) - 59.0 * v_star1(i,j,k) + 37.0 * v_star2(i,j,k) - 9.0 * v_star3(i,j,k) ) * dt / 24.0
        !    w0(i,j,k) = w(i,j,k) + ( 55.0 * w_star(i,j,k) - 59.0 * w_star1(i,j,k) + 37.0 * w_star2(i,j,k) - 9.0 * w_star3(i,j,k) ) * dt / 24.0
		
			!u_star4(i,j,k) = u_star3(i,j,k)
			!u_star3(i,j,k) = u_star2(i,j,k)
			!u_star2(i,j,k) = u_star1(i,j,k)		!3rd order
			!u_star1(i,j,k) = u_star(i,j,k)			!2nd order

			!v_star4(i,j,k) = v_star3(i,j,k)
			!v_star3(i,j,k) = v_star2(i,j,k)
			!v_star2(i,j,k) = v_star1(i,j,k)		!3rd order
			!v_star1(i,j,k) = v_star(i,j,k)			!2nd order

			!w_star4(i,j,k) = w_star3(i,j,k)
			!w_star3(i,j,k) = w_star2(i,j,k)
			!w_star2(i,j,k) = w_star1(i,j,k)		!3rd order
			!w_star1(i,j,k) = w_star(i,j,k)			!2nd order		

		! enddo; enddo; enddo
		! !$acc end parallel
		! !$OMP END PARALLEL DO
		
!
    !else if( istep == 5) then

		! !$OMP PARALLEL DO COLLAPSE(3) 
		! !$acc parallel loop independent collapse(3) gang vector 
		! do k=istart,iend 
		! do j=1,ny; do i=1,nx

        !    u0(i,j,k) = u(i,j,k) + ( 1901.0 * u_star(i,j,k) - 2774.0 * u_star1(i,j,k) + 2616.0 * u_star2(i,j,k) - 1274.0 * u_star3(i,j,k) + 251.0 * u_star4(i,j,k) ) * dt / 720.0
        !    v0(i,j,k) = v(i,j,k) + ( 1901.0 * v_star(i,j,k) - 2774.0 * v_star1(i,j,k) + 2616.0 * v_star2(i,j,k) - 1274.0 * v_star3(i,j,k) + 251.0 * v_star4(i,j,k) ) * dt / 720.0
        !    w0(i,j,k) = w(i,j,k) + ( 1901.0 * w_star(i,j,k) - 2774.0 * w_star1(i,j,k) + 2616.0 * w_star2(i,j,k) - 1274.0 * w_star3(i,j,k) + 251.0 * w_star4(i,j,k) ) * dt / 720.0

			!u_star4(i,j,k) = u_star3(i,j,k)
			!u_star3(i,j,k) = u_star2(i,j,k)
			!u_star2(i,j,k) = u_star1(i,j,k)		!3rd order
			!u_star1(i,j,k) = u_star(i,j,k)			!2nd order

			!v_star4(i,j,k) = v_star3(i,j,k)
			!v_star3(i,j,k) = v_star2(i,j,k)
			!v_star2(i,j,k) = v_star1(i,j,k)		!3rd order
			!v_star1(i,j,k) = v_star(i,j,k)			!2nd order

			!w_star4(i,j,k) = w_star3(i,j,k)
			!w_star3(i,j,k) = w_star2(i,j,k)
			!w_star2(i,j,k) = w_star1(i,j,k)		!3rd order
			!w_star1(i,j,k) = w_star(i,j,k)			!2nd order		

		! enddo; enddo; enddo
		! !$acc end parallel
		! !$OMP END PARALLEL DO	
  
    end if


!$acc end data	

	dt_2  = dt_1		!3rd order
	dt_1  = dt			!2nd order
	

end subroutine AdamsBashforth

