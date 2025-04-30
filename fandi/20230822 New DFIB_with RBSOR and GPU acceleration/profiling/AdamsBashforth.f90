! 22 Aug 2023 - FDS
subroutine AdamsBashforth()
    use variables
    implicit none


    !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps) 
    do k=istart,iend 
    do j=1,ny; do i=1,nx

        if( istep == 1) then

            u0(i,j,k) = u(i,j,k) + u_star(i,j,k) 
            v0(i,j,k) = v(i,j,k) + v_star(i,j,k)
            w0(i,j,k) = w(i,j,k) + w_star(i,j,k)
        
        else if( istep >= 2) then

            u0(i,j,k) = u(i,j,k) + ( 1.5d0 * u_star(i,j,k) - 0.5d0 * u_star1(i,j,k) )
            v0(i,j,k) = v(i,j,k) + ( 1.5d0 * v_star(i,j,k) - 0.5d0 * v_star1(i,j,k) )
            w0(i,j,k) = w(i,j,k) + ( 1.5d0 * w_star(i,j,k) - 0.5d0 * w_star1(i,j,k) )
            
        !else if( istep >= 3) then
!
		!	u0(i,j,k) = u(i,j,k) + ( 23.d0 * u_star(i,j,k) - 16.d0 * u_star1(i,j,k) + 5.d0*u_star2(i,j,k) ) / 12.d0 
        !    v0(i,j,k) = v(i,j,k) + ( 23.d0 * v_star(i,j,k) - 16.d0 * v_star1(i,j,k) + 5.d0*v_star2(i,j,k) ) / 12.d0 
        !    w0(i,j,k) = w(i,j,k) + ( 23.d0 * w_star(i,j,k) - 16.d0 * w_star1(i,j,k) + 5.d0*w_star2(i,j,k) ) / 12.d0 
!
        !else if( istep == 4) then
!
        !    u0(i,j,k) = u(i,j,k) + ( 55.0 * u_star(i,j,k) - 59.0 * u_star1(i,j,k) + 37.0 * u_star2(i,j,k) - 9.0 * u_star3(i,j,k) ) / 24.0
        !    v0(i,j,k) = v(i,j,k) + ( 55.0 * v_star(i,j,k) - 59.0 * v_star1(i,j,k) + 37.0 * v_star2(i,j,k) - 9.0 * v_star3(i,j,k) ) / 24.0
        !    w0(i,j,k) = w(i,j,k) + ( 55.0 * w_star(i,j,k) - 59.0 * w_star1(i,j,k) + 37.0 * w_star2(i,j,k) - 9.0 * w_star3(i,j,k) ) / 24.0
!
        !else if( istep >= 5) then
!
        !    u0(i,j,k) = u(i,j,k) + ( 1901.0 * u_star(i,j,k) - 2774.0 * u_star1(i,j,k) + 2616.0 * u_star2(i,j,k) - 1274.0 * u_star3(i,j,k) + 251.0 * u_star4(i,j,k) ) / 720.0
        !    v0(i,j,k) = v(i,j,k) + ( 1901.0 * v_star(i,j,k) - 2774.0 * v_star1(i,j,k) + 2616.0 * v_star2(i,j,k) - 1274.0 * v_star3(i,j,k) + 251.0 * v_star4(i,j,k) ) / 720.0
        !    w0(i,j,k) = w(i,j,k) + ( 1901.0 * w_star(i,j,k) - 2774.0 * w_star1(i,j,k) + 2616.0 * w_star2(i,j,k) - 1274.0 * w_star3(i,j,k) + 251.0 * w_star4(i,j,k) ) / 720.0
            
        end if


        !u_star4(i,j,k) = u_star3(i,j,k)
        !u_star3(i,j,k) = u_star2(i,j,k)
        !u_star2(i,j,k) = u_star1(i,j,k)
        u_star1(i,j,k) = u_star(i,j,k)


        !v_star4(i,j,k) = v_star3(i,j,k)
        !v_star3(i,j,k) = v_star2(i,j,k)
        !v_star2(i,j,k) = v_star1(i,j,k)
        v_star1(i,j,k) = v_star(i,j,k)


        !w_star4(i,j,k) = w_star3(i,j,k)
        !w_star3(i,j,k) = w_star2(i,j,k)
        !w_star2(i,j,k) = w_star1(i,j,k)
        w_star1(i,j,k) = w_star(i,j,k)

    enddo; enddo; enddo
    !$OMP END PARALLEL DO
	

end subroutine AdamsBashforth