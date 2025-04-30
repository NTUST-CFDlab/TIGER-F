subroutine initial_conditions()
use variables
implicit none

!---------------------------------------------------!
!         Initial conditions calculation            !
!---------------------------------------------------!

!$OMP PARALLEL DO PRIVATE(j,i)  
do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2

   !u(i,j,k) = 0.0
   !v(i,j,k) = 0.0
   w(i,j,k) = 1.0
   !u1(i,j,k) = u(i,j,k)
   !v1(i,j,k) = v(i,j,k)
   !w1(i,j,k) = w(i,j,k)
   !u2(i,j,k) = u(i,j,k)
   !v2(i,j,k) = v(i,j,k)
   !w2(i,j,k) = w(i,j,k)
   !p(i,j,k) = 0.0
   !u_star(i,j,k) = 0.0
   !v_star(i,j,k) = 0.0
   !w_star(i,j,k) = 0.0
   !ETA(i,j,k) = 0.0
   F_tavex(i,j,k) = 0.0
   F_tavey(i,j,k) = 0.0

end do; end do; end do
!$OMP END PARALLEL DO


do k=1,nz
!$OMP PARALLEL DO PRIVATE(i)  
do j=1,ny; do i=1,nx

   Den(i,j,k) = ( iDy(j) * iDz(k) / Dxs(i) + iDy(j) * iDz(k) / Dxs(i-1) &
                + iDx(i) * iDz(k) / Dys(j) + iDx(i) * iDz(k) / Dys(j-1) &
                + iDx(i) * iDy(j) / Dzs(k) + iDx(i) * iDy(j) / Dzs(k-1) )


   P_Den(i,j,k,1) = iDy(j) * iDz(k) / Dxs(i) 
   P_Den(i,j,k,2) = iDy(j) * iDz(k) / Dxs(i-1)
   P_Den(i,j,k,3) = iDx(i) * iDz(k) / Dys(j) 
   P_Den(i,j,k,4) = iDx(i) * iDz(k) / Dys(j-1)
   P_Den(i,j,k,5) = iDx(i) * iDy(j) / Dzs(k) 
   P_Den(i,j,k,6) = iDx(i) * iDy(j) / Dzs(k-1)

enddo; enddo
!$OMP END PARALLEL DO
enddo

open( 24,file = NACA_filename, form = 'FORMATTED' )
    do i = 1, poly
        read(24,*) iaz(i), iay(i)
    end do
close(24)


end subroutine initial_conditions



