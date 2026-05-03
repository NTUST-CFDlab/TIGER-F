! 22 Aug 2023 - FDS
subroutine calcul_new_velocity()
   use variables
   implicit none
	real*8 :: ETA_sub


   !-------------------------------------------------!
   !    Calculation of velocity field at t = dt*n+1  !
   !-------------------------------------------------!

   !!! In x direction !!!

  !$OMP PARALLEL
   !$OMP DO PRIVATE(i,j,ETA_sub,u_solid1) collapse(nclps)
   do k=istart,iend 
   do j=1,ny;do i=1,nx

      u1(i,j,k) = u0(i,j,k) - dt*(p(i+1,j,k)-p(i,j,k))/den_flu / Dxs(i)

		ETA_sub = 0.5d0*(ETA(i,j,k)+ETA(i+1,j,k))

      u_solid1 = u_solid

      u2(i,j,k) = ETA_sub * u_solid1 + (1.d0- ETA_sub) * u1(i,j,k)
      
      FX(i,j,k) = (u2(i,j,k) - u1(i,j,k)) * inv_dt

   end do; end do; end do
   !$OMP END DO


   !!! In y direction !!!


   !$OMP DO PRIVATE(i,j,ETA_sub,v_solid1) collapse(nclps)
   do k=istart,iend
   do j=1,ny;do i=1,nx

      v1(i,j,k) = v0(i,j,k) - dt*(p(i,j+1,k)-p(i,j,k))/den_flu / Dys(j)

		ETA_sub = 0.5d0*(ETA(i,j,k)+ETA(i,j+1,k))


		v_solid1= v_solid

	

      v2(i,j,k) = ETA_sub * v_solid1 + (1.d0- ETA_sub) * v1(i,j,k)

      FY(i,j,k) = ETA_sub * (v_solid1 - v1(i,j,k)) * inv_dt

   end do; end do; end do
   !$OMP END DO
   

   !!! In w direction !!!


   !$OMP DO PRIVATE(i,j,ETA_sub,w_solid1) collapse(nclps)
   do k=istart,iend 
   do j=1,ny;do i=1,nx

      w1(i,j,k) = w0(i,j,k) - dt*(p(i,j,k+1)-p(i,j,k))/den_flu  / Dzs(k)

		ETA_sub = 0.5d0*(ETA(i,j,k)+ETA(i,j,k+1))
		

		w_solid1= w_solid


   
      w2(i,j,k) = ETA_sub * w_solid1 + (1.d0- ETA_sub) * w1(i,j,k)
      
      FZ(i,j,k) = ETA_sub * (w_solid1 - w1(i,j,k)) * inv_dt

   end do;end do
   end do
   !$OMP END DO
  !$OMP END PARALLEL  


   ! !----------data collect among nodes----------!

   ! !Send my results back to the master
   ! if(myid>master)then
	  ! icount = igcount*(nx+4)*(ny+4)
      ! itag = 350
      ! call MPI_SEND( FX(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      ! itag = 360
      ! call MPI_SEND( FY(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      ! itag = 370
      ! call MPI_SEND( FZ(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
   ! end if

   ! if(myid==master)then
      ! do i = 1, (nproc-1)
            ! itag = 350
            ! call MPI_RECV( FX(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            ! itag = 360
            ! call MPI_RECV( FY(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            ! itag = 370
            ! call MPI_RECV( FZ(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )       
      ! end do
   ! end if
   ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   ! !Wait to receive results from each task
	  
   ! !----------data collect among nodes for filer----------!
 

end subroutine calcul_new_velocity


subroutine Updating_velocity()
   use variables
   implicit none
   !---------------------------------------------------------!
   !       loops to update the fluid velocity fields         !
   !---------------------------------------------------------!


   !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps) 
   do k=istart,iend  
   do j=1,ny;do i=1,nx

      u(i,j,k) = u2(i,j,k)    		!!! In x direction !!!
      v(i,j,k) = v2(i,j,k) 			!!! In y direction !!!
      w(i,j,k) = w2(i,j,k) 			!!! In z direction !!!

   end do; end do; end do
   !$OMP END PARALLEL DO
   
   
end subroutine Updating_velocity