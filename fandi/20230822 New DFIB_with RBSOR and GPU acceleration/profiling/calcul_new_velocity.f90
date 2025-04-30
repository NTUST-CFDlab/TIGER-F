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
   !$OMP DO PRIVATE(i,j,ETA_sub,u_solid) collapse(nclps)
   do k=istart,iend 
   do j=1,ny;do i=1,nx

      u1(i,j,k) = u0(i,j,k) - dt*(p(i+1,j,k)-p(i,j,k)) / Dxs(i)

		ETA_sub = 0.5d0*(ETA(i,j,k)+ETA(i+1,j,k))

      u_solid=0.d0

      u2(i,j,k) = ETA_sub * u_solid + (1.d0- ETA_sub) * u1(i,j,k)
      
      FX(i,j,k) = (u2(i,j,k) - u1(i,j,k)) / dt

   end do; end do; end do
   !$OMP END DO


   !!! In y direction !!!


   !$OMP DO PRIVATE(i,j,ETA_sub,v_solid) collapse(nclps)
   do k=istart,iend
   do j=1,ny;do i=1,nx

      v1(i,j,k) = v0(i,j,k) - dt*(p(i,j+1,k)-p(i,j,k)) / Dys(j)

		ETA_sub = 0.5d0*(ETA(i,j,k)+ETA(i,j+1,k))

	! Set the rotor + blade rotation

	   if( SQRT((Zs(k)-zc_t)**2 + (Ys(j)-yc_t)**2) .LT. (blade_r+dzSml)) then
			 v_solid= rotor_omega * (Zs(k)-zc) + blade_omega * (Zs(k)-zc_t) 	! Positive omega is CCW
	   else
			 v_solid= rotor_omega * (Zs(k)-zc)
	   endif


      v2(i,j,k) = ETA_sub * v_solid + (1.d0- ETA_sub) * v1(i,j,k)

      FY(i,j,k) = ETA_sub * (v_solid - v1(i,j,k)) / dt

   end do; end do; end do
   !$OMP END DO
   

   !!! In w direction !!!


   !$OMP DO PRIVATE(i,j,ETA_sub,w_solid) collapse(nclps)
   do k=istart,iend 
   do j=1,ny;do i=1,nx

      w1(i,j,k) = w0(i,j,k) - dt*(p(i,j,k+1)-p(i,j,k)) / Dzs(k)

		ETA_sub = 0.5d0*(ETA(i,j,k)+ETA(i,j,k+1))
		
	! Set the rotor + blade rotation

	   if( SQRT((Zs(k)-zc_t)**2 + (Ys(j)-yc_t)**2) .LT. (blade_r+dySml)) then
			 w_solid= -1.d0 * rotor_omega * (Ys(j)-yc) + (-1.d0 * blade_omega * (Ys(j)-yc_t))	! Positive omega is CCW
	   else
			 w_solid= -1.d0 * rotor_omega * (Ys(j)-yc)
	   endif

      
      w2(i,j,k) = ETA_sub * w_solid + (1.d0- ETA_sub) * w1(i,j,k)
      
      FZ(i,j,k) = ETA_sub * (w_solid - w1(i,j,k)) / dt

   end do;end do
   end do
   !$OMP END DO
  !$OMP END PARALLEL  


   !----------data collect among nodes----------!
   icount = igcount*(nx+4)*(ny+4)

   !Send my results back to the master
   if(myid>master)then

      itag = 350
      call MPI_SEND( FX(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      itag = 360
      call MPI_SEND( FY(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      itag = 370
      call MPI_SEND( FZ(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )


   end if


   !Wait to receive results from each task
   if(myid==master)then
      do i = 1, (nproc-1)

            itag = 350
            call MPI_RECV( FX(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 360
            call MPI_RECV( FY(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 370
            call MPI_RECV( FZ(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
               
      end do
   end if



   !icount = (nx+4)*(ny+4)*(nz+4)
   !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   !call MPI_BCAST( FX(-1,-1,-1), icount, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
   !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   !call MPI_BCAST( FY(-1,-1,-1), icount, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
   !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   !call MPI_BCAST( FZ(-1,-1,-1), icount, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
   !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   !----------data collect among nodes----------!
   

end subroutine calcul_new_velocity


subroutine Updating_velocity()
   use variables
   implicit none
   !---------------------------------------------------------!
   !    loops to update the main u and v velocity fields     !
   !---------------------------------------------------------!


   !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps) 
   do k=istart,iend  
   do j=1,ny;do i=1,nx

      u(i,j,k) = u2(i,j,k)    		!!! In x direction !!!
      v(i,j,k) = v2(i,j,k) 			!!! In y direction !!!
      w(i,j,k) = w2(i,j,k) 			!!! In z direction !!!

   end do; end do; end do
   !$OMP END PARALLEL DO
   
	!!! In x direction !!!
	
!   do k=istart,iend
!   !$OMP PARALLEL DO PRIVATE(i)  
!   do j=1,ny;do i=1,nx

!      u(i,j,k) = u2(i,j,k)    

!   end do;end do
!   !$OMP END PARALLEL DO
!   end do


   !!! In y direction !!!

!   do k=istart,iend
!   !$OMP PARALLEL DO PRIVATE(i)  
!   do j=1,ny;do i=1,nx

!      v(i,j,k) = v2(i,j,k)  

!   end do;end do
!   !$OMP END PARALLEL DO
!   end do


   !!! In w direction !!!

       
!   do k=istart,iend
!   !$OMP PARALLEL DO PRIVATE(i)  
!   do j=1,ny;do i=1,nx

!      w(i,j,k) = w2(i,j,k) 

!   end do;end do
!   !$OMP END PARALLEL DO
!   end do



   !call final_boundary_conditions()


   
   !----------------------------------------------------------------------------------------!
   !           Centering variables to get matrix (1:nx,1:ny,1:nz) for Tecplot               !
   !----------------------------------------------------------------------------------------!

   !!----------data transformation among nodes----------!
   !icount = (nx+4)*(ny+4)
!
   !itag = 380
   !call MPI_SENDRECV( w2(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
   !                   w2(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
   !            
   !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   !!----------data transformation among nodes----------!
!
!
!
!
	!do k=istart,iend
   !!$OMP PARALLEL DO PRIVATE(i)  
   !do j=1,ny;do i=1,nx
!
   !   uc(i,j,k) = 0.5d0*(u2(i,j,k)+u2(i-1,j,k))
   !   vc(i,j,k) = 0.5d0*(v2(i,j,k)+v2(i,j-1,k))
   !   wc(i,j,k) = 0.5d0*(w2(i,j,k)+w2(i,j,k-1))
!
   !end do;end do
   !!$OMP END PARALLEL DO
   !end do


end subroutine Updating_velocity