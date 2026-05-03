! 03 May 2026 by FDS

subroutine final_boundary_conditions()
use variables
implicit none


!---------------------------------------------------!
!         Boundary conditions calculation           !
!---------------------------------------------------!

!$acc data present(X,Y,Z,u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2), &
!$acc	u0(:,:,istart-2:iend+2),v0(:,:,istart-2:iend+2),w0(:,:,istart-2:iend+2))

!!!!!!!!Velocity boundary conditions!!!!!!!!

!West vertical wall i=0
    !$OMP PARALLEL
	!$OMP DO COLLAPSE(2)
	!$acc parallel
	!$acc loop independent collapse(2) gang vector
  	do k=istart-2,iend+2
	do j=-1,ny+2

  !-------------------------------------!
  	if(WestWall_u == 'NO_SLIP') then
    	u(0,j,k) = 0.d0
    	u(-1,j,k) = u(0,j,k)
  	else if(WestWall_u == 'NEUMANN') then
  	  	u(0,j,k) = u(1,j,k)
  	  	u(-1,j,k) = u(0,j,k)
  	else if(WestWall_u == 'DIRICHLET') then
  	  	u(0,j,k) = U_inf
  	  	u(-1,j,k) = u(0,j,k)
  	else if(WestWall_u == 'PERIODIC') then
  	  	u(0,j,k) = u(nx,j,k)
  	  	u(-1,j,k) = u(nx-1,j,k)
  	else if(WestWall_u == 'MOVING_REF') then
  	  	u(0,j,k) = -u_solid
  	  	u(-1,j,k) = u(0,j,k)
  	end if


  !-------------------------------------!
  
  !-------------------------------------!
  	if(WestWall_v == 'NO_SLIP') then
  	  	v(0,j,k) = 0.d0-v(1,j,k)
  	  	v(-1,j,k) = v(0,j,k)
  	else if(WestWall_v == 'NEUMANN') then
  	  	v(0,j,k) = v(1,j,k)
  	  	v(-1,j,k) = v(0,j,k)
  	else if(WestWall_v == 'DIRICHLET') then
  	  	v(0,j,k) = 2.d0*U_inf-v(1,j,k)
  	  	v(-1,j,k) = v(0,j,k)
  	else if(WestWall_v == 'PERIODIC') then
  	  	v(0,j,k) = v(nx,j,k)
  		v(-1,j,k) = v(nx-1,j,k)
  	end if


  	!-------------------------------------!

  	!-------------------------------------!
  	if(WestWall_w == 'NO_SLIP') then
  	  	w(0,j,k) = 0.d0-w(1,j,k)
  	  	w(-1,j,k) = w(0,j,k)
  	else if(WestWall_w == 'NEUMANN') then
  	  	w(0,j,k) = w(1,j,k)
  	  	w(-1,j,k) = w(0,j,k)
  	else if(WestWall_w == 'DIRICHLET') then
  	  	w(0,j,k) = 2.d0*U_inf-w(1,j,k)
  	  	w(-1,j,k) = w(0,j,k)
  	else if(WestWall_w == 'PERIODIC') then
  	  	w(0,j,k) = w(nx,j,k)
  	  	w(-1,j,k) = w(nx-1,j,k)
  	else if(WestWall_w == 'MOVING_REF') then
  	  	w(0,j,k) = -2.d0*w_solid-w(1,j,k)
  	  	w(-1,j,k) = w(0,j,k)
  	end if
  	!-------------------------------------!

	end do
  	end do
    !$OMP END DO

!East vertical wall
    !$OMP DO COLLAPSE(2)
	!$acc loop independent collapse(2) gang vector
  	do k=istart-2,iend+2
	do j=-1,ny+2

  !-------------------------------------!
  	if(EastWall_u == 'NO_SLIP') then
  	  	u(nx,j,k) = 0.d0
  	  	u(nx+1,j,k) = u(nx,j,k)
  	  	u(nx+2,j,k) = u(nx+1,j,k)
  	else if(EastWall_u == 'NEUMANN') then
  	  	u(nx,j,k) = u(nx-1,j,k)
  	  	u(nx+1,j,k) = u(nx,j,k)
  	  	u(nx+2,j,k) = u(nx+1,j,k)
  	else if(EastWall_u == 'DIRICHLET') then
  	  	u(nx,j,k) = U_inf
  	  	u(nx+1,j,k) = u(nx,j,k)
  	  	u(nx+2,j,k) = u(nx+1,j,k)
  	else if(EastWall_u == 'PERIODIC') then
  	  	!u(nx,j,k) = u()
  	  	u(nx+1,j,k) = u(1,j,k)
  	  	u(nx+2,j,k) = u(2,j,k)
  	else if(EastWall_u == 'MOVING_REF') then
  	  	u(nx,j,k) = -u_solid
  	  	u(nx+1,j,k) = u(nx,j,k)
  	  	u(nx+2,j,k) = u(nx+1,j,k)
  	end if


  !-------------------------------------!


  !-------------------------------------!
  	if(EastWall_v == 'NO_SLIP') then
   	v(nx+1,j,k) = 0.d0-v(nx,j,k)
   	v(nx+2,j,k) = v(nx+1,j,k)
  	else if(EastWall_v == 'NEUMANN') then
   	v(nx+1,j,k) = v(nx,j,k)
   	v(nx+2,j,k) = v(nx+1,j,k)
  	else if(EastWall_v == 'DIRICHLET') then
		v(nx+1,j,k) = 2.d0*U_inf-v(nx,j,k)
		v(nx+2,j,k) = v(nx+1,j,k)
  	else if(EastWall_v == 'PERIODIC') then
    	v(nx+1,j,k) = v(1,j,k)
    	v(nx+2,j,k) = v(2,j,k)
  	end if


  !-------------------------------------!

  !-------------------------------------!
  	if(EastWall_w == 'NO_SLIP') then
  	  	w(nx+1,j,k) = 0.d0-w(nx,j,k)
  	  	w(nx+2,j,k) = w(nx+1,j,k)
  	else if(EastWall_w == 'NEUMANN') then
  	  	w(nx+1,j,k) = w(nx,j,k)
  	  	w(nx+2,j,k) = w(nx+1,j,k)
  	else if(EastWall_w == 'DIRICHLET') then
  	  	w(nx+1,j,k) = 2.d0*U_inf-w(nx,j,k)
  	  	w(nx+2,j,k) = w(nx+1,j,k)
  	else if(EastWall_w == 'PERIODIC') then
  	  	w(nx+1,j,k) = w(1,j,k)
  	  	w(nx+2,j,k) = w(2,j,k)
  	else if(EastWall_w == 'MOVING_REF') then
  	  	w(nx+1,j,k) = -2.d0*w_solid-w(nx,j,k)
  	  	w(nx+2,j,k) = w(nx+1,j,k)
  	end if


  	!-------------------------------------!
  
	end do
	end do
	!$acc end parallel
    !$OMP END DO
	

!South horizontal wall
	!$OMP DO COLLAPSE(2)
	!$acc parallel
	!$acc loop independent collapse(2) gang vector
	do k=istart-2,iend+2
	do i=-1,nx+2

  	!-------------------------------------!
  	if(SouthWall_u == 'NO_SLIP') then
		u(i,0,k) = 0.d0-u(i,1,k)
		u(i,-1,k) = u(i,0,k)
  	else if(SouthWall_u == 'NEUMANN') then
		u(i,0,k) = u(i,1,k)
		u(i,-1,k) = u(i,0,k)
  	else if(SouthWall_u == 'DIRICHLET') then
		u(i,0,k) = 2.d0*U_inf-u(i,1,k)
		u(i,-1,k) = u(i,0,k)
  	else if(SouthWall_u == 'PERIODIC') then
  	  	u(i,0,k) = u(i,ny,k)
  	  	u(i,-1,k) = u(i,ny-1,k)	
  	else if(SouthWall_u == 'MOVING_REF') then
		u(i,0,k) = -2.d0*u_solid-u(i,1,k)
		u(i,-1,k) = u(i,0,k)
  	end if
  	!-------------------------------------!
  

  	!-------------------------------------!
  	if(SouthWall_v == 'NO_SLIP') then
		v(i,0,k) = 0.d0
	  	v(i,-1,k) = v(i,0,k)
  	else if(SouthWall_v == 'NEUMANN') then
  		v(i,0,k) = v(i,1,k)
	  	v(i,-1,k) = v(i,0,k)
  	else if(SouthWall_v == 'DIRICHLET') then
  		v(i,0,k) = U_inf
	  	v(i,-1,k) = v(i,0,k)
  	else if(SouthWall_v == 'PERIODIC') then
  	  	v(i,0,k) = v(i,ny,k)
  	  	v(i,-1,k) = v(i,ny-1,k)	
  	end if
  	!-------------------------------------!
  

  	!-------------------------------------!
  	if(SouthWall_w == 'NO_SLIP') then
  		w(i,0,k) = 0.d0-w(i,1,k)
	  	w(i,-1,k) = w(i,0,k)
  	else if(SouthWall_w == 'NEUMANN') then
  		w(i,0,k) = w(i,1,k)
	  	w(i,-1,k) = w(i,0,k)
  	else if(SouthWall_w == 'DIRICHLET') then
  		w(i,0,k) = 2.d0*U_inf-w(i,1,k)
	  	w(i,-1,k) = w(i,0,k)
  	else if(SouthWall_w == 'PERIODIC') then
  	  	w(i,0,k) = w(i,ny,k)
  	  	w(i,-1,k) = w(i,ny-1,k)	
  	else if(SouthWall_w == 'MOVING_REF') then
  		w(i,0,k) = -2.d0*w_solid-w(i,1,k)
	  	w(i,-1,k) = w(i,0,k)
  	end if
  	!-------------------------------------!

	end do
  	end do
	!$OMP END DO

!North horizontal wall j=ny
	!$OMP DO COLLAPSE(2)
	!$acc loop independent collapse(2) gang vector
  	do k=istart-2,iend+2
	do i=-1,nx+2

  	!-------------------------------------!
  	if(NorthWall_u == 'NO_SLIP') then
  	  	u(i,ny+1,k) = 0.d0-u(i,ny,k)
		u(i,ny+2,k) = u(i,ny+1,k)
  	else if(NorthWall_u == 'NEUMANN') then
  	  	u(i,ny+1,k) = u(i,ny,k)
		u(i,ny+2,k) = u(i,ny+1,k)
  	else if(NorthWall_u == 'DIRICHLET') then
  	  	u(i,ny+1,k) = 2.d0*U_inf-u(i,ny,k)
		u(i,ny+2,k) = u(i,ny+1,k)
  	else if(NorthWall_u == 'PERIODIC') then
  	  	u(i,ny+1,k) = u(i,1,k)
  	  	u(i,ny+2,k) = u(i,2,k)
  	else if(NorthWall_u == 'MOVING_REF') then
  	  	u(i,ny+1,k) = -2.d0*u_solid-u(i,ny,k)
		u(i,ny+2,k) = u(i,ny+1,k)
  	end if

  	!-------------------------------------!
  

  	!-------------------------------------!
  	if(NorthWall_v == 'NO_SLIP') then
  	  	v(i,ny,k) = 0.d0
	  	v(i,ny+1,k) = v(i,ny,k)
		v(i,ny+2,k) = v(i,ny+1,k)
  	else if(NorthWall_v == 'NEUMANN') then
  	  	v(i,ny,k) = v(i,ny-1,k)
	  	v(i,ny+1,k) = v(i,ny,k)
		v(i,ny+2,k) = v(i,ny+1,k)
  	else if(NorthWall_v == 'DIRICHLET') then
  	  	v(i,ny,k) = U_inf
	  	v(i,ny+1,k) = v(i,ny,k)
		v(i,ny+2,k) = v(i,ny+1,k)
  	else if(NorthWall_v == 'PERIODIC') then
  	  	v(i,ny,k) = v(i,0,k)
  	  	v(i,ny+1,k) = v(i,1,k)
  	  	v(i,ny+2,k) = v(i,2,k)
  	end if
	
  	!-------------------------------------!
  

  	!-------------------------------------!
  	if(NorthWall_w == 'NO_SLIP') then
  		w(i,ny+1,k) = 0.d0-w(i,ny,k)
	  	w(i,ny+2,k) = w(i,ny+1,k)
  	else if(NorthWall_w == 'NEUMANN') then
  		w(i,ny+1,k) = w(i,ny,k)
	  	w(i,ny+2,k) = w(i,ny+1,k)
  	else if(NorthWall_w == 'DIRICHLET') then
  		w(i,ny+1,k) = 2.d0*U_inf-w(i,ny,k)
	  	w(i,ny+2,k) = w(i,ny+1,k)
  	else if(NorthWall_w == 'PERIODIC') then
  	  	w(i,ny+1,k) = w(i,1,k)
  	  	w(i,ny+2,k) = w(i,2,k)
  	else if(NorthWall_w == 'MOVING_REF') then
  		w(i,ny+1,k) = -2.d0*w_solid-w(i,ny,k)
	  	w(i,ny+2,k) = w(i,ny+1,k)
  	end if

  	!-------------------------------------!
  
	end do
  	end do
    !$acc end parallel
	!$OMP END DO
    !$OMP END PARALLEL
	
!inlet
!Back horizontal wall k=0
    if(myid==0) then
		!$OMP PARALLEL DO COLLAPSE(2)
		!$acc parallel loop independent collapse(2) gang vector
		do j=-1,ny+2
		do i=-1,nx+2

		!-------------------------------------!
  	if(BackWall_u == 'NO_SLIP') then
  		u(i,j,0) =  0.d0-u(i,j,1)
		u(i,j,-1) = u(i,j,0)
  	else if(BackWall_u == 'NEUMANN') then
  		u(i,j,0) = u(i,j,1)
  	  	u(i,j,-1) = u(i,j,0)
	else if(BackWall_u == 'DIRICHLET') then
  		u(i,j,0) =  2.d0*U_inf-u(i,j,1)
  	  	u(i,j,-1) = u(i,j,0)
	else if(BackWall_u == 'PERIODIC') then
  	  	u(i,j,0) = u(i,j,nz)
  	  	u(i,j,-1) = u(i,j,nz-1)		
  	end if
		!-------------------------------------!
  

		!-------------------------------------!
  	if(BackWall_v == 'NO_SLIP') then
  		v(i,j,0) =  0.d0-v(i,j,1)
		v(i,j,-1) = v(i,j,0)
  	else if(BackWall_v == 'NEUMANN') then
  		v(i,j,0) = v(i,j,0)
  	  	v(i,j,-1) = v(i,j,0)
	else if(BackWall_v == 'DIRICHLET') then
  		v(i,j,0) =  2.d0*U_inf-v(i,j,1)
  	  	v(i,j,-1) = v(i,j,0)
	else if(BackWall_v == 'PERIODIC') then
  	  	v(i,j,0) = v(i,j,nz)
  		v(i,j,-1) = v(i,j,nz-1)	
	end if
		!-------------------------------------!
  

		!-------------------------------------!
  	if(BackWall_w == 'NO_SLIP') then
  		w(i,j,0) = 0.d0
		w(i,j,-1) = w(i,j,0)
  	else if(BackWall_w == 'NEUMANN') then
  		w(i,j,0) = w(i,j,1)
  	  	w(i,j,-1) = w(i,j,0)
	else if(BackWall_w == 'DIRICHLET') then
  		w(i,j,0) = U_inf
  	  	w(i,j,-1) = w(i,j,0)
	else if(BackWall_w == 'PERIODIC') then
  	  	w(i,j,0) = w(i,j,nz)
  	  	w(i,j,-1) = w(i,j,nz-1)	
	else if(BackWall_w == 'HALF_DIRICHLET') then
  		if (Y(j) .LT. y0) then
			w(i,j,0) = 0.d0
		else
			w(i,j,0) = U_inf
		end if
	  	w(i,j,-1) = w(i,j,0)
	end if
		!-------------------------------------!
  
		end do
		end do
		!$acc end parallel
		!$OMP END PARALLEL DO
	endif


!outlet
!Front horizontal wall
    if(myid==(nproc-1)) then
		!$OMP PARALLEL DO COLLAPSE(2)
		!$acc parallel loop independent collapse(2) gang vector
		do j=-1,ny+2
		do i=-1,nx+2

		!-------------------------------------!
  	if(FrontWall_u == 'NO_SLIP') then
  		u(i,j,nz+1) = 0.d0-u(i,j,nz)
		u(i,j,nz+2) = u(i,j,nz+1)
  	else if(FrontWall_u == 'NEUMANN') then
  		u(i,j,nz+1) = u(i,j,nz)
  	  	u(i,j,nz+2) = u(i,j,nz+1)
	else if(FrontWall_u == 'DIRICHLET') then
  		u(i,j,nz+1) = 2.d0*U_inf-u(i,j,nz)
  	  	u(i,j,nz+2) = u(i,j,nz+1)
	else if(FrontWall_u == 'PERIODIC') then
  	  	u(i,j,nz+1) = u(i,j,1)
  	  	u(i,j,nz+2) = u(i,j,2)	
	end if
		!-------------------------------------!
  

		!-------------------------------------!
  	if(FrontWall_v == 'NO_SLIP') then
  		v(i,j,nz+1) = 0.d0-v(i,j,nz)
		v(i,j,nz+2) = v(i,j,nz+1)
  	else if(FrontWall_v == 'NEUMANN') then
  		v(i,j,nz+1) = v(i,j,nz)
  	  	v(i,j,nz+2) = v(i,j,nz+1)
	else if(FrontWall_v == 'DIRICHLET') then
  		v(i,j,nz+1) = 2.d0*U_inf-v(i,j,nz)
  	  	v(i,j,nz+2) = v(i,j,nz+1)
	else if(FrontWall_v == 'PERIODIC') then
    	v(i,j,nz+1) = v(i,j,1)
    	v(i,j,nz+2) = v(i,j,2)	
	end if
		!-------------------------------------!


		!-------------------------------------!
  	if(FrontWall_w== 'NO_SLIP') then
  		w(i,j,nz) = 0.d0
	  	w(i,j,nz+1) = w(i,j,nz)
		w(i,j,nz+2) = w(i,j,nz+1)
  	else if(FrontWall_w == 'NEUMANN') then
  		w(i,j,nz) = w(i,j,nz-1)
  	  	w(i,j,nz+1) = w(i,j,nz)
		w(i,j,nz+2) = w(i,j,nz+1)
	else if(FrontWall_w == 'DIRICHLET') then
  		w(i,j,nz) = U_inf
  	  	w(i,j,nz+1) = w(i,j,nz)
		w(i,j,nz+2) = w(i,j,nz+1)
	else if(FrontWall_w == 'PERIODIC') then
  	  	w(i,j,nz+1) = w(i,j,1)
  	  	w(i,j,nz+2) = w(i,j,2)	
	end if
		!-------------------------------------!

		end do
		end do
		!$acc end parallel
		!$OMP END PARALLEL DO
	endif

    !$OMP PARALLEL DO COLLAPSE(3)
    !$acc parallel
	!$acc loop independent collapse(3) gang vector
   	do k=istart-2,iend+2 
  	do j=-1,ny+2
	do i=-1,nx+2

  	! u1(i,j,k)=u(i,j,k)
  	! v1(i,j,k)=v(i,j,k)
  	! w1(i,j,k)=w(i,j,k)

  	! u2(i,j,k)=u(i,j,k)
  	! v2(i,j,k)=v(i,j,k)
  	! w2(i,j,k)=w(i,j,k)

  	u0(i,j,k)=u(i,j,k)
  	v0(i,j,k)=v(i,j,k)
  	w0(i,j,k)=w(i,j,k)


	end do; end do; end do
    !$acc end parallel
    !$OMP END PARALLEL DO
!$acc end data

	call pressure_boundary_conditions()

end subroutine final_boundary_conditions



subroutine pressure_boundary_conditions()
use variables
implicit none

!!!!!!!Pressure boundary conditions!!!!!!!
!$acc data present(p(:,:,istart-2:iend+2))

!West vertical wall
   !$OMP PARALLEL
   !$OMP DO COLLAPSE(2)
   !$acc parallel
   !$acc loop independent collapse(2) gang vector   
  	do k=istart-2,iend+2 
	do j=-1,ny+2
	   if(WestWall_u == 'PERIODIC') then
		p(0,j,k)=p(nx,j,k)
		p(-1,j,k)=p(nx-1,j,k)
	   else
		p(0,j,k)=p(1,j,k)
		p(-1,j,k)=p(0,j,k)
	   endif
	end do
  	end do
	!$OMP END DO

!East vertical wall
   !$OMP DO COLLAPSE(2)
   !$acc loop independent collapse(2) gang vector 
  	do k=istart-2,iend+2
	do j=-1,ny+2
	   if(EastWall_u == 'PERIODIC') then
		p(nx+1,j,k)=p(1,j,k)
		p(nx+2,j,k)=p(2,j,k)
	   else
		p(nx+1,j,k)=p(nx,j,k)
		p(nx+2,j,k)=p(nx+1,j,k)
	   endif
	end do
  	end do
   !$acc end parallel
   !$OMP END DO

!South horizontal wall
   !$OMP DO COLLAPSE(2)
   !$acc parallel
   !$acc loop independent collapse(2) gang vector 
  	do k=istart-2,iend+2
	do i=-1,nx+2
	   if(SouthWall_v == 'PERIODIC') then
		p(i,0,k)=p(i,ny,k)
		p(i,-1,k)=p(i,ny-1,k)
	   else
		p(i,0,k)=p(i,1,k)
		p(i,-1,k)=p(i,0,k)
	   endif
	end do
  	end do
   !$OMP END DO
  
!North horizontal wall
   !$OMP DO COLLAPSE(2)
   !$acc loop independent collapse(2) gang vector 
  	do k=istart-2,iend+2
	do i=-1,nx+2
	   if(NorthWall_v == 'PERIODIC') then
		p(i,ny+1,k)=p(i,1,k)
		p(i,ny+2,k)=p(i,2,k)
	   else
		p(i,ny+1,k)=p(i,ny,k)
		p(i,ny+2,k)=p(i,ny+1,k)
	   endif
	end do
  	end do
  !$acc end parallel
  !$OMP END DO
  !$OMP END PARALLEL

!Back horizontal wall
	if(myid==0) then
		!$OMP PARALLEL DO COLLAPSE(2)
		!$acc parallel loop independent collapse(2) gang vector
			do j=-1,ny+2
			do i=-1,nx+2
			   if(BackWall_w == 'PERIODIC') then
				p(i,j,0)=p(i,j,nz)
				p(i,j,-1)=p(i,j,nz-1)
			   else
				p(i,j,0)=p(i,j,1)
				p(i,j,-1)=p(i,j,0)
			   endif				
			end do
			end do
		!$acc end parallel   
		!$OMP END PARALLEL DO
	endif

!Front horizontal wall
	if(myid==nproc-1) then
		!$OMP PARALLEL DO COLLAPSE(2)
		!$acc parallel loop independent collapse(2) gang vector
			do j=-1,ny+2
			do i=-1,nx+2
			   if(FrontWall_w == 'PERIODIC') then
				p(i,j,nz+1)=p(i,j,1)
				p(i,j,nz+2)=p(i,j,2)
			   else
			    p(i,j,nz+1)=p(i,j,nz)
			    p(i,j,nz+2)=p(i,j,nz+1)
			   endif
			end do
			end do
		!$acc end parallel
		!$OMP END PARALLEL DO
	endif
!$acc end data	

end subroutine pressure_boundary_conditions



subroutine u0_boundary_conditions()
use variables
implicit none

! ! Sending data to accelerator
!$acc data present(u0(:,:,istart-2:iend+2),v0(:,:,istart-2:iend+2),w0(:,:,istart-2:iend+2))

!West vertical wall i=0
  !$OMP PARALLEL
	!$OMP DO COLLAPSE(2)
	!$acc parallel
	!$acc loop independent collapse(2) gang vector
  	do k=istart-2,iend+2
	do j=-1,ny+2

  	!-------------------------------------!
  	if(WestWall_u == 'PERIODIC') then
  	  	u0(0,j,k) = u0(nx,j,k)
  	  	u0(-1,j,k) = u0(nx-1,j,k)
  	end if

  	!-------------------------------------!
  	if(WestWall_v == 'PERIODIC') then
  	  	v0(0,j,k) = v0(nx,j,k)
  		v0(-1,j,k) = v0(nx-1,j,k)
  	end if

  	!-------------------------------------!
  	if(WestWall_w == 'PERIODIC') then
  	  	w0(0,j,k) = w0(nx,j,k)
  	  	w0(-1,j,k) = w0(nx-1,j,k)
  	end if


	end do
  	end do
	!$OMP END DO


!East vertical wall
	!$OMP DO COLLAPSE(2)
	!$acc loop independent collapse(2) gang vector
  	do k=istart-2,iend+2
	do j=-1,ny+2

  	!-------------------------------------!
  	if(EastWall_u == 'PERIODIC') then
  	  	u0(nx+1,j,k) = u0(1,j,k)
  	  	u0(nx+2,j,k) = u0(2,j,k)
  	end if

  	!-------------------------------------!
  	if(EastWall_v == 'PERIODIC') then
    	v0(nx+1,j,k) = v0(1,j,k)
    	v0(nx+2,j,k) = v0(2,j,k)
  	end if

  	!-------------------------------------!
  	if(EastWall_w == 'PERIODIC') then
  	  	w0(nx+1,j,k) = w0(1,j,k)
  	  	w0(nx+2,j,k) = w0(2,j,k)
  	end if

	end do
  	end do
	!$acc end parallel
	!$OMP END DO
	
	
!South horizontal wall j=0
	!$OMP DO COLLAPSE(2)
	!$acc parallel
	!$acc loop independent collapse(2) gang vector
	do k=istart-2,iend+2
	do i=-1,nx+2

  	!-------------------------------------!
  	if(SouthWall_u == 'PERIODIC') then
  	  	u0(i,0,k) = u0(i,ny,k)
  	  	u0(i,-1,k) = u0(i,ny-1,k)
  	end if

  	!-------------------------------------!
  	if(SouthWall_v == 'PERIODIC') then
  	  	v0(i,0,k) = v0(i,ny,k)
  	  	v0(i,-1,k) = v0(i,ny-1,k)
  	end if

  	!-------------------------------------!
  	if(SouthWall_w == 'PERIODIC') then
  	  	w0(i,0,k) = w0(i,ny,k)
  	  	w0(i,-1,k) = w0(i,ny-1,k)
  	end if


	end do
  	end do
	!$OMP END DO

!North horizontal wall
	!$OMP DO COLLAPSE(2)
	!$acc loop independent collapse(2) gang vector
	do k=istart-2,iend+2
	do i=-1,nx+2

  	!-------------------------------------!
  	if(NorthWall_u == 'PERIODIC') then
  	  	u0(i,ny+1,k) = u0(i,1,k)
  	  	u0(i,ny+2,k) = u0(i,2,k)
  	end if

  	!-------------------------------------!
  	if(NorthWall_v == 'PERIODIC') then
  	  	v0(i,ny+1,k) = v0(i,1,k)
  	  	v0(i,ny+2,k) = v0(i,2,k)
  	end if

  	!-------------------------------------!
  	if(NorthWall_w == 'PERIODIC') then
  	  	w0(i,ny+1,k) = w0(i,1,k)
  	  	w0(i,ny+2,k) = w0(i,2,k)
  	end if

	end do
  	end do
	!$acc end parallel
	!$OMP END DO
  !$OMP END PARALLEL
  
!inlet
!Back horizontal wall k=0
    if(myid==0) then
	!$OMP PARALLEL DO COLLAPSE(2)
	!$acc parallel loop independent collapse(2) gang vector
  	do j=-1,ny+2
	do i=-1,nx+2

  	!-------------------------------------!
	if(BackWall_u == 'PERIODIC') then
  	  	u0(i,j,0) = u0(i,j,nz)
  	  	u0(i,j,-1) = u0(i,j,nz-1)		
  	end if
  	!-------------------------------------!
  

  	!-------------------------------------!
	if(BackWall_v == 'PERIODIC') then
  	  	v0(i,j,0) = v0(i,j,nz)
  		v0(i,j,-1) = v0(i,j,nz-1)	
	end if
  	!-------------------------------------!
  

  	!-------------------------------------!
	if(BackWall_w == 'PERIODIC') then
  	  	w0(i,j,0) = w0(i,j,nz)
  	  	w0(i,j,-1) = w0(i,j,nz-1)	
	end if
  	!-------------------------------------!
  
	end do
  	end do
	!$acc end parallel
	!$OMP END PARALLEL DO
	endif


!outlet
!Front horizontal wall
    if(myid==(nproc-1)) then
	!$OMP PARALLEL DO COLLAPSE(2)
	!$acc parallel loop independent collapse(2) gang vector
  	do j=-1,ny+2
	do i=-1,nx+2

  	!-------------------------------------!
	if(FrontWall_u == 'PERIODIC') then
  	  	u0(i,j,nz+1) = u0(i,j,1)
  	  	u0(i,j,nz+2) = u0(i,j,2)	
	end if
  	!-------------------------------------!
  

  	!-------------------------------------!
	if(FrontWall_v == 'PERIODIC') then
    	v0(i,j,nz+1) = v0(i,j,1)
    	v0(i,j,nz+2) = v0(i,j,2)	
	end if
  	!-------------------------------------!


  	!-------------------------------------!
	if(FrontWall_w == 'PERIODIC') then
  	  	w0(i,j,nz+1) = w0(i,j,1)
  	  	w0(i,j,nz+2) = w0(i,j,2)	
	end if
  	!-------------------------------------!

	end do
  	end do
	!$acc end parallel
	!$OMP END PARALLEL DO
	endif

!$acc end data	

end subroutine u0_boundary_conditions



subroutine theta_boundary_conditions()
use variables
implicit none

!!!!!!!Theta boundary conditions!!!!!!!
!$acc data present(theta(:,:,istart-2:iend+2))

!West vertical wall
   !$OMP PARALLEL
   !$OMP DO COLLAPSE(2)
   !$acc parallel
   !$acc loop independent collapse(2) gang vector   
  	do k=istart-2,iend+2 
	do j=-1,ny+2
	   if(WestWall == 'ADIABATIC') then
		theta(0,j,k)=theta(1,j,k)
		theta(-1,j,k)=theta(0,j,k)
	   else if(WestWall == 'ISOTHERMAL') then
		theta(0,j,k)=2.d0*Tiso_W-theta(1,j,k)
		theta(-1,j,k)=theta(0,j,k)
	   else if(WestWall == 'PERIODIC') then
		theta(0,j,k)=theta(nx,j,k)
		theta(-1,j,k)=theta(nx-1,j,k)
	   endif
	end do
  	end do
	!$OMP END DO

!East vertical wall
   !$OMP DO COLLAPSE(2)
   !$acc loop independent collapse(2) gang vector 
  	do k=istart-2,iend+2
	do j=-1,ny+2
	   if(EastWall == 'ADIABATIC') then
		theta(nx+1,j,k)=theta(nx,j,k)
		theta(nx+2,j,k)=theta(nx+1,j,k)
	   else if(EastWall == 'ISOTHERMAL') then
		theta(nx+1,j,k)=2.d0*Tiso_E-theta(nx,j,k)
		theta(nx+2,j,k)=theta(nx+1,j,k)
	   elseif(EastWall == 'PERIODIC') then
		theta(nx+1,j,k)=theta(1,j,k)
		theta(nx+2,j,k)=theta(2,j,k)
	   endif
	end do
  	end do
   !$acc end parallel
   !$OMP END DO

!South horizontal wall
   !$OMP DO COLLAPSE(2)
   !$acc parallel
   !$acc loop independent collapse(2) gang vector 
  	do k=istart-2,iend+2
	do i=-1,nx+2
	   if(SouthWall == 'ADIABATIC') then
		theta(i,0,k)=theta(i,1,k)
		theta(i,-1,k)=theta(i,0,k)
	   else if(SouthWall == 'ISOTHERMAL') then
		theta(i,0,k)=2.d0*Tiso_S-theta(i,1,k)
		theta(i,-1,k)=theta(i,0,k)
	   else if(SouthWall == 'PERIODIC') then
		theta(i,0,k)=theta(i,ny,k)
		theta(i,-1,k)=theta(i,ny-1,k)
	   endif
	end do
  	end do
   !$OMP END DO
  
!North horizontal wall
   !$OMP DO COLLAPSE(2)
   !$acc loop independent collapse(2) gang vector 
  	do k=istart-2,iend+2
	do i=-1,nx+2
	   if(NorthWall == 'ADIABATIC') then
		theta(i,ny+1,k)=theta(i,ny,k)
		theta(i,ny+2,k)=theta(i,ny+1,k)
	   else if(NorthWall == 'ISOTHERMAL') then
		theta(i,ny+1,k)=2.d0*Tiso_N-theta(i,ny,k)
		theta(i,ny+2,k)=theta(i,ny+1,k)
	   else if(NorthWall == 'PERIODIC') then
		theta(i,ny+1,k)=theta(i,1,k)
		theta(i,ny+2,k)=theta(i,2,k)
	   endif
	end do
  	end do
  !$acc end parallel
  !$OMP END DO
  !$OMP END PARALLEL

!Back horizontal wall
	if(myid==0) then
		!$OMP PARALLEL DO COLLAPSE(2)
		!$acc parallel loop independent collapse(2) gang vector
			do j=-1,ny+2
			do i=-1,nx+2
			   if(BackWall == 'ADIABATIC') then
				theta(i,j,0)=theta(i,j,1)
				theta(i,j,-1)=theta(i,j,0)
			   else if(BackWall == 'ISOTHERMAL') then
				theta(i,j,0)=2.d0*Tiso_B-theta(i,j,1)
				theta(i,j,-1)=theta(i,j,0)
			   else if(BackWall == 'PERIODIC') then
				theta(i,j,0)=theta(i,j,nz)
				theta(i,j,-1)=theta(i,j,nz-1)
			   endif				
			end do
			end do
		!$acc end parallel   
		!$OMP END PARALLEL DO
	endif

!Front horizontal wall
	if(myid==nproc-1) then
		!$OMP PARALLEL DO COLLAPSE(2)
		!$acc parallel loop independent collapse(2) gang vector
			do j=-1,ny+2
			do i=-1,nx+2
			   if(FrontWall == 'ADIABATIC') then
			    theta(i,j,nz+1)=theta(i,j,nz)
			    theta(i,j,nz+2)=theta(i,j,nz+1)
			   else if(FrontWall == 'ISOTHERMAL') then
			    theta(i,j,nz+1)=2.d0*Tiso_F-theta(i,j,nz)
			    theta(i,j,nz+2)=theta(i,j,nz+1)
			   else if(FrontWall == 'PERIODIC') then
				theta(i,j,nz+1)=theta(i,j,1)
				theta(i,j,nz+2)=theta(i,j,2)
			   endif
			end do
			end do
		!$acc end parallel
		!$OMP END PARALLEL DO
	endif
!$acc end data	

end subroutine theta_boundary_conditions



! subroutine theta0_boundary_conditions()
! use variables
! implicit none

! !!!!!!!Theta0 boundary conditions!!!!!!!
! !$acc data present(theta0(:,:,istart-2:iend+2))

! !West vertical wall
   ! !$OMP PARALLEL
   ! !$OMP DO COLLAPSE(2)
   ! !$acc parallel
   ! !$acc loop independent collapse(2) gang vector   
  	! do k=istart-2,iend+2 
	! do j=-1,ny+2
	   ! if(WestWall == 'PERIODIC') then
		! theta0(0,j,k)=theta0(nx,j,k)
		! theta0(-1,j,k)=theta0(nx-1,j,k)
	! end do
  	! end do
	! !$OMP END DO

! !East vertical wall
   ! !$OMP DO COLLAPSE(2)
   ! !$acc loop independent collapse(2) gang vector 
  	! do k=istart-2,iend+2
	! do j=-1,ny+2
	   ! if(EastWall == 'PERIODIC') then
		! theta0(nx+1,j,k)=theta0(1,j,k)
		! theta0(nx+2,j,k)=theta0(2,j,k)
	! end do
  	! end do
   ! !$acc end parallel
   ! !$OMP END DO

! !South horizontal wall
   ! !$OMP DO COLLAPSE(2)
   ! !$acc parallel
   ! !$acc loop independent collapse(2) gang vector 
  	! do k=istart-2,iend+2
	! do i=-1,nx+2
	   ! if(SouthWall == 'PERIODIC') then
		! theta0(i,0,k)=theta0(i,ny,k)
		! theta0(i,-1,k)=theta0(i,ny-1,k)
	! end do
  	! end do
   ! !$OMP END DO
  
! !North horizontal wall
   ! !$OMP DO COLLAPSE(2)
   ! !$acc loop independent collapse(2) gang vector 
  	! do k=istart-2,iend+2
	! do i=-1,nx+2
	   ! if(NorthWall == 'PERIODIC') then
		! theta0(i,ny+1,k)=theta0(i,1,k)
		! theta0(i,ny+2,k)=theta0(i,2,k)
	! end do
  	! end do
  ! !$acc end parallel
  ! !$OMP END DO
  ! !$OMP END PARALLEL

! !Back horizontal wall
	! if(myid==0) then
		! !$OMP PARALLEL DO COLLAPSE(2)
		! !$acc parallel loop independent collapse(2) gang vector
			! do j=-1,ny+2
			! do i=-1,nx+2
			   ! if(BackWall == 'PERIODIC') then
				! theta0(i,j,0)=theta0(i,j,nz)
				! theta0(i,j,-1)=theta0(i,j,nz-1)				
			! end do
			! end do
		! !$acc end parallel   
		! !$OMP END PARALLEL DO
	! endif

! !Front horizontal wall
	! if(myid==nproc-1) then
		! !$OMP PARALLEL DO COLLAPSE(2)
		! !$acc parallel loop independent collapse(2) gang vector
			! do j=-1,ny+2
			! do i=-1,nx+2
			   ! if(FrontWall == 'PERIODIC') then
				! theta0(i,j,nz+1)=theta0(i,j,1)
				! theta0(i,j,nz+2)=theta0(i,j,2)
			! end do
			! end do
		! !$acc end parallel
		! !$OMP END PARALLEL DO
	! endif
! !$acc end data	

! end subroutine theta0_boundary_conditions