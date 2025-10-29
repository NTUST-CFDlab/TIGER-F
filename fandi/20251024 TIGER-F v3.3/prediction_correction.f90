! 12 Jan 2025 - FDS

subroutine prediction_correction()
use variables
implicit none

!$acc data present(FX(:,:,istart-2:iend+2), FY(:,:,istart-2:iend+2), FZ(:,:,istart-2:iend+2), &
!$acc 	   		  FX1(:,:,istart-2:iend+2),FY1(:,:,istart-2:iend+2),FZ1(:,:,istart-2:iend+2))

    !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
    !$acc parallel loop independent collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx
        FX1(i,j,k) = FX(i,j,k)
        FY1(i,j,k) = FY(i,j,k)
        FZ1(i,j,k) = FZ(i,j,k)
    enddo; enddo; enddo
    !$OMP END PARALLEL DO
	!$acc end parallel


    do ccc=1,correction_stage
        call correction()
    enddo
    ccc=1


    !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
    !$acc parallel loop independent collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx
        FX(i,j,k) = FX1(i,j,k)
        FY(i,j,k) = FY1(i,j,k)
        FZ(i,j,k) = FZ1(i,j,k)
    enddo; enddo; enddo
    !$OMP END PARALLEL DO
	!$acc end parallel

!$acc end data

end subroutine prediction_correction




subroutine correction()
use variables
implicit none

!$acc data present(p(:,:,istart-2:iend+2),p_pre(:,:,istart-2:iend+2,:), &
!$acc	FX1(:,:,istart-2:iend+2),FY1(:,:,istart-2:iend+2),FZ1(:,:,istart-2:iend+2), &
!$acc    FX(:,:,istart-2:iend+2), FY(:,:,istart-2:iend+2), FZ(:,:,istart-2:iend+2), &
!$acc	 u0(:,:,istart-2:iend+2), v0(:,:,istart-2:iend+2), w0(:,:,istart-2:iend+2))


    !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
    !$acc parallel loop independent collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx

        u0(i,j,k) = u0(i,j,k) + FX(i,j,k)*dt
        v0(i,j,k) = v0(i,j,k) + FY(i,j,k)*dt
        w0(i,j,k) = w0(i,j,k) + FZ(i,j,k)*dt

    enddo; enddo; enddo
	!$acc end parallel
    !$OMP END PARALLEL DO

    call u0_boundary_conditions()


      !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
	  !$acc update self(w0(:,:,iend)) if(nproc>1)
	  !$acc update self(w0(:,:,istart-1)) if(nproc>1)
      icount = (nx+4)*(ny+4)
      itag = 301
      call MPI_SENDRECV( w0(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
                         w0(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	  !$acc update device(w0(:,:,istart-1)) if(nproc>1)     
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


	!$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
	!$acc parallel loop independent collapse(3) gang vector 
	do k=istart-2,iend+2; do j=-1,ny+2; do i=-1,nx+2
        p(i,j,k)=p_pre(i,j,k,ccc)
	end do; end do; enddo
	!$acc end parallel
	!$OMP END PARALLEL DO

      !------------------------------------------------------------------------------------------------------!
      if (pressure_solver == 1) then
       call RB_SOR()			! calculate pressure field using red-black SOR --> nx MUST be an even number
      elseif (pressure_solver == 2) then
	   call SOR()			 	! calculate pressure field using traditional SOR	  
      endif
!      call gauss_seidel()        ! calculate pressure field using traditional SOR (Old)
!      call BICG_stab()           ! calculate pressure field
      !------------------------------------------------------------------------------------------------------!

	!$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
	!$acc parallel loop independent collapse(3) gang vector 
	do k=istart-2,iend+2; do j=-1,ny+2; do i=-1,nx+2
        p_pre(i,j,k,ccc)=p(i,j,k)
	end do; end do; enddo
	!$acc end parallel
	!$OMP END PARALLEL DO


      ! if( time .GE. StartDynamic_time )then  ! Use in two-way simulation

		! !------------------------------------------------------------------------------------------------------!
		! ! Virtual force integrator
		! !		call virtualForceIntegrator_nima() 	! CD and CL only
				! call virtualForceTorqueIntegrator() ! CD, CL, CT, and CP
		! !		call virtualForceTorque_frozen()	! CD, CL, CT, and CP for frozen rotor simulation
		! !------------------------------------------------------------------------------------------------------!
		 
        ! !--------------------------------------------------------------------!
        ! call twoway_dynamic() 
        ! !--------------------------------------------------------------------!

      ! end if

      !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
	  !$acc update self(p(:,:,istart)) if(nproc>1)
	  !$acc update self(p(:,:,iend+1)) if(nproc>1)
      icount = (nx+4)*(ny+4)
      itag = 302
      call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	  !$acc update device(p(:,:,iend+1)) if(nproc>1)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !------------------------------------------------------------------------------------------------------!
    call calcul_new_velocity() ! update velocity field 
    !------------------------------------------------------------------------------------------------------!


    !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
    !$acc parallel loop independent collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx
        FX1(i,j,k) = FX(i,j,k) + FX1(i,j,k)
        FY1(i,j,k) = FY(i,j,k) + FY1(i,j,k)
        FZ1(i,j,k) = FZ(i,j,k) + FZ1(i,j,k)
    enddo; enddo; enddo
	!$acc end parallel
    !$OMP END PARALLEL DO

!$acc end data

end subroutine correction
