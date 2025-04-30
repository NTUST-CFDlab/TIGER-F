subroutine prediction_correction()
use variables
implicit none

    !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx
        FX1(i,j,k) = FX(i,j,k)
        FY1(i,j,k) = FY(i,j,k)
        FZ1(i,j,k) = FZ(i,j,k)
    enddo; enddo; enddo
    !$OMP END PARALLEL DO


    do ccc=1,correction_stage
        call correction()
    enddo
    ccc=1


    !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx
        FX(i,j,k) = FX1(i,j,k)
        FY(i,j,k) = FY1(i,j,k)
        FZ(i,j,k) = FZ1(i,j,k)
    enddo; enddo; enddo
    !$OMP END PARALLEL DO


end subroutine prediction_correction




subroutine correction()
use variables
implicit none

    !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx

        u0(i,j,k) = u0(i,j,k) + FX(i,j,k)*dt
        v0(i,j,k) = v0(i,j,k) + FY(i,j,k)*dt
        w0(i,j,k) = w0(i,j,k) + FZ(i,j,k)*dt

    enddo; enddo; enddo
    !$OMP END PARALLEL DO

    call u0_boundary_conditions()


    !----------data transformation among nodes----------!
    icount = (nx+4)*(ny+4)
        itag = 220
        call MPI_SENDRECV( w0(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
                           w0(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !----------data transformation among nodes----------!



    do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2
        p(i,j,k)=p_pre(i,j,k,ccc)
    end do; end do; enddo

      !------------------------------------------------------------------------------------------------------!
      if (pressure_solver == 1) then
       call SOR()				 ! calculate pressure field using traditional SOR (faster)
      elseif (pressure_solver == 2) then
	   call RB_SOR()			 ! calculate pressure field using red-black SOR --> nx MUST be an even number	  
      endif
!      call gauss_seidel()        ! calculate pressure field using traditional SOR (Old)
!      call BICG_stab()           ! calculate pressure field
      !------------------------------------------------------------------------------------------------------!

    do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2
        p_pre(i,j,k,ccc)=p(i,j,k)
    end do; end do; enddo


    ! !----------data transformation among nodes----------!  --> Moved to GS.f90
    ! icount = (nx+4)*(ny+4)

      ! itag = 230
      ! call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         ! p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! !----------data transformation among nodes----------!

    !------------------------------------------------------------------------------------------------------!
    call calcul_new_velocity() ! update velocity field 
    !------------------------------------------------------------------------------------------------------!

    !------------------------------------------------------------------------------------------------------!
    !call virtualForceIntegrator() ! COEFFICIENTS: CD and Center_CL,Distance,Velocity for solid  
    !call virtualForceIntegrator_nima()
    !------------------------------------------------------------------------------------------------------!


    !$OMP PARALLEL DO PRIVATE(j,i) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx
        FX1(i,j,k) = FX(i,j,k) + FX1(i,j,k)
        FY1(i,j,k) = FY(i,j,k) + FY1(i,j,k)
        FZ1(i,j,k) = FZ(i,j,k) + FZ1(i,j,k)
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

end subroutine correction
