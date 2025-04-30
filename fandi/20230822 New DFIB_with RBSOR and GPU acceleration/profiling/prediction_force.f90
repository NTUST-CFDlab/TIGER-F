subroutine Prediction_Force()
use variables
implicit none

    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(i)
    do j=1,ny; do i=1,nx
        FX1(i,j,k) = FX(i,j,k)
        FY1(i,j,k) = FY(i,j,k)
        FZ1(i,j,k) = FZ(i,j,k)
    enddo; enddo
    !$OMP END PARALLEL DO
    enddo
!
!
!
    do ccc=1,1
        call correction()
    enddo
    ccc=1


    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(i)
    do j=1,ny; do i=1,nx
        FX(i,j,k) = FX1(i,j,k)
        FY(i,j,k) = FY1(i,j,k)
        FZ(i,j,k) = FZ1(i,j,k)
    enddo; enddo
    !$OMP END PARALLEL DO
    enddo

   !----------data collect among nodes----------!
   icount = igcount*(nx+4)*(ny+4)

   !Send my results back to the master
   if(myid>master)then

      itag = 380
      call MPI_SEND( FX(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      itag = 390
      call MPI_SEND( FY(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      itag = 400
      call MPI_SEND( FZ(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )


   end if


   !Wait to receive results from each task
   if(myid==master)then
      do i = 1, (nproc-1)

            itag = 380
            call MPI_RECV( FX(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 390
            call MPI_RECV( FY(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 400
            call MPI_RECV( FZ(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
               
      end do
   end if

   !----------data collect among nodes----------!

    !call virtualForceIntegrator() ! COEFFICIENTS: CD and Center_CL,Distance,Velocity for solid  
    call virtualForceIntegrator_nima()


end subroutine Prediction_Force




subroutine correction()
use variables
implicit none

    do k=istart,iend
    !$OMP PARALLEL DO
    do j=1,ny; do i=1,nx

        u0(i,j,k) = u0(i,j,k) + FX(i,j,k)*dt
        v0(i,j,k) = v0(i,j,k) + FY(i,j,k)*dt
        w0(i,j,k) = w0(i,j,k) + FZ(i,j,k)*dt

    enddo; enddo
    !$OMP END PARALLEL DO
    enddo

    call u0_boundary_conditions()


    !----------data transformation among nodes----------!
    icount = (nx+4)*(ny+4)
        itag = 220
        call MPI_SENDRECV( w0(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
                           w0(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !----------data transformation among nodes----------!



    do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2

        p(i,j,k)=p_pre(i,j,k,ccc)

    end do; end do; enddo

    !------------------------------------------------------------------------------------------------------!
    call gauss_seidel()       ! calculate pressure field
    !------------------------------------------------------------------------------------------------------!

    do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2

        p_pre(i,j,k,ccc)=p(i,j,k)

    end do; end do; enddo


    !----------data transformation among nodes----------!
    icount = (nx+4)*(ny+4)

      itag = 230
      call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !----------data transformation among nodes----------!

    !------------------------------------------------------------------------------------------------------!
    call calcul_new_velocity() ! update velocity field 
    !------------------------------------------------------------------------------------------------------!

    !------------------------------------------------------------------------------------------------------!
    !call virtualForceIntegrator() ! COEFFICIENTS: CD and Center_CL,Distance,Velocity for solid  
    !call virtualForceIntegrator_nima()
    !------------------------------------------------------------------------------------------------------!


    do k=istart,iend
    !$OMP PARALLEL DO
    do j=1,ny; do i=1,nx
        FX1(i,j,k) = FX(i,j,k) + FX1(i,j,k)
        FY1(i,j,k) = FY(i,j,k) + FY1(i,j,k)
        FZ1(i,j,k) = FZ(i,j,k) + FZ1(i,j,k)
    enddo; enddo
    !$OMP END PARALLEL DO
    enddo

end subroutine correction



subroutine initial_Prediction()
use variables
implicit none

    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(i)
    do j=1,ny; do i=1,nx
        FX1(i,j,k) = FX(i,j,k)
        FY1(i,j,k) = FY(i,j,k)
        FZ1(i,j,k) = FZ(i,j,k)
    enddo; enddo
    !$OMP END PARALLEL DO
    enddo


end subroutine initial_Prediction
