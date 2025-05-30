! 12 Jan 2025 - FDS

subroutine RB_SOR() 
    use variables
    implicit none

! automatic omega upper & lower bounds, set according to ik from previous timestep
! the purpose is to avoid the use of high omega max once a stable # of P iteration is achieved
	if (ik > 500 .OR. istep <= 2) then	
		omega_max = 1.9d0
	elseif (ik > 50) then
		omega_max = 1.8d0
	else
		omega_max = 1.7d0		
	endif

	omega_min = 1.1d0
	d_omega = omega_max-omega_min


    ik=0
    pChangeMax = 1.d0
    pChangeMax_= 1.d0	

! Sending data to accelerator
!$acc data present(iDx,iDy,iDz,p(:,:,istart-2:iend+2), &
!$acc	u0(:,:,istart-2:iend+2),v0(:,:,istart-2:iend+2),w0(:,:,istart-2:iend+2), &
!$acc	P_Den(:,:,istart:iend,:),Den_inv(:,:,istart:iend)) &
!$acc create(mChange(:,:,istart:iend))

    !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
    !$acc parallel loop independent collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx

        mChange(i,j,k) =( ( u0(i,j,k) - u0(i-1,j,k) ) * iDy(j) * iDz(k) &
                        + ( v0(i,j,k) - v0(i,j-1,k) ) * iDx(i) * iDz(k) &
                        + ( w0(i,j,k) - w0(i,j,k-1) ) * iDx(i) * iDy(j) ) * den_flu * inv_dt

    enddo; enddo
	enddo
	!$acc end parallel
    !$OMP END PARALLEL DO


do while (pChangeMax_>zeta .AND. ik < itmax)

	! Setting dynamic omega response with respect to pChangeMax_
	if (pChangeMax_ < k_zeta_zeta) then	! omega decreases gradually from (k x zeta) to (zeta)
		omega = omega_min + d_omega * (LOG(pChangeMax_) - log_zeta)/log_k_zeta							! loglinear function (set parameters in variables module)
		! omega = omega_min + d_omega / (1.0d0 + EXP(-k_sigmoid * (LOG10(pChangeMax_) - x0_sigmoid)))	! Sigmoid function (set parameters in variables module)
	else
		omega = omega_max	! A constant of omega_max is used when pChangeMax_ > k_zeta*zeta
	endif
	
	
    ik=ik+1
    pChangeMax = 0.d0
    pChangeMax_= 0.d0

    !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
	!$acc update self(p(:,:,istart-1:istart)) if(nproc>1)
	!$acc update self(p(:,:,iend:iend+1)) if(nproc>1)
    icount = (nx+4)*(ny+4)

    itag = 201
    call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                       p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )

    itag = 202
    call MPI_SENDRECV( p(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
                       p(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	!$acc update device(p(:,:,iend+1)) if(nproc>1)
	!$acc update device(p(:,:,istart-1)) if(nproc>1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


! Red-Black ordering (RED)

    !$OMP PARALLEL
	!$OMP DO PRIVATE(i,j,pNEW) REDUCTION(max : pChangeMax) collapse(nclps)
	!$acc parallel vector_length(512)
	!$acc loop independent private(pNEW) reduction(max:pChangeMax) collapse(3) gang vector
    do k=istart,iend  
		do j=1,ny
            do i=1,nx
				!!$acc cache (p(i,j-8:j+8,k-8:k+8))			
				if (mod(i+j+k,2) == 1) then
                           ! pNEW = ( p(i+1,j,k) * P_Den(i,j,k,1) &
                                  ! + p(i-1,j,k) * P_Den(i,j,k,2) &
                           pNEW = ( p(MOD(i,nx)+1,j,k) * P_Den(i,j,k,1) & ! Periodic BC
                                  + p(MOD(i-2+nx,nx)+1,j,k) * P_Den(i,j,k,2) & ! Periodic BC
                                  + p(i,j+1,k) * P_Den(i,j,k,3) &
                                  + p(i,j-1,k) * P_Den(i,j,k,4) &
                                  + p(i,j,k+1) * P_Den(i,j,k,5) &
                                  + p(i,j,k-1) * P_Den(i,j,k,6) &
                                  - mChange(i,j,k) ) * Den_inv(i,j,k)	

				pChangeMax=MAX(pChangeMax,ABS( omega*(pNEW - p(i,j,k))))
				p(i,j,k) = (1.d0 - omega)*p(i,j,k) + omega*pNEW			
				endif
            enddo
        enddo
	enddo
	!$acc end parallel
    !$OMP END DO

! Red-Black ordering (BLACK)


	!$OMP DO PRIVATE(i,j,pNEW) REDUCTION(max : pChangeMax) collapse(nclps)
	!$acc parallel vector_length(512)
	!$acc loop independent private(pNEW) reduction(max:pChangeMax) collapse(3) gang vector
    do k=istart,iend  
		do j=1,ny
            do i=1,nx
				!!$acc cache (p(i,j-8:j+8,k-8:k+8))	
				if (mod(i+j+k,2) == 0) then
                           ! pNEW = ( p(i+1,j,k) * P_Den(i,j,k,1) &
                                  ! + p(i-1,j,k) * P_Den(i,j,k,2) &
                           pNEW = ( p(MOD(i,nx)+1,j,k) * P_Den(i,j,k,1) & ! Periodic BC
                                  + p(MOD(i-2+nx,nx)+1,j,k) * P_Den(i,j,k,2) & ! Periodic BC
                                  + p(i,j+1,k) * P_Den(i,j,k,3) &
                                  + p(i,j-1,k) * P_Den(i,j,k,4) &
                                  + p(i,j,k+1) * P_Den(i,j,k,5) &
                                  + p(i,j,k-1) * P_Den(i,j,k,6) &
                                  - mChange(i,j,k) ) * Den_inv(i,j,k)	

				pChangeMax=MAX(pChangeMax,ABS( omega*(pNEW - p(i,j,k))))
				p(i,j,k) = (1.d0 - omega)*p(i,j,k) + omega*pNEW	
				
				endif
            enddo       
        enddo
	enddo
	!$acc end parallel
    !$OMP END DO
	!$OMP END PARALLEL
	 
        call MPI_ALLREDUCE( pChangeMax, pChangeMax_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr )
	 
        if(myid==master .AND. istep <= 2)then
            write(71,*) REAL(pChangeMax_), ik, omega
        end if     
    
end do
!$acc end data



    if(myid==master)then
        write(71,*) 'Iterations GS',ccc,'=', ik , '  Convergence =', REAL(pChangeMax_)
    end if

    call pressure_boundary_conditions()


end subroutine RB_SOR




! Modified SOR 
subroutine SOR()
    use variables
    implicit none

! automatic omega upper & lower bounds, set according to ik from previous timestep
! the purpose is to avoid the use of high omega max once a stable #of P iteration is achieved
	if (ik > 500 .OR. istep <= 2) then	
		omega_max = 1.9d0
	elseif (ik > 50) then
		omega_max = 1.8d0
	else
		omega_max = 1.7d0		
	endif

	omega_min = 1.1d0
	d_omega = omega_max-omega_min


    ik=0
    pChangeMax = 1.d0
    pChangeMax_= 1.d0
	

    !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx

        mChange(i,j,k) =( ( u0(i,j,k) - u0(i-1,j,k) ) * iDy(j) * iDz(k) &
                        + ( v0(i,j,k) - v0(i,j-1,k) ) * iDx(i) * iDz(k) &
                        + ( w0(i,j,k) - w0(i,j,k-1) ) * iDx(i) * iDy(j) ) * den_flu * inv_dt

    enddo; enddo
	enddo
    !$OMP END PARALLEL DO



do while (pChangeMax_>zeta .AND. ik < itmax)

	! Setting dynamic omega response with respect to pChangeMax_
	if (pChangeMax_ < k_zeta_zeta) then	! omega decreases gradually from (k x zeta) to (zeta)
		omega = omega_min + d_omega * (LOG(pChangeMax_) - log_zeta)/log_k_zeta							! loglinear function (set parameters in variables module)
		! omega = omega_min + d_omega / (1.0d0 + EXP(-k_sigmoid * (LOG10(pChangeMax_) - x0_sigmoid)))	! Sigmoid function (set parameters in variables module)
	else
		omega = omega_max	! A constant of omega_max is used when pChangeMax_ > k_zeta*zeta
	endif
	
	
    ik=ik+1
    pChangeMax = 0.d0
    pChangeMax_= 0.d0

    !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
       
    icount = (nx+4)*(ny+4)

    itag = 201
    call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                       p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )

    itag = 202
    call MPI_SENDRECV( p(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
                       p(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    !$OMP PARALLEL DO PRIVATE(i,j,pNEW) REDUCTION(max : pChangeMAX) collapse(nclps)
    do k=istart,iend 		
        do j=1,ny
        
                    do i=2,nx-1
                    !do i=1,nx

                        !p_old(i,j,k)=p(i,j,k)

                        ! pNEW = ( p(MOD(i,nx)+1,j,k) * P_Den(i,j,k,1) & !Periodic boundary condition
                               ! + p(MOD(i-2+nx,nx)+1,j,k) * P_Den(i,j,k,2) & !Periodic boundary condition
                               ! + p(i,j+1,k) * P_Den(i,j,k,3) &
                               ! + p(i,j-1,k) * P_Den(i,j,k,4) &
                               ! + p(i,j,k+1) * P_Den(i,j,k,5) &
                               ! + p(i,j,k-1) * P_Den(i,j,k,6) &
                               ! - mChange(i,j,k) ) / Den(i,j,k)

                        pNEW = ( p(i+1,j,k) * P_Den(i,j,k,1) &
                               + p(i-1,j,k) * P_Den(i,j,k,2) &
                               + p(i,j+1,k) * P_Den(i,j,k,3) &
                               + p(i,j-1,k) * P_Den(i,j,k,4) &
                               + p(i,j,k+1) * P_Den(i,j,k,5) &
                               + p(i,j,k-1) * P_Den(i,j,k,6) &
                               - mChange(i,j,k) ) * Den_inv(i,j,k)

						pChangeMAX=MAX(pChangeMAX,ABS( omega*(pNew - p(i,j,k))))
						p(i,j,k) = (1.d0 - omega)*p(i,j,k) + omega*pNew
                    enddo


            !p_old(1,j,k)=p(1,j,k)

            pNEW = (  p(2,j,k)   * P_Den(1,j,k,1) &
                    + p(nx,j,k)  * P_Den(1,j,k,2) & ! Periodic boundary condition
                    + p(1,j+1,k) * P_Den(1,j,k,3) &
                    + p(1,j-1,k) * P_Den(1,j,k,4) &
                    + p(1,j,k+1) * P_Den(1,j,k,5) &
                    + p(1,j,k-1) * P_Den(1,j,k,6) &
                    - mChange(1,j,k) ) * Den_inv(1,j,k)
        
			pChangeMAX=MAX(pChangeMAX,ABS( omega*(pNew - p(1,j,k))))
			p(1,j,k) = (1.d0 - omega)*p(1,j,k) + omega*pNew
        
        
            !p_old(nx,j,k)=p(nx,j,k)

            pNEW = (  p(1,j,k)    * P_Den(nx,j,k,1) & ! Periodic boundary condition
                    + p(nx-1,j,k) * P_Den(nx,j,k,2) &
                    + p(nx,j+1,k) * P_Den(nx,j,k,3) &
                    + p(nx,j-1,k) * P_Den(nx,j,k,4) &
                    + p(nx,j,k+1) * P_Den(nx,j,k,5) &
                    + p(nx,j,k-1) * P_Den(nx,j,k,6) &
                    - mChange(nx,j,k) ) * Den_inv(nx,j,k)

			pChangeMAX=MAX(pChangeMAX,ABS( omega*(pNew - p(nx,j,k))))
			p(nx,j,k) = (1.d0 - omega)*p(nx,j,k) + omega*pNew
        
        enddo
	enddo
    !$OMP END PARALLEL DO


        call MPI_ALLREDUCE( pChangeMax, pChangeMax_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr )

        if(myid==master .AND. istep <= 2)then
            write(71,*) REAL(pChangeMax_), ik, omega
        end if
        
end do


    if(myid==master)then
        write(71,*) 'Iterations GS',ccc,'=', ik , '  Convergence =', REAL(pChangeMax_)
    end if

    call pressure_boundary_conditions()


end subroutine SOR






!*********Old SOR*********
! subroutine gauss_seidel()
    ! use variables
    ! implicit none



    ! ik=0
    ! pChangeMax = 1.d0
    ! pChangeMax_= 1.d0

    ! !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
    ! do k=istart,iend
    ! do j=1,ny; do i=1,nx

        ! mChange(i,j,k) =( ( u0(i,j,k) - u0(i-1,j,k) ) * iDy(j) * iDz(k) &
                        ! + ( v0(i,j,k) - v0(i,j-1,k) ) * iDx(i) * iDz(k) &
                        ! + ( w0(i,j,k) - w0(i,j,k-1) ) * iDx(i) * iDy(j) ) * den_flu * inv_dt

    ! enddo; enddo
	! enddo
    ! !$OMP END PARALLEL DO



! do while (pChangeMax_>zeta .AND. ik < itmax)

    ! ik=ik+1
    ! pChangeMax = 0.d0
    ! pChangeMax_= 0.d0

    ! !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
       
    ! icount = (nx+4)*(ny+4)

    ! itag = 201
    ! call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                       ! p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )

    ! itag = 202
    ! call MPI_SENDRECV( p(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
                       ! p(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    ! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



        ! !$OMP PARALLEL DO PRIVATE(i,j,pNEW) collapse(nclps)
        ! do k=istart,iend 
        ! do j=1,ny
        
                                ! do i=2,nx-1
                                ! !do i=1,nx

                                    ! p_old(i,j,k)=p(i,j,k)

                                        ! ! pNEW = ( p(MOD(i,nx)+1,j,k) * P_Den(i,j,k,1) & !Periodic boundary condition
                                              ! ! + p(MOD(i-2+nx,nx)+1,j,k) * P_Den(i,j,k,2) & !Periodic boundary condition
                                              ! ! + p(i,j+1,k) * P_Den(i,j,k,3) &
                                              ! ! + p(i,j-1,k) * P_Den(i,j,k,4) &
                                              ! ! + p(i,j,k+1) * P_Den(i,j,k,5) &
                                              ! ! + p(i,j,k-1) * P_Den(i,j,k,6) &
                                              ! ! - mChange(i,j,k) ) / Den(i,j,k)

                                        ! pNEW = ( p(i+1,j,k) * P_Den(i,j,k,1) &
                                               ! + p(i-1,j,k) * P_Den(i,j,k,2) &
                                               ! + p(i,j+1,k) * P_Den(i,j,k,3) &
                                               ! + p(i,j-1,k) * P_Den(i,j,k,4) &
                                               ! + p(i,j,k+1) * P_Den(i,j,k,5) &
                                               ! + p(i,j,k-1) * P_Den(i,j,k,6) &
                                               ! - mChange(i,j,k) ) * Den_inv(i,j,k)

                                    ! p(i,j,k) = p_old(i,j,k) + ( omega * (pNew - p_old(i,j,k)) )

                                ! enddo



            ! p_old(1,j,k)=p(1,j,k)

            ! pNEW = (  p(2,j,k)   * P_Den(1,j,k,1) &
                    ! + p(nx,j,k)  * P_Den(1,j,k,2) & ! Periodic boundary condition
                    ! + p(1,j+1,k) * P_Den(1,j,k,3) &
                    ! + p(1,j-1,k) * P_Den(1,j,k,4) &
                    ! + p(1,j,k+1) * P_Den(1,j,k,5) &
                    ! + p(1,j,k-1) * P_Den(1,j,k,6) &
                    ! - mChange(1,j,k) ) * Den_inv(1,j,k)
        
            ! p(1,j,k) = p_old(1,j,k) + ( omega * (pNew - p_old(1,j,k)) )
        
        
            ! p_old(nx,j,k)=p(nx,j,k)

            ! pNEW = (  p(1,j,k)    * P_Den(nx,j,k,1) & ! Periodic boundary condition
                    ! + p(nx-1,j,k) * P_Den(nx,j,k,2) &
                    ! + p(nx,j+1,k) * P_Den(nx,j,k,3) &
                    ! + p(nx,j-1,k) * P_Den(nx,j,k,4) &
                    ! + p(nx,j,k+1) * P_Den(nx,j,k,5) &
                    ! + p(nx,j,k-1) * P_Den(nx,j,k,6) &
                    ! - mChange(nx,j,k) ) * Den_inv(nx,j,k)

            ! p(nx,j,k) = p_old(nx,j,k) + ( omega * (pNew - p_old(nx,j,k)) )
        
        
        ! enddo
	! enddo
    ! !$OMP END PARALLEL DO

        ! do k=istart,iend; do j=1,ny; do i=1,nx
            ! pChangeMAX=MAX(pChangeMAX,ABS( p(i,j,k)-p_old(i,j,k) ) )
			! ! pChange=ABS( p(i,j,k)-p_old(i,j,k) )

            ! ! if (pChange > pChangeMax) then
                ! ! pChangeMax=pChange
            ! ! end if

        ! enddo; enddo; enddo

        ! call MPI_ALLREDUCE( pChangeMax, pChangeMax_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr )

        ! if(myid==master .AND. istep <= 2)then
            ! write(71,*) REAL(pChangeMax_), ik
        ! end if
        
! end do


    ! if(myid==master)then
        ! write(71,*) 'Iterations GS',ccc,'=', ik , '  Convergence =', REAL(pChangeMax_)
    ! end if

    ! call pressure_boundary_conditions()
    

! end subroutine gauss_seidel


