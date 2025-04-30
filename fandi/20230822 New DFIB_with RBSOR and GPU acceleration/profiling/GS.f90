! 22 Aug 2023 - FDS

! Red-black SOR
subroutine RB_SOR() 
    use variables
    implicit none


    ik=0
    pChangeMax = 1.0
    pChangeMax_= 1.0


    !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx

        mChange(i,j,k) =( ( u0(i,j,k) - u0(i-1,j,k) ) * iDy(j) * iDz(k) &
                        + ( v0(i,j,k) - v0(i,j-1,k) ) * iDx(i) * iDz(k) &
                        + ( w0(i,j,k) - w0(i,j,k-1) ) * iDx(i) * iDy(j) ) / dt

    enddo; enddo
	enddo
    !$OMP END PARALLEL DO


do while (pChangeMax_>zeta .AND. ik < itmax)

    ik=ik+1
    pChangeMax = 0.0
    pChangeMax_= 0.0

    !----------data transformation among nodes----------!
       
    icount = (nx+4)*(ny+4)

    itag = 250
    call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                       p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )

    itag = 260
    call MPI_SENDRECV( p(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
                       p(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !----------data transformation among nodes----------!


! Red-Black ordering

  !$OMP PARALLEL
	do rb=1,2

	!$OMP DO PRIVATE(i,j,pNEW) REDUCTION(max : pChangeMAX) collapse(nclps)
    do k=istart,iend 		
		do j=1,ny
            do i=1,nx
				if (mod(i+j+k,2) == rb - 1) then
				
                           pNEW = ( p(MOD(i,nx)+1,j,k) * P_Den(i,j,k,1) & !Periodic boundary condition
                                  + p(MOD(i-2+nx,nx)+1,j,k) * P_Den(i,j,k,2) & !Periodic boundary condition
                                  + p(i,j+1,k) * P_Den(i,j,k,3) &
                                  + p(i,j-1,k) * P_Den(i,j,k,4) &
                                  + p(i,j,k+1) * P_Den(i,j,k,5) &
                                  + p(i,j,k-1) * P_Den(i,j,k,6) &
                                  - mChange(i,j,k) ) * Den_inv(i,j,k)				

				pChangeMAX=DMAX1(pChangeMAX,DABS( omega*(pNEW - p(i,j,k))))
				p(i,j,k) = (1.d0 - omega)*p(i,j,k) + omega*pNEW
				  
				endif
    enddo; enddo
	enddo
    !$OMP END DO

	enddo
  !$OMP END PARALLEL


        call MPI_ALLREDUCE( pChangeMax, pChangeMax_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr )
	 
        if(myid==master .AND. istep <= 2)then
            write(71,*) REAL(pChangeMax_), ik
        end if     
    
end do


    if(myid==master)then
    write(71,*) 'Iterations GS',ccc,'=', ik , '  Convergence =', REAL(pChangeMax_)
    end if
	
    call pressure_boundary_conditions()
        

      !----------data transformation among nodes----------!
      icount = (nx+4)*(ny+4)

      itag = 270
      call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!   


end subroutine RB_SOR



subroutine RB_SOR_GPU() 
    use variables
	use mpi
    implicit none


    ik=0
    pChangeMax = 1.0

    !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx

        mChange(i,j,k) =( ( u0(i,j,k) - u0(i-1,j,k) ) * iDy(j) * iDz(k) &
                        + ( v0(i,j,k) - v0(i,j-1,k) ) * iDx(i) * iDz(k) &
                        + ( w0(i,j,k) - w0(i,j,k-1) ) * iDx(i) * iDy(j) ) / dt

    enddo; enddo
	enddo
    !$OMP END PARALLEL DO
	

    !----------data transfer to master----------!
       
      icount = igcount*(nx+4)*(ny+4)

      if(myid>master)then
         itag = 250
         call MPI_SEND( p(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 260
         call MPI_SEND( mChange(1,1,istart), igcount*(nx)*(ny), MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if
      !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !Wait to receive results from each task
      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount(i)*(nx+4)*(ny+4)
            itag = 250
            call MPI_RECV( p(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr ) 
			itag = 260
			call MPI_RECV( mChange(1,1,gstart(i)), gcount(i)*(nx)*(ny), MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr ) 			
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
    !----------data transfer to master----------!


if(myid==master)then

! Sending data to accelerator
!$acc data copy(p) copyin(P_Den,mChange,Den_inv) create(pNEW)
do while (pChangeMax>zeta .AND. ik < itmax)

    ik=ik+1
    pChangeMax = 0.0

! Red-Black ordering (RED)


	!$acc parallel
	!$acc loop independent private(pNEW) reduction(max:pChangeMAX) collapse(3)
    !do k=istart,iend 
    do k=1,nz 
        !!$acc loop private(pNEW) reduction(max:pChangeMAX)
		do j=1,ny
			!!$acc loop private(pNEW) reduction(max:pChangeMAX)
            do i=1,nx
				if (mod(i+j+k,2) == 0) then

                           pNEW = ( p(MOD(i,nx)+1,j,k) * P_Den(i,j,k,1) & ! Periodic BC
                                  + p(MOD(i-2+nx,nx)+1,j,k) * P_Den(i,j,k,2) & ! Periodic BC
                                  + p(i,j+1,k) * P_Den(i,j,k,3) &
                                  + p(i,j-1,k) * P_Den(i,j,k,4) &
                                  + p(i,j,k+1) * P_Den(i,j,k,5) &
                                  + p(i,j,k-1) * P_Den(i,j,k,6) &
                                  - mChange(i,j,k) ) * Den_inv(i,j,k)	

				pChangeMAX=DMAX1(pChangeMAX,DABS( omega*(pNEW - p(i,j,k))))
				p(i,j,k) = (1.d0 - omega)*p(i,j,k) + omega*pNEW			
				endif
            enddo
        enddo
	enddo
	!$acc end parallel


! Red-Black ordering (BLACK)


	!$acc parallel
	!$acc loop independent private(pNEW) reduction(max:pChangeMAX) collapse(3)
    !do k=istart,iend 
    do k=1,nz 
        !!$acc loop private(pNEW) reduction(max:pChangeMAX)
		do j=1,ny
			!!$acc loop private(pNEW) reduction(max:pChangeMAX)
            !do i=2,nx-1
            do i=1,nx
				if (mod(i+j+k,2) == 1) then

                           pNEW = ( p(MOD(i,nx)+1,j,k) * P_Den(i,j,k,1) & ! Periodic BC
                                  + p(MOD(i-2+nx,nx)+1,j,k) * P_Den(i,j,k,2) & ! Periodic BC
                                  + p(i,j+1,k) * P_Den(i,j,k,3) &
                                  + p(i,j-1,k) * P_Den(i,j,k,4) &
                                  + p(i,j,k+1) * P_Den(i,j,k,5) &
                                  + p(i,j,k-1) * P_Den(i,j,k,6) &
                                  - mChange(i,j,k) ) * Den_inv(i,j,k)	

				pChangeMAX=DMAX1(pChangeMAX,DABS( omega*(pNEW - p(i,j,k))))
				p(i,j,k) = (1.d0 - omega)*p(i,j,k) + omega*pNEW	
				
				endif
            enddo       
        enddo
	enddo
	!$acc end parallel

	 
        if(istep <= 2)then
            write(71,*) REAL(pChangeMax), ik
        end if
        
        
    
end do
!$acc end data

    write(71,*) 'Iterations GS',ccc,'=', ik , '  Convergence =', REAL(pChangeMax)

    call pressure_boundary_conditions()
        
endif

	  

      !----------data transfer from master to all nodes----------!

      if(myid==master) then
         do i = 1, (nproc-1)
            icount = (gcount(i)+1)*(nx+4)*(ny+4)
            itag = 270
            call MPI_SEND( p(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, ierr ) 		
         end do
	  else
		icount = (igcount+1)*(nx+4)*(ny+4)
		itag = 270
        call MPI_RECV( p(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, status, ierr )
      end if
       
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !Wait to receive results from each task

      
      !----------data transfer from master to all nodes----------!
	  

      ! icount = (nz+4)*(nx+4)*(ny+4)	  
      ! call MPI_BCAST ( p, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)

end subroutine RB_SOR_GPU



! Modified SOR 
subroutine SOR()
    use variables
    implicit none


    ik=0
    pChangeMax = 1.0
    pChangeMax_= 1.0

    !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx

        mChange(i,j,k) =( ( u0(i,j,k) - u0(i-1,j,k) ) * iDy(j) * iDz(k) &
                        + ( v0(i,j,k) - v0(i,j-1,k) ) * iDx(i) * iDz(k) &
                        + ( w0(i,j,k) - w0(i,j,k-1) ) * iDx(i) * iDy(j) ) / dt

    enddo; enddo
	enddo
    !$OMP END PARALLEL DO



do while (pChangeMax_>zeta .AND. ik < itmax)

    ik=ik+1
    pChangeMax = 0.0
    pChangeMax_= 0.0

    !----------data transformation among nodes----------!
       
    icount = (nx+4)*(ny+4)

    itag = 250
    call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                       p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )

    itag = 260
    call MPI_SENDRECV( p(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
                       p(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !----------data transformation among nodes----------!


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

						pChangeMAX=DMAX1(pChangeMAX,DABS( omega*(pNew - p(i,j,k))))
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
        
			pChangeMAX=DMAX1(pChangeMAX,DABS( omega*(pNew - p(1,j,k))))
			p(1,j,k) = (1.d0 - omega)*p(1,j,k) + omega*pNew
        
        
            !p_old(nx,j,k)=p(nx,j,k)

            pNEW = (  p(1,j,k)    * P_Den(nx,j,k,1) & ! Periodic boundary condition
                    + p(nx-1,j,k) * P_Den(nx,j,k,2) &
                    + p(nx,j+1,k) * P_Den(nx,j,k,3) &
                    + p(nx,j-1,k) * P_Den(nx,j,k,4) &
                    + p(nx,j,k+1) * P_Den(nx,j,k,5) &
                    + p(nx,j,k-1) * P_Den(nx,j,k,6) &
                    - mChange(nx,j,k) ) * Den_inv(nx,j,k)

			pChangeMAX=DMAX1(pChangeMAX,DABS( omega*(pNew - p(nx,j,k))))
			p(nx,j,k) = (1.d0 - omega)*p(nx,j,k) + omega*pNew
        
        enddo
	enddo
    !$OMP END PARALLEL DO


        call MPI_ALLREDUCE( pChangeMax, pChangeMax_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr )

        if(myid==master .AND. istep <= 2)then
            write(71,*) REAL(pChangeMax_), ik
        end if
        
end do


    if(myid==master)then
        write(71,*) 'Iterations GS',ccc,'=', ik , '  Convergence =', REAL(pChangeMax_)
    end if

    call pressure_boundary_conditions()
        

      !----------data transformation among nodes----------!
      icount = (nx+4)*(ny+4)

      itag = 270
      call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!  

end subroutine SOR






!*********Old SOR*********
subroutine gauss_seidel()
    use variables
    implicit none



    ik=0
    pChangeMax = 1.0
    pChangeMax_= 1.0

    !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
    do k=istart,iend
    do j=1,ny; do i=1,nx

        mChange(i,j,k) =( ( u0(i,j,k) - u0(i-1,j,k) ) * iDy(j) * iDz(k) &
                        + ( v0(i,j,k) - v0(i,j-1,k) ) * iDx(i) * iDz(k) &
                        + ( w0(i,j,k) - w0(i,j,k-1) ) * iDx(i) * iDy(j) ) / dt

    enddo; enddo
	enddo
    !$OMP END PARALLEL DO



do while (pChangeMax_>zeta .AND. ik < itmax)

    ik=ik+1
    pChangeMax = 0.0
    pChangeMax_= 0.0

    !----------data transformation among nodes----------!
       
    icount = (nx+4)*(ny+4)

    itag = 250
    call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                       p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )

    itag = 260
    call MPI_SENDRECV( p(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
                       p(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !----------data transformation among nodes----------!



        !$OMP PARALLEL DO PRIVATE(i,j,pNEW) collapse(nclps)
        do k=istart,iend 
        do j=1,ny
        
                                do i=2,nx-1
                                !do i=1,nx

                                    p_old(i,j,k)=p(i,j,k)

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

                                    p(i,j,k) = p_old(i,j,k) + ( omega * (pNew - p_old(i,j,k)) )

                                enddo



            p_old(1,j,k)=p(1,j,k)

            pNEW = (  p(2,j,k)   * P_Den(1,j,k,1) &
                    + p(nx,j,k)  * P_Den(1,j,k,2) & ! Periodic boundary condition
                    + p(1,j+1,k) * P_Den(1,j,k,3) &
                    + p(1,j-1,k) * P_Den(1,j,k,4) &
                    + p(1,j,k+1) * P_Den(1,j,k,5) &
                    + p(1,j,k-1) * P_Den(1,j,k,6) &
                    - mChange(1,j,k) ) * Den_inv(1,j,k)
        
            p(1,j,k) = p_old(1,j,k) + ( omega * (pNew - p_old(1,j,k)) )
        
        
            p_old(nx,j,k)=p(nx,j,k)

            pNEW = (  p(1,j,k)    * P_Den(nx,j,k,1) & ! Periodic boundary condition
                    + p(nx-1,j,k) * P_Den(nx,j,k,2) &
                    + p(nx,j+1,k) * P_Den(nx,j,k,3) &
                    + p(nx,j-1,k) * P_Den(nx,j,k,4) &
                    + p(nx,j,k+1) * P_Den(nx,j,k,5) &
                    + p(nx,j,k-1) * P_Den(nx,j,k,6) &
                    - mChange(nx,j,k) ) * Den_inv(nx,j,k)

            p(nx,j,k) = p_old(nx,j,k) + ( omega * (pNew - p_old(nx,j,k)) )
        
        
        enddo
	enddo
    !$OMP END PARALLEL DO

        do k=istart,iend; do j=1,ny; do i=1,nx
            pChangeMAX=DMAX1(pChangeMAX,DABS( p(i,j,k)-p_old(i,j,k) ) )
			! pChange=DABS( p(i,j,k)-p_old(i,j,k) )

            ! if (pChange > pChangeMax) then
                ! pChangeMax=pChange
            ! end if

        enddo; enddo; enddo

        call MPI_ALLREDUCE( pChangeMax, pChangeMax_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr )

        if(myid==master .AND. istep <= 2)then
            write(71,*) REAL(pChangeMax_), ik
        end if
        
end do


    if(myid==master)then
        write(71,*) 'Iterations GS',ccc,'=', ik , '  Convergence =', REAL(pChangeMax_)
    end if

    call pressure_boundary_conditions()
        

      !----------data transformation among nodes----------!
      icount = (nx+4)*(ny+4)

      itag = 270
      call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!  
	  

end subroutine gauss_seidel


