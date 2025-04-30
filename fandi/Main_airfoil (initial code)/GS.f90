subroutine gauss_seidel() 
    use variables
    implicit none



    ik=0
    pChangeMax = 1.0
    pChangeMax_= 1.0

    !$OMP PARALLEL DO PRIVATE(i,j)  
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

     
        !$OMP PARALLEL DO PRIVATE(i,j,pNEW)  
        do k=istart,iend
        do j=1,ny
        
                                do i=2,nx-1
                                !do i=1,nx

                                    p_old(i,j,k)=p(i,j,k)

                                        !pNEW = ( p(MOD(i,nx)+1,j,k) * P_Den(i,j,k,1) &
                                        !       + p(MOD(i-2+nx,nx)+1,j,k) * P_Den(i,j,k,2) &
                                        !       + p(i,j+1,k) * P_Den(i,j,k,3) &
                                        !       + p(i,j-1,k) * P_Den(i,j,k,4) &
                                        !       + p(i,j,k+1) * P_Den(i,j,k,5) &
                                        !       + p(i,j,k-1) * P_Den(i,j,k,6) &
                                        !       - mChange(i,j,k) ) / Den(i,j,k)

                                        pNEW = ( p(i+1,j,k) * P_Den(i,j,k,1) &
                                               + p(i-1,j,k) * P_Den(i,j,k,2) &
                                               + p(i,j+1,k) * P_Den(i,j,k,3) &
                                               + p(i,j-1,k) * P_Den(i,j,k,4) &
                                               + p(i,j,k+1) * P_Den(i,j,k,5) &
                                               + p(i,j,k-1) * P_Den(i,j,k,6) &
                                               - mChange(i,j,k) ) / Den(i,j,k)

                                    p(i,j,k) = p_old(i,j,k) + ( omega * (pNew - p_old(i,j,k)) )

                                enddo



            p_old(1,j,k)=p(1,j,k)

            pNEW = (  p(2,j,k)   * P_Den(1,j,k,1) &
                    + p(nx,j,k)  * P_Den(1,j,k,2) &
                    + p(1,j+1,k) * P_Den(1,j,k,3) &
                    + p(1,j-1,k) * P_Den(1,j,k,4) &
                    + p(1,j,k+1) * P_Den(1,j,k,5) &
                    + p(1,j,k-1) * P_Den(1,j,k,6) &
                    - mChange(1,j,k) ) / Den(1,j,k)
        
            p(1,j,k) = p_old(1,j,k) + ( omega * (pNew - p_old(1,j,k)) )
        
        
            p_old(nx,j,k)=p(nx,j,k)

            pNEW = (  p(1,j,k)    * P_Den(nx,j,k,1) &
                    + p(nx-1,j,k) * P_Den(nx,j,k,2) &
                    + p(nx,j+1,k) * P_Den(nx,j,k,3) &
                    + p(nx,j-1,k) * P_Den(nx,j,k,4) &
                    + p(nx,j,k+1) * P_Den(nx,j,k,5) &
                    + p(nx,j,k-1) * P_Den(nx,j,k,6) &
                    - mChange(nx,j,k) ) / Den(nx,j,k)

            p(nx,j,k) = p_old(nx,j,k) + ( omega * (pNew - p_old(nx,j,k)) )
        
        
        enddo
        enddo
        !$OMP END PARALLEL DO



        do k=istart,iend; do j=1,ny; do i=1,nx
            pChange=ABS( p(i,j,k)-p_old(i,j,k) )

            if (pChange > pChangeMax) then
                pChangeMax=pChange
            end if

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
        


end subroutine gauss_seidel
