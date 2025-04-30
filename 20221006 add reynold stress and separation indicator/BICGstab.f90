subroutine BICG_stab()
    use variables
    implicit none


    call write_matrix_b()


    !$OMP PARALLEL DO PRIVATE(c)
        do i=1,s
            r0(i)=bm(i)
            c=(i-1)*7+1
            do while ( INT(aa(1,c)) == i )    
                r0(i) = r0(i) - aa(3,c) * xm( INT(aa(2,c)) )
                c=c+1
            enddo
            rm(i)=r0(i)
        enddo
    !$OMP END PARALLEL DO



    !$OMP PARALLEL DO 
        do i=1,s;   p0(i)=rm(i)
                    ap(i)=rm(i)  ; enddo
    !$OMP END PARALLEL DO

    alpha=1.
    om=1.
    gama=1.


    ik=0
    pChangeMax_= 1.0

    !-----------------------------CONJUGATE LOOP-------------------------------------
    do while (pChangeMax_>zeta .AND. ik < itmax)  
        ik=ik+1


        sub=0.d0
        do k=istart,iend
        !$OMP PARALLEL DO REDUCTION(+:sub)
        do i=(k-1)*nxy+1,k*nxy  
                xmp(i)=xm(i)
                sub=sub+r0(i)*rm(i)
        enddo
        !$OMP END PARALLEL DO
        enddo

        !----------------------------------ALLGATHER matrix(i)----------------------------------
        call MPI_ALLREDUCE(sub, gama, nproc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !----------------------------------ALLGATHER matrix(i)----------------------------------
        gamap=gama

        call create_ap()

        sub=0.d0
        do k=istart,iend
        !$OMP PARALLEL DO REDUCTION(+:sub)
        do i=(k-1)*nxy+1,k*nxy  
                sub=sub+ap(i)*r0(i)
        enddo
        !$OMP END PARALLEL DO
        enddo

        !----------------------------------ALLGATHER matrix(i)----------------------------------
        call MPI_ALLREDUCE(sub, apr, nproc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !----------------------------------ALLGATHER matrix(i)----------------------------------

        alpha = gama/apr

        do k=istart,iend
        !$OMP PARALLEL DO 
        do i=(k-1)*nxy+1,k*nxy
                ss(i)=rm(i)-alpha*ap(i)
        enddo
        !$OMP END PARALLEL DO
        enddo

        !----------------------------------ALLGATHER matrix(i)----------------------------------
        if(myid>master)then
            call MPI_SEND( ss((istart-1)*nxy+1), igcount*nxy, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
        end if

        if(myid==master)then
            do i = 1, (nproc-1)
                call MPI_RECV( ss((gstart(i)-1)*nxy+1), igcount*nxy, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            end do
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST( ss(1), s, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !----------------------------------ALLGATHER matrix(i)----------------------------------
        
        call create_as()

        sub=0.d0
        do k=istart,iend
        !$OMP PARALLEL DO REDUCTION(+:sub)
        do i=(k-1)*nxy+1,k*nxy  
                sub=sub+as(i)*ss(i)
        enddo
        !$OMP END PARALLEL DO
        enddo

        !----------------------------------ALLGATHER matrix(i)----------------------------------
        call MPI_ALLREDUCE(sub, ass, nproc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !----------------------------------ALLGATHER matrix(i)----------------------------------


        sub=0.d0
        do k=istart,iend
        !$OMP PARALLEL DO REDUCTION(+:sub)
        do i=(k-1)*nxy+1,k*nxy  
                sub=sub+as(i)*as(i)
        enddo
        !$OMP END PARALLEL DO
        enddo

        !----------------------------------ALLGATHER matrix(i)----------------------------------
        call MPI_ALLREDUCE(sub, asas, nproc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !----------------------------------ALLGATHER matrix(i)----------------------------------

        om=ass/asas

        do k=istart,iend
        !$OMP PARALLEL DO 
        do i=(k-1)*nxy+1,k*nxy
                xm(i)=xm(i)+alpha*p0(i)+om*ss(i)
        enddo
        !$OMP END PARALLEL DO
        enddo

        !----------------------------------ALLGATHER matrix(i)----------------------------------
        if(myid>master)then
            call MPI_SEND( xm((istart-1)*nxy+1), igcount*nxy, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
        end if

        if(myid==master)then
            do i = 1, (nproc-1)
                call MPI_RECV( xm((gstart(i)-1)*nxy+1), igcount*nxy, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            end do
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST( xm(1), s, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !----------------------------------ALLGATHER matrix(i)----------------------------------
        

        do k=istart,iend
        !$OMP PARALLEL DO 
        do i=(k-1)*nxy+1,k*nxy
                rm(i)=om*as(i)
                ss(i)=rm(i)
        enddo
        !$OMP END PARALLEL DO
        enddo

        sub=0.d0
        do k=istart,iend
        !$OMP PARALLEL DO REDUCTION(+:sub)
        do i=(k-1)*nxy+1,k*nxy  
                sub=sub+r0(i)*rm(i)
        enddo
        !$OMP END PARALLEL DO
        enddo

        !----------------------------------ALLGATHER matrix(i)----------------------------------
        call MPI_ALLREDUCE(sub, gama, nproc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !----------------------------------ALLGATHER matrix(i)----------------------------------
        
       
        beta=gama*alpha/gamap/om

        do k=istart,iend
        !$OMP PARALLEL DO 
        do i=(k-1)*nxy+1,k*nxy
                p0(i)=rm(i)+beta*(p0(i)-om*ap(i))
        enddo
        !$OMP END PARALLEL DO
        enddo

        !----------------------------------ALLGATHER matrix(i)----------------------------------
        if(myid>master)then
            call MPI_SEND( p0((istart-1)*nxy+1), igcount*nxy, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
        end if

        if(myid==master)then
            do i = 1, (nproc-1)
                call MPI_RECV( p0((gstart(i)-1)*nxy+1), igcount*nxy, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            end do
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST( p0(1), s, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !----------------------------------ALLGATHER matrix(i)----------------------------------
        


        pChangeMax=0.
        do k=istart,iend
        do i=(k-1)*nxy+1,k*nxy
            if ( ABS(xm(i)-xmp(i)) > pChangeMax ) then
                pChangeMax=ABS(xm(i)-xmp(i))
            endif
        enddo; enddo
        
        call MPI_ALLREDUCE( pChangeMax, pChangeMax_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        if(myid==master .AND. istep <= 20)then
                write(*,*) REAL(pChangeMax_), ik
        end if

    enddo


    !----------------------------------ALLGATHER matrix(i)----------------------------------
        if(myid>master)then
            call MPI_SEND( xm((istart-1)*nxy+1), igcount*nxy, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
        end if

        if(myid==master)then
            do i = 1, (nproc-1)
                call MPI_RECV( xm((gstart(i)-1)*nxy+1), igcount*nxy, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            end do
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST( xm(1), s, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !----------------------------------ALLGATHER matrix(i)----------------------------------

    call update_pressure()

    call pressure_boundary_conditions()
    
end subroutine





subroutine write_matrix_a()
    use variables
    implicit none
    real*8:: ano,aso,aea,awe,afr,aba
    integer:: n, cso
    aa=0.

    !$OMP PARALLEL DO PRIVATE(c,i,j,k,ano,aso,aea,awe,afr,aba,cso)
        do n=1,s
            c=(n-1)*7+1

            ano=1.d0
            aso=1.d0
            aea=1.d0
            awe=1.d0
            afr=1.d0
            aba=1.d0

            i=MOD(INT(n-1),nx)+1
            j=MOD(INT((n-1)/nx),ny)+1
            k=INT((n-1)/nxy)+1
            
            if ( k==1 ) then;             aba=0.d0
            else                  
                aa(1,c)= n
                aa(2,c)= n-nxy
                aa(3,c)= -1.d0*iDx(i) * iDy(j) / Dzs(k-1)
                c=c+1
            endif

            if ( j==1 ) then;             aso=0.d0
            else                  
                aa(1,c)= n
                aa(2,c)= n-nx
                aa(3,c)= -1.d0*iDx(i) * iDz(k) / Dys(j-1)
                c=c+1
            endif

            if ( i==1 ) then;             !awe=0
                aa(1,c)= n
                aa(2,c)= n+nx-1
                aa(3,c)= -1.d0*iDy(j) * iDz(k) / Dxs(i-1)
                c=c+1
            else             
                aa(1,c)= n
                aa(2,c)= n-1
                aa(3,c)= -1.d0*iDy(j) * iDz(k) / Dxs(i-1)
                c=c+1
            endif

            cso=c
            c=c+1

            if ( i==nx ) then;            !aea=0
                aa(1,c)= n
                aa(2,c)= n-nx+1
                aa(3,c)= -1.d0*iDy(j) * iDz(k) / Dxs(i)
                c=c+1
            else                   
                aa(1,c)= n
                aa(2,c)= n+1
                aa(3,c)= -1.d0*iDy(j) * iDz(k) / Dxs(i)
                c=c+1
            endif

            if ( j==ny ) then;            ano=0.d0
            else                 
                aa(1,c)= n
                aa(2,c)= n+nx
                aa(3,c)= -1.d0*iDx(i) * iDz(k) / Dys(j)
                c=c+1
            endif

            if ( k==nz ) then;           afr=0.d0
            else                  
                aa(1,c)= n
                aa(2,c)= n+nxy
                aa(3,c)= -1.d0*iDx(i) * iDy(j) / Dzs(k)
            endif


            aa(1,cso)= n
            aa(2,cso)= n
            aa(3,cso)= 1.d0*( aea * iDy(j) * iDz(k) / Dxs(i)   &
                            + awe * iDy(j) * iDz(k) / Dxs(i-1) &
                            + ano * iDx(i) * iDz(k) / Dys(j)   &
                            + aso * iDx(i) * iDz(k) / Dys(j-1) &
                            + afr * iDx(i) * iDy(j) / Dzs(k)   &
                            + aba * iDx(i) * iDy(j) / Dzs(k-1) )

        enddo
    !$OMP END PARALLEL DO

    !write(*,*) aa
    
end subroutine




subroutine write_matrix_b()
    use variables
    implicit none

    !$OMP PARALLEL DO PRIVATE(j,i)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                bm((k-1)*nxy+(j-1)*nx+i)=-1.d0*( ( u0(i,j,k) - u0(i-1,j,k) ) * iDy(j) * iDz(k) &
                                               + ( v0(i,j,k) - v0(i,j-1,k) ) * iDx(i) * iDz(k) &
                                               + ( w0(i,j,k) - w0(i,j,k-1) ) * iDx(i) * iDy(j) ) /dt
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

end subroutine





subroutine create_ap()
    use variables
    implicit none

    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(c)
    do i=(k-1)*nxy+1,k*nxy
        ap(i)=0.0d0
        c=(i-1)*7+1
        do while ( INT(aa(1,c)) == i )    
            ap(i)=ap(i) + aa(3,c) * p0( INT(aa(2,c)) )
            c=c+1
        enddo 
    enddo
    !$OMP END PARALLEL DO
    enddo

end subroutine





subroutine create_as()
    use variables
    implicit none

    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(c)
    do i=(k-1)*nxy+1,k*nxy
        as(i)=0.0d0
        c=(i-1)*7+1
        do while ( INT(aa(1,c)) == i )    
            as(i)=as(i) + aa(3,c) * ss( INT(aa(2,c)) )
            c=c+1
        enddo 
    enddo
    !$OMP END PARALLEL DO
    enddo

end subroutine





subroutine update_pressure()
    use variables
    implicit none

    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(i)
    do j=1,ny; do i=1,nx
            p(i,j,k)=xm((k-1)*nx*ny+(j-1)*nx+i)
    enddo; enddo
    !$OMP END PARALLEL DO
    enddo


    if(myid==master)then
        write(*,*) 'Iterations BICG =',ik,'     pChangeMax = ',REAL(pChangeMax_)
    end if


    do k=istart,iend
    !$OMP PARALLEL DO PRIVATE(i)
    do j=1,ny; do i=1,nx

        pre(i,j,k) = p(i,j,k)

    enddo; enddo
    !$OMP END PARALLEL DO
    enddo

end subroutine
