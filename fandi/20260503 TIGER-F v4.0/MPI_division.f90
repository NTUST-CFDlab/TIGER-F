subroutine Mpi_division()
use variables
implicit none

allocate(    gstart(0:nproc-1),     gend(0:nproc-1),     gcount(0:nproc-1))
allocate(gstart_vos(0:nproc-1), gend_vos(0:nproc-1), gcount_vos(0:nproc-1))

    Zdv = (nz) / nproc
    Zr  = (nz) - Zdv * nproc 
    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !i = myid
    do i=0,(nproc-1)

        if(i < Zr) then
            gstart(i) = 1 + i * (Zdv+1)
            gend(i) = gstart(i) + Zdv
        else
            gstart(i) = 1 + i * Zdv + Zr
            gend(i) = gstart(i) + Zdv - 1
        end if
        
        gcount(i) = gend(i) - gstart(i) + 1

    end do
    
end subroutine Mpi_division




subroutine Mpi_division_unbalance2gpu()
use variables 
implicit none

allocate(    gstart(0:nproc-1),     gend(0:nproc-1),     gcount(0:nproc-1))
allocate(gstart_vos(0:nproc-1), gend_vos(0:nproc-1), gcount_vos(0:nproc-1))

Zdv = 0.55*nz			! proportion of the first process (master)

	if(nproc==1) then
		gstart(0) = 1
		gend(0)  = nz
        gcount(0) = nz
	elseif(nproc==2) then
		gstart(0) = 1
		gend(0)  = Zdv
		gstart(1) = Zdv+1
		gend(1)  = nz
        gcount(0) = gend(0) - gstart(0) + 1
        gcount(1) = gend(1) - gstart(1) + 1
    end if
        

end subroutine Mpi_division_unbalance2gpu
