subroutine check_steady()
use variables
implicit none


        VelocityDifference = 0

        do k=1,nz; do j=1,ny; do i=1,nx
            if(abs( w(i,j,k)- last_velocity(i,j,k) ) > VelocityDifference ) then
                VelocityDifference =  abs( w(i,j,k)- last_velocity(i,j,k) )
            end if
        end do; end do; end do

        if(myid==master)then
            open (16,file='velocity_difference.dat',position='append')
            write(16,*) REAL(time),REAL(VelocityDifference)
            close(16)
        end if


    do k=1,nz; do j=1,ny; do i=1,nx
        last_velocity(i,j,k) = w(i,j,k)
    end do; end do; end do

end subroutine check_steady
