! 22 Aug 2023 - FDS
subroutine gridder_equal()
    use variables
    implicit none

    
    !-----------------unequal grid intervals----------------!
    do i=1,nx+3
        if(i == 1) then
            X(i) = 0.0
            X(i-1) = X(i) - dx
            X(i-2) = X(i-1) - dx
        else
            X(i) = X(i-1) + dx
        end if 
    end do

    do j=1,ny+3
        if(j == 1) then
            Y(j) = 0.0
            Y(j-1) = Y(j) - dy
            Y(j-2) = Y(j-1) - dy
        else
            Y(j) = Y(j-1) + dy
        end if 
    end do

    do k=1,nz+3
        if(k == 1) then
            Z(k) = 0.0
            Z(k-1) = Z(k) - dz
            Z(k-2) = Z(k-1) - dz
        else
            Z(k) = Z(k-1) + dz
        end if 
    end do
    !-----------------unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=1,nx-1

        iDx(i) = ( X(i+1) - X(i) )
        Dxs(i) = ( X(i+2) - X(i) ) *0.5d0

    end do

    do j=1,ny-1

        iDy(j) = ( Y(j+1) - Y(j) )
        Dys(j) = ( Y(j+2) - Y(j) ) *0.5d0

    end do

    do k=1,nz-1

        iDz(k) = ( Z(k+1) - Z(k) )
        Dzs(k) = ( Z(k+2) - Z(k) ) *0.5d0

    end do



    !Ghost boundary grid lengths
    iDx(0) = iDx(1)
    iDx(-1) = iDx(1)
    iDx(nx) = X(nx+1) - X(nx)
    iDx(nx+1) = iDx(nx)
    iDx(nx+2) = iDx(nx)

    Dxs(0) = Dxs(1)
    Dxs(-1) = Dxs(1)
    Dxs(nx) = Dxs(nx-1)
    Dxs(nx+1) = Dxs(nx-1)
    Dxs(nx+2) = Dxs(nx-1)


    iDy(0) = iDy(1)
    iDy(-1) = iDy(1)
    iDy(ny) = Y(ny+1) - Y(ny)
    iDy(ny+1) = iDy(ny)
    iDy(ny+2) = iDy(ny)

    Dys(0) = Dys(1)
    Dys(-1) = Dys(1)
    Dys(ny) = Dys(ny-1)
    Dys(ny+1) = Dys(ny-1)
    Dys(ny+2) = Dys(ny-1)


    iDz(0) = iDz(1)
    iDz(-1) = iDz(1)
    iDz(nz) = Z(nz+1) - Z(nz)
    iDz(nz+1) = iDz(nz)
    iDz(nz+2) = iDz(nz)

    Dzs(0) = Dzs(1)
    Dzs(-1) = Dzs(1)
    Dzs(nz) = Dzs(nz-1)
    Dzs(nz+1) = Dzs(nz-1)
    Dzs(nz+2) = Dzs(nz-1)

    !Modifying the index of X, Y and Z arrays to represent the actual grid
    !do i=1,nx+1
    !    Xa(i) = X(i)
    !end do
!
    !do j=1,ny+1
    !    Ya(j) = Y(j)
    !end do
!
    !do k=1,nz+1
    !    Za(k) = Z(k)
    !end do

    !Defining the midpoint values of the grids
	do i=1,nx
        Xs(i) = 0.5d0 * ( X(i+1) + X(i) )
		if ((X(i) .LE. x_start) .AND. (X(i+1) .GT. x_start) ) then 	! filer boundary
			i_start = i
		elseif (x_end .GE. lx) then
			i_end = nx
		elseif ((X(i) .LT. x_end) .AND. (X(i+1) .GE. x_end) ) then
			i_end = i
		endif
    end do

	do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
		if ((Y(j) .LE. y_start) .AND. (Y(j+1) .GT. y_start) ) then 	! filer boundary
			j_start = j
		elseif (y_end .GE. ly) then
			j_end = ny
		elseif ((Y(j) .LT. y_end) .AND. (Y(j+1) .GE. y_end) ) then
			j_end = j
		endif
    end do
	
	do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
		if ((Z(k) .LE. z_start) .AND. (Z(k+1) .GT. z_start) ) then 	! filer boundary
			k_start = k
		elseif (z_end .GE. lz) then
			k_end = nz
		elseif ((Z(k) .LT. z_end) .AND. (Z(k+1) .GE. z_end) ) then
			k_end = k
		endif
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
        Yout(i,j,k) = Ys(j)
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo


    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output3D/'))
        open (unit=1,form='unformatted',file=TRIM('output3D/')//'3d_mesh.x')
        write(1) nblocks
        write(1) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1

        write(1)    (((Xout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Yout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Zout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end)

        close(1)
    end if
	
    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output2D/'))
        open (unit=2,form='unformatted',file=TRIM('output2D/')//'2d_mesh.x')
        write(2) nblocks
        write(2) 1, ny, nz

        write(2)    (((Xout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,1),j=1,ny),k=1,nz)

        close(2)
    end if


end subroutine gridder_equal

!-------------------------MODIFIED-FANDI-------------------------------------!

subroutine gridder_sin4() ! Grid spacing method proposed by Kuyper
use variables
implicit none
real*8 :: nratio

real*8, parameter :: lyLrg1 = GridderYc - lySml*0.5d0
real*8, parameter :: lyLrg2 = ly - lyLrg1 - lySml

integer, parameter :: nyLrg1 = nyLow !INT( (ny-nySml)/(ly-lySml)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nySml

real*8, parameter :: dyLrg1 = lyLrg1/nyLrg1*1.d0
real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: lzLrg1 = GridderZc - lzSml*0.5d0
real*8, parameter :: lzLrg2 = lz - lzLrg1 - lzSml

integer, parameter :: nzLrg1 = nzUps !INT( (nz-nzSml)/(lz-lzSml)*lzLrg1 )  !1 
integer, parameter :: nzLrg2 = nz-nzLrg1-nzSml                     !89

real*8, parameter :: dzLrg1 = lzLrg1/nzLrg1*1.d0
real*8, parameter :: dzLrg2 = lzLrg2/nzLrg2*1.d0

real*8, parameter :: thetaa = 4.d0*ATAN(1.d0)
real*8, parameter :: tune = 1.d0



!----------------Unequal grid intervals----------------!
    do i=1,nx+3
        if(i == 1) then
            X(i) = 0.0
            X(i-1) = X(i) - dx
            X(i-2) = X(i-1) - dx
        else
            X(i) = X(i-1) + dx
        end if 
    end do


    do j=1,ny+1

        if ( j==1 ) then
            Y(1) = 0.d0

        elseif( j <= nyLrg1+1 ) then                                                                                !------------------------Sin dyLrg1-----------------------!

            nratio = (j-1.d0)/nyLrg1*1.d0 
            Y(j) = (lyLrg1-dySml*nyLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dySml*(j-1.d0)

        elseif( j > nyLrg1+1 .AND. j <= nyLrg1+nySml+1 ) then                                                       !------------------------Uniform dySml-----------------------!

            Y(j) = Y(j-1) + dySml

        elseif( j > nyLrg1+nySml+1 .AND. j <= ny+1 ) then                                                           !------------------------Sin dyLrg2-----------------------!

            nratio = (j-nyLrg1-nySml-1.d0)/nyLrg2*1.d0 
            Y(j) = lyLrg1+lySml + (lyLrg2-dySml*nyLrg2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + dySml*(j-nyLrg1-nySml-1.d0)
        endif

    end do

    Y(0)=Y(1)-(Y(2)-Y(1))
    Y(-1)=Y(0)-(Y(1)-Y(0))
    Y(ny+2)=Y(ny+1)+(Y(ny+1)-Y(ny))
    Y(ny+3)=Y(ny+2)+(Y(ny+2)-Y(ny+1))

    dy=Y(2)-Y(1)

    do k=1,nz+1
        
        if ( k==1 ) then
            Z(k) = 0.d0

        elseif( k <= nzLrg1 ) then                                                                                  !------------------------Sin dzLrg1-----------------------!
        
            nratio = (k-1.d0)/nzLrg1*1.d0 
            Z(k) = (lzLrg1-dzSml*nzLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dzSml*(k-1.d0)

        elseif( k > nzLrg1 .AND. k <= nzLrg1+nzSml+1 ) then                                                         !------------------------Uniform dzSml-----------------------!
            Z(k) = Z(k-1) + dzSml

        elseif( k > nzLrg1+nzSml+1 .AND. k <= nz+1 ) then                                                           !------------------------Sin dzLrg2-----------------------!

            nratio = (k-nzLrg1-nzSml-1.d0)/nzLrg2*1.d0 
            Z(k) =  lzLrg1+lzSml + (lzLrg2-dzSml*nzLrg2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + dzSml*(k-nzLrg1-nzSml-1.d0)

        endif

    end do

    Z(0)=Z(1)-(Z(nz+1)-Z(nz))
    Z(-1)=Z(0)-(Z(nz)-Z(nz-1))
    Z(nz+2)=Z(nz+1)+(Z(2)-Z(1))
    Z(nz+3)=Z(nz+2)+(Z(3)-Z(2))

    dz=Z(2)-Z(1)

    !----------------Unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=-1,nx+1

        iDx(i) = ( X(i+1) - X(i) )
        Dxs(i) = ( X(i+2) - X(i) ) *0.5d0

    end do

    do j=-1,ny+1

        iDy(j) = ( Y(j+1) - Y(j) )
        Dys(j) = ( Y(j+2) - Y(j) ) *0.5d0

    end do

    do k=-1,nz+1

        iDz(k) = ( Z(k+1) - Z(k) )
        Dzs(k) = ( Z(k+2) - Z(k) ) *0.5d0

    end do



    !Ghost boundary grid lengths
    iDx(nx+2) = iDx(nx)
    Dxs(nx+2) = Dxs(nx-1)

    iDy(ny+2) = iDy(ny)
    Dys(ny+2) = Dys(ny-1)

    iDz(nz+2) = iDz(nz)
    Dzs(nz+2) = Dzs(nz-1)


    !Defining the midpoint values of the grids
	do i=1,nx
        Xs(i) = 0.5d0 * ( X(i+1) + X(i) )
		if ((X(i) .LE. x_start) .AND. (X(i+1) .GT. x_start) ) then 	! filer boundary
			i_start = i
		elseif (x_end .GE. lx) then
			i_end = nx
		elseif ((X(i) .LT. x_end) .AND. (X(i+1) .GE. x_end) ) then
			i_end = i
		endif
    end do

	do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
		if ((Y(j) .LE. y_start) .AND. (Y(j+1) .GT. y_start) ) then 	! filer boundary
			j_start = j
		elseif (y_end .GE. ly) then
			j_end = ny
		elseif ((Y(j) .LT. y_end) .AND. (Y(j+1) .GE. y_end) ) then
			j_end = j
		endif
    end do
	
	do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
		if ((Z(k) .LE. z_start) .AND. (Z(k+1) .GT. z_start) ) then 	! filer boundary
			k_start = k
		elseif (z_end .GE. lz) then
			k_end = nz
		elseif ((Z(k) .LT. z_end) .AND. (Z(k+1) .GE. z_end) ) then
			k_end = k
		endif
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
        Yout(i,j,k) = Ys(j)
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo


    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output3D/'))
        open (unit=1,form='unformatted',file=TRIM('output3D/')//'3d_mesh.x')
        write(1) nblocks
        write(1) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1

        write(1)    (((Xout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Yout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Zout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end)

        close(1)
    end if
	
    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output2D/'))
        open (unit=2,form='unformatted',file=TRIM('output2D/')//'2d_mesh.x')
        write(2) nblocks
        write(2) 1, ny, nz

        write(2)    (((Xout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,1),j=1,ny),k=1,nz)

        close(2)
    end if

end subroutine gridder_sin4

!-------------------------------------------------------------------------------------------!

subroutine gridder_sin4_3D() ! Grid spacing method proposed by Kuyper-3D
use variables
implicit none
real*8 :: nratio

real*8, parameter :: lxLrg1 = xc - lxSml*0.5d0
real*8, parameter :: lxLrg2 = lx - lxLrg1 - lxSml

integer, parameter :: nxLrg1 = INT( (nx-nxSml)/(lx-lxSml)*lxLrg1 )
integer, parameter :: nxLrg2 = nx-nxLrg1-nxSml

real*8, parameter :: dxLrg1 = lxLrg1/nxLrg1*1.d0
real*8, parameter :: dxLrg2 = lxLrg2/nxLrg2*1.d0

real*8, parameter :: lyLrg1 = GridderYc - lySml*0.5d0
real*8, parameter :: lyLrg2 = ly - lyLrg1 - lySml

integer, parameter :: nyLrg1 = INT( (ny-nySml)/(ly-lySml)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nySml

real*8, parameter :: dyLrg1 = lyLrg1/nyLrg1*1.d0
real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: lzLrg1 = GridderZc - lzSml*0.5d0
real*8, parameter :: lzLrg2 = lz - lzLrg1 - lzSml

integer, parameter :: nzLrg1 = nzUps !INT( (nz-nzSml)/(lz-lzSml)*lzLrg1 )  !1 
integer, parameter :: nzLrg2 = nz-nzLrg1-nzSml                     !89

real*8, parameter :: dzLrg1 = lzLrg1/nzLrg1*1.d0
real*8, parameter :: dzLrg2 = lzLrg2/nzLrg2*1.d0

real*8, parameter :: thetaa = 4.d0*ATAN(1.d0)
real*8, parameter :: tune = 1.d0



!----------------Unequal grid intervals----------------!
    do i=1,nx+1

        if ( i==1 ) then
            X(1) = 0.d0

        elseif( i <= nxLrg1+1 ) then                                                                                !------------------------Sin dxLrg1-----------------------!

            nratio = (i-1.d0)/nxLrg1*1.d0 
            X(i) = (lxLrg1-dxSml*nxLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dxSml*(i-1.d0)

        elseif( i > nxLrg1+1 .AND. i <= nxLrg1+nxSml+1 ) then                                                       !------------------------Uniform dxSml-----------------------!

            X(i) = X(i-1) + dxSml

        elseif( i > nxLrg1+nxSml+1 .AND. i <= nx+1 ) then                                                           !------------------------Sin dxLrg2-----------------------!

            nratio = (i-nxLrg1-nxSml-1.d0)/nxLrg2*1.d0 
            X(i) = lxLrg1+lxSml + (lxLrg2-dxSml*nxLrg2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + dxSml*(i-nxLrg1-nxSml-1.d0)
        endif

    end do

    X(0)=X(1)-(X(2)-X(1))
    X(-1)=X(0)-(X(1)-X(0))
    X(nx+2)=X(nx+1)+(X(nx+1)-X(nx))
    X(nx+3)=X(nx+2)+(X(nx+2)-X(nx+1))

    dx=X(2)-X(1)


    do j=1,ny+1

        if ( j==1 ) then
            Y(1) = 0.d0

        elseif( j <= nyLrg1+1 ) then                                                                                !------------------------Sin dyLrg1-----------------------!

            nratio = (j-1.d0)/nyLrg1*1.d0 
            Y(j) = (lyLrg1-dySml*nyLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dySml*(j-1.d0)

        elseif( j > nyLrg1+1 .AND. j <= nyLrg1+nySml+1 ) then                                                       !------------------------Uniform dySml-----------------------!

            Y(j) = Y(j-1) + dySml

        elseif( j > nyLrg1+nySml+1 .AND. j <= ny+1 ) then                                                           !------------------------Sin dyLrg2-----------------------!

            nratio = (j-nyLrg1-nySml-1.d0)/nyLrg2*1.d0 
            Y(j) = lyLrg1+lySml + (lyLrg2-dySml*nyLrg2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + dySml*(j-nyLrg1-nySml-1.d0)
        endif

    end do

    Y(0)=Y(1)-(Y(2)-Y(1))
    Y(-1)=Y(0)-(Y(1)-Y(0))
    Y(ny+2)=Y(ny+1)+(Y(ny+1)-Y(ny))
    Y(ny+3)=Y(ny+2)+(Y(ny+2)-Y(ny+1))

    dy=Y(2)-Y(1)

    do k=1,nz+1
        
        if ( k==1 ) then
            Z(k) = 0.d0

        elseif( k <= nzLrg1 ) then                                                                                  !------------------------Sin dzLrg1-----------------------!
        
            nratio = (k-1.d0)/nzLrg1*1.d0 
            Z(k) = (lzLrg1-dzSml*nzLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dzSml*(k-1.d0)

        elseif( k > nzLrg1 .AND. k <= nzLrg1+nzSml+1 ) then                                                         !------------------------Uniform dzSml-----------------------!
            Z(k) = Z(k-1) + dzSml

        elseif( k > nzLrg1+nzSml+1 .AND. k <= nz+1 ) then                                                           !------------------------Sin dzLrg2-----------------------!

            nratio = (k-nzLrg1-nzSml-1.d0)/nzLrg2*1.d0 
            Z(k) =  lzLrg1+lzSml + (lzLrg2-dzSml*nzLrg2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + dzSml*(k-nzLrg1-nzSml-1.d0)

        endif

    end do

    Z(0)=Z(1)-(Z(nz+1)-Z(nz))
    Z(-1)=Z(0)-(Z(nz)-Z(nz-1))
    Z(nz+2)=Z(nz+1)+(Z(2)-Z(1))
    Z(nz+3)=Z(nz+2)+(Z(3)-Z(2))

    dz=Z(2)-Z(1)

    !----------------Unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=-1,nx+1

        iDx(i) = ( X(i+1) - X(i) )
        Dxs(i) = ( X(i+2) - X(i) ) *0.5d0

    end do

    do j=-1,ny+1

        iDy(j) = ( Y(j+1) - Y(j) )
        Dys(j) = ( Y(j+2) - Y(j) ) *0.5d0

    end do

    do k=-1,nz+1

        iDz(k) = ( Z(k+1) - Z(k) )
        Dzs(k) = ( Z(k+2) - Z(k) ) *0.5d0

    end do



    !Ghost boundary grid lengths
    iDx(nx+2) = iDx(nx)
    Dxs(nx+2) = Dxs(nx-1)

    iDy(ny+2) = iDy(ny)
    Dys(ny+2) = Dys(ny-1)

    iDz(nz+2) = iDz(nz)
    Dzs(nz+2) = Dzs(nz-1)


    !Defining the midpoint values of the grids
	do i=1,nx
        Xs(i) = 0.5d0 * ( X(i+1) + X(i) )
		if ((X(i) .LE. x_start) .AND. (X(i+1) .GT. x_start) ) then 	! filer boundary
			i_start = i
		elseif (x_end .GE. lx) then
			i_end = nx
		elseif ((X(i) .LT. x_end) .AND. (X(i+1) .GE. x_end) ) then
			i_end = i
		endif
    end do

	do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
		if ((Y(j) .LE. y_start) .AND. (Y(j+1) .GT. y_start) ) then 	! filer boundary
			j_start = j
		elseif (y_end .GE. ly) then
			j_end = ny
		elseif ((Y(j) .LT. y_end) .AND. (Y(j+1) .GE. y_end) ) then
			j_end = j
		endif
    end do
	
	do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
		if ((Z(k) .LE. z_start) .AND. (Z(k+1) .GT. z_start) ) then 	! filer boundary
			k_start = k
		elseif (z_end .GE. lz) then
			k_end = nz
		elseif ((Z(k) .LT. z_end) .AND. (Z(k+1) .GE. z_end) ) then
			k_end = k
		endif
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
        Yout(i,j,k) = Ys(j)
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo


    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output3D/'))
        open (unit=1,form='unformatted',file=TRIM('output3D/')//'3d_mesh.x')
        write(1) nblocks
        write(1) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1

        write(1)    (((Xout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Yout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Zout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end)

        close(1)
    end if
	
    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output2D/'))
        open (unit=2,form='unformatted',file=TRIM('output2D/')//'2d_mesh.x')
        write(2) nblocks
        write(2) 1, ny, nz

        write(2)    (((Xout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,1),j=1,ny),k=1,nz)

        close(2)
    end if


end subroutine gridder_sin4_3D


subroutine gridder_sin5() ! Double uniform & sine grid spacing method : Uniform - SIN - Uniform in y and z
use variables
implicit none
real*8 :: nratio

real*8, parameter :: lyLrg1 = GridderYc - lyMid*0.5d0
real*8, parameter :: lyLrg2 = ly - lyLrg1 - lyMid

integer, parameter :: nyLrg1 = INT( (ny-nyMid-nySml)/(ly-lyMid)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nyMid-nySml

real*8, parameter :: dyLrg1 = lyLrg1/nyLrg1*1.d0
real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: lyMid1 = GridderYc - lyLrg1 - lySml*0.5d0
real*8, parameter :: lyMid2 = lyMid - lyMid1 - lySml

integer, parameter :: nyMid1 = INT( nyMid/(lyMid-lySml)*lyMid1 )
integer, parameter :: nyMid2 = nyMid-nyMid1

real*8, parameter :: dyMid1 = lyMid1/nyMid1*1.d0
real*8, parameter :: dyMid2 = lyMid2/nyMid2*1.d0

real*8, parameter :: lzLrg1 = GridderZc - lzSml*0.5d0
real*8, parameter :: lzLrg2 = lz - lzLrg1 - lzMid

integer, parameter :: nzLrg1 = nzUps !INT( (nz-nzMid-nzSml)/(lz-lzMid)*lzLrg1 )  !1 
integer, parameter :: nzLrg2 = nz-nzLrg1-nzMid-nzSml                     !89

real*8, parameter :: dzLrg1 = lzLrg1/nzLrg1*1.d0
real*8, parameter :: dzLrg2 = lzLrg2/nzLrg2*1.d0

real*8, parameter :: lzMid2 = lzMid - lzSml

integer, parameter :: nzMid2 = nzMid

real*8, parameter :: dzMid2 = lzMid2/nzMid2*1.d0

real*8, parameter :: thetaa = 4.d0*ATAN(1.d0)
real*8, parameter :: tune = 1.d0



!----------------Unequal grid intervals----------------!
    do i=1,nx+3
        if(i == 1) then
            X(i) = 0.0
            X(i-1) = X(i) - dx
            X(i-2) = X(i-1) - dx
        else
            X(i) = X(i-1) + dx
        end if 
    end do


    do j=1,ny+1

        if ( j==1 ) then
            Y(1) = 0.d0

        elseif( j <= nyLrg1+1 ) then                                                                                 !------------------------Uniform dyLrg1-----------------------!

            Y(j) = Y(j-1) + dyLrg1

        elseif( j > nyLrg1+1 .AND. j <= nyLrg1+nyMid1+1 ) then                                                       !------------------------Sin dyMid1-------------------!

            nratio = (j-nyLrg1-1.d0)/nyMid1*1.d0 
            Y(j) = lyLrg1+ (lyMid1-dySml*nyMid1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dySml*(j-nyLrg1-1.d0)

        elseif( j > nyLrg1+nyMid1+1 .AND. j <= nyLrg1+nyMid1+nySml+1 ) then                                          !------------------------Uniform dySml--------------------!

            Y(j) = Y(j-1) + dySml

        elseif( j > nyLrg1+nyMid1+nySml+1 .AND. j <= nyLrg1+nyMid1+nySml+nyMid2+1 ) then                             !------------------------Sin dyMid2-------------------!

            nratio = (j-nyLrg1-nyMid1-nySml-1.d0)/nyMid2*1.d0 
            Y(j) = lyLrg1+lyMid1+lySml + (lyMid2-dySml*nyMid2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + &
                   dySml*(j-nyLrg1-nyMid1-nySml-1.d0)

        elseif( j > nyLrg1+nyMid1+nySml+nyMid2+1 .AND. j <= ny+1 ) then                                              !------------------------Uniform dyLrg2-----------------------!

            Y(j) = Y(j-1) + dyLrg2
        endif

    end do

    Y(0)=Y(1)-(Y(2)-Y(1))
    Y(-1)=Y(0)-(Y(1)-Y(0))
    Y(ny+2)=Y(ny+1)+(Y(ny+1)-Y(ny))
    Y(ny+3)=Y(ny+2)+(Y(ny+2)-Y(ny+1))

    dy=Y(2)-Y(1)

    do k=1,nz+1
        
        if ( k==1 ) then
            Z(k) = 0.d0

        elseif( k <= nzLrg1+1 ) then                                                                                !------------------------Sin dzLrg1-----------------------!
        
            nratio = (k-1.d0)/nzLrg1*1.d0 
            Z(k) = (lzLrg1-dzMid*nzLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dzMid*(k-1.d0)

        elseif( k > nzLrg1+1 .AND. k <= nzLrg1+nzSml+1 ) then                                                       !------------------------Uniform dzSml--------------------!
            Z(k) = Z(k-1) + dzSml

        elseif( k > nzLrg1+nzSml+1 .AND. k <= nzLrg1+nzSml+nzMid2+1 ) then                                          !------------------------Sin dzMid2-------------------!
            nratio = (k-nzLrg1-nzSml-1.d0)/nzMid2*1.d0 
            Z(k) =  lzLrg1+lzSml + (lzMid2-dzSml*nzMid2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + &
                    dzSml*(k-nzLrg1-nzSml-1.d0)

        elseif( k > nzLrg1+nzSml+nzMid2+1  .AND. k <= nz+1 ) then                                                   !------------------------Uniform dzLrg2-----------------------!

            Z(k) = Z(k-1) + dzLrg2

        endif

    end do

    Z(0)=Z(1)-(Z(nz+1)-Z(nz))
    Z(-1)=Z(0)-(Z(nz)-Z(nz-1))
    Z(nz+2)=Z(nz+1)+(Z(2)-Z(1))
    Z(nz+3)=Z(nz+2)+(Z(3)-Z(2))

    dz=Z(2)-Z(1)

    !----------------Unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=-1,nx+1

        iDx(i) = ( X(i+1) - X(i) )
        Dxs(i) = ( X(i+2) - X(i) ) *0.5d0

    end do

    do j=-1,ny+1

        iDy(j) = ( Y(j+1) - Y(j) )
        Dys(j) = ( Y(j+2) - Y(j) ) *0.5d0

    end do

    do k=-1,nz+1

        iDz(k) = ( Z(k+1) - Z(k) )
        Dzs(k) = ( Z(k+2) - Z(k) ) *0.5d0

    end do



    !Ghost boundary grid lengths
    iDx(nx+2) = iDx(nx)
    Dxs(nx+2) = Dxs(nx-1)

    iDy(ny+2) = iDy(ny)
    Dys(ny+2) = Dys(ny-1)

    iDz(nz+2) = iDz(nz)
    Dzs(nz+2) = Dzs(nz-1)


    !Defining the midpoint values of the grids
	do i=1,nx
        Xs(i) = 0.5d0 * ( X(i+1) + X(i) )
		if ((X(i) .LE. x_start) .AND. (X(i+1) .GT. x_start) ) then 	! filer boundary
			i_start = i
		elseif (x_end .GE. lx) then
			i_end = nx
		elseif ((X(i) .LT. x_end) .AND. (X(i+1) .GE. x_end) ) then
			i_end = i
		endif
    end do

	do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
		if ((Y(j) .LE. y_start) .AND. (Y(j+1) .GT. y_start) ) then 	! filer boundary
			j_start = j
		elseif (y_end .GE. ly) then
			j_end = ny
		elseif ((Y(j) .LT. y_end) .AND. (Y(j+1) .GE. y_end) ) then
			j_end = j
		endif
    end do
	
	do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
		if ((Z(k) .LE. z_start) .AND. (Z(k+1) .GT. z_start) ) then 	! filer boundary
			k_start = k
		elseif (z_end .GE. lz) then
			k_end = nz
		elseif ((Z(k) .LT. z_end) .AND. (Z(k+1) .GE. z_end) ) then
			k_end = k
		endif
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
        Yout(i,j,k) = Ys(j)
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo


    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output3D/'))
        open (unit=1,form='unformatted',file=TRIM('output3D/')//'3d_mesh.x')
        write(1) nblocks
        write(1) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1

        write(1)    (((Xout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Yout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Zout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end)

        close(1)
    end if
	
    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output2D/'))
        open (unit=2,form='unformatted',file=TRIM('output2D/')//'2d_mesh.x')
        write(2) nblocks
        write(2) 1, ny, nz

        write(2)    (((Xout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,1),j=1,ny),k=1,nz)

        close(2)
    end if

end subroutine gridder_sin5

subroutine gridder_sin5_3D() ! Double uniform & sine grid spacing method : Uniform - SIN - Uniform in x,y, and z
use variables
implicit none
real*8 :: nratio

real*8, parameter :: lxLrg1 = xc - lxMid*0.5d0
real*8, parameter :: lxLrg2 = lx - lxLrg1 - lxMid

integer, parameter :: nxLrg1 = INT( (nx-nxMid-nxSml)/(lx-lxMid)*lxLrg1 )
integer, parameter :: nxLrg2 = nx-nxLrg1-nxMid-nxSml

real*8, parameter :: dxLrg1 = lxLrg1/nxLrg1*1.d0
real*8, parameter :: dxLrg2 = lxLrg2/nxLrg2*1.d0

real*8, parameter :: lxMid1 = xc - lxLrg1 - lxSml*0.5d0
real*8, parameter :: lxMid2 = lxMid - lxMid1 - lxSml

integer, parameter :: nxMid1 = INT( nxMid/(lxMid-lxSml)*lxMid1 )
integer, parameter :: nxMid2 = nxMid-nxMid1

real*8, parameter :: dxMid1 = lxMid1/nxMid1*1.d0
real*8, parameter :: dxMid2 = lxMid2/nxMid2*1.d0

real*8, parameter :: lyLrg1 = GridderYc - lyMid*0.5d0
real*8, parameter :: lyLrg2 = ly - lyLrg1 - lyMid

integer, parameter :: nyLrg1 = INT( (ny-nyMid-nySml)/(ly-lyMid)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nyMid-nySml

real*8, parameter :: dyLrg1 = lyLrg1/nyLrg1*1.d0
real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: lyMid1 = GridderYc - lyLrg1 - lySml*0.5d0
real*8, parameter :: lyMid2 = lyMid - lyMid1 - lySml

integer, parameter :: nyMid1 = INT( nyMid/(lyMid-lySml)*lyMid1 )
integer, parameter :: nyMid2 = nyMid-nyMid1

real*8, parameter :: dyMid1 = lyMid1/nyMid1*1.d0
real*8, parameter :: dyMid2 = lyMid2/nyMid2*1.d0

real*8, parameter :: lzLrg1 = GridderZc - lzSml*0.5d0
real*8, parameter :: lzLrg2 = lz - lzLrg1 - lzMid

integer, parameter :: nzLrg1 = nzUps !INT( (nz-nzMid-nzSml)/(lz-lzMid)*lzLrg1 )  !1 
integer, parameter :: nzLrg2 = nz-nzLrg1-nzMid-nzSml                     !89

real*8, parameter :: dzLrg1 = lzLrg1/nzLrg1*1.d0
real*8, parameter :: dzLrg2 = lzLrg2/nzLrg2*1.d0

real*8, parameter :: lzMid2 = lzMid - lzSml

integer, parameter :: nzMid2 = nzMid

real*8, parameter :: dzMid2 = lzMid2/nzMid2*1.d0

real*8, parameter :: thetaa = 4.d0*ATAN(1.d0)
real*8, parameter :: tune = 1.d0



!----------------Unequal grid intervals----------------!
    do i=1,nx+1

        if ( i==1 ) then
            X(1) = 0.d0

        elseif( i <= nxLrg1+1 ) then                                                                                 !------------------------Uniform dxLrg1-----------------------!

            X(i) = X(i-1) + dxLrg1

        elseif( i > nxLrg1+1 .AND. i <= nxLrg1+nxMid1+1 ) then                                                       !------------------------Sin dxMid1-------------------!

            nratio = (i-nxLrg1-1.d0)/nxMid1*1.d0 
            X(i) = lxLrg1+ (lxMid1-dxSml*nxMid1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dxSml*(i-nxLrg1-1.d0)

        elseif( i > nxLrg1+nxMid1+1 .AND. i <= nxLrg1+nxMid1+nxSml+1 ) then                                          !------------------------Uniform dxSml--------------------!

            X(i) = X(i-1) + dxSml

        elseif( i > nxLrg1+nxMid1+nxSml+1 .AND. i <= nxLrg1+nxMid1+nxSml+nxMid2+1 ) then                             !------------------------Sin dxMid2-------------------!

            nratio = (i-nxLrg1-nxMid1-nxSml-1.d0)/nxMid2*1.d0 
            X(i) = lxLrg1+lxMid1+lxSml + (lxMid2-dxSml*nxMid2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + &
                   dxSml*(i-nxLrg1-nxMid1-nxSml-1.d0)

        elseif( i > nxLrg1+nxMid1+nxSml+nxMid2+1 .AND. i <= nx+1 ) then                                              !------------------------Uniform dxLrg2-----------------------!

            X(i) = X(i-1) + dxLrg2
        endif

    end do

    X(0)=X(1)-(X(2)-X(1))
    X(-1)=X(0)-(X(1)-X(0))
    X(nx+2)=X(nx+1)+(X(nx+1)-X(nx)) 
    X(nx+3)=X(nx+2)+(X(nx+2)-X(nx+1))

    dx=X(2)-X(1)


    do j=1,ny+1

        if ( j==1 ) then
            Y(1) = 0.d0

        elseif( j <= nyLrg1+1 ) then                                                                                 !------------------------Uniform dyLrg1-----------------------!

            Y(j) = Y(j-1) + dyLrg1

        elseif( j > nyLrg1+1 .AND. j <= nyLrg1+nyMid1+1 ) then                                                       !------------------------Sin dyMid1-------------------!

            nratio = (j-nyLrg1-1.d0)/nyMid1*1.d0 
            Y(j) = lyLrg1+ (lyMid1-dySml*nyMid1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dySml*(j-nyLrg1-1.d0)

        elseif( j > nyLrg1+nyMid1+1 .AND. j <= nyLrg1+nyMid1+nySml+1 ) then                                          !------------------------Uniform dySml--------------------!

            Y(j) = Y(j-1) + dySml

        elseif( j > nyLrg1+nyMid1+nySml+1 .AND. j <= nyLrg1+nyMid1+nySml+nyMid2+1 ) then                             !------------------------Sin dyMid2-------------------!

            nratio = (j-nyLrg1-nyMid1-nySml-1.d0)/nyMid2*1.d0 
            Y(j) = lyLrg1+lyMid1+lySml + (lyMid2-dySml*nyMid2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + &
                   dySml*(j-nyLrg1-nyMid1-nySml-1.d0)

        elseif( j > nyLrg1+nyMid1+nySml+nyMid2+1 .AND. j <= ny+1 ) then                                              !------------------------Uniform dyLrg2-----------------------!

            Y(j) = Y(j-1) + dyLrg2
        endif

    end do

    Y(0)=Y(1)-(Y(2)-Y(1))
    Y(-1)=Y(0)-(Y(1)-Y(0))
    Y(ny+2)=Y(ny+1)+(Y(ny+1)-Y(ny))
    Y(ny+3)=Y(ny+2)+(Y(ny+2)-Y(ny+1))

    dy=Y(2)-Y(1)

    do k=1,nz+1
        
        if ( k==1 ) then
            Z(k) = 0.d0

        elseif( k <= nzLrg1+1 ) then                                                                                !------------------------Sin dzLrg1-----------------------!
        
            nratio = (k-1.d0)/nzLrg1*1.d0 
            Z(k) = (lzLrg1-dzMid*nzLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dzMid*(k-1.d0)

        elseif( k > nzLrg1+1 .AND. k <= nzLrg1+nzSml+1 ) then                                                       !------------------------Uniform dzSml--------------------!
            Z(k) = Z(k-1) + dzSml

        elseif( k > nzLrg1+nzSml+1 .AND. k <= nzLrg1+nzSml+nzMid2+1 ) then                                          !------------------------Sin dzMid2-------------------!
            nratio = (k-nzLrg1-nzSml-1.d0)/nzMid2*1.d0 
            Z(k) =  lzLrg1+lzSml + (lzMid2-dzSml*nzMid2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + &
                    dzSml*(k-nzLrg1-nzSml-1.d0)

        elseif( k > nzLrg1+nzSml+nzMid2+1  .AND. k <= nz+1 ) then                                                   !------------------------Uniform dzLrg2-----------------------!

            Z(k) = Z(k-1) + dzLrg2

        endif

    end do

    Z(0)=Z(1)-(Z(nz+1)-Z(nz))
    Z(-1)=Z(0)-(Z(nz)-Z(nz-1))
    Z(nz+2)=Z(nz+1)+(Z(2)-Z(1))
    Z(nz+3)=Z(nz+2)+(Z(3)-Z(2))

    dz=Z(2)-Z(1)

    !----------------Unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=-1,nx+1

        iDx(i) = ( X(i+1) - X(i) )
        Dxs(i) = ( X(i+2) - X(i) ) *0.5d0

    end do

    do j=-1,ny+1

        iDy(j) = ( Y(j+1) - Y(j) )
        Dys(j) = ( Y(j+2) - Y(j) ) *0.5d0

    end do

    do k=-1,nz+1

        iDz(k) = ( Z(k+1) - Z(k) )
        Dzs(k) = ( Z(k+2) - Z(k) ) *0.5d0

    end do



    !Ghost boundary grid lengths
    iDx(nx+2) = iDx(nx)
    Dxs(nx+2) = Dxs(nx-1)

    iDy(ny+2) = iDy(ny)
    Dys(ny+2) = Dys(ny-1)

    iDz(nz+2) = iDz(nz)
    Dzs(nz+2) = Dzs(nz-1)


    !Defining the midpoint values of the grids
	do i=1,nx
        Xs(i) = 0.5d0 * ( X(i+1) + X(i) )
		if ((X(i) .LE. x_start) .AND. (X(i+1) .GT. x_start) ) then 	! filer boundary
			i_start = i
		elseif (x_end .GE. lx) then
			i_end = nx
		elseif ((X(i) .LT. x_end) .AND. (X(i+1) .GE. x_end) ) then
			i_end = i
		endif
    end do

	do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
		if ((Y(j) .LE. y_start) .AND. (Y(j+1) .GT. y_start) ) then 	! filer boundary
			j_start = j
		elseif (y_end .GE. ly) then
			j_end = ny
		elseif ((Y(j) .LT. y_end) .AND. (Y(j+1) .GE. y_end) ) then
			j_end = j
		endif
    end do
	
	do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
		if ((Z(k) .LE. z_start) .AND. (Z(k+1) .GT. z_start) ) then 	! filer boundary
			k_start = k
		elseif (z_end .GE. lz) then
			k_end = nz
		elseif ((Z(k) .LT. z_end) .AND. (Z(k+1) .GE. z_end) ) then
			k_end = k
		endif
    end do
	

    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
        Yout(i,j,k) = Ys(j)
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo


    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output3D/'))
        open (unit=1,form='unformatted',file=TRIM('output3D/')//'3d_mesh.x')
        write(1) nblocks
        write(1) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1

        write(1)    (((Xout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Yout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Zout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end)

        close(1)
    end if
	
    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output2D/'))
        open (unit=2,form='unformatted',file=TRIM('output2D/')//'2d_mesh.x')
        write(2) nblocks
        write(2) 1, ny, nz

        write(2)    (((Xout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,1),j=1,ny),k=1,nz)

        close(2)
    end if


end subroutine gridder_sin5_3D

subroutine gridder_ground() ! Double uniform & sine grid spacing method : Uniform - SIN - Uniform in x,y, and z for object lies on the bottom wall
use variables
implicit none
real*8 :: nratio

real*8, parameter :: lxLrg1 = xc - lxMid*0.5d0
real*8, parameter :: lxLrg2 = lx - lxLrg1 - lxMid

integer, parameter :: nxLrg1 = INT( (nx-nxMid-nxSml)/(lx-lxMid)*lxLrg1 )
integer, parameter :: nxLrg2 = nx-nxLrg1-nxMid-nxSml

real*8, parameter :: dxLrg1 = lxLrg1/nxLrg1*1.d0
real*8, parameter :: dxLrg2 = lxLrg2/nxLrg2*1.d0

real*8, parameter :: lxMid1 = xc - lxLrg1 - lxSml*0.5d0
real*8, parameter :: lxMid2 = lxMid - lxMid1 - lxSml

integer, parameter :: nxMid1 = INT( nxMid/(lxMid-lxSml)*lxMid1 )
integer, parameter :: nxMid2 = nxMid-nxMid1

real*8, parameter :: dxMid1 = lxMid1/nxMid1*1.d0
real*8, parameter :: dxMid2 = lxMid2/nxMid2*1.d0

real*8, parameter :: lyLrg2 = ly - lyMid

integer, parameter :: nyLrg2 = ny-nyMid-nySml

real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: lyMid2 = lyMid - lySml

integer, parameter :: nyMid2 = nyMid

real*8, parameter :: dyMid2 = lyMid2/nyMid2*1.d0

real*8, parameter :: lzLrg1 = GridderZc - lzSml*0.5d0
real*8, parameter :: lzLrg2 = lz - lzLrg1 - lzMid

integer, parameter :: nzLrg1 = nzUps !INT( (nz-nzMid-nzSml)/(lz-lzMid)*lzLrg1 )  !1 
integer, parameter :: nzLrg2 = nz-nzLrg1-nzMid-nzSml                     !89

real*8, parameter :: dzLrg1 = lzLrg1/nzLrg1*1.d0
real*8, parameter :: dzLrg2 = lzLrg2/nzLrg2*1.d0

real*8, parameter :: lzMid2 = lzMid - lzSml

integer, parameter :: nzMid2 = nzMid

real*8, parameter :: dzMid2 = lzMid2/nzMid2*1.d0

real*8, parameter :: thetaa = 4.d0*ATAN(1.d0)
real*8, parameter :: tune = 1.d0



!----------------Unequal grid intervals----------------!
    do i=1,nx+1

        if ( i==1 ) then
            X(1) = 0.d0

        elseif( i <= nxLrg1+1 ) then                                                                                 !------------------------Uniform dxLrg1-----------------------!

            X(i) = X(i-1) + dxLrg1

        elseif( i > nxLrg1+1 .AND. i <= nxLrg1+nxMid1+1 ) then                                                       !------------------------Sin dxMid1-------------------!

            nratio = (i-nxLrg1-1.d0)/nxMid1*1.d0 
            X(i) = lxLrg1+ (lxMid1-dxSml*nxMid1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dxSml*(i-nxLrg1-1.d0)

        elseif( i > nxLrg1+nxMid1+1 .AND. i <= nxLrg1+nxMid1+nxSml+1 ) then                                          !------------------------Uniform dxSml--------------------!

            X(i) = X(i-1) + dxSml

        elseif( i > nxLrg1+nxMid1+nxSml+1 .AND. i <= nxLrg1+nxMid1+nxSml+nxMid2+1 ) then                             !------------------------Sin dxMid2-------------------!

            nratio = (i-nxLrg1-nxMid1-nxSml-1.d0)/nxMid2*1.d0 
            X(i) = lxLrg1+lxMid1+lxSml + (lxMid2-dxSml*nxMid2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + &
                   dxSml*(i-nxLrg1-nxMid1-nxSml-1.d0)

        elseif( i > nxLrg1+nxMid1+nxSml+nxMid2+1 .AND. i <= nx+1 ) then                                              !------------------------Uniform dxLrg2-----------------------!

            X(i) = X(i-1) + dxLrg2
        endif

    end do

    X(0)=X(1)-(X(2)-X(1))
    X(-1)=X(0)-(X(1)-X(0))
    X(nx+2)=X(nx+1)+(X(nx+1)-X(nx)) 
    X(nx+3)=X(nx+2)+(X(nx+2)-X(nx+1))

    dx=X(2)-X(1)


    do j=1,ny+1

        if ( j==1 ) then
            Y(1) = 0.d0

        elseif( j <= nySml+1 ) then                                                                                 !------------------------Uniform dySml1-----------------------!

            Y(j) = Y(j-1) + dySml

        elseif( j > nySml+1 .AND. j <= nySml+nyMid2+1 ) then                                                        !------------------------Sin dyMid2-------------------!

            nratio = (j-nySml-1.d0)/nyMid2*1.d0 
            Y(j) = lySml + (lyMid2-dySml*nyMid2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + &
                   dySml*(j-nySml-1.d0)

        elseif( j > nySml+nyMid2+1 .AND. j <= ny+1 ) then                                                           !------------------------Uniform dyLrg2-----------------------!

            Y(j) = Y(j-1) + dyLrg2
        endif

    end do

    Y(0)=Y(1)-(Y(2)-Y(1))
    Y(-1)=Y(0)-(Y(1)-Y(0))
    Y(ny+2)=Y(ny+1)+(Y(ny+1)-Y(ny))
    Y(ny+3)=Y(ny+2)+(Y(ny+2)-Y(ny+1))

    dy=Y(2)-Y(1)

    do k=1,nz+1
        
        if ( k==1 ) then
            Z(k) = 0.d0

        elseif( k <= nzLrg1+1 ) then                                                                                !------------------------Sin dzLrg1-----------------------!
        
            nratio = (k-1.d0)/nzLrg1*1.d0 
            Z(k) = (lzLrg1-dzMid*nzLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dzMid*(k-1.d0)

        elseif( k > nzLrg1+1 .AND. k <= nzLrg1+nzSml+1 ) then                                                       !------------------------Uniform dzSml--------------------!
            Z(k) = Z(k-1) + dzSml

        elseif( k > nzLrg1+nzSml+1 .AND. k <= nzLrg1+nzSml+nzMid2+1 ) then                                          !------------------------Sin dzMid2-------------------!
            nratio = (k-nzLrg1-nzSml-1.d0)/nzMid2*1.d0 
            Z(k) =  lzLrg1+lzSml + (lzMid2-dzSml*nzMid2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + &
                    dzSml*(k-nzLrg1-nzSml-1.d0)

        elseif( k > nzLrg1+nzSml+nzMid2+1  .AND. k <= nz+1 ) then                                                   !------------------------Uniform dzLrg2-----------------------!

            Z(k) = Z(k-1) + dzLrg2

        endif

    end do

    Z(0)=Z(1)-(Z(nz+1)-Z(nz))
    Z(-1)=Z(0)-(Z(nz)-Z(nz-1))
    Z(nz+2)=Z(nz+1)+(Z(2)-Z(1))
    Z(nz+3)=Z(nz+2)+(Z(3)-Z(2))

    dz=Z(2)-Z(1)

    !----------------Unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=-1,nx+1

        iDx(i) = ( X(i+1) - X(i) )
        Dxs(i) = ( X(i+2) - X(i) ) *0.5d0

    end do

    do j=-1,ny+1

        iDy(j) = ( Y(j+1) - Y(j) )
        Dys(j) = ( Y(j+2) - Y(j) ) *0.5d0

    end do

    do k=-1,nz+1

        iDz(k) = ( Z(k+1) - Z(k) )
        Dzs(k) = ( Z(k+2) - Z(k) ) *0.5d0

    end do



    !Ghost boundary grid lengths
    iDx(nx+2) = iDx(nx)
    Dxs(nx+2) = Dxs(nx-1)

    iDy(ny+2) = iDy(ny)
    Dys(ny+2) = Dys(ny-1)

    iDz(nz+2) = iDz(nz)
    Dzs(nz+2) = Dzs(nz-1)


    !Defining the midpoint values of the grids
	do i=1,nx
        Xs(i) = 0.5d0 * ( X(i+1) + X(i) )
		if ((X(i) .LE. x_start) .AND. (X(i+1) .GT. x_start) ) then 	! filer boundary
			i_start = i
		elseif (x_end .GE. lx) then
			i_end = nx
		elseif ((X(i) .LT. x_end) .AND. (X(i+1) .GE. x_end) ) then
			i_end = i
		endif
    end do

	do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
		if ((Y(j) .LE. y_start) .AND. (Y(j+1) .GT. y_start) ) then 	! filer boundary
			j_start = j
		elseif (y_end .GE. ly) then
			j_end = ny
		elseif ((Y(j) .LT. y_end) .AND. (Y(j+1) .GE. y_end) ) then
			j_end = j
		endif
    end do
	
	do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
		if ((Z(k) .LE. z_start) .AND. (Z(k+1) .GT. z_start) ) then 	! filer boundary
			k_start = k
		elseif (z_end .GE. lz) then
			k_end = nz
		elseif ((Z(k) .LT. z_end) .AND. (Z(k+1) .GE. z_end) ) then
			k_end = k
		endif
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
        Yout(i,j,k) = Ys(j)
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo


    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output3D/'))
        open (unit=1,form='unformatted',file=TRIM('output3D/')//'3d_mesh.x')
        write(1) nblocks
        write(1) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1

        write(1)    (((Xout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Yout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Zout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end)

        close(1)
    end if
	
    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output2D/'))
        open (unit=2,form='unformatted',file=TRIM('output2D/')//'2d_mesh.x')
        write(2) nblocks
        write(2) 1, ny, nz

        write(2)    (((Xout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,1),j=1,ny),k=1,nz)

        close(2)
    end if


end subroutine gridder_ground


! subroutine gridder_read_txt()
    ! use variables
    ! implicit none


! !---------------Load grid file in a single X Y Z format----------
	
	! ! OPEN(UNIT=3,FILE='grid.txt')
        ! ! do i=1,nx+1 ; do j=1,ny+1 ; do k=1,nz+1
           ! ! READ(3,*)X(i),Y(j),Z(k)
        ! ! enddo; enddo; enddo
	! ! CLOSE(3)

! !---------------Load grid file in separated X,Y,Z format----------

	! OPEN(UNIT=3,FILE='X.txt')
        ! do i=1,nx+1
           ! READ(3,*)X(i)
        ! enddo
	! CLOSE(3)
	
	! OPEN(UNIT=5,FILE='Y.txt')
        ! do j=1,ny+1
           ! READ(5,*)Y(j)
        ! enddo
	! CLOSE(5)
	
	! OPEN(UNIT=7,FILE='Z.txt')
        ! do k=1,nz+1
           ! READ(7,*)Z(k)
        ! enddo
	! CLOSE(7)
! !-----------------------------------------------------------------

    ! X(0)=X(1)-(X(2)-X(1))
    ! X(-1)=X(0)-(X(1)-X(0))
    ! X(nx+2)=X(nx+1)+(X(nx+1)-X(nx))
    ! X(nx+3)=X(nx+2)+(X(nx+2)-X(nx+1))

    ! Y(0)=Y(1)-(Y(2)-Y(1))
    ! Y(-1)=Y(0)-(Y(1)-Y(0))
    ! Y(ny+2)=Y(ny+1)+(Y(ny+1)-Y(ny))
    ! Y(ny+3)=Y(ny+2)+(Y(ny+2)-Y(ny+1))

    ! Z(0)=Z(1)-(Z(nz+1)-Z(nz))
    ! Z(-1)=Z(0)-(Z(nz)-Z(nz-1))
    ! Z(nz+2)=Z(nz+1)+(Z(2)-Z(1))
    ! Z(nz+3)=Z(nz+2)+(Z(3)-Z(2))


    ! !Define each of the directional grid lengths
    ! do i=1,nx-1

        ! iDx(i) = ( X(i+1) - X(i) )
        ! Dxs(i) = ( X(i+2) - X(i) ) *0.5d0

    ! end do

    ! do j=1,ny-1

        ! iDy(j) = ( Y(j+1) - Y(j) )
        ! Dys(j) = ( Y(j+2) - Y(j) ) *0.5d0

    ! end do

    ! do k=1,nz-1

        ! iDz(k) = ( Z(k+1) - Z(k) )
        ! Dzs(k) = ( Z(k+2) - Z(k) ) *0.5d0

    ! end do



    ! !Ghost boundary grid lengths
    ! iDx(0) = iDx(1)
    ! iDx(-1) = iDx(1)
    ! iDx(nx) = X(nx+1) - X(nx)
    ! iDx(nx+1) = iDx(nx)
    ! iDx(nx+2) = iDx(nx)

    ! Dxs(0) = Dxs(1)
    ! Dxs(-1) = Dxs(1)
    ! Dxs(nx) = Dxs(nx-1)
    ! Dxs(nx+1) = Dxs(nx-1)
    ! Dxs(nx+2) = Dxs(nx-1)


    ! iDy(0) = iDy(1)
    ! iDy(-1) = iDy(1)
    ! iDy(ny) = Y(ny+1) - Y(ny)
    ! iDy(ny+1) = iDy(ny)
    ! iDy(ny+2) = iDy(ny)

    ! Dys(0) = Dys(1)
    ! Dys(-1) = Dys(1)
    ! Dys(ny) = Dys(ny-1)
    ! Dys(ny+1) = Dys(ny-1)
    ! Dys(ny+2) = Dys(ny-1)


    ! iDz(0) = iDz(1)
    ! iDz(-1) = iDz(1)
    ! iDz(nz) = Z(nz+1) - Z(nz)
    ! iDz(nz+1) = iDz(nz)
    ! iDz(nz+2) = iDz(nz)

    ! Dzs(0) = Dzs(1)
    ! Dzs(-1) = Dzs(1)
    ! Dzs(nz) = Dzs(nz-1)
    ! Dzs(nz+1) = Dzs(nz-1)
    ! Dzs(nz+2) = Dzs(nz-1)

    ! ! dx = MAXVAL(iDx) 
    ! ! dy = MAXVAL(iDy)
    ! ! dz = MAXVAL(iDz)
	  
	! dxSml = MINVAL(iDx)
    ! dySml = MINVAL(iDy)
    ! dzSml = MINVAL(iDz)
	
	

    ! !Defining the midpoint values of the grids
	! do i=1,nx
        ! Xs(i) = 0.5d0 * ( X(i+1) + X(i) )
		! if ((X(i) .LE. x_start) .AND. (X(i+1) .GT. x_start) ) then 	! filer boundary
			! i_start = i
		! elseif (x_end .GE. lx) then
			! i_end = nx
		! elseif ((X(i) .LT. x_end) .AND. (X(i+1) .GE. x_end) ) then
			! i_end = i
		! endif
    ! end do

	! do j=1,ny
        ! Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
		! if ((Y(j) .LE. y_start) .AND. (Y(j+1) .GT. y_start) ) then 	! filer boundary
			! j_start = j
		! elseif (y_end .GE. ly) then
			! j_end = ny
		! elseif ((Y(j) .LT. y_end) .AND. (Y(j+1) .GE. y_end) ) then
			! j_end = j
		! endif
    ! end do
	
	! do k=1,nz
        ! Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
		! if ((Z(k) .LE. z_start) .AND. (Z(k+1) .GT. z_start) ) then 	! filer boundary
			! k_start = k
		! elseif (z_end .GE. lz) then
			! k_end = nz
		! elseif ((Z(k) .LT. z_end) .AND. (Z(k+1) .GE. z_end) ) then
			! k_end = k
		! endif
    ! end do


    ! !Output values of the grids
    ! do k=1,nz; do j=1,ny; do i=1,nx
        ! Xout(i,j,k) = Xs(i)
        ! Yout(i,j,k) = Ys(j)
        ! Zout(i,j,k) = Zs(k)
    ! enddo; enddo; enddo


    ! if(myid==master)then
	    ! CALL SYSTEM('mkdir -p ' // TRIM('output3D/'))
        ! open (unit=1,form='unformatted',file=TRIM('output3D/')//'3d_mesh.x')
        ! write(1) nblocks
        ! write(1) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1

        ! write(1)    (((Xout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    ! (((Yout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    ! (((Zout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end)

        ! close(1)
    ! end if
	
    ! if(myid==master)then
	    ! CALL SYSTEM('mkdir -p ' // TRIM('output2D/'))
        ! open (unit=2,form='unformatted',file=TRIM('output2D/')//'2d_mesh.x')
        ! write(2) nblocks
        ! write(2) 1, ny, nz

        ! write(2)    (((Xout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    ! (((Yout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    ! (((Zout(i,j,k),i=1,1),j=1,ny),k=1,nz)

        ! close(2)
    ! end if

	! ! !------------------- Write formatted grid output -----------------------
    ! ! if(myid==master)then
	    ! ! CALL SYSTEM('mkdir -p ' // TRIM('output3D_formatted/'))
        ! ! open (unit=4,file=TRIM('output3D_formatted/')//'3d_mesh.DAT')

        ! ! do i=1,nx ; do j=1,ny ; do k=1,nz
           ! ! write(4,*)Xs(i),Ys(j),Zs(k)
        ! ! enddo; enddo; enddo

        ! ! close(4)
    ! ! end if
	! ! !------------------- Write formatted grid output -----------------------


! end subroutine gridder_read_txt

!------------------------------------------------------------------------------------------------------------------

subroutine gridder_fall()
    use variables
    implicit none
real*8 :: nratio

real*8, parameter :: lxLrg1 = xc - lxSml*0.5d0
real*8, parameter :: lxLrg2 = lx - lxLrg1 - lxSml

integer, parameter :: nxLrg1 = nxLow !INT( (ny-nySml)/(ly-lySml)*lyLrg1 )
integer, parameter :: nxLrg2 = nx-nxLrg1-nxSml

real*8, parameter :: dxLrg1 = lxLrg1/nxLrg1*1.d0
real*8, parameter :: dxLrg2 = lxLrg2/nxLrg2*1.d0

real*8, parameter :: lyLrg1 = GridderYc - lySml*0.5d0
real*8, parameter :: lyLrg2 = ly - lyLrg1 - lySml

integer, parameter :: nyLrg1 = nyLow !INT( (ny-nySml)/(ly-lySml)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nySml

real*8, parameter :: dyLrg1 = lyLrg1/nyLrg1*1.d0
real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: thetaa = 4.d0*ATAN(1.d0)
real*8, parameter :: tune = 1.d0
    
    !-----------------unequal grid intervals----------------!
    do i=1,nx+1

        if ( i==1 ) then
            X(1) = 0.d0

        elseif( i <= nxLrg1+1 ) then                                                                                !------------------------Sin dxLrg1-----------------------!

            nratio = (i-1.d0)/nxLrg1*1.d0 
            X(i) = (lxLrg1-dxSml*nxLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dxSml*(i-1.d0)

        elseif( i > nxLrg1+1 .AND. i <= nxLrg1+nxSml+1 ) then                                                       !------------------------Uniform dxSml-----------------------!

            X(i) = X(i-1) + dxSml

        elseif( i > nxLrg1+nxSml+1 .AND. i <= nx+1 ) then                                                           !------------------------Sin dxLrg2-----------------------!

            nratio = (i-nxLrg1-nxSml-1.d0)/nxLrg2*1.d0 
            X(i) = lxLrg1+lxSml + (lxLrg2-dxSml*nxLrg2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + dxSml*(i-nxLrg1-nxSml-1.d0)
        endif

    end do

    X(0)=X(1)-(X(2)-X(1))
    X(-1)=X(0)-(X(1)-X(0))
    X(nx+2)=X(nx+1)+(X(nx+1)-X(nx))
    X(nx+3)=X(nx+2)+(X(nx+2)-X(nx+1))

    dx=X(2)-X(1)

     do j=1,ny+1

        if ( j==1 ) then
            Y(1) = 0.d0

        elseif( j <= nyLrg1+1 ) then                                                                                !------------------------Sin dyLrg1-----------------------!

            nratio = (j-1.d0)/nyLrg1*1.d0 
            Y(j) = (lyLrg1-dySml*nyLrg1)*(nratio+tune/thetaa*SIN(nratio*thetaa)) + dySml*(j-1.d0)

        elseif( j > nyLrg1+1 .AND. j <= nyLrg1+nySml+1 ) then                                                       !------------------------Uniform dySml-----------------------!

            Y(j) = Y(j-1) + dySml

        elseif( j > nyLrg1+nySml+1 .AND. j <= ny+1 ) then                                                           !------------------------Sin dyLrg2-----------------------!

            nratio = (j-nyLrg1-nySml-1.d0)/nyLrg2*1.d0 
            Y(j) = lyLrg1+lySml + (lyLrg2-dySml*nyLrg2)*(nratio-tune/thetaa*SIN(nratio*thetaa)) + dySml*(j-nyLrg1-nySml-1.d0)
        endif

    end do

    Y(0)=Y(1)-(Y(2)-Y(1))
    Y(-1)=Y(0)-(Y(1)-Y(0))
    Y(ny+2)=Y(ny+1)+(Y(ny+1)-Y(ny))
    Y(ny+3)=Y(ny+2)+(Y(ny+2)-Y(ny+1))

    dy=Y(2)-Y(1)

    do k=1,nz+3
        if(k == 1) then
            Z(k) = 0.0
            Z(k-1) = Z(k) - dz
            Z(k-2) = Z(k-1) - dz
        else
            Z(k) = Z(k-1) + dz
        end if 
    end do
    !-----------------unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=-1,nx+1

        iDx(i) = ( X(i+1) - X(i) )
        Dxs(i) = ( X(i+2) - X(i) ) *0.5d0

    end do

    do j=-1,ny+1

        iDy(j) = ( Y(j+1) - Y(j) )
        Dys(j) = ( Y(j+2) - Y(j) ) *0.5d0

    end do

    do k=1,nz-1

        iDz(k) = ( Z(k+1) - Z(k) )
        Dzs(k) = ( Z(k+2) - Z(k) ) *0.5d0

    end do



    !Ghost boundary grid lengths
    iDx(nx+2) = iDx(nx)
    Dxs(nx+2) = Dxs(nx-1)


    ! iDy(0) = iDy(1)
    ! iDy(-1) = iDy(1)
    ! iDy(ny) = Y(ny+1) - Y(ny)
    ! iDy(ny+1) = iDy(ny)
    ! iDy(ny+2) = iDy(ny)

    ! Dys(0) = Dys(1)
    ! Dys(-1) = Dys(1)
    ! Dys(ny) = Dys(ny-1)
    ! Dys(ny+1) = Dys(ny-1)
    ! Dys(ny+2) = Dys(ny-1)
	
	iDy(ny+2) = iDy(ny)
    Dys(ny+2) = Dys(ny-1)


    iDz(0) = iDz(1)
    iDz(-1) = iDz(1)
    iDz(nz) = Z(nz+1) - Z(nz)
    iDz(nz+1) = iDz(nz)
    iDz(nz+2) = iDz(nz)

    Dzs(0) = Dzs(1)
    Dzs(-1) = Dzs(1)
    Dzs(nz) = Dzs(nz-1)
    Dzs(nz+1) = Dzs(nz-1)
    Dzs(nz+2) = Dzs(nz-1)

    !Modifying the index of X, Y and Z arrays to represent the actual grid
    !do i=1,nx+1
    !    Xa(i) = X(i)
    !end do
!
    !do j=1,ny+1
    !    Ya(j) = Y(j)
    !end do
!
    !do k=1,nz+1
    !    Za(k) = Z(k)
    !end do

    !Defining the midpoint values of the grids
	do i=1,nx
        Xs(i) = 0.5d0 * ( X(i+1) + X(i) )
		if ((X(i) .LE. x_start) .AND. (X(i+1) .GT. x_start) ) then 	! filer boundary
			i_start = i
		elseif (x_end .GE. lx) then
			i_end = nx
		elseif ((X(i) .LT. x_end) .AND. (X(i+1) .GE. x_end) ) then
			i_end = i
		endif
    end do

	do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
		if ((Y(j) .LE. y_start) .AND. (Y(j+1) .GT. y_start) ) then 	! filer boundary
			j_start = j
		elseif (y_end .GE. ly) then
			j_end = ny
		elseif ((Y(j) .LT. y_end) .AND. (Y(j+1) .GE. y_end) ) then
			j_end = j
		endif
    end do
	
	do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
		if ((Z(k) .LE. z_start) .AND. (Z(k+1) .GT. z_start) ) then 	! filer boundary
			k_start = k
		elseif (z_end .GE. lz) then
			k_end = nz
		elseif ((Z(k) .LT. z_end) .AND. (Z(k+1) .GE. z_end) ) then
			k_end = k
		endif
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
        Yout(i,j,k) = Ys(j)
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo


    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output3D/'))
        open (unit=1,form='unformatted',file=TRIM('output3D/')//'3d_mesh.x')
        write(1) nblocks
        write(1) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1

        write(1)    (((Xout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Yout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end), &
                    (((Zout(i,j,k),i=i_start,i_end),j=j_start,j_end),k=k_start,k_end)

        close(1)
    end if
	
    if(myid==master)then
	    CALL SYSTEM('mkdir -p ' // TRIM('output2D/'))
        open (unit=2,form='unformatted',file=TRIM('output2D/')//'2d_mesh.x')
        write(2) nblocks
        write(2) 1, ny, nz

        write(2)    (((Xout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,1),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,1),j=1,ny),k=1,nz)

        close(2)
    end if


end subroutine gridder_fall