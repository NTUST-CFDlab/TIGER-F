subroutine gridder_Unequal()
!--------------------------ORIGINAL------------------------------
!use variables
!implicit none
!integer :: ncount, ncounts

!------------Following are calculated using the above variables------------!
!real*8 ,parameter :: ySBgn = GridderYc 
!real*8 ,parameter :: ySEnd = GridderYc + lySml
!real*8 ,parameter :: yMBgn = GridderYc 
!real*8 ,parameter :: yMEnd = GridderYc + lyMid

!real*8 ,parameter :: zSBgn = GridderZc
!real*8 ,parameter :: zSEnd = GridderZc + lzSml
!real*8 ,parameter :: zMBgn = GridderZc
!real*8 ,parameter :: zMEnd = GridderZc + lzMid


!real*8 :: yNextLrgValue
!real*8 :: yNextMidValue
!real*8 :: yNextSmlValue

!real*8 :: zNextLrgValue
!real*8 :: zNextMidValue
!real*8 :: zNextSmlValue
!------------Following are calculated using the above variables------------!


!----------------Unequal grid intervals----------------!
!    do i=1,nx+3
!        if(i == 1) then
!            X(i) = 0.0
!            X(i-1) = X(i) - dx
!            X(i-2) = X(i-1) - dx
!        else
!            X(i) = X(i-1) + dx
!        end if 
!    end do


!    do j=1,ny+1
!        yNextLrgValue =  Y(j-1) + dy   
!        yNextMidValue =  Y(j-1) + dyMid   
!        yNextSmlValue =  Y(j-1) + dySml   

!        if (j==1) then
!            Y(j) = 0.0
!            ncounts=0
!            ncount=0
!        elseif(yNextSmlValue > ySBgn .AND. ncounts < nySml) then
!            Y(j) = yNextSmlValue
!            ncounts=ncounts+1
!        elseif(yNextMidValue > yMBgn .AND. ncount < nyMid-nySml) then
!            Y(j) = yNextMidValue
!            ncount=ncount+1
!        else
!            Y(j) = yNextLrgValue
!        end if
!    end do

   

!    do k=1,nz+1
!        zNextLrgValue =  Z(k-1) + dz   
!        zNextMidValue =  Z(k-1) + dzMid   
!        zNextSmlValue =  Z(k-1) + dzSml   

!        if (k==1) then
!            Z(k) = 0.0
!            ncounts=0
!            ncount=0
!        elseif(zNextSmlValue > zSBgn .AND. ncounts < nzSml) then
!            Z(k) = zNextSmlValue
!            ncounts=ncounts+1
!        elseif(zNextMidValue > zMBgn .AND. ncount < nzMid-nzSml) then
!            Z(k) = zNextMidValue
!            ncount=ncount+1
!        else
!            Z(k) = zNextLrgValue
!        end if
!    end do
	

!-----------------------MODIFIED-FANDI-------------------------------------!	
use variables
implicit none
integer :: ncount

!------------Following are calculated using the above variables------------!
real*8 ,parameter :: lyLrg1 = GridderYc - 0.5*lyMid 
real*8, parameter :: lyLrg2 = ly - lyLrg1 - lyMid

integer, parameter :: nyLrg1 = INT( (ny-nyMid-nySml)/(ly-lyMid)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nyMid

real*8, parameter :: dyLrg = ( ly-lyMid ) / ( ny - nyMid - nySml )

real*8 ,parameter :: lyMid1 = 0.5*(lyMid-lySml)
real*8, parameter :: lyMid2 = lyMid-lyMid1

integer, parameter :: nyMid1 = INT( 0.5*nyMid )
integer, parameter :: nyMid2 = nyMid-nyMid1

real*8 ,parameter :: lzLrg1 = GridderZc - 0.5*lzMid 
real*8, parameter :: lzLrg2 = lz - lzLrg1 - lzMid

integer, parameter :: nzLrg1 = INT( (nz-nzMid-nzSml)/(lz-lzMid)*lzLrg1 )
integer, parameter :: nzLrg2 = nz-nzLrg1-nzMid

real*8, parameter :: dzLrg = ( lz-lzMid ) / ( nz - nzMid - nzSml )

real*8 ,parameter :: lzMid1 = 0.5*(lzMid-lzSml)
real*8, parameter :: lzMid2 = lzMid-lzMid1

integer, parameter :: nzMid1 = INT( 0.5*nzMid )
integer, parameter :: nzMid2 = nzMid-nzMid1

real*8 :: yNextStep
real*8 :: yNextLrgValue
real*8 :: yNextMidValue
real*8 :: yNextSmlValue

real*8 :: zNextStep
real*8 :: zNextLrgValue
real*8 :: zNextMidValue
real*8 :: zNextSmlValue
!------------Following are calculated using the above variables------------!


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
        yNextLrgValue =  Y(j-1) + dyLrg   
        yNextMidValue =  Y(j-1) + dyMid   
        yNextSmlValue =  Y(j-1) + dySml
		yNextStep = Y(j-1) + 0.1*dySml

        if (j==1) then
            Y(j) = 0.0
            ncount=0
		elseif(yNextStep <= lyLrg1-dyMid .AND. ncount <= nyLrg1) then
			Y(j) = yNextLrgValue
            ncount=ncount+1			
		elseif(yNextStep > lyLrg1-dyMid .AND. ncount <= nyLrg1+nyMid1) then
			Y(j) = yNextMidValue
            ncount=ncount+1	
        elseif(yNextStep > lyLrg1+lyMid1-dySml .AND. ncount <= nyLrg1+nyMid1+nySml) then
            Y(j) = yNextSmlValue
            ncount=ncount+1	
        elseif(yNextStep > lyLrg1+lyMid1+lySml .AND. ncount <= nyLrg1+nyMid1+nySml+nyMid2) then
            Y(j) = yNextMidValue
            ncount=ncount+1
        else
            Y(j) = yNextLrgValue
        end if
    end do
!--------------------------ORIGINAL------------------------------
   

    do k=1,nz+1
        zNextLrgValue =  Z(k-1) + dzLrg   
        zNextMidValue =  Z(k-1) + dzMid   
        zNextSmlValue =  Z(k-1) + dzSml  
		zNextStep = Z(k-1) + 0.1*dzSml		

        if (k==1) then
            Z(k) = 0.0
            ncount=0
		elseif(zNextStep <= lzLrg1-dzMid .AND. ncount <= nzLrg1) then
			Z(k) = zNextLrgValue
            ncount=ncount+1			
		elseif(zNextStep > lzLrg1-dzMid .AND. ncount <= nzLrg1+nzMid1) then
			Z(k) = zNextMidValue
            ncount=ncount+1	
        elseif(zNextStep > lzLrg1+lzMid1-dzSml .AND. ncount <= nzLrg1+nzMid1+nzSml) then
            Z(k) = zNextSmlValue
            ncount=ncount+1	
        elseif(zNextStep > lzLrg1+lzMid1+lzSml .AND. ncount <= nzLrg1+nzMid1+nzSml+nzMid2) then
            Z(k) = zNextMidValue
            ncount=ncount+1
        else
            Z(k) = zNextLrgValue
        end if
    end do
    !----------------Unequal grid intervals----------------!

    !Define each of the directional grid lengths
    do i=1,nx-1

        iDx(i) = ( X(i+1) - X(i) )                      ! Cell length(x)
        Dxs(i) = ( X(i+2) - X(i) ) *0.5d0               ! Distance between node centers (E to P or P to W)

    end do

    do j=1,ny-1

        iDy(j) = ( Y(j+1) - Y(j) )                      ! Cell length(y)
        Dys(j) = ( Y(j+2) - Y(j) ) *0.5d0               ! Distance between node centers (N to P or P to S)

    end do

    do k=1,nz-1

        iDz(k) = ( Z(k+1) - Z(k) )                      ! Cell length(z)
        Dzs(k) = ( Z(k+2) - Z(k) ) *0.5d0               ! Distance between node centers (B to P or P to F)

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


    !Defining the midpoint values of the grids
    do i=1,nx
        Xs(i) = 0.5d0 * ( X(i+1) + X(i) )
    end do

    do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
    end do

    do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
    end do



    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Yout(i,j,k) = Ys(j)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo

    if(myid==master)then
        open (unit=1,form='unformatted',file='mesh.x')
        write(1) nblocks
        write(1) nx, ny, nz

        write(1)    (((Xout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,nx),j=1,ny),k=1,nz)

        close(1)
    end if


end subroutine gridder_Unequal


subroutine gridder_sin()
use variables
implicit none
integer :: ncount
real*8 :: thetaa
real*8, parameter :: lyMid1 = (lyMid-lySml)*0.3d0
real*8, parameter :: lyMid2 = (lyMid-lySml)*0.7d0

integer, parameter :: nyMid1 = INT( (nyMid-nySml)*0.3d0 )
integer, parameter :: nyMid2 = nyMid-nySml-nyMid1

real*8, parameter :: lyLrg1 = GridderYc - lySml*0.5d0 - lyMid1
real*8, parameter :: lyLrg2 = ly - lyLrg1 - lyMid

integer, parameter :: nyLrg1 = INT( (ny-nyMid)/(ly-lyMid)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nyMid

real*8, parameter :: dyLrg1 = lyLrg1/nyLrg1*1.d0
real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: zLrg1 = GridderZc
real*8, parameter :: zLrg2 = lz - GridderZc - lzSml

integer, parameter :: nzLrg1 = INT( (nz-nzMid)/(lz-lzMid)*zLrg1 )  !1 
integer, parameter :: nzLrg2 = nz-nzLrg1-nzMid                  !  89

real*8, parameter :: tune = 0.2d0



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

            thetaa = 0.5d0*PI*(j-1.d0)/nyLrg1*1.d0 
            Y(j) = (lyLrg1-tune*dyLrg1*nyLrg1)*sin(thetaa) + tune*dyLrg1*(j-1.d0)

        elseif( j > nyLrg1+1 .AND. j <= nyLrg1+nyMid1+1 ) then                                                      !------------------------Sin dyMid1-----------------------!

            thetaa = 0.5d0*PI*(j-nyLrg1-1.d0)/nyMid1*1.d0 
            Y(j) = lyLrg1 + (lyMid1-dySml*nyMid1)*sin(thetaa) + dySml*(j-nyLrg1-1.d0)

        elseif( j > nyLrg1+nyMid1+1 .AND. j <= nyLrg1+nyMid1+nySml+1 ) then                                         !------------------------Uniform dySml-----------------------!

            Y(j) = lyLrg1+lyMid1+ dySml*(j-nyLrg1-nyMid1-1)

        elseif( j > nyLrg1+nyMid1+nySml+1 .AND. j <= nyLrg1+nyMid+1 ) then                                          !------------------------Sin dyMid2-----------------------!

            thetaa = -0.5d0*PI + 0.5d0*PI*(j-nyLrg1-nyMid1-nySml-1.d0)/nyMid2*1.d0 
            Y(j) = GridderYc + lySml*0.5d0 + (lyMid2-dySml*nyMid2)*(1.d0+sin(thetaa)) + dySml*(j-nyLrg1-nyMid1-nySml-1.d0)

        elseif( j > nyLrg1+nyMid+1 .AND. j <= ny+1 ) then                                                           !------------------------Sin dyLrg-----------------------!

            thetaa = -0.5d0*PI + 0.5d0*PI*(j-nyLrg1-nyMid-1.d0)/nyLrg2*1.d0 
            Y(j) = GridderYc + lySml*0.5d0 + lyMid2 + (lyLrg2-tune*dyLrg2*nyLrg2)*(1.d0+sin(thetaa)) + tune*dyLrg2*(j-nyLrg1-nyMid-1.d0)

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
        
            thetaa = 0.5d0*PI*(k-1.d0)/nzLrg1*1.d0 
            Z(k) = (zLrg1-dzSml*nzLrg1)*sin(thetaa) + dzSml*(k-1.d0)

        elseif( k > nzLrg1 .AND. k <= nzLrg1+nzSml+1 ) then                                                         !------------------------Uniform dzSml-----------------------!
            Z(k) = zLrg1 + dzSml*(k-nzLrg1-1.d0)

        elseif( k > nzLrg1+nzSml+1 .AND. k <= nzLrg1+nzMid+1 ) then                                                 !------------------------Sin dzMid-----------------------!

            thetaa = -0.5d0*PI + 0.5d0*PI*(k-nzLrg1-nzsml-1.d0)/(nzMid-nzSml)*1.d0 
            Z(k) = zLrg1 + lzSml + ((lzMid-lzSml)-dzsml*(nzMid-nzSml))*(1.d0+sin(thetaa)) + dzSml*(k-nzLrg1-nzSml-1.d0)

        elseif( k > nzLrg1+nzMid+1 .AND. k <= nz+1 ) then                                                           !------------------------Sin dzLrg2-----------------------!

            thetaa = -0.5d0*PI + 0.5d0*PI*(k-nzLrg1-nzMid-1.d0)/nzLrg2*1.d0 
            Z(k) = zLrg1 + lzMid + (zLrg2-dzMid*nzLrg2)*(1.d0+sin(thetaa)) + dzMid*(k-nzLrg1-nzMid-1.d0)

        endif

    end do

    Z(0)=Z(1)-(Z(2)-Z(1))
    Z(-1)=Z(0)-(Z(1)-Z(0))
    Z(nz+2)=Z(nz+1)+(Z(ny+1)-Z(ny))
    Z(nz+3)=Z(nz+2)+(Z(ny+2)-Z(ny+1))

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
    end do

    do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
    end do

    do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Yout(i,j,k) = Ys(j)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo

    if(myid==master)then
        open (unit=1,form='unformatted',file='mesh.x')
        write(1) nblocks
        write(1) nx, ny, nz

        write(1)    (((Xout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,nx),j=1,ny),k=1,nz)

        close(1)
    end if


end subroutine gridder_sin



subroutine gridder_sin2()
use variables
implicit none
integer :: ncount
real*8 :: thetaa
real*8, parameter :: lyLrg1 = GridderYc - lySml*0.5d0
real*8, parameter :: lyLrg2 = ly - (lyLrg1 + lyMid)

real*8, parameter :: lyMid1 = lyMid-lySml

integer, parameter :: nyLrg1 = INT( (ny-nyMid)/(ly-lyMid)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nyMid

real*8, parameter :: dyLrg1 = lyLrg1/nyLrg1*1.d0
real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: zLrg1 = GridderZc
real*8, parameter :: zLrg2 = lz - GridderZc - lzSml

integer, parameter :: nzLrg1 = INT( (nz-nzMid)/(lz-lzMid)*zLrg1 )
integer, parameter :: nzLrg2 = nz-nzLrg1-nzMid

real*8, parameter :: tune = 0.1d0



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

            thetaa = 0.5d0*PI*(j-1.d0)/nyLrg1*1.d0 
            Y(j) = (lyLrg1-dySml*nyLrg1)*SIN(thetaa) + dySml*(j-1.d0)

        elseif( j > nyLrg1+1 .AND. j <= nyLrg1+nySml+1 ) then                                                       !------------------------Uniform dySml-----------------------!

            Y(j) = lyLrg1+ dySml*(j-nyLrg1-1)

        elseif( j > nyLrg1+nySml+1 .AND. j <= nyLrg1+nyMid+1 ) then                                                 !------------------------Sin dyMid-----------------------!

            thetaa = -0.5d0*PI + 0.5d0*PI*(j-nyLrg1-nySml-1.d0)/(nyMid-nySml)*1.d0 
            Y(j) = GridderYc + lySml*0.5d0 + (lyMid1-dySml*(nyMid-nySml))*(1.d0+SIN(thetaa)) + dySml*(j-nyLrg1-nySml-1.d0)

        elseif( j > nyLrg1+nyMid+1 .AND. j <= ny+1 ) then                                                           !------------------------Sin dyLrg2-----------------------!

            thetaa = -0.5d0*PI + 0.5d0*PI*(j-nyLrg1-nyMid-1.d0)/nyLrg2*1.d0 
            Y(j) = lyLrg1 + lyMid + (lyLrg2-tune*dyLrg2*nyLrg2)*(1.d0+SIN(thetaa)) + tune*dyLrg2*(j-nyLrg1-nyMid-1.d0)

        endif

    end do


    Y(0)=Y(1)-(Y(ny+1)-Y(ny))
    Y(-1)=Y(0)-(Y(ny)-Y(ny-1))
    Y(ny+2)=Y(ny+1)+(Y(2)-Y(1))
    Y(ny+3)=Y(ny+2)+(Y(3)-Y(2))

    dy=Y(2)-Y(1)



    do k=1,nz+1
        
        if ( k==1 ) then
            Z(k) = 0.d0
        elseif( k <= nzLrg1 ) then                                                                                 !------------------------Sin dzLrg1-----------------------!
        
            thetaa = 0.5d0*PI*(k-1.d0)/nzLrg1*1.d0 
            Z(k) = (zLrg1-dzMid*nzLrg1)*SIN(thetaa) + dzMid*(k-1.d0)

        elseif( k > nzLrg1 .AND. k <= nzLrg1+nzSml+1 ) then                                                        !------------------------Uniform dzSml-----------------------!
            Z(k) = zLrg1 + dzSml*(k-nzLrg1-1.d0)

        elseif( k > nzLrg1+nzSml+1 .AND. k <= nzLrg1+nzMid+1 ) then                                                !------------------------Uniform dzMid-----------------------!
            Z(k) = zLrg1 + lzSml + dzMid*(k-nzLrg1-nzSml-1.d0)

        elseif( k > nzLrg1+nzMid+1 .AND. k <= nz+1 ) then                                                          !------------------------Sin dyLrg2-----------------------!

            thetaa = -0.5d0*PI + 0.5d0*PI*(k-nzLrg1-nzMid-1.d0)/nzLrg2*1.d0 
            Z(k) = zLrg1 + lzMid + (zLrg2-dzMid*nzLrg2)*(1+SIN(thetaa)) + dzMid*(k-nzLrg1-nzMid-1.d0)

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
    end do

    do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
    end do

    do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Yout(i,j,k) = Ys(j)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo

    if(myid==master)then
        open (unit=1,form='unformatted',file='mesh.x')
        write(1) nblocks
        write(1) nx, ny, nz

        write(1)    (((Xout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,nx),j=1,ny),k=1,nz)

        close(1)
    end if


end subroutine gridder_sin2





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
    end do

    do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
    end do

    do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
    end do

    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Yout(i,j,k) = Ys(j)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo

    if(myid==master)then
        open (unit=1,form='unformatted',file='mesh.x')
        write(1) nblocks
        write(1) nx, ny, nz

        write(1)  (((Xout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                  (((Yout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                  (((Zout(i,j,k),i=1,nx),j=1,ny),k=1,nz)

        

        close(1)
    end if


end subroutine gridder_equal

!-----------------------MODIFIED-FANDI-------------------------------------!

subroutine gridder_sin3() ! Simple SIN function from 0 to pi/2
use variables
implicit none
integer :: ncount
real*8 :: thetaa

real*8, parameter :: lyLrg1 = GridderYc - lySml*0.5d0
real*8, parameter :: lyLrg2 = ly - lyLrg1 - lySml

integer, parameter :: nyLrg1 = INT( (ny-nySml)/(ly-lySml)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nySml

real*8, parameter :: dyLrg1 = lyLrg1/nyLrg1*1.d0
real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: lzLrg1 = GridderZc - lzSml*0.5d0
real*8, parameter :: lzLrg2 = lz - lzLrg1 - lzSml

integer, parameter :: nzLrg1 = INT( (nz-nzSml)/(lz-lzSml)*lzLrg1 )  !1 
integer, parameter :: nzLrg2 = nz-nzLrg1-nzSml                     !89

real*8, parameter :: dzLrg1 = lzLrg1/nzLrg1*1.d0
real*8, parameter :: dzLrg2 = lzLrg2/nzLrg2*1.d0

! real*8, parameter :: tune = 0.2d0



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

            thetaa = 0.5d0*PI*(j-1.d0)/nyLrg1*1.d0 
            Y(j) = (lyLrg1-dySml*nyLrg1)*SIN(thetaa) + dySml*(j-1.d0)

        elseif( j > nyLrg1+1 .AND. j <= nyLrg1+nySml+1 ) then                                                       !------------------------Uniform dySml-----------------------!

            Y(j) = Y(j-1) + dySml

        elseif( j > nyLrg1+nySml+1 .AND. j <= ny+1 ) then                                                           !------------------------Sin dyLrg2-----------------------!

            thetaa = -0.5d0*PI + 0.5d0*PI*(j-nyLrg1-nySml-1.d0)/nyLrg2*1.d0 
			Y(j) = lyLrg1+lySml + (lyLrg2-dySml*nyLrg2)*(1.d0+SIN(thetaa)) + dySml*(j-nyLrg1-nySml-1.d0)
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
        
            thetaa = 0.5d0*PI*(k-1.d0)/nzLrg1*1.d0 
            Z(k) = (lzLrg1-dzSml*nzLrg1)*SIN(thetaa) + dzSml*(k-1.d0)

        elseif( k > nzLrg1 .AND. k <= nzLrg1+nzSml+1 ) then                                                         !------------------------Uniform dzSml-----------------------!
            Z(k) = Z(k-1) + dzSml

        elseif( k > nzLrg1+nzSml+1 .AND. k <= nz+1 ) then                                                           !------------------------Sin dzLrg2-----------------------!

            thetaa = -0.5d0*PI + 0.5d0*PI*(k-nzLrg1-nzSml-1.d0)/nzLrg2*1.d0 
			Z(k) =  lzLrg1+lzSml + (lzLrg2-dzSml*nzLrg2)*(1.d0+SIN(thetaa)) + dzSml*(k-nzLrg1-nzSml-1.d0)

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
    end do

    do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
    end do

    do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Yout(i,j,k) = Ys(j)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo

    if(myid==master)then
        open (unit=1,form='unformatted',file='mesh.x')
        write(1) nblocks
        write(1) nx, ny, nz

        write(1)    (((Xout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,nx),j=1,ny),k=1,nz)

        close(1)
    end if


end subroutine gridder_sin3


!-----------------------MODIFIED-FANDI-------------------------------------!

subroutine gridder_sin4() ! Grid spacing method proposed by Kuyper
use variables
implicit none
integer :: ncount
real*8 :: nratio

real*8, parameter :: lyLrg1 = GridderYc - lySml*0.5d0
real*8, parameter :: lyLrg2 = ly - lyLrg1 - lySml

integer, parameter :: nyLrg1 = INT( (ny-nySml)/(ly-lySml)*lyLrg1 )
integer, parameter :: nyLrg2 = ny-nyLrg1-nySml

real*8, parameter :: dyLrg1 = lyLrg1/nyLrg1*1.d0
real*8, parameter :: dyLrg2 = lyLrg2/nyLrg2*1.d0

real*8, parameter :: lzLrg1 = GridderZc - lzSml*0.5d0
real*8, parameter :: lzLrg2 = lz - lzLrg1 - lzSml

integer, parameter :: nzLrg1 = INT( (nz-nzSml)/(lz-lzSml)*lzLrg1 )  !1 
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
    end do

    do j=1,ny
        Ys(j) = 0.5d0 * ( Y(j+1) + Y(j) )
    end do

    do k=1,nz
        Zs(k) = 0.5d0 * ( Z(k+1) + Z(k) )
    end do


    !Output values of the grids
    do k=1,nz; do j=1,ny; do i=1,nx
        Xout(i,j,k) = Xs(i)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Yout(i,j,k) = Ys(j)
    enddo; enddo; enddo

    do k=1,nz; do j=1,ny; do i=1,nx
        Zout(i,j,k) = Zs(k)
    enddo; enddo; enddo

    if(myid==master)then
        open (unit=1,form='unformatted',file='mesh.x')
        write(1) nblocks
        write(1) nx, ny, nz

        write(1)    (((Xout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Yout(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
                    (((Zout(i,j,k),i=1,nx),j=1,ny),k=1,nz)

        close(1)
    end if


end subroutine gridder_sin4