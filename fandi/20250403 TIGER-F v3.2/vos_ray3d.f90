! 31 Oct 2024 - FDS

subroutine read_stl()
use variables
use mpi
implicit none
integer                             :: n
character(len=20)                   :: useless
real*4                              :: tun



    !-------------------------------------------------READ STL FILE------------------------------------------------!

    open(9,file=stl_file)
    read(9,*)                                                   ! solid name
    do n=1,nf
        read(9,*) useless                                       ! seven lines to be read until reach to endsolid <name>
        if(useless=='endsolid')then                             ! determine file length (numbers of rows)
            nf=n-1
            exit 
        end if
        read(9,*)
        read(9,*) 
        read(9,*)                                               ! three facet normal + three vertecies + loop
        read(9,*) 
        read(9,*)
        read(9,*)
    end do 
    close(9)
    
        if(myid==master)then
                write(*,*)'ASCII, NF=',nf
        endif

    allocate(iFN(3,nf))
    allocate(iPx(3*nf))                                          ! Size of arrays determined for facet normal & vertices
    allocate(iPy(3*nf))
    allocate(iPz(3*nf))

    !--------------------------------------Adjust heading direction-----------------------------------

    open(9,file=stl_file)
    read(9,*)
    do n=1,nf
        read(9,*) useless, useless, iFN(1,n), iFN(2,n), iFN(3,n)                  ! facet normal nx ny nz (or swapped axis)
        read(9,*)                                                              ! outer loop
        read(9,*) useless, iPx((n-1)*3+1), iPy((n-1)*3+1), iPz((n-1)*3+1)         ! vertex v1x v1y v1z (or swapped axis)
        read(9,*) useless, iPx((n-1)*3+2), iPy((n-1)*3+2), iPz((n-1)*3+2)         ! vertex v2x v2y v2z (or swapped axis)
        read(9,*) useless, iPx((n-1)*3+3), iPy((n-1)*3+3), iPz((n-1)*3+3)         ! vertex v3x v3y v3z (or swapped axis)
        read(9,*)                                                              ! endloop
        read(9,*)                                                              ! endfacet
    enddo
    close(9)

    !!$OMP PARALLEL DO
	!!$acc parallel loop independent
    !do n=1,nf
    !    FN(3,n)=-FN(3,n)					! to mirror the solid in z direction ( xy plane )
    !    Pz((n-1)*3+1)=-Pz((n-1)*3+1)
    !    Pz((n-1)*3+2)=-Pz((n-1)*3+2)
    !    Pz((n-1)*3+3)=-Pz((n-1)*3+3)
    !enddo
	!!$acc end parallel
    !!$OMP END PARALLEL DO



    tun = (MAXVAL(iPz) - MINVAL(iPz))/ z_length3D


	!$OMP PARALLEL DO 
	!$acc parallel device_type(host)
	!$acc loop
    do n=1,nf*3
        iPx(n) = iPx(n)/tun                      ! to scale the solid dimensions 
        iPy(n) = iPy(n)/tun
        iPz(n) = iPz(n)/tun	
    enddo
	!$acc end parallel
    !$OMP END PARALLEL DO


end subroutine read_stl




subroutine vos_ray3d()				! For GPU, but need improvement
use variables
implicit none
integer ,parameter                  :: bub=60
real*4                              :: side_a,side_b,side_c,side_s,range_Z,range_Y,facet_area,facet_avg
real*4,dimension(1:nf*3)			:: Px,Py,Pz
real*4,dimension(1:3,1:nf)			:: FN
integer                             :: index1,index2
real*8  ,dimension(-1:nx+2,-1:ny+2,-1:nz+2)         :: ETA_sub
real*4                              :: Xbgn, Xend, Ybgn, Yend, Zbgn, Zend, Xsub, Ysub, Zsub
real*4                              :: c1, c2, c3, bubble
real*4  ,dimension(1:bub)           :: ray_inter
integer                             :: n, nc, bubi, bubj, nsy, nsz, off, nn
integer ,dimension(:,:), allocatable:: nk
integer ,dimension(:,:,:), allocatable :: n_



!	Enter solid transformation here    (running on CPU to save GPU VRAM. feel free to move to GPU if neeed)

	!$OMP PARALLEL DO 
	!$acc parallel device_type(host)
	!$acc loop
    do n=1,nf
		!!! Update FN
		!!! a function of f(iPz,iPy,iPx,iFN,AOA,x0_t,y0_t,z0_t)
		
		FN(1,n) = iFN(1,n)
		FN(2,n) = iFN(2,n)
		FN(3,n) = iFN(3,n)
		
    end do
	!$acc end parallel
    !$OMP END PARALLEL DO
	
	!$OMP PARALLEL DO 
	!$acc parallel device_type(host)
	!$acc loop
    do n=1,nf*3
		!!! Update Px,Py,Pz
		!!! a function of f(iPz,iPy,iPx,iFN,AOA,x0_t,y0_t,z0_t)
		
		Px(n) = iPx(n) + x0_t
		Py(n) = iPy(n) + y0_t
		Pz(n) = iPz(n) + z0_t
		
    end do
	!$acc end parallel
    !$OMP END PARALLEL DO



!$acc data present(Zs,Ys,Xs,iDz,iDy,iDx,ETA) copyin(FN,Px,Py,Pz) create(ray_inter,ETA_sub)



! updating bounds after the solid transformation 

	!$acc kernels
		Xbgn = MINVAL(Px)
		Xend = MAXVAL(Px)
		Ybgn = MINVAL(Py)
		Yend = MAXVAL(Py)
		Zbgn = MINVAL(Pz)
		Zend = MAXVAL(Pz)
	!$acc end kernels



!--------- Get  Border of solid body in the computational domain using mesh points & Define area applied VOS ..........
!................ Bounding BOX  ......................

    off=0
    do i= 1,nx
        if ( Xs(i)>Xbgn .AND. off==0 )then
            iBgnVOS=i-1
            off=1
        elseif ( Xs(i)>Xend .AND. off==1 )then
            iEndVOS=i
            exit
        endif
    enddo
    off=0
    do j= 1,ny
        if ( Ys(j)>Ybgn .AND. off==0 )then
            jBgnVOS=j-1
            off=1
        elseif ( Ys(j)>Yend .AND. off==1 )then
            jEndVOS=j
            exit
        endif
    enddo
    off=0
    do k= 0,nz
        if ( Zs(k)>Zbgn .AND. off==0 )then
            kBgnVOS=k-1
            off=1
        elseif ( Zs(k)>Zend .AND. off==1 )then
            kEndVOS=k
            exit
        endif
    enddo
 



    !-----------------------MPI DIVISION-------------------------!
    Zdv_vos = (kEndVOS-kBgnVOS+1) / nproc
    Zr_vos  = (kEndVOS-kBgnVOS+1) - Zdv_vos * nproc 
    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !i = myid
    do i=0,(nproc-1)

        if(i < Zr_vos) then
            gstart_vos(i) = kBgnVOS + i * (Zdv_vos+1)
            gend0_vos(i) = gstart_vos(i) + Zdv_vos
        else
            gstart_vos(i) = kBgnVOS + i * Zdv_vos + Zr_vos
            gend0_vos(i) = gstart_vos(i) + Zdv_vos - 1
        end if
        
        gcount_vos(i) = gend0_vos(i) - gstart_vos(i) + 1
        gend_vos(i) = gcount_vos(i) + 2

    end do

    !----------for nz vos----------!
    istart_vos = gstart_vos(myid)  !
    iend_vos = gend0_vos(myid)     !
    igcount_vos = gcount_vos(myid) !
    !----------for nz vos----------!

    !-----------------------MPI DIVISION-------------------------!

	allocate(nk(jbgnvos:jendvos,kbgnvos:kendvos) )
	allocate(n_(jbgnvos:jendvos,kbgnvos:kendvos,1:nf) )


!$acc data create(nk(jbgnvos:jendvos,istart_vos:iend_vos),n_(jbgnvos:jendvos,istart_vos:iend_vos,1:nf))


    !$OMP PARALLEL DO PRIVATE(j,nn,n)
	!$acc parallel vector_length(64)
	!$acc loop independent private(nn,n) collapse(2) gang vector
	do k=istart_vos,iend_vos
    do j=jBgnVOS,jEndVOS
		nn = 0
		!$acc loop independent reduction(+:nn)
		!!$acc loop seq
		do n=1,nf
			if( ( MAX(Pz((n-1)*3+1),Pz((n-1)*3+2),Pz((n-1)*3+3)) .GE. Zs(k)-iDz(k)/2.) .AND. (MIN(Pz((n-1)*3+1),Pz((n-1)*3+2),Pz((n-1)*3+3)) .LE. Zs(k)+iDz(k)/2.) .AND. &
				( MAX(Py((n-1)*3+1),Py((n-1)*3+2),Py((n-1)*3+3)) .GE. Ys(j)-iDy(j)/2.) .AND. (MIN(Py((n-1)*3+1),Py((n-1)*3+2),Py((n-1)*3+3)) .LE. Ys(j)+iDy(j)/2.) ) then
				nn = nn+1
				n_(j,k,nn) = n
			endif
		enddo
		nk(j,k) = nn
	enddo; enddo
	!$acc end parallel
    !$OMP END PARALLEL DO	


! !------------------------------------  VOS using Ray0asting & Subgrid -----------------------------------------

    ! !$OMP PARALLEL DO PRIVATE(nsy, nsz, n, nc, bubi, bubj, bubble, i, j, ETA_sub, ray_inter, Zsub, Ysub, c1, c2, c3, ETAs)
	! !$acc parallel vector_length(64)
	! !$acc loop independent private(nsy, nsz, n, nc, bubi, bubj, bubble, i, ETA_sub, ray_inter, Zsub, Ysub, c1, c2, c3, ETAs) collapse(2) gang vector
    ! do k=istart_vos,iend_vos
    ! do j=jBgnVOS,jEndVOS
    ! ETA_sub=0.0
	! !!$acc loop independent private(n, nc, bubi, bubj, bubble, i, ray_inter, Zsub, Ysub, c1, c2, c3, ETAs) reduction(+:ETA_sub) collapse(2) vector
	! !$acc loop seq
    ! do nsz=1,nSubGrids_3d
	! !$acc loop seq
	! do nsy=1,nSubGrids_3d         ! Each point (k,j) in the bounding box is divided into subgrids
        ! ray_inter=Xend+5.0                         ! ray starts from Xend + 5.0  ( from spanwise direction)
        ! nc=1

        ! Zsub=Zs(k)-iDz(k)/2.+(nsz-0.5)*iDz(k)/nSubGrids_3d
        ! Ysub=Ys(j)-iDy(j)/2.+(nsy-0.5)*iDy(j)/nSubGrids_3d
        
! !        !$OMP PARALLEL DO PRIVATE(c1,c2,c3) SHARED(nc,ray_inter)
		! !$acc loop seq
        ! do n=1,nk(j,k)                                                     ! Only facets with interiors intersected by the ray are included.			 
				! nn = n_(j,k,n)
                ! if ( (Zsub .GE. MIN(Pz((nn-1)*3+1),Pz((nn-1)*3+2),Pz((nn-1)*3+3)) .AND. Zsub .LE. MAX(Pz((nn-1)*3+1),Pz((nn-1)*3+2),Pz((nn-1)*3+3))) .AND.  &
					 ! (Ysub .GE. MIN(Py((nn-1)*3+1),Py((nn-1)*3+2),Py((nn-1)*3+3)) .AND. Ysub .LE. MAX(Py((nn-1)*3+1),Py((nn-1)*3+2),Py((nn-1)*3+3))) ) then
                        ! ! Intersection test using the line and point equation (cross product) (Position of a Point Relative to a Line or side of triangle)
                        ! c1 = (Py((nn-1)*3+2) - Py((nn-1)*3+1))*( Pz((nn-1)*3+1) - Zsub ) - (Pz((nn-1)*3+2) - Pz((nn-1)*3+1))*( Py((nn-1)*3+1) - Ysub )
                        ! c2 = (Py((nn-1)*3+3) - Py((nn-1)*3+2))*( Pz((nn-1)*3+2) - Zsub ) - (Pz((nn-1)*3+3) - Pz((nn-1)*3+2))*( Py((nn-1)*3+2) - Ysub )
                        ! c3 = (Py((nn-1)*3+1) - Py((nn-1)*3+3))*( Pz((nn-1)*3+3) - Zsub ) - (Pz((nn-1)*3+1) - Pz((nn-1)*3+3))*( Py((nn-1)*3+3) - Ysub )
                        ! ! Finding the ray or x-coordinate of the intersection point (ray_inter)
                        ! if (( c1 > 0 .AND. c2 > 0 .AND. c3 > 0) .OR. (c1 < 0 .AND. c2 < 0 .AND. c3 < 0 )) then
! !                                !$OMP CRITICAL
                                ! ray_inter(nc)=( FN(1,nn)* (Px((nn-1)*3+1) + Px((nn-1)*3+2) + Px((nn-1)*3+3) )/3. + &
												! FN(2,nn)*((Py((nn-1)*3+1) + Py((nn-1)*3+2) + Py((nn-1)*3+3) )/3.-Ysub ) + &
												! FN(3,nn)*((Pz((nn-1)*3+1) + Pz((nn-1)*3+2) + Pz((nn-1)*3+3) )/3.-Zsub ) ) /FN(1,nn)
                                ! nc=nc+1
! !                                !$OMP END CRITICAL
                        ! endif
                ! endif
        ! enddo
! !        !$OMP END PARALLEL DO

		! !!$acc loop independent private(bubj,bubble,ray_inter)
		! !$acc loop seq
        ! do bubi=1,bub-1
		! !$acc loop seq
		! do bubj=bubi+1,bub
            ! if ( ray_inter(bubi) > ray_inter(bubj) )then
                ! bubble=ray_inter(bubi)
                ! ray_inter(bubi)=ray_inter(bubj)
                ! ray_inter(bubj)=bubble
            ! endif
        ! enddo; enddo
		

        ! nc=1
        ! ETAs=0.0
		! !$acc loop seq
        ! do i=iBgnVOS,iEndVOS
            ! if ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) .AND. abs(Xs(i)-ray_inter(nc+1)) < (iDx(i)/2.) ) then
                ! ETA_sub(i)=ETA_sub(i) + ABS( ( ray_inter(nc+1)-ray_inter(nc) )/iDx(i) -ETAs) 
                ! nc=nc+2
            ! elseif ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) ) then
                ! ETA_sub(i)=ETA_sub(i) + ABS( ( ray_inter(nc)-Xs(i) )/iDx(i) +ETAs -0.5 ) 
                ! ETAs=ABS(ETAs-1.0)
                ! nc=nc+1
            ! else 
                ! ETA_sub(i)=ETA_sub(i)+ETAs
            ! endif		
        ! enddo
		
    ! enddo ; enddo

        ! !$OMP PARALLEL DO
		! !$acc loop independent
        ! do i=iBgnVOS,iEndVOS
            ! ETA(i,j,k)=ETA_sub(i)/nSubGrids_3d/nSubGrids_3d
        ! enddo
        ! !$OMP END PARALLEL DO
    ! enddo

    ! enddo
	! !$acc end parallel
    ! !$OMP END PARALLEL DO



!------------------------------------  VOS using Ray0asting & Subgrid -----------------------------------------

ETA_sub = 0.

    !$OMP PARALLEL DO PRIVATE(nsy, nsz, n, nn, nc, bubi, bubj, bubble, i, j, ray_inter, Zsub, Ysub, c1, c2, c3, ETAs) collapse(2)
	!$acc parallel vector_length(32)
	!$acc loop independent private(ray_inter) collapse(2) gang vector
    do k=istart_vos,iend_vos
    do j=jBgnVOS,jEndVOS
	!$acc loop independent private(ray_inter) collapse(2)
	!!$acc loop seq
    do nsz=1,nSubGrids_3d
	!!$acc loop seq
	do nsy=1,nSubGrids_3d         ! Each point (k,j) in the bounding box is divided into subgrids
        ray_inter=Xend+5.                         ! ray starts from Xend + 5.0  ( from spanwise direction)
        nc=1

        Zsub=Z(k)+(nsz-0.5)*iDz(k)/nSubGrids_3d
        Ysub=Y(j)+(nsy-0.5)*iDy(j)/nSubGrids_3d
        
!        !$OMP PARALLEL DO PRIVATE(c1,c2,c3) SHARED(nc,ray_inter)
		!$acc loop independent private(ray_inter)
		!!$acc loop seq
        do n=1,nk(j,k)                                                     ! Only facets with interiors intersected by the ray are included.			 
				nn = n_(j,k,n)
                if ( (Zsub .GE. MIN(Pz((nn-1)*3+1),Pz((nn-1)*3+2),Pz((nn-1)*3+3)) .AND. Zsub .LE. MAX(Pz((nn-1)*3+1),Pz((nn-1)*3+2),Pz((nn-1)*3+3))) .AND.  &
					 (Ysub .GE. MIN(Py((nn-1)*3+1),Py((nn-1)*3+2),Py((nn-1)*3+3)) .AND. Ysub .LE. MAX(Py((nn-1)*3+1),Py((nn-1)*3+2),Py((nn-1)*3+3))) ) then
                        ! Intersection test using the line and point equation (cross product) (Position of a Point Relative to a Line or side of triangle)
                        c1 = (Py((nn-1)*3+2) - Py((nn-1)*3+1))*( Pz((nn-1)*3+1) - Zsub ) - (Pz((nn-1)*3+2) - Pz((nn-1)*3+1))*( Py((nn-1)*3+1) - Ysub )
                        c2 = (Py((nn-1)*3+3) - Py((nn-1)*3+2))*( Pz((nn-1)*3+2) - Zsub ) - (Pz((nn-1)*3+3) - Pz((nn-1)*3+2))*( Py((nn-1)*3+2) - Ysub )
                        c3 = (Py((nn-1)*3+1) - Py((nn-1)*3+3))*( Pz((nn-1)*3+3) - Zsub ) - (Pz((nn-1)*3+1) - Pz((nn-1)*3+3))*( Py((nn-1)*3+3) - Ysub )
                        ! Finding the ray or x-coordinate of the intersection point (ray_inter)
                        if (( c1 > 0 .AND. c2 > 0 .AND. c3 > 0) .OR. (c1 < 0 .AND. c2 < 0 .AND. c3 < 0 )) then
!                                !$OMP CRITICAL
                                ray_inter(nc)=( FN(1,nn)* (Px((nn-1)*3+1) + Px((nn-1)*3+2) + Px((nn-1)*3+3) )/3. + &
												FN(2,nn)*((Py((nn-1)*3+1) + Py((nn-1)*3+2) + Py((nn-1)*3+3) )/3.-Ysub ) + &
												FN(3,nn)*((Pz((nn-1)*3+1) + Pz((nn-1)*3+2) + Pz((nn-1)*3+3) )/3.-Zsub ) ) /FN(1,nn)
                                nc=nc+1
!                                !$OMP END CRITICAL
                        endif
                endif
        enddo
!        !$OMP END PARALLEL DO

		!$acc loop independent private(ray_inter)
		!!$acc loop seq
        do bubi=1,bub-1
		!$acc loop seq
		do bubj=bubi+1,bub
            if ( ray_inter(bubi) > ray_inter(bubj) )then
                bubble=ray_inter(bubi)
                ray_inter(bubi)=ray_inter(bubj)
                ray_inter(bubj)=bubble
            endif
        enddo; enddo
		

        nc=1
        ETAs=0.0
		!$acc loop seq
        do i=iBgnVOS,iEndVOS
            if ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) .AND. abs(Xs(i)-ray_inter(nc+1)) < (iDx(i)/2.) ) then
                ETA_sub(i,j,k)=ETA_sub(i,j,k) + ABS( ( ray_inter(nc+1)-ray_inter(nc) )/iDx(i) -ETAs) 
                nc=nc+2
            elseif ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) ) then
                ETA_sub(i,j,k)=ETA_sub(i,j,k) + ABS( ( ray_inter(nc)-Xs(i) )/iDx(i) +ETAs -0.5 ) 
                ETAs=ABS(ETAs-1.)
                nc=nc+1
            else 
                ETA_sub(i,j,k)=ETA_sub(i,j,k)+ETAs
            endif		
        enddo
		
    enddo ; enddo

        !$OMP PARALLEL DO
		!$acc loop independent
        do i=iBgnVOS,iEndVOS
            ETA(i,j,k)=ETA_sub(i,j,k)/(nSubGrids_3d*nSubGrids_3d)
        enddo
        !$OMP END PARALLEL DO
    enddo

    enddo
	!$acc end parallel
    !$OMP END PARALLEL DO


      !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
	  !$acc update self(ETA(:,:,istart_vos:iend_vos)) if(nproc>1)
      icount = igcount_vos*(nx+4)*(ny+4)
      !Send my results back to the master
      if(myid>master)then
         itag = 401
         call MPI_SEND( ETA(-1,-1,istart_vos), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      !Wait to receive results from each task
      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount_vos(i)*(nx+4)*(ny+4)
            itag = 401
            call MPI_RECV( ETA(-1,-1,gstart_vos(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )    
         end do
      end if
	  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


	!>>>>>>>>> data transformation from master to all nodes <<<<<<<<<<<<<<<<
	icount= (nz+4)*(nx+4)*(ny+4)
	call MPI_BCAST ( ETA, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	!$acc update device(ETA) if(nproc>1)
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!$acc end data
!$acc end data


    if(myid==master)then
        open (62,file='geometry_info.dat',position='append')
        write(62,*)'Geometry_file = ',stl_file
        write(62,*)'                   '
        write(62,*)'Total number of facets = ',nf
        write(62,*)'                   '
        write(62,*)'z length of VOS = ',Zend-Zbgn
        write(62,*)'y length of VOS = ',Yend-Ybgn
        write(62,*)'                   '
    endif


    if(myid==master)then
        write(*,*)'*******Ray casting SUCCESSFUL*******'
    endif

end subroutine vos_ray3d



! subroutine vos_ray3d()			! This is a previous algorithm with too many variables uploaded into the GPU (maybe faster for CPU)
! use variables
! use mpi
! implicit none
! integer ,parameter                  :: bub=60
! real*4                              :: side_a,side_b,side_c,side_s,range_Z,range_Y,facet_area,facet_avg
! real*4                              :: x0_stl, y0_stl, z0_stl
! real*8  ,dimension(1:10000)         :: ETA_sub
! real*4                              :: Xbgn, Xend, Ybgn, Yend, Zbgn, Zend, tun, Xsub, Ysub, Zsub
! real*4                              :: c1, c2, c3, bubble
! real*4  ,dimension(1:bub)           :: ray_inter
! integer                             :: n, nc, bubi, bubj, nsy, nsz, off, nn
! integer ,dimension(:,:), allocatable:: nk
! integer ,dimension(:,:,:), allocatable :: n_


    ! !----------------------------MPI----------------------------!

    ! integer                      :: Zdv_vos, Zr_vos

    ! integer, dimension(1:1024)   :: gstart_vos ,gend_vos, gen_vos, gcount_vos

    ! integer                      :: iend_vos, istart_vos, igcount_vos

    ! !----------------------------MPI----------------------------!


    ! Xbgn = MINVAL(Px)                                           !  max & min limit of vertex in x, y and z
    ! Xend = MAXVAL(Px)
    ! x0_stl=0.5*( Xbgn+Xend )
    ! Ybgn = MINVAL(Py)
    ! Yend = MAXVAL(Py)
    ! y0_stl=0.5*( Ybgn+Yend )
    ! Zbgn = MINVAL(Pz)
    ! Zend = MAXVAL(Pz)
    ! z0_stl=0.5*( Zbgn+Zend )
    ! tun = Zend - Zbgn
	


! !$acc data present(Zs,Ys,Xs,iDz,iDy,iDx,ETA) copyin(Px,Py,Pz,M0,FN,Vy,Vz) create(ray_inter,ETA_sub)


    ! !$OMP PARALLEL DO
	! !$acc parallel loop independent gang vector
    ! do n=1,nf*3
        ! Px(n) = ( Px(n)-x0_stl )/tun + x0                       ! to scale the solid dimensions 
        ! Py(n) = ( Py(n)-y0_stl )/tun + y0
        ! Pz(n) = ( Pz(n)-z0_stl )/tun + z0
    ! enddo
	! !$acc end parallel
    ! !$OMP END PARALLEL DO


    ! !---------------------------------Calculate the average triangular:cell area ratio------------------------
    ! ! facet_area=0.
    ! ! !$OMP PARALLEL DO PRIVATE(side_a,side_b,side_c,side_s) REDUCTION(+:facet_area)
	! ! !$acc parallel loop independent private(side_a,side_b,side_c,side_s) reduction(max : facet_area)
	
    ! ! do n=1,nf
        ! ! side_a=SQRT((Px((n-1)*3+1)-Px((n-1)*3+2))**2+(Py((n-1)*3+1)-Py((n-1)*3+2))**2+(Pz((n-1)*3+1)-Pz((n-1)*3+2))**2)
        ! ! side_b=SQRT((Px((n-1)*3+1)-Px((n-1)*3+3))**2+(Py((n-1)*3+1)-Py((n-1)*3+3))**2+(Pz((n-1)*3+1)-Pz((n-1)*3+3))**2)
        ! ! side_c=SQRT((Px((n-1)*3+2)-Px((n-1)*3+3))**2+(Py((n-1)*3+2)-Py((n-1)*3+3))**2+(Pz((n-1)*3+2)-Pz((n-1)*3+3))**2)
        ! ! side_s=(side_a+side_b+side_c)/2.
        ! ! facet_area=facet_area+SQRT(side_s*(side_s-side_a)*(side_s-side_b)*(side_s-side_c))
    ! ! enddo
	! ! !$acc end parallel
    ! ! !$OMP END PARALLEL DO

        ! ! facet_avg=facet_area/(nf*1.)


! ! updated after scale the solid dimensions 
	! !$acc update self(Px,Py,Pz)
    ! Xbgn = MINVAL(Px)
    ! Xend = MAXVAL(Px)
    ! Ybgn = MINVAL(Py)
    ! Yend = MAXVAL(Py)
    ! Zbgn = MINVAL(Pz)
    ! Zend = MAXVAL(Pz)

! !--------- Get  Border of solid body in the computational domain using mesh points & Define area applied VOS ..........
! !................ Bounding BOX  ......................

    ! off=0
    ! do i= 1,nx
        ! if ( Xs(i)>Xbgn .AND. off==0 )then
            ! iBgnVOS=i-1
            ! off=1
        ! elseif ( Xs(i)>Xend .AND. off==1 )then
            ! iEndVOS=i
            ! off=2
        ! endif
    ! enddo
    ! off=0
    ! do j= 1,ny
        ! if ( Ys(j)>Ybgn .AND. off==0 )then
            ! jBgnVOS=j-1
            ! off=1
        ! elseif ( Ys(j)>Yend .AND. off==1 )then
            ! jEndVOS=j
            ! off=2
        ! endif
    ! enddo
    ! off=0
    ! do k= 0,nz
        ! if ( Zs(k)>Zbgn .AND. off==0 )then
            ! kBgnVOS=k-1
            ! off=1
        ! elseif ( Zs(k)>Zend .AND. off==1 )then
            ! kEndVOS=k
            ! off=2
        ! endif
    ! enddo
 



    ! !-----------------------MPI DIVISION-------------------------!
    ! Zdv_vos = (kEndVOS-kBgnVOS+1) / nproc
    ! Zr_vos  = (kEndVOS-kBgnVOS+1) - Zdv_vos * nproc 
    ! !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! !i = myid
    ! do i=0,(nproc-1)

        ! if(i < Zr_vos) then
            ! gstart_vos(i) = kBgnVOS + i * (Zdv_vos+1)
            ! gen_vos(i) = gstart_vos(i) + Zdv_vos
        ! else
            ! gstart_vos(i) = kBgnVOS + i * Zdv_vos + Zr_vos
            ! gen_vos(i) = gstart_vos(i) + Zdv_vos - 1
        ! end if
        
        ! gcount_vos(i) = gen_vos(i) - gstart_vos(i) + 1
        ! gend_vos(i) = gcount_vos(i) + 2

    ! end do

    ! !----------for nz vos----------!
    ! istart_vos = gstart_vos(myid)  !
    ! iend_vos = gen_vos(myid)     !
    ! igcount_vos = gcount_vos(myid) !
    ! !----------for nz vos----------!

    ! !-----------------------MPI DIVISION-------------------------!

	! allocate(nk(jbgnvos:jendvos,kbgnvos:kendvos) )
	! allocate(n_(jbgnvos:jendvos,kbgnvos:kendvos,1:nf) )


! !$acc data create(nk(jbgnvos:jendvos,istart_vos:iend_vos),n_(jbgnvos:jendvos,istart_vos:iend_vos,1:nf))


    ! !---------------------  facet normal position and side length using Vertex value ---------------------------

    ! !$OMP PARALLEL DO
	! !$acc parallel loop independent gang vector
    ! do n=1,nf
        ! M0(1,n)=( Px((n-1)*3+1) + Px((n-1)*3+2) + Px((n-1)*3+3) )/3.   !  centroid facet normal position
        ! M0(2,n)=( Py((n-1)*3+1) + Py((n-1)*3+2) + Py((n-1)*3+3) )/3.
        ! M0(3,n)=( Pz((n-1)*3+1) + Pz((n-1)*3+2) + Pz((n-1)*3+3) )/3.
        ! Vy((n-1)*3+1)=Py((n-1)*3+2) - Py((n-1)*3+1)					! coordinates a side (line)  of triangle in x,y, z 
        ! Vy((n-1)*3+2)=Py((n-1)*3+3) - Py((n-1)*3+2)
        ! Vy((n-1)*3+3)=Py((n-1)*3+1) - Py((n-1)*3+3)
        ! Vz((n-1)*3+1)=Pz((n-1)*3+2) - Pz((n-1)*3+1)
        ! Vz((n-1)*3+2)=Pz((n-1)*3+3) - Pz((n-1)*3+2)
        ! Vz((n-1)*3+3)=Pz((n-1)*3+1) - Pz((n-1)*3+3)
    ! enddo
	! !$acc end parallel
    ! !$OMP END PARALLEL DO

    ! !$OMP PARALLEL DO PRIVATE(j,nn,n)
	! !$acc parallel 
	! !$acc loop independent private(nn,n) collapse(2) gang vector
	! do k=istart_vos,iend_vos
    ! do j=jBgnVOS,jEndVOS
		! nn = 0
		! !$acc loop independent reduction(+:nn)
		! !!$acc loop seq
		! do n=1,nf
			! if( ( MAX(Pz((n-1)*3+1),Pz((n-1)*3+2),Pz((n-1)*3+3)) .GE. Zs(k)-iDz(k)/2.) .AND. (MIN(Pz((n-1)*3+1),Pz((n-1)*3+2),Pz((n-1)*3+3)) .LE. Zs(k)+iDz(k)/2.) .AND. &
				! ( MAX(Py((n-1)*3+1),Py((n-1)*3+2),Py((n-1)*3+3)) .GE. Ys(j)-iDy(j)/2.) .AND. (MIN(Py((n-1)*3+1),Py((n-1)*3+2),Py((n-1)*3+3)) .LE. Ys(j)+iDy(j)/2.) ) then
				! nn = nn+1
				! n_(j,k,nn) = n
			! endif
		! enddo
		! nk(j,k) = nn
	! enddo; enddo
	! !$acc end parallel
    ! !$OMP END PARALLEL DO	


    ! if(myid==master)then
        ! open (62,file='geometry_info.dat',position='append')
        ! write(62,*)'Geometry_file = ',stl_file
        ! write(62,*)'                   '
        ! write(62,*)'Total number of facets = ',nf
        ! write(62,*)'Facet average area = ',facet_avg
        ! write(62,*)'                   '
        ! write(62,*)'z length of VOS = ',Zend-Zbgn
        ! write(62,*)'y length of VOS = ',Yend-Ybgn
        ! write(62,*)'                   '
    ! endif


    ! if(myid==master)then
        ! write(*,*) 'Ray0asting progress (%)'
    ! endif

! !------------------------------------  VOS using Ray0asting & Subgrid -----------------------------------------

    ! !$OMP PARALLEL DO PRIVATE(nsy, nsz, n, nc, bubi, bubj, i, j, ETA_sub, ray_inter, Zsub, Ysub, c1, c2, c3, ETAs)
	! !$acc parallel
	! !$acc loop independent private(nsy, nsz, n, nc, bubi, bubj, i, ETA_sub, ray_inter, Zsub, Ysub, c1, c2, c3, ETAs) tile(4,4) gang vector
    ! do k=istart_vos,iend_vos
    ! do j=jBgnVOS,jEndVOS
    ! ETA_sub=0.0
	! !$acc loop independent private(nc, bubi, bubj, i, ray_inter, Zsub, Ysub, c1, c2, c3, ETAs) reduction(+:ETA_sub)
	! !!$acc loop seq
    ! do nsz=1,nSubGrids_3d
	! !!$acc loop seq
	! do nsy=1,nSubGrids_3d         ! Each point (k,j) in the bounding box is divided into subgrids
        ! ray_inter=Xend+5.0                         ! ray starts from Xend + 5.0  ( from spanwise direction)
        ! nc=1

        ! Zsub=Zs(k)-iDz(k)/2.+(nsz-0.5)*iDz(k)/nSubGrids_3d
        ! Ysub=Ys(j)-iDy(j)/2.+(nsy-0.5)*iDy(j)/nSubGrids_3d
        
! !        !$OMP PARALLEL DO PRIVATE(c1,c2,c3) SHARED(nc,ray_inter)
		! !$acc loop seq
        ! do n=1,nk(j,k)                                                     ! Only facets with interiors intersected by the ray are included.			 
                ! if ( (Zsub .GE. MIN(Pz((n_(j,k,n)-1)*3+1),Pz((n_(j,k,n)-1)*3+2),Pz((n_(j,k,n)-1)*3+3)) .AND. Zsub .LE. MAX(Pz((n_(j,k,n)-1)*3+1),Pz((n_(j,k,n)-1)*3+2),Pz((n_(j,k,n)-1)*3+3))) .AND.  &
					 ! (Ysub .GE. MIN(Py((n_(j,k,n)-1)*3+1),Py((n_(j,k,n)-1)*3+2),Py((n_(j,k,n)-1)*3+3)) .AND. Ysub .LE. MAX(Py((n_(j,k,n)-1)*3+1),Py((n_(j,k,n)-1)*3+2),Py((n_(j,k,n)-1)*3+3))) ) then
                        ! ! Intersection test using the line and point equation (cross product) (Position of a Point Relative to a Line or side of triangle)
                        ! c1 = Vy( (n_(j,k,n)-1)*3+1 )*( Pz((n_(j,k,n)-1)*3+1) - Zsub ) - Vz( (n_(j,k,n)-1)*3+1 )*( Py((n_(j,k,n)-1)*3+1) - Ysub )
                        ! c2 = Vy( (n_(j,k,n)-1)*3+2 )*( Pz((n_(j,k,n)-1)*3+2) - Zsub ) - Vz( (n_(j,k,n)-1)*3+2 )*( Py((n_(j,k,n)-1)*3+2) - Ysub )
                        ! c3 = Vy( (n_(j,k,n)-1)*3+3 )*( Pz((n_(j,k,n)-1)*3+3) - Zsub ) - Vz( (n_(j,k,n)-1)*3+3 )*( Py((n_(j,k,n)-1)*3+3) - Ysub )
                        ! ! Finding the ray or x-coordinate of the intersection point (ray_inter)
                        ! if (( c1 > 0 .AND. c2 > 0 .AND. c3 > 0) .OR. (c1 < 0 .AND. c2 < 0 .AND. c3 < 0 )) then
! !                                !$OMP CRITICAL
                                ! ray_inter(nc)=( FN(1,n_(j,k,n))*M0(1,n_(j,k,n)) + FN(2,n_(j,k,n))*( M0(2,n_(j,k,n))-Ysub ) + FN(3,n_(j,k,n))*( M0(3,n_(j,k,n))-Zsub ) ) /FN(1,n_(j,k,n))
                                ! nc=nc+1
! !                                !$OMP END CRITICAL
                        ! endif
                ! endif
        ! enddo
! !        !$OMP END PARALLEL DO

		! !$acc loop independent private(bubi, bubj, ray_inter) 
		! !!$acc loop seq
        ! do bubi=1,bub-1; do bubj=bubi+1,bub
            ! if ( ray_inter(bubi) > ray_inter(bubj) )then
                ! bubble=ray_inter(bubi)
                ! ray_inter(bubi)=ray_inter(bubj)
                ! ray_inter(bubj)=bubble
            ! endif
        ! enddo; enddo
		

        ! nc=1
        ! ETAs=0.0
		! !$acc loop seq
        ! do i=iBgnVOS,iEndVOS
            ! if ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) .AND. abs(Xs(i)-ray_inter(nc+1)) < (iDx(i)/2.) ) then
                ! ETA_sub(i)=ETA_sub(i) + ABS( ( ray_inter(nc+1)-ray_inter(nc) )/iDx(i) -ETAs) 
                ! nc=nc+2
            ! elseif ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) ) then
                ! ETA_sub(i)=ETA_sub(i) + ABS( ( ray_inter(nc)-Xs(i) )/iDx(i) +ETAs -0.5 ) 
                ! ETAs=ABS(ETAs-1.0)
                ! nc=nc+1
            ! else 
                ! ETA_sub(i)=ETA_sub(i)+ETAs
            ! endif		
        ! enddo
		
    ! enddo ; enddo

        ! !$OMP PARALLEL DO
		! !$acc loop independent
        ! do i=iBgnVOS,iEndVOS
            ! ETA(i,j,k)=ETA_sub(i)/nSubGrids_3d/nSubGrids_3d
        ! enddo
        ! !$OMP END PARALLEL DO
    ! enddo

    ! ! if(myid==master)then
		! ! write(*,'(F6.2)') (k-istart_vos)*100./(iend_vos-istart_vos)*1.
    ! ! endif
    ! enddo
	! !$acc end parallel
    ! !$OMP END PARALLEL DO


      ! !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
	  ! !$acc update self(ETA(:,:,istart_vos:iend_vos)) if(nproc>1)
      ! icount = igcount_vos*(nx+4)*(ny+4)
      ! !Send my results back to the master
      ! if(myid>master)then
         ! itag = 401
         ! call MPI_SEND( ETA(-1,-1,istart_vos), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      ! end if
! !      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      ! !Wait to receive results from each task
      ! if(myid==master)then
         ! do i = 1, (nproc-1)
            ! icount = gcount_vos(i)*(nx+4)*(ny+4)
            ! itag = 401
            ! call MPI_RECV( ETA(-1,-1,gstart_vos(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )    
         ! end do
      ! end if
	  ! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


	! !>>>>>>>>> data transformation from master to all nodes <<<<<<<<<<<<<<<<
	! icount= (nz+4)*(nx+4)*(ny+4)
	! call MPI_BCAST ( ETA, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
	! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	! !$acc update device(ETA) if(nproc>1)
	! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! !$acc end data
! !$acc end data

    ! if(myid==master)then
        ! write(*,*)'*******Ray casting SUCCESSFUL*******'
    ! endif

! end subroutine vos_ray3d






! subroutine vos_ray3d()				! For CPU-MPI
! use variables
! implicit none
! integer ,parameter                  :: bub=60
! real*4                              :: side_a,side_b,side_c,side_s,range_Z,range_Y,facet_area,facet_avg
! real*4,dimension(1:nf*3)			:: Px,Py,Pz
! real*4,dimension(1:3,1:nf)			:: FN
! real*8  ,dimension(1:10000)         :: ETA_sub
! real*4                              :: Xbgn, Xend, Ybgn, Yend, Zbgn, Zend, Xsub, Ysub, Zsub
! real*4                              :: c1, c2, c3, bubble
! real*4  ,dimension(1:bub)           :: ray_inter
! integer                             :: n, nc, bubi, bubj, nsy, nsz, off, nn



! !	Enter solid transformation here  

	! !$OMP PARALLEL DO 
    ! do n=1,nf
		! !!! Update FN
		! !!! a function of f(iPz,iPy,iPx,iFN,AOA,x0_t,y0_t,z0_t)
		
		! FN(1,n) = iFN(1,n)
		! FN(2,n) = iFN(2,n)
		! FN(3,n) = iFN(3,n)
		
    ! end do
    ! !$OMP END PARALLEL DO
	
	! !$OMP PARALLEL DO 
    ! do n=1,nf*3
		! !!! Update Px,Py,Pz
		! !!! a function of f(iPz,iPy,iPx,iFN,AOA,x0_t,y0_t,z0_t)
		
		! Px(n) = iPx(n) + x0_t
		! Py(n) = iPy(n) + y0_t
		! Pz(n) = iPz(n) + z0_t
		
    ! end do
    ! !$OMP END PARALLEL DO




! ! updating bounds after the solid transformation 

		! Xbgn = MINVAL(Px)
		! Xend = MAXVAL(Px)
		! Ybgn = MINVAL(Py)
		! Yend = MAXVAL(Py)
		! Zbgn = MINVAL(Pz)
		! Zend = MAXVAL(Pz)




! !--------- Get  Border of solid body in the computational domain using mesh points & Define area applied VOS ..........
! !................ Bounding BOX  ......................

    ! off=0
    ! do i= 1,nx
        ! if ( Xs(i)>Xbgn .AND. off==0 )then
            ! iBgnVOS=i-1
            ! off=1
        ! elseif ( Xs(i)>Xend .AND. off==1 )then
            ! iEndVOS=i
            ! exit
        ! endif
    ! enddo
    ! off=0
    ! do j= 1,ny
        ! if ( Ys(j)>Ybgn .AND. off==0 )then
            ! jBgnVOS=j-1
            ! off=1
        ! elseif ( Ys(j)>Yend .AND. off==1 )then
            ! jEndVOS=j
            ! exit
        ! endif
    ! enddo
    ! off=0
    ! do k= 0,nz
        ! if ( Zs(k)>Zbgn .AND. off==0 )then
            ! kBgnVOS=k-1
            ! off=1
        ! elseif ( Zs(k)>Zend .AND. off==1 )then
            ! kEndVOS=k
            ! exit
        ! endif
    ! enddo
 



    ! !-----------------------MPI DIVISION-------------------------!
    ! Zdv_vos = (kEndVOS-kBgnVOS+1) / nproc
    ! Zr_vos  = (kEndVOS-kBgnVOS+1) - Zdv_vos * nproc 
    ! !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! !i = myid
    ! do i=0,(nproc-1)

        ! if(i < Zr_vos) then
            ! gstart_vos(i) = kBgnVOS + i * (Zdv_vos+1)
            ! gend0_vos(i) = gstart_vos(i) + Zdv_vos
        ! else
            ! gstart_vos(i) = kBgnVOS + i * Zdv_vos + Zr_vos
            ! gend0_vos(i) = gstart_vos(i) + Zdv_vos - 1
        ! end if
        
        ! gcount_vos(i) = gend0_vos(i) - gstart_vos(i) + 1
        ! gend_vos(i) = gcount_vos(i) + 2

    ! end do

    ! !----------for nz vos----------!
    ! istart_vos = gstart_vos(myid)  !
    ! iend_vos = gend0_vos(myid)     !
    ! igcount_vos = gcount_vos(myid) !
    ! !----------for nz vos----------!

    ! !-----------------------MPI DIVISION-------------------------!


! !------------------------------------  VOS using RayCasting & Subgrid -----------------------------------------

    ! do k=istart_vos,iend_vos
    ! !$OMP PARALLEL DO PRIVATE(nsy, nsz, nc, bubi, bubj, i, ETA_sub, ray_inter, Zsub, Ysub, c1, c2, c3, ETAs)  
    ! do j=jBgnVOS,jEndVOS
    ! ETA_sub=0.0
    ! do nsz=1,nSubGrids_3d; do nsy=1,nSubGrids_3d         ! Each point (k,j) in the bounding box is divided into subgrids
        ! ray_inter=Xend+5.0                         ! ray starts from Xend + 5.0  ( from spanwise direction)
        ! nc=1

        ! Zsub=Z(k)+(nsz-0.5)*iDz(k)/nSubGrids_3d
        ! Ysub=Y(j)+(nsy-0.5)*iDy(j)/nSubGrids_3d
        
! !        !$OMP PARALLEL DO PRIVATE(c1,c2,c3) SHARED(nc,ray_inter)	
        ! do n=1,nf                                                     ! Only facets with centroid position within the possible range in z & y are included
                ! if ( (Zsub .GE. MIN(Pz((n-1)*3+1),Pz((n-1)*3+2),Pz((n-1)*3+3)) .AND. Zsub .LE. MAX(Pz((n-1)*3+1),Pz((n-1)*3+2),Pz((n-1)*3+3))) .AND.  &
					 ! (Ysub .GE. MIN(Py((n-1)*3+1),Py((n-1)*3+2),Py((n-1)*3+3)) .AND. Ysub .LE. MAX(Py((n-1)*3+1),Py((n-1)*3+2),Py((n-1)*3+3))) ) then
                        ! ! Intersection test using the line and point equation (cross product) (Position of a Point Relative to a Line or side of triangle)
                        ! c1 = (Py((n-1)*3+2) - Py((n-1)*3+1))*( Pz((n-1)*3+1) - Zsub ) - (Pz((n-1)*3+2) - Pz((n-1)*3+1))*( Py((n-1)*3+1) - Ysub )
                        ! c2 = (Py((n-1)*3+3) - Py((n-1)*3+2))*( Pz((n-1)*3+2) - Zsub ) - (Pz((n-1)*3+3) - Pz((n-1)*3+2))*( Py((n-1)*3+2) - Ysub )
                        ! c3 = (Py((n-1)*3+1) - Py((n-1)*3+3))*( Pz((n-1)*3+3) - Zsub ) - (Pz((n-1)*3+1) - Pz((n-1)*3+3))*( Py((n-1)*3+3) - Ysub )
                        ! ! Finding the ray or x-coordinate of the intersection point (ray_inter)
                        ! if (( c1 > 0 .AND. c2 > 0 .AND. c3 > 0) .OR. (c1 < 0 .AND. c2 < 0 .AND. c3 < 0 )) then
! !                                !$OMP CRITICAL
                                ! ray_inter(nc)=( FN(1,n)* (Px((n-1)*3+1) + Px((n-1)*3+2) + Px((n-1)*3+3) )/3. + &
												! FN(2,n)*((Py((n-1)*3+1) + Py((n-1)*3+2) + Py((n-1)*3+3) )/3.-Ysub ) + &
												! FN(3,n)*((Pz((n-1)*3+1) + Pz((n-1)*3+2) + Pz((n-1)*3+3) )/3.-Zsub ) ) /FN(1,n)
                                ! nc=nc+1
! !                                !$OMP END CRITICAL
                        ! endif
                ! endif
        ! enddo
! !        !$OMP END PARALLEL DO


        ! do bubi=1,bub-1
		! do bubj=bubi+1,bub
            ! if ( ray_inter(bubi) > ray_inter(bubj) )then
                ! bubble=ray_inter(bubi)
                ! ray_inter(bubi)=ray_inter(bubj)
                ! ray_inter(bubj)=bubble
            ! endif
        ! enddo; enddo
		

        ! nc=1
        ! ETAs=0.0
        ! do i=iBgnVOS,iEndVOS
            ! if ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) .AND. abs(Xs(i)-ray_inter(nc+1)) < (iDx(i)/2.) ) then
                ! ETA_sub(i)=ETA_sub(i) + ABS( ( ray_inter(nc+1)-ray_inter(nc) )/iDx(i) -ETAs) 
                ! nc=nc+2
            ! elseif ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) ) then
                ! ETA_sub(i)=ETA_sub(i) + ABS( ( ray_inter(nc)-Xs(i) )/iDx(i) +ETAs -0.5 ) 
                ! ETAs=ABS(ETAs-1.0)
                ! nc=nc+1
            ! else 
                ! ETA_sub(i)=ETA_sub(i)+ETAs
            ! endif
        ! enddo 
    ! enddo ; enddo


        ! do i=iBgnVOS,iEndVOS
            ! ETA(i,j,k)=ETA_sub(i,j,k)/(nSubGrids_3d*nSubGrids_3d)
        ! enddo

    ! enddo
    ! !$OMP END PARALLEL DO
    ! enddo


      ! !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
  
      ! icount = igcount_vos*(nx+4)*(ny+4)
      ! !Send my results back to the master
      ! if(myid>master)then
         ! itag = 401
         ! call MPI_SEND( ETA(-1,-1,istart_vos), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      ! end if
! !      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      ! !Wait to receive results from each task
      ! if(myid==master)then
         ! do i = 1, (nproc-1)
            ! icount = gcount_vos(i)*(nx+4)*(ny+4)
            ! itag = 401
            ! call MPI_RECV( ETA(-1,-1,gstart_vos(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )    
         ! end do
      ! end if
	  ! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! !>>>>>>>>> data transformation from master to all nodes <<<<<<<<<<<<<<<<
	! icount= (nz+4)*(nx+4)*(ny+4)
	! call MPI_BCAST ( ETA, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
	! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	

    ! if(myid==master)then
        ! open (62,file='geometry_info.dat',position='append')
        ! write(62,*)'Geometry_file = ',stl_file
        ! write(62,*)'                   '
        ! write(62,*)'Total number of facets = ',nf
        ! write(62,*)'                   '
        ! write(62,*)'z length of VOS = ',Zend-Zbgn
        ! write(62,*)'y length of VOS = ',Yend-Ybgn
        ! write(62,*)'                   '
    ! endif


    ! if(myid==master)then
        ! write(*,*)'*******Ray casting SUCCESSFUL*******'
    ! endif

! end subroutine vos_ray3d