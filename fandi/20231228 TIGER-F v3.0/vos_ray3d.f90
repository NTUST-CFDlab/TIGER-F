! 28 Jun 2023 - FDS

subroutine vos_ray3d()
use variables
use mpi
implicit none
integer*4                           :: nf=99999999
integer ,parameter                  :: bub=60
real*8  ,allocatable                :: FN(:,:),M(:,:)
real*8  ,allocatable                :: Px(:), Py(:), Pz(:), Vx(:), Vy(:), Vz(:)
real*8                              :: side_a,side_b,side_c,side_s,range_Z,range_Y,facet_area,facet_avg
real*8                              :: Xc_stl, Yc_stl, Zc_stl
real*8  ,dimension(1:10000)         :: ETA_sub
real*8                              :: Xbng, Xend, Ybng, Yend, Zbng, Zend, tun, Xsub, Ysub, Zsub
real*8                              :: c1, c2, c3, bubble
real*8  ,dimension(1:bub)           :: ray_inter
integer                             :: n, nc, bubi, bubj, nsx, nsy, nsz, off
character(len=20)                   :: useless
integer                             :: iBgnVOS, iEndVOS, jBgnVOS, jEndVOS, kBgnVOS, kEndVOS
character*80                        :: header


    !----------------------------MPI----------------------------!

    integer                      :: Zdv_vos, Zr_vos

    integer, dimension(1:1024)   :: gstart_vos ,gend_vos, gend0_vos, gcount_vos

    integer                      :: iend_vos, istart_vos, igcount_vos

    !----------------------------MPI----------------------------!


    totalstarttime = MPI_WTIME()


    !-------------------------------------------------READ STL FILE------------------------------------------------!
    !-------------------------------------------------READ STL FILE------------------------------------------------!
    !-------------------------------------------------READ STL FILE------------------------------------------------!
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

    allocate(FN(3,nf))
    allocate(M(3,nf))
    allocate(Px(3*nf))                                          ! Size of arrays determined for facet normal & vertecies
    allocate(Py(3*nf))
    allocate(Pz(3*nf))
    allocate(Vx(3*nf))
    allocate(Vy(3*nf))
    allocate(Vz(3*nf))

    !--------------------------------------Adjust heading direction-----------------------------------

    open(9,file=stl_file)
    read(9,*)
    do n=1,nf
        read(9,*) useless, useless, FN(3,n), FN(2,n), FN(1,n)                  ! facet normal nx ny nz (or swapped axis)
        read(9,*)                                                              ! outer loop
        read(9,*) useless, Pz((n-1)*3+1), Py((n-1)*3+1), Px((n-1)*3+1)         ! vertex v1x v1y v1z (or swapped axis)
        read(9,*) useless, Pz((n-1)*3+2), Py((n-1)*3+2), Px((n-1)*3+2)         ! vertex v2x v2y v2z (or swapped axis)
        read(9,*) useless, Pz((n-1)*3+3), Py((n-1)*3+3), Px((n-1)*3+3)         ! vertex v3x v3y v3z (or swapped axis)
        read(9,*)                                                              ! endloop
        read(9,*)                                                              ! endfacet
    enddo
    close(9)

    !!$OMP PARALLEL DO
    !do n=1,nf
    !    FN(3,n)=-FN(3,n)					! to mirror the solid in z direction ( xy plane )
    !    Pz((n-1)*3+1)=-Pz((n-1)*3+1)
    !    Pz((n-1)*3+2)=-Pz((n-1)*3+2)
    !    Pz((n-1)*3+3)=-Pz((n-1)*3+3)
    !enddo
    !!$OMP END PARALLEL DO


    Xbng = MINVAL(Px)                                           !  max & min limit of vertex in x, y and z
    Xend = MAXVAL(Px)
    Xc_stl=0.5*( Xbng+Xend )
    Ybng = MINVAL(Py)
    Yend = MAXVAL(Py)
    Yc_stl=0.5*( Ybng+Yend )
    Zbng = MINVAL(Pz)
    Zend = MAXVAL(Pz)
    Zc_stl=0.5*( Zbng+Zend )
    tun = Zend - Zbng
	

    !$OMP PARALLEL DO
    do n=1,nf*3
        Px(n) = ( Px(n)-Xc_stl )/tun + xc                       ! to scale the solid dimensions 
        Py(n) = ( Py(n)-Yc_stl )/tun + yc
        Pz(n) = ( Pz(n)-Zc_stl )/tun + zc
    enddo
    !$OMP END PARALLEL DO


    !---------------------------------Calculate the average triangular:cell area ratio------------------------
    facet_area=0.d0
    !$OMP PARALLEL DO PRIVATE(side_a,side_b,side_c,side_s) REDUCTION(+:facet_area)
    do n=1,nf
        side_a=SQRT((Px((n-1)*3+1)-Px((n-1)*3+2))**2+(Py((n-1)*3+1)-Py((n-1)*3+2))**2+(Pz((n-1)*3+1)-Pz((n-1)*3+2))**2)
        side_b=SQRT((Px((n-1)*3+1)-Px((n-1)*3+3))**2+(Py((n-1)*3+1)-Py((n-1)*3+3))**2+(Pz((n-1)*3+1)-Pz((n-1)*3+3))**2)
        side_c=SQRT((Px((n-1)*3+2)-Px((n-1)*3+3))**2+(Py((n-1)*3+2)-Py((n-1)*3+3))**2+(Pz((n-1)*3+2)-Pz((n-1)*3+3))**2)
        side_s=(side_a+side_b+side_c)/2.d0
        facet_area=facet_area+SQRT(side_s*(side_s-side_a)*(side_s-side_b)*(side_s-side_c))
    enddo
    !$OMP END PARALLEL DO

        facet_avg=facet_area/(nf*1.d0)


! updated after scale the solid dimensions 

    Xbng = MINVAL(Px)
    Xend = MAXVAL(Px)
    Ybng = MINVAL(Py)
    Yend = MAXVAL(Py)
    Zbng = MINVAL(Pz)
    Zend = MAXVAL(Pz)

!--------- Get  Border of solid body in the computational domain using mesh points & Define area applied VOS ..........
!................ Bounding BOX  ......................

    off=0
    do i= 1,nx
        if ( Xs(i)>Xbng .AND. off==0 )then
            iBgnVOS=i-1
            off=1
        elseif ( Xs(i)>Xend .AND. off==1 )then
            iEndVOS=i
            off=2
        endif
    enddo
    off=0
    do j= 1,ny
        if ( Ys(j)>Ybng .AND. off==0 )then
            jBgnVOS=j-1
            off=1
        elseif ( Ys(j)>Yend .AND. off==1 )then
            jEndVOS=j
            off=2
        endif
    enddo
    off=0
    do k= 0,nz
        if ( Zs(k)>Zbng .AND. off==0 )then
            kBgnVOS=k-1
            off=1
        elseif ( Zs(k)>Zend .AND. off==1 )then
            kEndVOS=k
            off=2
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


    !---------------------  facet normal position and side length using Vertex value ---------------------------

    !$OMP PARALLEL DO
    do n=1,nf
        M(1,n)=( Px((n-1)*3+1) + Px((n-1)*3+2) + Px((n-1)*3+3) )/3.   !  centroid facet normal position
        M(2,n)=( Py((n-1)*3+1) + Py((n-1)*3+2) + Py((n-1)*3+3) )/3.
        M(3,n)=( Pz((n-1)*3+1) + Pz((n-1)*3+2) + Pz((n-1)*3+3) )/3.
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO
    do n=1,nf
        Vx((n-1)*3+1)=Px((n-1)*3+2) - Px((n-1)*3+1)                   ! coordinates a side (line)  of triangle in x,y, z 
        Vx((n-1)*3+2)=Px((n-1)*3+3) - Px((n-1)*3+2)
        Vx((n-1)*3+3)=Px((n-1)*3+1) - Px((n-1)*3+3)
        Vy((n-1)*3+1)=Py((n-1)*3+2) - Py((n-1)*3+1)
        Vy((n-1)*3+2)=Py((n-1)*3+3) - Py((n-1)*3+2)
        Vy((n-1)*3+3)=Py((n-1)*3+1) - Py((n-1)*3+3)
        Vz((n-1)*3+1)=Pz((n-1)*3+2) - Pz((n-1)*3+1)
        Vz((n-1)*3+2)=Pz((n-1)*3+3) - Pz((n-1)*3+2)
        Vz((n-1)*3+3)=Pz((n-1)*3+1) - Pz((n-1)*3+3)
    enddo
    !$OMP END PARALLEL DO

    !---------------------------------The longest centroid from vertex in z direction------------------------
    range_Z=0.d0
        do n=1,nf
                if     (ABS(Pz((n-1)*3+1)-M(3,n)) > range_Z) then
                        range_Z = ABS(Pz((n-1)*3+1)-M(3,n))
                elseif (ABS(Pz((n-1)*3+2)-M(3,n)) > range_Z) then
                        range_Z = ABS(Pz((n-1)*3+2)-M(3,n))
                elseif (ABS(Pz((n-1)*3+3)-M(3,n)) > range_Z) then
                        range_Z = ABS(Pz((n-1)*3+3)-M(3,n))
                endif
    enddo

    !---------------------------------The longest centroid from vertex in y direction------------------------
    range_Y=0.d0
        do n=1,nf
                if     (ABS(Py((n-1)*3+1)-M(2,n)) > range_Y) then
                        range_Y = ABS(Py((n-1)*3+1)-M(2,n))
                elseif (ABS(Py((n-1)*3+2)-M(2,n)) > range_Y) then
                        range_Y = ABS(Py((n-1)*3+2)-M(2,n))
                elseif (ABS(Py((n-1)*3+3)-M(2,n)) > range_Y) then
                        range_Y = ABS(Py((n-1)*3+3)-M(2,n))
                endif
    enddo

    if(myid==master)then
        open (62,file='geometry_info.dat',position='append')
        write(62,*)'Geometry_file = ',stl_file
        write(62,*)'                   '
        write(62,*)'Total number of facets = ',nf
        write(62,*)'Facet average area = ',facet_avg
        write(62,*)'                   '
        write(62,*)'z length of VOS = ',Zend-Zbng
        write(62,*)'y length of VOS = ',Yend-Ybng
        write(62,*)'Max distance of facets from a subgrid to be included in the intersection test'
        write(62,*)'	In z direction = ',range_Z
        write(62,*)'	In y direction = ',range_Y
        write(62,*)'                   '
    endif


    if(myid==master)then
        write(*,*) 'Raycasting progress (%)'
    endif

!------------------------------------  VOS using Raycasting & Subgrid -----------------------------------------
    do k=istart_vos,iend_vos
    !$OMP PARALLEL DO PRIVATE(nsy, nsz, nc, bubi, bubj, i, ETA_sub, ray_inter, Zsub, Ysub, c1, c2, c3, ETAs)  
    do j=jBgnVOS,jEndVOS
    ETA_sub=0.0
    do nsz=1,nSubGrids_3d; do nsy=1,nSubGrids_3d         ! Each point (k,j) in the bounding box is divided into subgrids
        ray_inter=Xend+5.0                         ! ray starts from Xend + 5.0  ( from spanwise direction)
        nc=1

        Zsub=Zs(k)-iDz(k)/2.+(nsz-0.5)*iDz(k)/nSubGrids_3d
        Ysub=Ys(j)-iDy(j)/2.+(nsy-0.5)*iDy(j)/nSubGrids_3d
        
!        !$OMP PARALLEL DO PRIVATE(c1,c2,c3) SHARED(nc,ray_inter)	
        do n=1,nf                                                     ! Only facets with centroid position within the possible range in z & y are included
                if ((ABS(M(3,n)-Zsub) .LE. range_Z) .AND. (ABS(M(2,n)-Ysub) .LE. range_Y)) then
                        ! Intersection test using the line and point equation (cross product) (Position of a Point Relative to a Line or side of triangle)
                        c1 = Vy( (n-1)*3+1 )*( Pz((n-1)*3+1) - Zsub ) - Vz( (n-1)*3+1 )*( Py((n-1)*3+1) - Ysub )
                        c2 = Vy( (n-1)*3+2 )*( Pz((n-1)*3+2) - Zsub ) - Vz( (n-1)*3+2 )*( Py((n-1)*3+2) - Ysub )
                        c3 = Vy( (n-1)*3+3 )*( Pz((n-1)*3+3) - Zsub ) - Vz( (n-1)*3+3 )*( Py((n-1)*3+3) - Ysub )
                        ! Finding the ray or x-coordinate of the intersection point (ray_inter)
                        if (( c1 > 0 .AND. c2 > 0 .AND. c3 > 0) .OR. (c1 < 0 .AND. c2 < 0 .AND. c3 < 0 )) then
!                                !$OMP CRITICAL
                                ray_inter(nc)=( FN(1,n)*M(1,n) + FN(2,n)*( M(2,n)-Ysub ) + FN(3,n)*( M(3,n)-Zsub ) ) /FN(1,n)
                                nc=nc+1
!                                !$OMP END CRITICAL
                        endif
                endif
        enddo
!        !$OMP END PARALLEL DO


        do bubi=1,bub-1; do bubj=bubi+1,bub
            if ( ray_inter(bubi) > ray_inter(bubj) )then
                bubble=ray_inter(bubi)
                ray_inter(bubi)=ray_inter(bubj)
                ray_inter(bubj)=bubble
            endif
        enddo; enddo
		

        nc=1
        ETAs=0.0
        do i=iBgnVOS,iEndVOS
            if ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) .AND. abs(Xs(i)-ray_inter(nc+1)) < (iDx(i)/2.) ) then
                ETA_sub(i)=ETA_sub(i) + ABS( ( ray_inter(nc+1)-ray_inter(nc) )/iDx(i) -ETAs) 
                nc=nc+2
            elseif ( abs(Xs(i)-ray_inter(nc)) < (iDx(i)/2.) ) then
                ETA_sub(i)=ETA_sub(i) + ABS( ( ray_inter(nc)-Xs(i) )/iDx(i) +ETAs -0.5 ) 
                ETAs=ABS(ETAs-1.0)
                nc=nc+1
            else 
                ETA_sub(i)=ETA_sub(i)+ETAs
            endif
        enddo 
    enddo ; enddo

!        !$OMP PARALLEL DO
        do i=iBgnVOS,iEndVOS
            ETA(i,j,k)=0.0
            ETA(i,j,k)=ETA_sub(i)/nSubGrids_3d/nSubGrids_3d
        enddo
!        !$OMP END PARALLEL DO
    enddo
    !$OMP END PARALLEL DO
    if(myid==master)then
		write(*,'(F6.2)') (k-istart_vos)*100.d0/(iend_vos-istart_vos)*1.d0
    endif
    enddo


      !----------data collect among nodes----------!
  
      icount = igcount_vos*(nx+4)*(ny+4)
      !Send my results back to the master
      if(myid>master)then
         itag = 410
         call MPI_SEND( ETA(-1,-1,istart_vos), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      !Wait to receive results from each task
      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount_vos(i)*(nx+4)*(ny+4)
            itag = 410
            call MPI_RECV( ETA(-1,-1,gstart_vos(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )    
         end do
      end if


	!----------data transformation from master to all nodes----------!
	icount= (nz+4)*(nx+4)*(ny+4)
	call MPI_BCAST ( ETA, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	
    
    if(myid==master)then
        write(*,*)'*******Ray casting SUCCESSFUL*******'
    endif

end
