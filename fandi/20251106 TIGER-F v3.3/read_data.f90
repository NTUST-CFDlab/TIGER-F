! 12 Jan 2025 - FDS

subroutine read_data()
   use variables
   implicit none


   !*************** Toggle actions ****************
  	if (dimensionality==0) then
		Re = fixed_Re
		nu= (U_inf*L_ch)/Re
	else if(dimensionality==1) then
		nu = fixed_nu
		Re= (U_inf*L_ch)/nu
	endif

  	if (LES_mode=='ON') then
		LES = 1
	else if(LES_mode=='OFF') then
		LES = 0
	endif	
	

   !***************** Dynamic SOR *****************	
   k_zeta_zeta	= k_zeta*zeta						 ! Dynamic omega parameter
   log_zeta 	= LOG(zeta)							 ! Dynamic omega parameter
   x0_sigmoid 	= LOG10(0.5d0*zeta*(k_zeta + 1.d0))	 ! The inflection point (midpoint) of sigmoid


   !*********** Gridder initial processes *********
   if(Gridder=='non-uniform-sin4')then
      !Small interval
      dySml = lySml/nySml
      dzSml = lzSml/nzSml

      dx = lx / (nx*1.d0)     
	  
   ! else if(Gridder=='non-uniform-sin4-3D')then
      ! !Small interval
      ! dySml = lySml/nySml
      ! dzSml = lzSml/nzSml
      ! dxSml = lxSml/nxSml

   else if(Gridder=='uniform')then

      dx = lx / nx 
      dy = ly / ny
      dz = lz / nz
	  
	  dxSml = dx
      dySml = dy
      dzSml = dz 

   ! else if(Gridder=='ground-3D')then
      ! !Small interval
      ! dxMid = (lxMid-lxSml)/nxMid
      ! dyMid = (lyMid-lySml)/nyMid
      ! dzMid = (lzMid-lzSml)/nzMid

      ! dxSml = lxSml/nxSml
      ! dySml = lySml/nySml
      ! dzSml = lzSml/nzSml 
	  
   ! else if(Gridder=='non-uniform-sin5')then
      ! !Small interval
      ! dyMid = (lyMid-lySml)/nyMid
      ! dzMid = (lzMid-lzSml)/nzMid

      ! dySml = lySml/nySml
      ! dzSml = lzSml/nzSml

      ! dx = lx / (nx*1.d0)   

   ! else if(Gridder=='non-uniform-sin5-3D')then
      ! !Small interval
      ! dxMid = (lxMid-lxSml)/nxMid
      ! dyMid = (lyMid-lySml)/nyMid
      ! dzMid = (lzMid-lzSml)/nzMid

      ! dxSml = lxSml/nxSml
      ! dySml = lySml/nySml
      ! dzSml = lzSml/nzSml  	  

   else if(Gridder=='non-uniform-fall')then
      !Small interval
      dx = lx / nx 
      dy = ly / ny
      dz = lz / nz
	  
	  dxSml = dx
      dySml = lySml/nySml
      dzSml = dz  

   else if(Gridder=='non-uniform-fall-wall')then
      !Small interval  
	  dxSml = lxSml/nxSml
      dySml = lySml/nySml

      dzSml = lzSml/nzSml  

   ! else if(Gridder=='gridder_flatplate')then
      ! !Small interval

      ! dySml = lySml/nySml
      ! dzSml = lzSml/nzSml

      ! dx = lx / (nx*1.d0)


!   else if(Gridder=='non-uniform-sin2')then
      !Small interval
!      dySml = lySml/nySml
!      dzSml = lzSml/nzSml

!      dyMid = (lyMid-lySml)/(nyMid-nySml)
!      dzMid = (lzMid-lzSml)/(nzMid-nzSml)
      !Large interval
!      dy = ( ly-lyMid ) / ( ny - nyMid )
!      dz = ( lz-lzMid ) / ( nz - nzMid )

!      dx = lx / (nx*1.d0) 
  
!   else if(Gridder=='non-uniform-sin3')then
      !Small interval
!      dySml = lySml/nySml
!      dzSml = lzSml/nzSml

      !Large interval
!      dy = ( ly-lySml ) / ( ny - nySml )
!      dz = ( lz-lzSml ) / ( nz - nzSml )

!      dx = lx / (nx*1.d0)  
  

   end if

end subroutine read_data



subroutine read_backupfile()
use variables
implicit none

open (28,file='DYNAMIC.bak',form='unformatted',status='old', iostat=ierr)
if (ierr /= 0) then
	open (28,file='STATIC.bak',form='unformatted',status='old', iostat=ierr)
    if (ierr /= 0) then
	  if(myid==master)then
      write(*,*) 'Error: could not open file'
	  endif
      stop
	else
		inputfile = 'STATIC.bak'
    endif
else
	inputfile = 'DYNAMIC.bak'
endif

read(28) inblocks
read(28) inx, iny, inz
read(28) temp, temp, temp, temp
read(28) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
read(28) lastsaved_time,dt
read(28) totalFX_,totalFY_,totalFZ_
read(28) totalTorqx,totalTorqy,totalTorqz,totalTorq   
read(28) u_solid,v_solid,w_solid,rotor_omega
read(28) rotate_sx,rotate_sy,rotate_sz
read(28) x0_t,y0_t,z0_t,AOA
! read(28) totalFX1_,totalFY1_,totalFZ1_	! tandem cylinder
! read(28) totalFX2_,totalFY2_,totalFZ2_	! tandem cylinder
! read(28) u_solid_1,v_solid_1,w_solid_1	! tandem cylinder
! read(28) u_solid_2,v_solid_2,w_solid_2	! tandem cylinder
! read(28) x0_t1,y0_t1,z0_t1				! tandem cylinder
! read(28) x0_t2,y0_t2,z0_t2   			! tandem cylinder
close(28)


   do k=1,nz; do j=1,ny; do i=1,nx  

      p(i,j,k) = Qout(i,j,k,1)
      p_no(i,j,k) = Qout(i,j,k,1)
      p_pre(i,j,k,1) = Qout(i,j,k,1)   

   end do; end do; end do  
   

   do k=1,nz; do j=1,ny; do i=1,nx  

      u(i,j,k) = Qout(i,j,k,2)
      v(i,j,k) = Qout(i,j,k,3)
      w(i,j,k) = Qout(i,j,k,4)  

   enddo; enddo; enddo 


   temp=0.

end subroutine read_backupfile


subroutine read_backup_Avg()
use variables
implicit none

open (29,file='RUNNING_AVG.bak',form='unformatted',status='old', iostat=ierr)
if (ierr /= 0) then
	if(myid==master)then
      write(*,*) 'Error: could not open file'
      stop
    endif
else
	inputfile = 'RUNNING_AVG.bak'
endif

read(29) inblocks
read(29) inx, iny, inz
read(29) temp, temp, temp, temp
read(29) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
close(29)


   do k=1,nz; do j=1,ny; do i=1,nx  

      p_TimeAvg(i,j,k) = Qout(i,j,k,1)
      u_TimeAvg(i,j,k) = Qout(i,j,k,2)
      v_TimeAvg(i,j,k) = Qout(i,j,k,3)
      w_TimeAvg(i,j,k) = Qout(i,j,k,4)
			TKE(i,j,k) = Qout(i,j,k,5)  	  

   enddo; enddo; enddo 


   temp=0.

end subroutine read_backup_Avg


subroutine read_backup_Avg_full1()
use variables
implicit none

open (30,file='RUNNING_AVG1.bak',form='unformatted',status='old', iostat=ierr)
if (ierr /= 0) then
	if(myid==master)then
      write(*,*) 'Error: could not open file'
      stop
    endif
else
	inputfile = 'RUNNING_AVG1.bak'
endif

read(30) inblocks
read(30) inx, iny, inz
read(30) temp, temp, temp, temp
read(30) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
close(30)


   do k=1,nz; do j=1,ny; do i=1,nx  

      p_TimeAvg(i,j,k) = Qout(i,j,k,1)
      u_TimeAvg(i,j,k) = Qout(i,j,k,2)
      v_TimeAvg(i,j,k) = Qout(i,j,k,3)
      w_TimeAvg(i,j,k) = Qout(i,j,k,4)
		  RS_ww(i,j,k) = Qout(i,j,k,5)  	  

   enddo; enddo; enddo 


   temp=0.

end subroutine read_backup_Avg_full1



subroutine read_backup_Avg_full2()
use variables
implicit none

open (31,file='RUNNING_AVG2.bak',form='unformatted',status='old', iostat=ierr)
if (ierr /= 0) then
	if(myid==master)then
      write(*,*) 'Error: could not open file'
      stop
    endif
else
	inputfile = 'RUNNING_AVG2.bak'
endif

read(31) inblocks
read(31) inx, iny, inz
read(31) temp, temp, temp, temp
read(31) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
close(31)


   do k=1,nz; do j=1,ny; do i=1,nx  

		RS_uu(i,j,k) = Qout(i,j,k,1)
		RS_uv(i,j,k) = Qout(i,j,k,2)
		RS_uw(i,j,k) = Qout(i,j,k,3)
		RS_vv(i,j,k) = Qout(i,j,k,4)
		RS_vw(i,j,k) = Qout(i,j,k,5)  	  

   enddo; enddo; enddo 


   temp=0.

end subroutine read_backup_Avg_full2