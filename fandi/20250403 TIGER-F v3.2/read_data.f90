! 12 Jan 2025 - FDS

subroutine read_data()
   use variables
   implicit none

   nstep  	 	= NINT((total_time-initial_time)/dt) 	 ! number of timesteps for the simulation (total time-initial data time)
   isto3d       = NINT(isto3d_int/dt)          		 	 ! filer3d data storing steps interval
   istoVF		= NINT(istoVF_int/dt)					 ! Virtual force data storing steps interval
   istoAvg    	= NINT(istoAvg_int/dt)          		 ! TKE data storing steps interval
   isto2d       = NINT(isto2d_int/dt)          		 	 ! filer2d data storing steps interval
   istocp       = NINT(istocp_int/dt)          		 	 ! cp data storing steps interval
   istea        = NINT(istea_int/dt)         		     ! steps interval to check steadiness
   ibackup	 	= NINT(backup_int/dt)          	     	 ! backup data stored every 'ibackup' steps

   k_zeta_zeta	= k_zeta*zeta							 ! Dynamic omega parameter
   log_zeta 	= LOG(zeta)							 ! Dynamic omega parameter
   x0_sigmoid 	= LOG10(0.5d0*zeta*(k_zeta + 1.d0))	 ! The inflection point (midpoint) of sigmoid

	
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
read(28) lastsaved_time
read(28) totalFX_,totalFY_,totalFZ_
read(28) totalTorqx,totalTorqy,totalTorqz,totalTorq   
read(28) u_solid,v_solid,w_solid,rotor_omega
read(28) rotate_sx,rotate_sy,rotate_sz
read(28) x0_t,y0_t,z0_t,AOA
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


subroutine read_backup_uvwp_Avg()
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
read(29) istep2
close(29)


   do k=1,nz; do j=1,ny; do i=1,nx  

      p_TimeAvg(i,j,k) = Qout(i,j,k,1)
      u_TimeAvg(i,j,k) = Qout(i,j,k,2)
      v_TimeAvg(i,j,k) = Qout(i,j,k,3)
      w_TimeAvg(i,j,k) = Qout(i,j,k,4)
			TKE(i,j,k) = Qout(i,j,k,5)  	  

   enddo; enddo; enddo 


   temp=0.

end subroutine read_backup_uvwp_Avg