! 22 Aug 2023 - FDS
subroutine read_data()
   use variables
   implicit none

   nstep  	 	= NINT((total_time-initial_time)/dt) 	 ! number of timesteps for the simulation (total time-initial data time)
   isto3d       = NINT(isto3d_int/dt)          		 	 ! filer3d data storing steps interval
   isto2d       = NINT(isto2d_int/dt)          		 	 ! filer2d data storing steps interval
   istocp       = NINT(istocp_int/dt)          		 	 ! filer2d data storing steps interval
   istea        = NINT(istea_int/dt)         		     ! steps interval to check steadiness
   ibackup	 	= NINT(backup_int/dt)          	     	 ! backup data stored every 'ibackup' steps
    
   if(Gridder=='non-uniform-sin4')then
      !Small interval
      dySml = lySml/nySml
      dzSml = lzSml/nzSml

      dx = lx / (nx*1.d0)     
	  
   else if(Gridder=='non-uniform-sin4-3D')then
      !Small interval
      dySml = lySml/nySml
      dzSml = lzSml/nzSml
      dxSml = lxSml/nxSml

   else if(Gridder=='uniform')then

      dx = lx / nx 
      dy = ly / ny
      dz = lz / nz
	  
	  dxSml = dx
      dySml = dy
      dzSml = dz 

   else if(Gridder=='ground-3D')then
      !Small interval
      dxMid = (lxMid-lxSml)/nxMid
      dyMid = (lyMid-lySml)/nyMid
      dzMid = (lzMid-lzSml)/nzMid

      dxSml = lxSml/nxSml
      dySml = lySml/nySml
      dzSml = lzSml/nzSml 
	  
   else if(Gridder=='non-uniform-sin5')then
      !Small interval
      dyMid = (lyMid-lySml)/nyMid
      dzMid = (lzMid-lzSml)/nzMid

      dySml = lySml/nySml
      dzSml = lzSml/nzSml

      dx = lx / (nx*1.d0)   

   else if(Gridder=='non-uniform-sin5-3D')then
      !Small interval
      dxMid = (lxMid-lxSml)/nxMid
      dyMid = (lyMid-lySml)/nyMid
      dzMid = (lzMid-lzSml)/nzMid

      dxSml = lxSml/nxSml
      dySml = lySml/nySml
      dzSml = lzSml/nzSml  	  


!   else if(Gridder=='non-uniform')then
      !Small interval
!      dySml = lySml/nySml
!      dzSml = lzSml/nzSml

      !Middle intervel
!      dyMid = (lyMid-lySml)/(nyMid)
!      dzMid = (lzMid-lzSml)/(nzMid)

      !Large interval
!      dy = ( ly-lyMid ) / ( ny - nyMid )
!      dz = ( lz-lzMid ) / ( nz - nzMid )

!      dx = lx / nx 

!   else if(Gridder=='non-uniform-sin')then
      !Small interval
!      dySml = lySml/nySml
!      dzSml = lzSml/nzSml

!      dyMid = (lyMid-lySml)/(nyMid-nySml)
!      dzMid = (lzMid-lzSml)/(nzMid-nzSml)
      !Large interval
!      dy = ( ly-lySml ) / ( ny - nySml )
!      dz = ( lz-lzMid ) / ( nz - nzMid )

!      dx = lx / (nx*1.d0) 

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

open (18,file='OOXX.bak',form='unformatted',status='old', iostat=ierr)
if (ierr /= 0) then
	open (18,file='Static0.bak',form='unformatted',status='old', iostat=ierr)
    if (ierr /= 0) then
	  if(myid==master)then
      write(*,*) 'Error: could not open file'
	  endif
      stop
	else
		inputfile = 'Static0.bak'
    endif
else
	inputfile = 'OOXX.bak'
endif

read(18) inblocks
read(18) inx, iny, inz
read(18) temp, temp, temp, temp
read(18) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
read(18) lastsaved_time
read(18) u_solid,v_solid,w_solid,rotor_omega,AOA
read(18) xc_t,yc_t,zc_t
read(18) totalFX_,totalFY_,totalFZ_,totalTorq
close(18)


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


