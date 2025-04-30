! 22 Aug 2023 - FDS
subroutine reading_data()
   use variables
   implicit none

   nstep  	 	= NINT((total_time-initial_time)/dt) 	 ! number of timesteps for the simulation (total time-initial data time)
   isto3d       = NINT(isto3d_int/dt)          		 	 ! filer3d data storing steps interval
   isto2d       = NINT(isto2d_int/dt)          		 	 ! filer2d data storing steps interval
   istocp       = NINT(istocp_int/dt)          		 	 ! filer2d data storing steps interval
   istea        = NINT(istea_int/dt)         		     ! steps interval to check steadiness
   ibackup	 	= NINT(backup_int/dt)          	     	 ! backup data stored every 'ibackup' steps

   inv_rho = 1.d0/den_flu
   inv_dt  = 1.d0/dt

    
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

      nu= (U_inf*L_ch)/Re

      !-------------Damping Wall Function (stationary)-------------------!

      Cf=  (2*log10(Re)-0.65)**(-2.3)
      !Cf = 0.026/ (Re)**(1/7)
      tau = 0.5*Cf*den_flu *U_inf**2
      ufric = sqrt( tau/den_flu)
      Yplus =  (thickness*ufric)/nu
      fwall = (1.0 - exp(-Yplus/25.d0))
	  
      !-------------Damping Wall Function (dynamic)-------------------!

      Cf_d=  (2*log10(Re*(1.d0+DABS(blade_alpha)))-0.65)**(-2.3)
      !Cf = 0.026/ (Re)**(1/7)
      tau_d = 0.5*Cf_d*den_flu *(U_inf*(1.d0+DABS(blade_alpha)))**2
      ufric_d = sqrt( tau_d/den_flu)
      Yplus_d =  (thickness*ufric_d)/nu
      fwall_d = (1.0 - exp(-Yplus_d/25.d0))



end subroutine reading_data




subroutine reading_variables()
use variables
implicit none

open (18,file='OOXX.Q',form='unformatted',status='old', iostat=ierr)
if (ierr /= 0) then
	open (18,file='Static0.Q',form='unformatted',status='old', iostat=ierr)
    if (ierr /= 0) then
	  if(myid==master)then
      write(*,*) 'Error: could not open file'
	  endif
      stop
	else
		inputfile = 'Static0.Q'
    endif
else
	inputfile = 'OOXX.Q'
endif

read(18) inblocks
read(18) inx, iny, inz
read(18) temp, temp, temp, temp
read(18) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
read(18) lastsaved_time
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

end subroutine reading_variables


