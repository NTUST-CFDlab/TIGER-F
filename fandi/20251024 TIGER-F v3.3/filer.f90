! 12 Jan 2025 - FDS

subroutine filer_final()
use variables
implicit none
   
   open (11,file='information.dat',position='append')
   write(11,*) ' ' 
   write(11,*) 'final simulation time = ',REAL(time),'sec'
   write(11,*) 'Last backup file time = ',REAL(lastbackuptime),'sec'
   write(11,*) ' '    
   write(11,*) '==================Simulation cost time======================'   
   write(11,*) 'total iterations = ',istep
   write(11,*) 'total cost time = ',REAL(totalcosttime),'sec'
   write(11,*) '============================================================'
   write(11,*) ' ' 
   close(11)

   close(21)
   close(31)
  
end subroutine filer_final

subroutine filerInfo()
use variables
implicit none
   real*8 :: xx(1:2), yy(1:2), zz(1:2)
   real*8 :: xratio, yratio, zratio
   real*8 :: xMax, yMax, zMax, xpos, ypos, zpos
   real*8 :: lyLrg1, lzLrg1, nyLrg1, nzLrg1
   
   open (1,file='information.dat',position='append')

   write(1,*) ' '
   write(1,*) ' '
   write(1,*) '!-----------------------------------------------------------!'
   write(1,*) '!  Finite Volume Method by projection method with mpi       !'
   write(1,*) '!-----------------------------------------------------------!'
   write(1,*) ' '
   write(1,*) '============================================================='
   write(1,*) 'number of processor(MPI)  = ',nproc
   write(1,*) 'number of threads(OpenMP) = ',nthreads
   write(1,*) 'number of grid = ',nx*ny*nz/1000000.d0,'M'
!   write(1,*) 'AOA = ',REAL(AOA)
   write(1,*) ' '
   write(1,*) 'x-component length = ',REAL(lx)
   write(1,*) 'y-component length = ',REAL(ly)
   write(1,*) 'z-component length = ',REAL(lz)
   write(1,*) ' '
   write(1,*) 'nx = ',nx
   write(1,*) 'ny = ',ny  
   write(1,*) 'nz = ',nz
   write(1,*) ' '
   write(1,*) '=========================gridder============================='
!   if(Gridder=='non-uniform')then
!   write(1,*) '========================non-uniform=========================='
!   write(1,*) 'dx = ',dx
!   write(1,*) 'Large dy = ',dy
!   write(1,*) 'Large dz = ',dz
!   write(1,*) ' '
!   write(1,*) 'nySml = ',nySml
!   write(1,*) 'nzSml = ',nzSml
!   write(1,*) ' '
!   write(1,*) 'lyMid = ',lyMid
!   write(1,*) 'lzMid = ',lzMid
!   write(1,*) ' '
!   write(1,*) 'Middle dy = ',dyMid
!   write(1,*) 'Middle dz = ',dzMid
!   write(1,*) ' '
!   write(1,*) 'Small dy = ',dySml
!   write(1,*) 'Small dz = ',dzSml
!   write(1,*) ' '
!   write(1,*) 'Small grid y0 = ',Griddery0
!   write(1,*) 'Small grid z0 = ',Gridderz0
!   write(1,*) ' '
!   write(1,*) 'lySml  = ',lySml
!   write(1,*) 'lzSml  = ',lzSml
!   write(1,*) '============================================================='
!   write(1,*) ' '
!   end if
!      if(Gridder=='non-uniform-sin2')then
!   write(1,*) '=====================non-uniform-sin2========================'
!   write(1,*) 'dx = ',REAL(dx)
!   write(1,*) 'Large dy = ',REAL(dy)
!   write(1,*) 'Large dz = ',REAL(dz)
!   write(1,*) ' '
!   write(1,*) 'Middle dy = ',REAL(dyMid)
!   write(1,*) 'Middle dz = ',REAL(dzMid)
!   write(1,*) ' '
!   write(1,*) 'Small dy = ',REAL(dySml)
!   write(1,*) 'Small dz = ',REAL(dzSml)
!   write(1,*) ' '
!   write(1,*) 'Small grid y0 = ',REAL(Griddery0)
!   write(1,*) 'Small grid z0 = ',REAL(Gridderz0)
!   write(1,*) ' '
!   write(1,*) 'lyMid = ',REAL(lyMid)
!   write(1,*) 'nyMid = ',nyMid
!   write(1,*) ' '
!   write(1,*) 'lySml = ',REAL(lySml)
!   write(1,*) 'nySml = ',nySml
!   write(1,*) ' '
!   write(1,*) 'lzMid = ',REAL(lzMid)
!   write(1,*) 'nzMid = ',nzMid
!   write(1,*) ' '
!   write(1,*) 'lzSml = ',REAL(lzSml)
!   write(1,*) 'nzSml = ',nzSml
!   write(1,*) '============================================================='
!   write(1,*) ' '
!   end if
!      if(Gridder=='non-uniform-sin')then
!   write(1,*) '======================non-uniform-sin========================'
!   write(1,*) 'Small grid y0 = ',REAL(Griddery0)
!   write(1,*) 'Small grid z0 = ',REAL(Gridderz0)
!   write(1,*) ' '
!   write(1,*) 'dx = ',REAL(dx)
!   write(1,*) 'Large dy = ',REAL(dy)
!   write(1,*) 'Large dz = ',REAL(dz)
!   write(1,*) ' '
!   write(1,*) 'Middle dy = ',REAL(dyMid)
!   write(1,*) 'Middle dz = ',REAL(dzMid)
!   write(1,*) ' '
!   write(1,*) 'Small dy = ',REAL(dySml)
!   write(1,*) 'Small dz = ',REAL(dzSml)
!   write(1,*) ' '
!   write(1,*) 'lyMid = ',REAL(lyMid)
!   write(1,*) 'nyMid = ',nyMid
!   write(1,*) ' '
!   write(1,*) 'lySml = ',REAL(lySml)
!   write(1,*) 'nySml = ',nySml
!   write(1,*) ' '
!   write(1,*) 'lzMid = ',REAL(lzMid)
!   write(1,*) 'nzMid = ',nzMid
!   write(1,*) ' '
!   write(1,*) 'lzSml = ',REAL(lzSml)
!   write(1,*) 'nzSml = ',nzSml
!   write(1,*) '============================================================='
!   write(1,*) ' '
!      end if
!      if(Gridder=='non-uniform-sin3')then
!   write(1,*) '=====================non-uniform-sin3========================'
!   write(1,*) 'dx = ',REAL(dx)
!   write(1,*) 'Large dy = ',REAL(dy)
!   write(1,*) 'Large dz = ',REAL(dz)
!   write(1,*) ' '
!   write(1,*) 'Small dy = ',REAL(dySml)
!   write(1,*) 'Small dz = ',REAL(dzSml)
!   write(1,*) ' '
!   write(1,*) 'lySml = ',REAL(lySml)
!   write(1,*) 'nySml = ',nySml
!   write(1,*) ' '
!   write(1,*) 'lzSml = ',REAL(lzSml)
!   write(1,*) 'nzSml = ',nzSml
!   write(1,*) ' '
                ! Grid ratio check
!        xMax=0.
!        do i=1,nx-1
!                xx(1) = iDx(i)
!                xx(2) = iDx(i+1)
!                xratio = maxval(xx)/minval(xx)*1.d0
!                if (xratio > xMax) then
!                        xMax = xratio
!                end if
!        end do

!        zMax=0.
!        do i=1,nz-1
!                zz(1) = iDz(i)
!                zz(2) = iDz(i+1)
!                zratio = maxval(zz)/minval(zz)*1.d0
!                if (zratio > zMax) then
!                        zMax = zratio
!                end if
!        end do

!        yMax=0.
!        do i=1,ny-1
!                yy(1) = iDy(i)
!                yy(2) = iDy(i+1)
!                yratio = maxval(yy)/minval(yy)*1.d0
!                if (yratio > yMax) then
!                        yMax = yratio
!                end if
!        end do

!   write(1,*) 'largest z ratio:',zMax,' at Z =',zpos
!   write(1,*) 'largest y ratio:',yMax,' at Y =',ypos
!   write(1,*) 'largest x ratio:',xMax   
!   write(1,*) '============================================================='
!   write(1,*) ' '
!   end if
      if(Gridder=='non-uniform-sin4')then 
   write(1,*) '=====================non-uniform-sin4========================'
   write(1,*) 'dx = ',REAL(dx)
   write(1,*) ' '
   write(1,*) 'Small dy = ',REAL(dySml)
   write(1,*) 'Small dz = ',REAL(dzSml)
   write(1,*) ' '
   write(1,*) 'lySml = ',REAL(lySml)
   write(1,*) 'nySml = ',nySml
   write(1,*) ' '
   write(1,*) 'lzSml = ',REAL(lzSml)
   write(1,*) 'nzSml = ',nzSml
   write(1,*) ' '
                ! Grid ratio check
        xMax=0.
        do i=1,nx-1
                xx(1) = iDx(i)
                xx(2) = iDx(i+1)
                xratio = maxval(xx)/minval(xx)*1.d0
                if (xratio > xMax) then
                        xMax = xratio
                end if
        end do

        zMax=0.
        do i=1,nz-1
                zz(1) = iDz(i)
                zz(2) = iDz(i+1)
                zratio = maxval(zz)/minval(zz)*1.d0
                if (zratio > zMax) then
                        zMax = zratio
                        zpos = Z(i)
                end if
        end do

        yMax=0.
        do i=1,ny-1
                yy(1) = iDy(i)
                yy(2) = iDy(i+1)
                yratio = maxval(yy)/minval(yy)*1.d0
                if (yratio > yMax) then
                        yMax = yratio
                        ypos = Y(i)
                end if
        end do

   write(1,*) 'largest z ratio:',zMax,' at Z =',zpos
   write(1,*) 'largest y ratio:',yMax,' at Y =',ypos
   write(1,*) 'largest x ratio:',xMax  
   write(1,*) '============================================================='
   write(1,*) ' '
   end if  
      ! if(Gridder=='non-uniform-sin4-3D')then 
   ! write(1,*) '====================non-uniform-sin4-3D======================'
   ! write(1,*) 'Small dx = ',REAL(dxSml)
   ! write(1,*) 'Small dy = ',REAL(dySml)
   ! write(1,*) 'Small dz = ',REAL(dzSml)
   ! write(1,*) ' '
   ! write(1,*) 'lxSml = ',REAL(lxSml)
   ! write(1,*) 'nxSml = ',nxSml
   ! write(1,*) ' '
   ! write(1,*) 'lySml = ',REAL(lySml)
   ! write(1,*) 'nySml = ',nySml
   ! write(1,*) ' '
   ! write(1,*) 'lzSml = ',REAL(lzSml)
   ! write(1,*) 'nzSml = ',nzSml
   ! write(1,*) ' '
                ! ! Grid ratio check
        ! xMax=0.
        ! do i=1,nx-1
                ! xx(1) = iDx(i)
                ! xx(2) = iDx(i+1)
                ! xratio = maxval(xx)/minval(xx)*1.d0
                ! if (xratio > xMax) then
                        ! xMax = xratio
                        ! xpos = X(i)
                ! end if
        ! end do

        ! zMax=0.
        ! do i=1,nz-1
                ! zz(1) = iDz(i)
                ! zz(2) = iDz(i+1)
                ! zratio = maxval(zz)/minval(zz)*1.d0
                ! if (zratio > zMax) then
                        ! zMax = zratio
                        ! zpos = Z(i)
                ! end if
        ! end do

        ! yMax=0.
        ! do i=1,ny-1
                ! yy(1) = iDy(i)
                ! yy(2) = iDy(i+1)
                ! yratio = maxval(yy)/minval(yy)*1.d0
                ! if (yratio > yMax) then
                        ! yMax = yratio
                        ! ypos = Y(i)
                ! end if
        ! end do


   ! write(1,*) 'largest x ratio:',xMax,' at X =',xpos
   ! write(1,*) 'largest y ratio:',yMax,' at Y =',ypos
   ! write(1,*) 'largest z ratio:',zMax,' at Z =',zpos 
   ! write(1,*) '============================================================='
   ! write(1,*) ' '
   ! end if  
   
      ! if(Gridder=='non-uniform-sin5')then 
   ! write(1,*) '=====================non-uniform-sin5========================'
   ! write(1,*) 'dx = ',REAL(dx)
   ! write(1,*) 'Mid dy = ',REAL(dyMid)
   ! write(1,*) 'Mid dz = ',REAL(dzMid)
   ! write(1,*) ' '
   ! write(1,*) 'Small dy = ',REAL(dySml)
   ! write(1,*) 'Small dz = ',REAL(dzSml)
   ! write(1,*) ' '
   ! write(1,*) 'lySml = ',REAL(lySml)
   ! write(1,*) 'nySml = ',nySml
   ! write(1,*) ' '
   ! write(1,*) 'lzSml = ',REAL(lzSml)
   ! write(1,*) 'nzSml = ',nzSml
   ! write(1,*) ' '
                ! ! Grid ratio check
        ! xMax=0.
        ! do i=1,nx-1
                ! xx(1) = iDx(i)
                ! xx(2) = iDx(i+1)
                ! xratio = maxval(xx)/minval(xx)*1.d0
                ! if (xratio > xMax) then
                        ! xMax = xratio
                ! end if
        ! end do

        ! zMax=0.
        ! do i=1,nz-1
                ! zz(1) = iDz(i)
                ! zz(2) = iDz(i+1)
                ! zratio = maxval(zz)/minval(zz)*1.d0
                ! if (zratio > zMax) then
                        ! zMax = zratio
                        ! zpos = Z(i)
                ! end if
        ! end do

        ! yMax=0.
        ! do i=1,ny-1
                ! yy(1) = iDy(i)
                ! yy(2) = iDy(i+1)
                ! yratio = maxval(yy)/minval(yy)*1.d0
                ! if (yratio > yMax) then
                        ! yMax = yratio
                        ! ypos = Y(i)
                ! end if
        ! end do


   ! write(1,*) 'largest z ratio:',zMax,' at Z =',zpos
   ! write(1,*) 'largest y ratio:',yMax,' at Y =',ypos
   ! write(1,*) 'largest x ratio:',xMax  
   ! write(1,*) '============================================================='
   ! write(1,*) ' '
   ! end if  

      ! if(Gridder=='non-uniform-sin5-3D')then 
   ! write(1,*) '=====================non-uniform-sin5-3D====================='
   ! write(1,*) 'Small dx = ',REAL(dxSml)
   ! write(1,*) 'Small dy = ',REAL(dySml)
   ! write(1,*) 'Small dz = ',REAL(dzSml)
   ! write(1,*) ' '
   ! write(1,*) 'lxSml = ',REAL(lxSml)
   ! write(1,*) 'nxSml = ',nxSml
   ! write(1,*) ' '
   ! write(1,*) 'lySml = ',REAL(lySml)
   ! write(1,*) 'nySml = ',nySml
   ! write(1,*) ' '
   ! write(1,*) 'lzSml = ',REAL(lzSml)
   ! write(1,*) 'nzSml = ',nzSml
   ! write(1,*) ' '
                ! ! Grid ratio check
        ! xMax=0.
        ! do i=1,nx-1
                ! xx(1) = iDx(i)
                ! xx(2) = iDx(i+1)
                ! xratio = maxval(xx)/minval(xx)*1.d0
                ! if (xratio > xMax) then
                        ! xMax = xratio
                        ! xpos = X(i)
                ! end if
        ! end do

        ! zMax=0.
        ! do i=1,nz-1
                ! zz(1) = iDz(i)
                ! zz(2) = iDz(i+1)
                ! zratio = maxval(zz)/minval(zz)*1.d0
                ! if (zratio > zMax) then
                        ! zMax = zratio
                        ! zpos = Z(i)
                ! end if
        ! end do

        ! yMax=0.
        ! do i=1,ny-1
                ! yy(1) = iDy(i)
                ! yy(2) = iDy(i+1)
                ! yratio = maxval(yy)/minval(yy)*1.d0
                ! if (yratio > yMax) then
                        ! yMax = yratio
                        ! ypos = Y(i)
                ! end if
        ! end do


   ! write(1,*) 'largest x ratio:',xMax,' at X =',xpos
   ! write(1,*) 'largest y ratio:',yMax,' at Y =',ypos
   ! write(1,*) 'largest z ratio:',zMax,' at Z =',zpos 
   ! write(1,*) '============================================================='
   ! write(1,*) ' '
   ! end if  

      ! if(Gridder=='ground-3D')then 
   ! write(1,*) '=========================ground-3D==========================='
   ! write(1,*) 'Small dx = ',REAL(dxSml)
   ! write(1,*) 'Small dy = ',REAL(dySml)
   ! write(1,*) 'Small dz = ',REAL(dzSml)
   ! write(1,*) ' '
   ! write(1,*) 'lxSml = ',REAL(lxSml)
   ! write(1,*) 'nxSml = ',nxSml
   ! write(1,*) ' '
   ! write(1,*) 'lySml = ',REAL(lySml)
   ! write(1,*) 'nySml = ',nySml
   ! write(1,*) ' '
   ! write(1,*) 'lzSml = ',REAL(lzSml)
   ! write(1,*) 'nzSml = ',nzSml
   ! write(1,*) ' '
                ! ! Grid ratio check
        ! xMax=0.
        ! do i=1,nx-1
                ! xx(1) = iDx(i)
                ! xx(2) = iDx(i+1)
                ! xratio = maxval(xx)/minval(xx)*1.d0
                ! if (xratio > xMax) then
                        ! xMax = xratio
                        ! xpos = X(i)
                ! end if
        ! end do

        ! zMax=0.
        ! do i=nzLrg1+nzSml,nz-1
                ! zz(1) = iDz(i)
                ! zz(2) = iDz(i+1)
                ! zratio = maxval(zz)/minval(zz)*1.d0
                ! if (zratio > zMax) then
                        ! zMax = zratio
                        ! zpos = Z(i)
                ! end if
        ! end do

        ! yMax=0.
        ! do i=1,ny-1
                ! yy(1) = iDy(i)
                ! yy(2) = iDy(i+1)
                ! yratio = maxval(yy)/minval(yy)*1.d0
                ! if (yratio > yMax) then
                        ! yMax = yratio
                        ! ypos = Y(i)
                ! end if
        ! end do


   ! write(1,*) 'largest x ratio:',xMax,' at X =',xpos
   ! write(1,*) 'largest y ratio:',yMax,' at Y =',ypos
   ! write(1,*) 'largest z ratio:',zMax,' at Z =',zpos 
   ! write(1,*) '============================================================='
   ! write(1,*) ' '
   ! end if  

      if(Gridder=='non-uniform-fall' .OR. Gridder=='non-uniform-fall-wall')then 
   write(1,*) '=========================ground-3D==========================='
   write(1,*) 'Small dx = ',REAL(dxSml)
   write(1,*) 'Small dy = ',REAL(dySml)
   write(1,*) 'Small dz = ',REAL(dzSml)
   write(1,*) ' '
   write(1,*) 'lxSml = ',REAL(lxSml)
   write(1,*) 'nxSml = ',nxSml
   write(1,*) ' '
   write(1,*) 'lySml = ',REAL(lySml)
   write(1,*) 'nySml = ',nySml
   write(1,*) ' '
   write(1,*) 'lzSml = ',REAL(lzSml)
   write(1,*) 'nzSml = ',nzSml
   write(1,*) ' '
                ! Grid ratio check
        xMax=0.
        do i=1,nx-1
                xx(1) = iDx(i)
                xx(2) = iDx(i+1)
                xratio = maxval(xx)/minval(xx)*1.d0
                if (xratio > xMax) then
                        xMax = xratio
                        xpos = X(i)
                end if
        end do

        zMax=0.
        do i=1,nz-1
                zz(1) = iDz(i)
                zz(2) = iDz(i+1)
                zratio = maxval(zz)/minval(zz)*1.d0
                if (zratio > zMax) then
                        zMax = zratio
                        zpos = Z(i)
                end if
        end do

        yMax=0.
        do i=1,ny-1
                yy(1) = iDy(i)
                yy(2) = iDy(i+1)
                yratio = maxval(yy)/minval(yy)*1.d0
                if (yratio > yMax) then
                        yMax = yratio
                        ypos = Y(i)
                end if
        end do


   write(1,*) 'largest x ratio:',xMax,' at X =',xpos
   write(1,*) 'largest y ratio:',yMax,' at Y =',ypos
   write(1,*) 'largest z ratio:',zMax,' at Z =',zpos 
   write(1,*) '============================================================='
   write(1,*) ' '
   end if  

      if(Gridder=='read_txt')then 
   write(1,*) '=====================read_txt========================'
   write(1,*) 'dx = ',REAL(dx)
   write(1,*) ' '
   write(1,*) 'Small dy = ',REAL(dySml)
   write(1,*) 'Small dz = ',REAL(dzSml)
   write(1,*) ' '
   write(1,*) 'lySml = ',REAL(lySml)
   write(1,*) 'nySml = ',nySml
   write(1,*) ' '
   write(1,*) 'lzSml = ',REAL(lzSml)
   write(1,*) 'nzSml = ',nzSml
   write(1,*) ' '
                ! Grid ratio check
        xMax=0.
        do i=1,nx-1
                xx(1) = iDx(i)
                xx(2) = iDx(i+1)
                xratio = maxval(xx)/minval(xx)*1.d0
                if (xratio > xMax) then
                        xMax = xratio
                end if
        end do

        zMax=0.
        do i=1,nz-1
                zz(1) = iDz(i)
                zz(2) = iDz(i+1)
                zratio = maxval(zz)/minval(zz)*1.d0
                if (zratio > zMax) then
                        zMax = zratio
                        zpos = Z(i)
                end if
        end do

        yMax=0.
        do i=1,ny-1
                yy(1) = iDy(i)
                yy(2) = iDy(i+1)
                yratio = maxval(yy)/minval(yy)*1.d0
                if (yratio > yMax) then
                        yMax = yratio
                        ypos = Y(i)
                end if
        end do

   write(1,*) 'largest z ratio:',zMax,' at Z =',zpos
   write(1,*) 'largest y ratio:',yMax,' at Y =',ypos
   write(1,*) 'largest x ratio:',xMax  
   write(1,*) '============================================================='
   write(1,*) ' '
   end if  

      if(Gridder=='uniform')then
   write(1,*) '==========================uniform============================'
   write(1,*) 'dx = ',dx
   write(1,*) 'dy = ',dy
   write(1,*) 'dz = ',dz
   write(1,*) '============================================================='
   write(1,*) ' '
   endif
   
   
   write(1,*) '=========================VoS Info============================'
   if (VOS_by == 2) then; write(1,*) 'VOS created by ray0asting2D'
						  write(1,'(A)') ' Input DAT file = '//dat_file
   else if(VOS_by == 3) then; write(1,*) 'VOS created by ray0asting3D'
						  write(1,'(A)') ' Input STL file = '//stl_file
   else if(VOS_by == 1) then; write(1,*) 'VOS created by curve equation'
   endif
   write(1,*) ' '   
   write(1,*) 'DFIB x0 = ',REAL(x0)
   write(1,*) 'DFIB y0 = ',REAL(y0)
   write(1,*) 'DFIB z0 = ',REAL(z0)
	if (select_ref_area == 1) then						! Reference area in virtual force calculation
		write(1,*) 'Reference area:Z_castingarea =', REAL(Z_castingarea)
    elseif (select_ref_area == 2) then		
		write(1,*) 'Reference area:Y_castingarea =', REAL(Y_castingarea)
    elseif (select_ref_area == 3) then		
		write(1,*) 'Reference area:X_castingarea =', REAL(X_castingarea)
    elseif (select_ref_area == 4) then		
		write(1,*) 'Reference area:user defined =', REAL(ref_area)
	endif
   write(1,*)'Casting_volume = ', REAL(solid_volume)
   write(1,*)'Conversion cost = ',REAL(vosfinaltime-vosstarttime),'sec'
   write(1,*) '============================================================='
   write(1,*) ' ' 
   write(1,*) '====================simulation properties====================='
   write(1,*)'Time-stepping strategy :'
   if (time_stepping_by=='DT') then
		write(1,*) '	constant dt = ',REAL(constant_dt)
   else if(time_stepping_by=='CFL') then
		write(1,*) '	constant CFL = ',REAL(constant_cfl)
   endif
   write(1,*) ' '
   if (dimensionality==0) then
		write(1,*) 'Non-Dimensional simulation'
		write(1,*) '	Re = ',REAL(fixed_Re)
   elseif (dimensionality==1) then
		write(1,*) 'Dimensional simulation'
		write(1,*) '	nu = ',REAL(fixed_nu)
		write(1,*) '	Fluid density = ',REAL(den_flu)
		write(1,*) '	Free stream velocity = ', REAL(U_inf) ,' m/s '   
   endif
   write(1,*) ' ' 
   if (solid_motion == 1) then 
   write(1,*) 'SOLID MOTION: 1-way coupling', ', start time = ', REAL(StartDynamic_time)
   write(1,*) ' '
   write(1,*) '===================dynamic wall parameters==================='
   write(1,*) 'CFL x = ', REAL(U_inf*(1.d0+ABS(blade_alpha))*dt / dzsml),&
          '  , CFL y = ', REAL(U_inf*(1.d0+ABS(blade_alpha))*dt / dySml)
!   write(1,*) 'Cf  = ', REAL(Cf) ,'  '
!   write(1,*) 'Tau = ', REAL(tau) ,'  '
!   write(1,*) 'Y plus wall distance estimation  = ', REAL(Yplus) ,'  '   
!   write(1,*) 'Van Driest  damping function  = ', REAL(fwall) ,'  ' 
   write(1,*) '============================================================='
   write(1,*) ' '
   else if (solid_motion == 2) then
   write(1,*) 'SOLID MOTION: 2-way coupling', ', start time = ', REAL(StartDynamic_time)
   write(1,*) ' '
   write(1,*) '===================dynamic wall parameters==================='
   write(1,*) 'CFL x = ', REAL(U_inf*(1.d0+ABS(blade_alpha))*dt / dzsml),&
          '  , CFL y = ', REAL(U_inf*(1.d0+ABS(blade_alpha))*dt / dySml)
!   write(1,*) 'Cf  = ', REAL(Cf) ,'  '
!   write(1,*) 'Tau = ', REAL(tau) ,'  '
!   write(1,*) 'Y plus wall distance estimation  = ', REAL(Yplus) ,'  '   
!   write(1,*) 'Van Driest  damping function  = ', REAL(fwall) ,'  ' 
   write(1,*) '============================================================='
   write(1,*) ' '
   else
   write(1,*) 'SOLID MOTION: STATIONARY'
   write(1,*) ' '
   write(1,*) '====================static wall parameters==================='
!   write(1,*) 'Moment of intertia = ', REAL(I_moment) ,' kg.m2 '
   write(1,*) 'CFL x = ', REAL(U_inf*dt / dzsml), '  , CFL y = ', REAL(U_inf*dt / dySml)
!   write(1,*) 'Fo = ', dt/Re/dySml/dzSml
   write(1,*) 'Cf  = ', REAL(Cf) ,'  '
!   write(1,*) 'Tau = ', REAL(tau) ,'  '
!   write(1,*) 'friction velocity  = ', REAL(ufric) ,'  '
!   write(1,*) 'Y plus wall distance estimation  = ', REAL(Yplus) ,'  '   
!   write(1,*) 'Van Driest  damping function  = ', REAL(fwall) ,'  ' 
   write(1,*) '============================================================='
   write(1,*) ' '
   end if
   write(1,*) '====================Simulation time=========================='
   if (resume=='ON') then; write(1,*) 'resume mode = ON'
						  write(1,'(A)') ' Input file = '//inputfile
   else if(resume=='OFF') then; write(1,*) 'resume mode = OFF'
   end if
   write(1,*) ' '
   write(1,*) 'Initial time   = ',REAL(initial_time),'sec'
   write(1,*) 'Max time step  = ',REAL(total_time),'sec'
   write(1,*) '============================================================='
   write(1,*) ' '
   write(1,*) '======================Output files==========================='
   if (filer_cp=='ON') then
   write(1,*) '>>>>> Cp plot'
   write(1,*) 'start at ',REAL(coeff_start),'sec,  every ',REAL(filerCp_int),'sec'
   else if(filer_cp=='OFF') then; write(1,*) '>>>>> Cp plot: OFF'
   end if
   write(1,*) ' '
   if (filer3d=='ON') then
   write(1,*) '>>>>> filer3D'
   write(1,*) 'start at ',REAL(startfiler3d_time),'sec,  every ',REAL(filer3d_int),'sec'
   write(1,*) 'Domain clipping: '
   write(1,*) '	X = ',x_start,'<--->',x_end
   write(1,*) '	Y = ',y_start,'<--->',y_end
   write(1,*) '	Z = ',z_start,'<--->',z_end
   write(1,*) ' '
   write(1,*) '	RHO   = P'
   write(1,*) '	RHO-U = u'
   write(1,*) '	RHO-V = v'
   write(1,*) '	RHO-W = w'
   write(1,*) '	E     = ETA'   
   else if(filer3d=='OFF') then; write(1,*) '>>>>> filer3D: OFF'
   end if
   write(1,*) ' '
   if (filer2d=='AVG') then
   write(1,*) '>>>>> filer2D Spanwise avg '
   write(1,*) 'start at ',REAL(startfiler2d_time),'sec,  every ',REAL(filer2d_int),'sec'
   write(1,*) ' '
   write(1,*) '	RHO   = P'
   write(1,*) '	RHO-U = u'
   write(1,*) '	RHO-V = v'
   write(1,*) '	RHO-W = w'
   write(1,*) '	E     = ETA' 
   else if (filer2d=='MID') then
   write(1,*) '>>>>> filer2D Spanwise Midplane slice'
   write(1,*) 'start at ',REAL(startfiler2d_time),'sec,  every ',REAL(filer2d_int),'sec'
   write(1,*) ' '
   write(1,*) '	RHO   = P'
   write(1,*) '	RHO-U = u'
   write(1,*) '	RHO-V = v'
   write(1,*) '	RHO-W = w'
   write(1,*) '	E     = ETA' 
   else if(filer2d=='OFF') then; write(1,*) '>>>>> filer2D: OFF'
   end if
   write(1,*) ' '
   if (filer_virtual_force=='ON') then
   write(1,*) '>>>>> filer virtual force'
   write(1,*) 'start at ',REAL(startfilerVF_time),'sec,  every ',REAL(filerVF_int),'sec'
   write(1,*) 'Domain clipping: '
   write(1,*) '	X = ',x_start,'<--->',x_end
   write(1,*) '	Y = ',y_start,'<--->',y_end
   write(1,*) '	Z = ',z_start,'<--->',z_end
   write(1,*) ' '
   write(1,*) '	RHO   = P'
   write(1,*) '	RHO-U = solid FX'
   write(1,*) '	RHO-V = solid FY'
   write(1,*) '	RHO-W = solid FZ'
   write(1,*) '	E     = ETA'   
   else if(filer_virtual_force=='OFF') then; write(1,*) '>>>>> filer virtual force: OFF'
   end if
   write(1,*) ' '
   if (filer_running_avg=='STANDARD' .OR. filer_running_avg=='ON') then
   write(1,*) '>>>>> filer P,u,v,w,TKE time averaged'
   write(1,*) 'start at ',REAL(startfilerAvg_time),'sec,  every ',REAL(filerAvg_int),'sec'
   write(1,*) 'Domain clipping: '
   write(1,*) '	X = ',x_start,'<--->',x_end
   write(1,*) '	Y = ',y_start,'<--->',y_end
   write(1,*) '	Z = ',z_start,'<--->',z_end
   write(1,*) ' '
   write(1,*) '	RHO   = P (avg)'
   write(1,*) '	RHO-U = u (avg)'
   write(1,*) '	RHO-V = v (avg'
   write(1,*) '	RHO-W = w (avg)'
   write(1,*) '	E     = TKE (avg)'   
   else if (filer_running_avg=='FULL') then
   write(1,*) '>>>>> filer P,u,v,w,Reynolds stresses time averaged'
   write(1,*) 'start at ',REAL(startfilerAvg_time),'sec,  every ',REAL(filerAvg_int),'sec'
   write(1,*) 'Domain clipping: '
   write(1,*) '	X = ',x_start,'<--->',x_end
   write(1,*) '	Y = ',y_start,'<--->',y_end
   write(1,*) '	Z = ',z_start,'<--->',z_end
   write(1,*) ' '
   write(1,*) 'Avg1_                    Avg2_'
   write(1,*) '	RHO   = P (avg)          RHO   = Rs_uu (avg)'
   write(1,*) '	RHO-U = u (avg)          RHO-U = Rs_uv (avg)'
   write(1,*) '	RHO-V = v (avg)          RHO-V = Rs_uw (avg)'
   write(1,*) '	RHO-W = w (avg)          RHO-W = Rs_vv (avg)'
   write(1,*) '	E     = RS_ww (avg)      E     = Rs_vw (avg)' 
   else if(filer_running_avg=='OFF') then; write(1,*) '>>>>> filer P,u,v,w,TKE time averaged: OFF'
   end if
   write(1,*) '============================================================='
   ! write(1,*) ' '
   ! write(1,*) '============================================================='
   ! if (DBD == 1) then; write(1,*) 'DBD actuator on'
   ! else if(DBD == 0) then; write(1,*) 'DBD actuator off'
   ! end if
   ! write(1,*) '============================================================='
   write(1,*) ' '
   write(1,*) '============================================================='
   if (LES_mode=='ON') then; write(1,*) 'LES mode on'
   else if(LES_mode=='OFF') then; write(1,*) 'LES mode off'
   end if
   write(1,*) 'LES Cs = ',REAL(Cs)
!   write(1,*) 'Reduced Frequency k = ',REAL(reduce_frequency_k)
   write(1,*) '============================================================='
   write(1,*) ' '
   write(1,*) '============================================================='
   if (correction_stage == 0) then
	write(1,*) 'Prediction-correction = No correction stage'
   else
	write(1,*) 'Prediction-correction =',correction_stage,' correction stage'
   end if
   write(1,*) '============================================================='
   write(1,*) ' '
   write(1,*) '============================================================='
   write(1,*) 'p Residual =',REAL(zeta)
!   write(1,*) 'velocity Residual =',REAL(zeta_vel)
   write(1,*) '============================================================='
   close(1)



   open (21,file='cd_time.dat',position='append')
   write(21,*) ' TITLE     = "" '
   write(21,*) ' VARIABLES = t*,C<sub>D</sub>,C<sub>L</sub> '
   !write(21,*) ' VARIABLES = t*,AOA(<math>0</math>),C<sub>D</sub>,C<sub>L</sub> '
  

   open (31,file='ct_time.dat',position='append')
   write(31,*) ' TITLE     = "" '
   write(31,*) ' VARIABLES = t*,AOA,C<sub>T</sub>,C<sub>P</sub>,T<sub>gen</sub>,P<sub>gen</sub>'


   open (51,file='cp.dat',position='append')
   write(51,*) ' TITLE     = "" '
   write(51,*) ' VARIABLES = x/c,C<sub>p</sub>,Tan_vel'



end subroutine filerInfo



subroutine filerProcess()
use variables
implicit none

   write(21,*) REAL(time), REAL(cDrag), REAL(cLift)
   !write(21,*) REAL(time), REAL(AOA), REAL(cDrag), REAL(cLift)

   write(31,'(F12.7,3X,F14.6,4(3X,F12.7))') REAL(time), REAL(AOA), REAL(cTorq), REAL(cPower), REAL(T_gen), REAL(P_gen)
   
end subroutine filerProcess



subroutine filerProcess_cp()
use variables
implicit none
real*8 :: x_c, p_avg, tan_vel_avg, w_mid, v_mid, dot
integer :: cell_shift


!-----------Cp for cylinder of D = 1 with theta from 0 to 180 (CW from x-) & pressure reference at inlet (1,1)--------------
!------------------------------------ also calculate near boundary velocity parallel to surface-----------------------------
!write(51,*) 'ZONE solutiontime=',  REAL(time)
!  !write(51,*) 'ZONE T = "p at <greek>h</greek>=0.5 , p<sub><math>%</math></sub>=(1,1)"'
!     do k=1,nz
!      p_avg=0.d0
!      tan_vel_avg=0.d0
!      do j=1,ny
!         if (ETA(1,j,k)>=0.5d0 .AND. ETA(1,j+1,k)<0.5d0) then
!            if (ABS(Zs(k)-z0) .GT. 0.5d0) cy0le
!            x_c=DACOS(-2.d0*(Zs(k)-z0))*180.d0/PI 
!            ratio=(ETA(1,j+1,k)-0.5d0)/(ETA(1,j+1,k)-ETA(1,j,k))
!            do i=1,nx
!               p_avg=p_avg+p(i,j,k)*ratio+p(i,j+1,k)*(1.d0-ratio)-p(i,1,1)
!                        w_mid=w(i,j+1,k)*ratio+w(i,j+2,k)*(1.d0-ratio)
!                        v_mid=v(i,j+1,k)*ratio+v(i,j+2,k)*(1.d0-ratio)   
!                        tan_vel_avg=tan_vel_avg+(w_mid*SIN(x_c*PI/180.d0)+v_mid*COS(x_c*PI/180.d0))
!            enddo
!            p_avg=p_avg/(nx*1.d0)
!            tan_vel_avg=tan_vel_avg/(nx*1.d0)
!            write(51,*) REAL(x_c), REAL(2.d0*p_avg), REAL(tan_vel_avg)
!            exit
!        endif
!      enddo
!   enddo


!-----------Cp for cylinder of D = 1 with theta from 0 to 360 (CCW from x+) & pressure reference at inlet (1,1)-------------
!------------------------------------ also calculate near boundary velocity parallel to surface-----------------------------
write(51,*) 'ZONE solutiontime=',  REAL(time)
  !write(51,*) 'ZONE T = "p at <greek>h</greek>=0.5 , p<sub><math>%</math></sub>=(1,1)"'

cell_shift = 1

    do k=kEndVOS,kBgnVOS,-1											! for 0<theta<180 (upper half)
      p_avg=0.d0
      tan_vel_avg=0.d0
	  do j=jBgnVOS,jEndVOS
         if (ETA(1,j,k)>=0.5d0 .AND. ETA(1,j+1,k)<0.5d0) then
            !if (ABS(Zs(k)-z0) .GT. 0.5d0) cycle			
			x_c=ACOS(2.d0*(Zs(k)-z0))*180.d0/PI 
            ratio=(ETA(1,j,k)-0.5d0)/(ETA(1,j,k)-ETA(1,j+1,k))  
			do i=1,nx
               p_avg=p_avg+p(i,j+1+cell_shift,k)*ratio+p(i,j+cell_shift,k)*(1.d0-ratio)-p(i,1,1)
                        w_mid=0.5d0*(w(i,j+1+cell_shift,k)+w(i,j+1+cell_shift,k-1))*ratio+0.5d0*(w(i,j+cell_shift,k)+w(i,j+cell_shift,k-1))*(1.d0-ratio)
                        v_mid=0.5d0*(v(i,j+cell_shift,k)  +v(i,j+1+cell_shift,k))  *ratio+0.5d0*(v(i,j+cell_shift,k)+v(i,j-1+cell_shift,k))*(1.d0-ratio)     
                        tan_vel_avg=tan_vel_avg+(-w_mid*SIN(x_c*PI/180.d0)+v_mid*COS(x_c*PI/180.d0))
            enddo
            p_avg=p_avg/(nx*1.d0)
            tan_vel_avg=tan_vel_avg/(nx*1.d0)
            write(51,*) REAL(x_c), REAL(2.d0*p_avg), REAL(tan_vel_avg)
            exit
         endif
      enddo
    enddo
	
    do k=kBgnVOS,kEndVOS												! for 180<theta<360 (lower half)
      p_avg=0.d0
      tan_vel_avg=0.d0
      do j=jBgnVOS,jEndVOS
         if (ETA(1,j,k)>=0.5d0 .AND. ETA(1,j-1,k)<0.5d0) then
            !if (ABS(Zs(k)-z0) .GT. 0.5d0) cycle	
			x_c=ACOS(-2.d0*(Zs(k)-z0))*180.d0/PI 
            ratio=(ETA(1,j,k)-0.5d0)/(ETA(1,j,k)-ETA(1,j-1,k))      
			do i=1,nx
               p_avg=p_avg+p(i,j-1-cell_shift,k)*ratio+p(i,j-cell_shift,k)*(1.d0-ratio)-p(i,1,1)
                        w_mid=0.5d0*(w(i,j-1-cell_shift,k)+w(i,j-1-cell_shift,k-1))*ratio+0.5d0*(w(i,j-cell_shift,k)+w(i,j-cell_shift,k-1))*(1.d0-ratio)
                        v_mid=0.5d0*(v(i,j-1-cell_shift,k)+v(i,j-2-cell_shift,k))  *ratio+0.5d0*(v(i,j-cell_shift,k)+v(i,j-1-cell_shift,k))*(1.d0-ratio)  
                        tan_vel_avg=tan_vel_avg+(w_mid*SIN(x_c*PI/180.d0)-v_mid*COS(x_c*PI/180.d0))
            enddo
            p_avg=p_avg/(nx*1.d0)
            tan_vel_avg=tan_vel_avg/(nx*1.d0)
            write(51,*) REAL(180.d0+x_c), REAL(2.d0*p_avg), REAL(tan_vel_avg)
            exit
         endif
      enddo
    enddo


!---------------------------------------- Airfoil  ---------------------------------------
   !p_ref=0.d0
   !do j=1,ny; do i=1,nx
   !   p_ref=p_ref+p(i,j,nz)*iDx(i)*iDy(j)
   !enddo; enddo
   !p_ref=p_ref/(lx*ly)

   !write(51,*) 'ZONE solutiontime=',REAL(time)
   !do k=1,nz
   !   p_avg=0.d0
   !   do j=1,ny
         
   !       if (ETA(5,j,k)>=0.5d0 .AND. ETA(5,j+1,k)<0.5d0) then
         
   !         ratio=(ETA(1,j+1,k)-0.5d0)/(ETA(1,j+1,k)-ETA(1,j,k))

   !         x_c= COS((AOA)*PI/180.d0)*(Zs(k)-z0-0.25d0) - SIN((AOA)*PI/180.d0)*(Ys(j)*ratio+Ys(j+1)*(1.d0-ratio)-y0) + 0.25d0  

   !         do i=1,nx
   !            p_avg=p_avg+p(i,j,k)*ratio+p(i,j+1,k)*(1.d0-ratio)
   !         enddo

   !         p_avg=p_avg/(nx*1.d0)-p_ref

   !         write(51,*) REAL(x_c), REAL(2.d0*p_avg)
                  
   !         exit

   !      endif
   !   enddo
   !enddo

   !do k=nz,1,-1
   !   p_avg=0.d0
   !   do j=1,ny
         
   !      if (ETA(5,j,k)<=0.5d0 .AND. ETA(5,j+1,k)>0.5d0) then
         
   !         ratio=(ETA(1,j+1,k)-0.5d0)/(ETA(1,j+1,k)-ETA(1,j,k))

   !         x_c= COS((AOA)*PI/180.d0)*(Zs(k)-z0-0.25d0) - SIN((AOA)*PI/180.d0)*(Ys(j)*ratio+Ys(j+1)*(1.d0-ratio)-y0) + 0.25d0  

   !         do i=1,nx
   !            p_avg=p_avg+p(i,j,k)*ratio+p(i,j+1,k)*(1.d0-ratio)
   !         enddo

   !         p_avg=p_avg/(nx*1.d0)-p_ref

   !         write(51,*) REAL(x_c), REAL(2.d0*p_avg)
                  
   !         exit

   !      endif
   !   enddo
   !enddo

end subroutine filerProcess_cp



! subroutine filerProber()
! use variables
! implicit none
    ! integer 					 	:: n
    ! real*8							:: trilinear_interp
    ! real*8, dimension(1:2,1:2,1:2)  :: u_s,v_s,w_s,p_s

! do n=1,n_points
   ! if (near_pt(3*n) .GE. istart .AND. near_pt(3*n) .LE. iend) then
	! do k=1,2 ; do j=1,2 ; do i=1,2
	! u_s(i,j,k) = 0.5d0*( u(near_pt(3*n-2)  +(i-1),near_pt(3*n-1)+(j-1),near_pt(3*n)+(k-1))+ &
						 ! u(near_pt(3*n-2)-1+(i-1),near_pt(3*n-1)+(j-1),near_pt(3*n)+(k-1)) )

	! v_s(i,j,k) = 0.5d0*( v(near_pt(3*n-2)+(i-1),near_pt(3*n-1)  +(j-1),near_pt(3*n)+(k-1))+ &
						 ! v(near_pt(3*n-2)+(i-1),near_pt(3*n-1)-1+(j-1),near_pt(3*n)+(k-1)) )

	! w_s(i,j,k) = 0.5d0*( w(near_pt(3*n-2)+(i-1),near_pt(3*n-1)+(j-1),near_pt(3*n)  +(k-1))+ &
						 ! w(near_pt(3*n-2)+(i-1),near_pt(3*n-1)+(j-1),near_pt(3*n)-1+(k-1)) )

	! p_s(i,j,k) = p(near_pt(3*n-2)+(i-1),near_pt(3*n-1)+(j-1),near_pt(3*n)+(k-1))

	! enddo; enddo; enddo

	! probe_u(n) = trilinear_interp(probe(3*n-2), probe(3*n-1), probe(3*n), Xs(near_pt(3*n-2)), Xs(near_pt(3*n-2)+1), &
				! Ys(near_pt(3*n-1)), Ys(near_pt(3*n-1)+1), Zs(near_pt(3*n)), Zs(near_pt(3*n)+1), u_s)

	! probe_v(n) = trilinear_interp(probe(3*n-2), probe(3*n-1), probe(3*n), Xs(near_pt(3*n-2)), Xs(near_pt(3*n-2)+1), &
				! Ys(near_pt(3*n-1)), Ys(near_pt(3*n-1)+1), Zs(near_pt(3*n)), Zs(near_pt(3*n)+1), v_s)

	! probe_w(n) = trilinear_interp(probe(3*n-2), probe(3*n-1), probe(3*n), Xs(near_pt(3*n-2)), Xs(near_pt(3*n-2)+1), &
				! Ys(near_pt(3*n-1)), Ys(near_pt(3*n-1)+1), Zs(near_pt(3*n)), Zs(near_pt(3*n)+1), w_s)

	! probe_p(n) = trilinear_interp(probe(3*n-2), probe(3*n-1), probe(3*n), Xs(near_pt(3*n-2)), Xs(near_pt(3*n-2)+1), &
				! Ys(near_pt(3*n-1)), Ys(near_pt(3*n-1)+1), Zs(near_pt(3*n)), Zs(near_pt(3*n)+1), p_s)

	  ! ! Slave processes write directly to file
		! write(filename,'(A,I3.3)')'probe_',n
		! fileformat = '.dat'

		! open (74+n,file=TRIM(filename)//fileformat, status='old', position='append')	
	
		! write(74+n,'(F9.4,4(2x,F10.5))') time,probe_u(n),probe_v(n),probe_w(n),probe_p(n)
		! close(74+n)
   ! endif
        ! call MPI_Barrier(MPI_COMM_WORLD, ierr)

! enddo
! end subroutine filerProber



! function trilinear_interp(x, y, z, x0, x1, y0, y1, z0, z1, f) result(val)
  ! implicit none
  ! real*8, intent(in) :: x, y, z
  ! real*8, intent(in) :: x0, x1, y0, y1, z0, z1
  ! real*8, intent(in) :: f(2,2,2)
  ! real*8 :: val
  ! real*8 :: tx, ty, tz

  ! tx = (x - x0) / (x1 - x0)
  ! ty = (y - y0) / (y1 - y0)
  ! tz = (z - z0) / (z1 - z0)

  ! val = f(1,1,1)*(1.d0-tx)*(1.d0-ty)*(1.d0-tz) + &
        ! f(2,1,1)*tx       *(1.d0-ty)*(1.d0-tz) + &
        ! f(1,2,1)*(1.d0-tx)*ty       *(1.d0-tz) + &
        ! f(2,2,1)*tx       *ty       *(1.d0-tz) + &
        ! f(1,1,2)*(1.d0-tx)*(1.d0-ty)*tz        + &
        ! f(2,1,2)*tx       *(1.d0-ty)*tz        + &
        ! f(1,2,2)*(1.d0-tx)*ty       *tz        + &
        ! f(2,2,2)*tx       *ty       *tz
! end function trilinear_interp



!**************************************************
! Running average
!**************************************************


subroutine running_avg()
use variables
implicit none
real*8  :: uc_,vc_,wc_, tnew, told, FlucS_u, FlucS_v, FlucS_w, FlucS_Sum
!real*8 ,dimension(1:nx,1:ny,1:nz)  :: TKE_u, TKE_v, TKE_w, u_filter_avg, v_filter_avg, w_filter_avg

!$acc data present(TKE(:,:,istart-2:iend+2), &
!$acc	u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2), p(:,:,istart-2:iend+2), &
!$acc	u_TimeAvg(:,:,istart-2:iend+2),v_TimeAvg(:,:,istart-2:iend+2),w_TimeAvg(:,:,istart-2:iend+2),p_TimeAvg(:,:,istart-2:iend+2))


   tnew = 1.d0/(time-startfilerAvg_time)
   told = time-dt-startfilerAvg_time
   !$OMP PARALLEL DO PRIVATE(i,j,uc_,vc_,wc_,FlucS_u,FlucS_v,FlucS_w,FlucS_Sum) collapse(nclps)
   !$acc parallel loop independent private(uc_,vc_,wc_,FlucS_u,FlucS_v,FlucS_w,FlucS_Sum) collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx
      uc_ = 0.5d0*( u(i,j,k)+u(i-1,j,k) )
      vc_ = 0.5d0*( v(i,j,k)+v(i,j-1,k) )
      wc_ = 0.5d0*( w(i,j,k)+w(i,j,k-1) )

	  u_TimeAvg(i,j,k) = (u_TimeAvg(i,j,k) * told + uc_ * dt) * tnew
	  v_TimeAvg(i,j,k) = (v_TimeAvg(i,j,k) * told + vc_ * dt) * tnew
	  w_TimeAvg(i,j,k) = (w_TimeAvg(i,j,k) * told + wc_ * dt) * tnew
	  p_TimeAvg(i,j,k) = (p_TimeAvg(i,j,k) * told + p(i,j,k) * dt) * tnew

	  FlucS_u = (uc_ - u_TimeAvg(i,j,k))*(uc_ - u_TimeAvg(i,j,k))
	  FlucS_v = (vc_ - v_TimeAvg(i,j,k))*(vc_ - v_TimeAvg(i,j,k))
	  FlucS_w = (wc_ - w_TimeAvg(i,j,k))*(wc_ - w_TimeAvg(i,j,k))
	  FlucS_Sum = 0.5d0 * (FlucS_u + FlucS_v + FlucS_w)
	  TKE(i,j,k)= (TKE(i,j,k) * told + FlucS_Sum * dt) * tnew
	  

   enddo; enddo; enddo
   !$acc end parallel
   !$OMP END PARALLEL DO

!$acc end data
  
end subroutine running_avg



subroutine running_avg_full()
use variables
implicit none

real*8  :: uc_,vc_,wc_, tnew, told

!$acc data present(u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2),p(:,:,istart-2:iend+2), &
!$acc	u_TimeAvg(:,:,istart-2:iend+2),v_TimeAvg(:,:,istart-2:iend+2),w_TimeAvg(:,:,istart-2:iend+2), &
!$acc	p_TimeAvg(:,:,istart-2:iend+2),RS_uu(:,:,istart-2:iend+2),RS_uv(:,:,istart-2:iend+2),RS_uw(:,:,istart-2:iend+2), &
!$acc	RS_vv(:,:,istart-2:iend+2),RS_vw(:,:,istart-2:iend+2),RS_ww(:,:,istart-2:iend+2))


   tnew = 1.d0/(time-startfilerAvg_time)
   told = time-dt-startfilerAvg_time
   !$OMP PARALLEL DO PRIVATE(i,j,uc_,vc_,wc_) collapse(nclps)
   !$acc parallel loop independent private(uc_,vc_,wc_) collapse(3) gang vector
    do k=istart,iend
    do j=1,ny; do i=1,nx
      uc_ = 0.5d0*( u(i,j,k)+u(i-1,j,k) )
      vc_ = 0.5d0*( v(i,j,k)+v(i,j-1,k) )
      wc_ = 0.5d0*( w(i,j,k)+w(i,j,k-1) )

	  u_TimeAvg(i,j,k) = (u_TimeAvg(i,j,k) * told + uc_ * dt) * tnew
	  v_TimeAvg(i,j,k) = (v_TimeAvg(i,j,k) * told + vc_ * dt) * tnew
	  w_TimeAvg(i,j,k) = (w_TimeAvg(i,j,k) * told + wc_ * dt) * tnew
	  p_TimeAvg(i,j,k) = (p_TimeAvg(i,j,k) * told + p(i,j,k) * dt) * tnew
	  
	  ! Formula:	<u>_{n+1} = <u>_{n} * (n/(n+1)) + u_{n} * (n/(n+1)^2)
	  RS_uu(i,j,k) = RS_uu(i,j,k)*told*tnew + (u_TimeAvg(i,j,k) - uc_)*(u_TimeAvg(i,j,k) - uc_)*dt*tnew
	  RS_uv(i,j,k) = RS_uv(i,j,k)*told*tnew + (u_TimeAvg(i,j,k) - uc_)*(v_TimeAvg(i,j,k) - vc_)*dt*tnew
	  RS_uw(i,j,k) = RS_uw(i,j,k)*told*tnew + (u_TimeAvg(i,j,k) - uc_)*(w_TimeAvg(i,j,k) - wc_)*dt*tnew
	  RS_vv(i,j,k) = RS_vv(i,j,k)*told*tnew + (v_TimeAvg(i,j,k) - vc_)*(v_TimeAvg(i,j,k) - vc_)*dt*tnew
	  RS_vw(i,j,k) = RS_vw(i,j,k)*told*tnew + (v_TimeAvg(i,j,k) - vc_)*(w_TimeAvg(i,j,k) - wc_)*dt*tnew
	  RS_ww(i,j,k) = RS_ww(i,j,k)*told*tnew + (w_TimeAvg(i,j,k) - wc_)*(w_TimeAvg(i,j,k) - wc_)*dt*tnew
	  

   enddo; enddo; enddo
   !$acc end parallel
   !$OMP END PARALLEL DO

!$acc end data
  
end subroutine running_avg_full


!**************************************************
! Output filer
!**************************************************


subroutine filereachtime3d()
use variables
implicit none

   !$OMP PARALLEL
   !$OMP DO collapse(3)
   !$acc parallel
   !$acc loop independent collapse(3)
   do k=k_start,k_end; do j=j_start,j_end; do i=i_start,i_end
      uc(i,j,k) = 0.5d0*( u(i,j,k)+u(i-1,j,k) )
      vc(i,j,k) = 0.5d0*( v(i,j,k)+v(i,j-1,k) )
      wc(i,j,k) = 0.5d0*( w(i,j,k)+w(i,j,k-1) )
   enddo; enddo; enddo
   !$OMP END DO

   !$OMP DO collapse(3)
   !$acc loop independent collapse(3)
   do k=k_start,k_end; do j=j_start,j_end; do i=i_start,i_end
      Qout(i,j,k,1) = p(i,j,k)
      Qout(i,j,k,2) = uc(i,j,k)
      Qout(i,j,k,3) = vc(i,j,k)
      Qout(i,j,k,4) = wc(i,j,k)
      Qout(i,j,k,5) = ETA(i,j,k)
   enddo; enddo; enddo
   !$acc end parallel
   !$OMP END DO   
   !$OMP END PARALLEL


   write(filename,'(A,I6.6)')'3d_',NINT(time/filer3d_int)
   fileformat = '.q'

   open (17,file=TRIM('output3D/')//TRIM(filename)//fileformat,form='unformatted')
   write(17) nblocks
   write(17) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1
   write(17) temp, temp, temp, REAL(time)
   write(17) ( ( ( ( Qout(i,j,k,h), i = i_start, i_end), j = j_start, j_end), k = k_start, k_end), h = 1, 5 )

   close(17)

end subroutine filereachtime3d

subroutine filereachtime2d_Avg() ! use this for 2D spanwise (Z axis) averaged
use variables
implicit none

   !$OMP PARALLEL
   !$OMP DO PRIVATE(i,j) collapse(nclps)
   !$acc parallel loop independent collapse(3)
   do k=1,nz; do j=1,ny; do i=1,nx
      uc(i,j,k) = 0.5d0*( u(i,j,k)+u(i-1,j,k) )
      vc(i,j,k) = 0.5d0*( v(i,j,k)+v(i,j-1,k) )
      wc(i,j,k) = 0.5d0*( w(i,j,k)+w(i,j,k-1) )
   enddo; enddo; enddo
   !$acc end parallel
   !$OMP END DO
   

   !$OMP DO PRIVATE(i,j) reduction(+: p_spansum,uc_spansum,vc_spansum,wc_spansum) collapse(nclps)
   !$acc parallel
   !$acc loop independent reduction(+: p_spansum,uc_spansum,vc_spansum,wc_spansum)
   do k=1,nz
   !$acc loop seq
   do j=1,ny
	p_spansum = 0.d0; uc_spansum = 0.d0; vc_spansum = 0.d0; wc_spansum = 0.d0
   !$acc loop seq
   do i=1,nx
	   p_spansum =  p_spansum +  p(i,j,k)
	  uc_spansum = uc_spansum + uc(i,j,k) 
	  vc_spansum = vc_spansum + vc(i,j,k)
	  wc_spansum = wc_spansum + wc(i,j,k)
   enddo
	 p_spanavg(1,j,k) =  p_spansum/nx*1.d0
    uc_spanavg(1,j,k) = uc_spansum/nx*1.d0
    vc_spanavg(1,j,k) = vc_spansum/nx*1.d0
	wc_spanavg(1,j,k) = wc_spansum/nx*1.d0
   enddo; enddo
   !$acc end parallel
   !$OMP END DO
   !$OMP END PARALLEL
   

   !$OMP PARALLEL DO PRIVATE(j) collapse(nclps)
   !$acc parallel loop independent collapse(2)
   do k=1,nz; do j=1,ny
      Qout(1,j,k,1) =  p_spanavg(1,j,k)
      Qout(1,j,k,2) = uc_spanavg(1,j,k)
      Qout(1,j,k,3) = vc_spanavg(1,j,k)
      Qout(1,j,k,4) = wc_spanavg(1,j,k)
      Qout(1,j,k,5) =        ETA(1,j,k)
   enddo; enddo
   !$acc end parallel   
    !$OMP END PARALLEL DO


   write(filename,'(A,I6.6)')'2d_',NINT(time/filer2d_int)
   fileformat = '.q'
   
   open (20,file=TRIM('output2D/')//TRIM(filename)//fileformat,form='unformatted')
   write(20) nblocks
   write(20) 1, ny, nz
   write(20) temp, temp, temp, REAL(time)
   write(20) ( ( ( ( Qout(i,j,k,h), i = 1, 1), j = 1, ny), k = 1, nz), h = 1, 5 )

   close(20)

end subroutine filereachtime2d_Avg


subroutine filereachtime2d_Mid()	! use this for 2D slice at the middle of Z plane
use variables
implicit none

i = INT((nx+1)/2.)

   !$OMP PARALLEL
   !$OMP DO PRIVATE(j) collapse(nclps) 
   !$acc parallel loop independent collapse(2)
   do k=1,nz; do j=1,ny
      uc(i,j,k) = 0.5d0*( u(i,j,k)+u(i-1,j,k) )
      vc(i,j,k) = 0.5d0*( v(i,j,k)+v(i,j-1,k) )
      wc(i,j,k) = 0.5d0*( w(i,j,k)+w(i,j,k-1) )
   enddo; enddo
   !$acc end parallel
   !$OMP END DO
   

   !$OMP DO PRIVATE(j) collapse(nclps) 
   !$acc parallel loop independent collapse(2)
   do k=1,nz; do j=1,ny
      Qout(1,j,k,1) = p(i,j,k)
      Qout(1,j,k,2) = uc(i,j,k)
      Qout(1,j,k,3) = vc(i,j,k)
      Qout(1,j,k,4) = wc(i,j,k)
      Qout(1,j,k,5) = ETA(i,j,k)
   enddo; enddo
   !$acc end parallel   
   !$OMP END DO
   !$OMP END PARALLEL


   write(filename,'(A,I6.6)')'2d_',NINT(time/filer2d_int)
   fileformat = '.q'
   
   open (20,file=TRIM('output2D/')//TRIM(filename)//fileformat,form='unformatted')
   write(20) nblocks
   write(20) 1, ny, nz
   write(20) temp, temp, temp, REAL(time)
   write(20) ( ( ( ( Qout(i,j,k,h), i = 1, 1), j = 1, ny), k = 1, nz), h = 1, 5 )

   close(20)

end subroutine filereachtime2d_Mid


subroutine filereachtime2d_wall()	! use this for 2D slice at the middle of Z plane
use variables
implicit none

i = INT((nx+1)/2.)
  
   !$OMP PARALLEL
   !$OMP DO PRIVATE(j) collapse(nclps) 
   !$acc parallel loop independent collapse(2)
   do k=1,nz; do j=1,ny
      Qout(1,j,k,1) = dist_print(i,j,k)
      Qout(1,j,k,2) = yplus_print(i,j,k)
      Qout(1,j,k,3) = uplus_print(i,j,k)
      Qout(1,j,k,4) = nut(i,j,k)
      Qout(1,j,k,5) = ETA(i,j,k)
   enddo; enddo
   !$acc end parallel   
   !$OMP END DO
   !$OMP END PARALLEL


   write(filename,'(A,I6.6)')'wall_',NINT(time/filer2d_int)
   fileformat = '.q'
   
   open (22,file=TRIM('output2D/')//TRIM(filename)//fileformat,form='unformatted')
   write(22) nblocks
   write(22) 1, ny, nz
   write(22) temp, temp, temp, REAL(time)
   write(22) ( ( ( ( Qout(i,j,k,h), i = 1, 1), j = 1, ny), k = 1, nz), h = 1, 5 )

   close(22)

end subroutine filereachtime2d_wall


subroutine filereachtime3d_VF()	! filer for virtual force field
use variables
implicit none

    !$OMP PARALLEL
    !$OMP DO collapse(3)
	!$acc parallel loop independent collapse(3)
	do k=k_start,k_end; do j=j_start,j_end; do i=i_start,i_end
      Qout(i,j,k,1) = p(i,j,k)
      Qout(i,j,k,2) = -FX(i,j,k)
      Qout(i,j,k,3) = -FY(i,j,k)
      Qout(i,j,k,4) = -FZ(i,j,k)
      Qout(i,j,k,5) = ETA(i,j,k)
	enddo; enddo; enddo
	!$acc end parallel
	!$OMP END DO   
	!$OMP END PARALLEL

	write(filename,'(A,I6.6)')'VF_',NINT(time/filerVF_int)
	fileformat = '.q'
   
	open (59,file=TRIM('output3D/')//TRIM(filename)//fileformat,form='unformatted')
	write(59) nblocks
	write(59) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1
	write(59) temp, temp, temp, REAL(time)
	write(59) ( ( ( ( Qout(i,j,k,h), i = i_start, i_end), j = j_start, j_end), k = k_start, k_end), h = 1, 5 )

	close(59)

end subroutine filereachtime3d_VF


subroutine filer3d_avg()
use variables
implicit none

    !$OMP PARALLEL
    !$OMP DO collapse(3)
	!$acc parallel loop independent collapse(3)
	do k=k_start,k_end; do j=j_start,j_end; do i=i_start,i_end
		Qout(i,j,k,1) = p_TimeAvg(i,j,k)
		Qout(i,j,k,2) = u_TimeAvg(i,j,k)
		Qout(i,j,k,3) = v_TimeAvg(i,j,k)
		Qout(i,j,k,4) = w_TimeAvg(i,j,k)
		Qout(i,j,k,5) = TKE(i,j,k)
	enddo; enddo; enddo
	!$acc end parallel
	!$OMP END DO   
	!$OMP END PARALLEL

	write(filename,'(A,I4.4)')'Avg_',NINT(time/filerAvg_int)
	fileformat = '.q'
   
	open (60,file=TRIM('output3D/')//TRIM(filename)//fileformat,form='unformatted')
	write(60) nblocks
	write(60) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1
	write(60) temp, temp, temp, REAL(time)
	write(60) ( ( ( ( Qout(i,j,k,h), i = i_start, i_end), j = j_start, j_end), k = k_start, k_end), h = 1, 5 )

	close(60)

end subroutine filer3d_avg



subroutine filer3d_avg_full1()
use variables
implicit none

    !$OMP PARALLEL
    !$OMP DO collapse(3)
	!$acc parallel loop independent collapse(3)
	do k=k_start,k_end; do j=j_start,j_end; do i=i_start,i_end
		Qout(i,j,k,1) = p_TimeAvg(i,j,k)
		Qout(i,j,k,2) = u_TimeAvg(i,j,k)
		Qout(i,j,k,3) = v_TimeAvg(i,j,k)
		Qout(i,j,k,4) = w_TimeAvg(i,j,k)
		Qout(i,j,k,5) = RS_ww(i,j,k)
	enddo; enddo; enddo
	!$acc end parallel
	!$OMP END DO   
	!$OMP END PARALLEL

	write(filename,'(A,I4.4)')'Avg1_',NINT(time/filerAvg_int)
	fileformat = '.q'
   
	open (60,file=TRIM('output3D/')//TRIM(filename)//fileformat,form='unformatted')
	write(60) nblocks
	write(60) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1
	write(60) temp, temp, temp, REAL(time)
	write(60) ( ( ( ( Qout(i,j,k,h), i = i_start, i_end), j = j_start, j_end), k = k_start, k_end), h = 1, 5 )

	close(60)

end subroutine filer3d_avg_full1


subroutine filer3d_avg_full2()
use variables
implicit none

    !$OMP PARALLEL
    !$OMP DO collapse(3)
	!$acc parallel loop independent collapse(3)
	do k=k_start,k_end; do j=j_start,j_end; do i=i_start,i_end
		Qout(i,j,k,1) = RS_uu(i,j,k)
		Qout(i,j,k,2) = RS_uv(i,j,k)
		Qout(i,j,k,3) = RS_uw(i,j,k)
		Qout(i,j,k,4) = RS_vv(i,j,k)
		Qout(i,j,k,5) = RS_vw(i,j,k)
	enddo; enddo; enddo
	!$acc end parallel
	!$OMP END DO   
	!$OMP END PARALLEL

	write(filename,'(A,I4.4)')'Avg2_',NINT(time/filerAvg_int)
	fileformat = '.q'
   
	open (60,file=TRIM('output3D/')//TRIM(filename)//fileformat,form='unformatted')
	write(60) nblocks
	write(60) i_end-i_start+1, j_end-j_start+1, k_end-k_start+1
	write(60) temp, temp, temp, REAL(time)
	write(60) ( ( ( ( Qout(i,j,k,h), i = i_start, i_end), j = j_start, j_end), k = k_start, k_end), h = 1, 5 )

	close(60)

end subroutine filer3d_avg_full2


!***************************************************
! Backup filer (make consistent with read_data.f90)
!***************************************************


subroutine backupfile()
use variables
implicit none

   !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)  
   !$acc parallel loop independent collapse(3) 
   do k=1,nz; do j=1,ny; do i=1,nx
      Qout(i,j,k,1) = p(i,j,k)
      Qout(i,j,k,2) = u(i,j,k)
      Qout(i,j,k,3) = v(i,j,k)
      Qout(i,j,k,4) = w(i,j,k)
      Qout(i,j,k,5) = ETA(i,j,k)
   enddo; enddo; enddo
   !$acc end parallel   
   !$OMP END PARALLEL DO


   if(solid_motion == 0 .OR. time < StartDynamic_time+dt) then
      open (unit=18,form='unformatted',file='STATIC.bak')
   else 
      open (unit=18,form='unformatted',file='DYNAMIC.bak')
   endif
   write(18) nblocks
   write(18) nx, ny, nz
   write(18) temp, temp, temp, REAL(time)
   write(18) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
   write(18) time,dt
   write(18) totalFX_,totalFY_,totalFZ_
   write(18) totalTorqx,totalTorqy,totalTorqz,totalTorq   
   write(18) u_solid,v_solid,w_solid,rotor_omega
   write(18) rotate_sx,rotate_sy,rotate_sz
   write(18) x0_t,y0_t,z0_t,AOA
   close(18)

end subroutine backupfile


subroutine backupAvg()
use variables
implicit none

   !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps) 
   !$acc parallel loop independent collapse(3) 
   do k=1,nz; do j=1,ny; do i=1,nx
      Qout(i,j,k,1) = p_TimeAvg(i,j,k)
      Qout(i,j,k,2) = u_TimeAvg(i,j,k)
      Qout(i,j,k,3) = v_TimeAvg(i,j,k)
      Qout(i,j,k,4) = w_TimeAvg(i,j,k)
      Qout(i,j,k,5) = TKE(i,j,k)
   enddo; enddo; enddo
   !$acc end parallel   
   !$OMP END PARALLEL DO


   open (unit=19,form='unformatted',file='RUNNING_AVG.bak')
   write(19) nblocks
   write(19) nx, ny, nz
   write(19) temp, temp, temp, REAL(time)
   write(19) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
   close(19)

end subroutine backupAvg


subroutine backup_avg_full1()
use variables
implicit none

   !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps) 
   !$acc parallel loop independent collapse(3) 
   do k=1,nz; do j=1,ny; do i=1,nx
      Qout(i,j,k,1) = p_TimeAvg(i,j,k)
      Qout(i,j,k,2) = u_TimeAvg(i,j,k)
      Qout(i,j,k,3) = v_TimeAvg(i,j,k)
      Qout(i,j,k,4) = w_TimeAvg(i,j,k)
      Qout(i,j,k,5) = RS_ww(i,j,k)
   enddo; enddo; enddo
   !$acc end parallel   
   !$OMP END PARALLEL DO


   open (unit=23,form='unformatted',file='RUNNING_AVG1.bak')
   write(23) nblocks
   write(23) nx, ny, nz
   write(23) temp, temp, temp, REAL(time)
   write(23) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
   close(23)

end subroutine backup_avg_full1

subroutine backup_avg_full2()
use variables
implicit none

   !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps) 
   !$acc parallel loop independent collapse(3) 
   do k=1,nz; do j=1,ny; do i=1,nx
	  Qout(i,j,k,1) = RS_uu(i,j,k)
	  Qout(i,j,k,2) = RS_uv(i,j,k)
	  Qout(i,j,k,3) = RS_uw(i,j,k)
	  Qout(i,j,k,4) = RS_vv(i,j,k)
	  Qout(i,j,k,5) = RS_vw(i,j,k)
   enddo; enddo; enddo
   !$acc end parallel   
   !$OMP END PARALLEL DO


   open (unit=24,form='unformatted',file='RUNNING_AVG2.bak')
   write(24) nblocks
   write(24) nx, ny, nz
   write(24) temp, temp, temp, REAL(time)
   write(24) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
   close(24)

end subroutine backup_avg_full2


!subroutine filer_bodyforce()
!use variables
!implicit none

!   !$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)  
!   do k=1,nz; do j=1,ny; do i=1,nx
!      Qout(i,j,k,1)=edelta(i,j,k)
!      Qout(i,j,k,2)=0.d0
!      Qout(i,j,k,3)=F_tavey(i,j,k)
!      Qout(i,j,k,4)=F_tavex(i,j,k)
!      Qout(i,j,k,5)=ETA(i,j,k)
!   enddo; enddo; enddo
!   !$OMP END PARALLEL DO


   !open (unit=18,form='unformatted',file='PLASMA.Q')
   !write(filename,'(I4.4)')num
   !fileformat = 'plasma.q'

   !open (18,file=TRIM(filename)//fileformat,position='append',form='unformatted')
!   open (unit=18,form='unformatted',file='PLASMA.q')
!   write(18) nblocks
!   write(18) nx, ny, nz
!   write(18) temp, temp, temp, REAL(time)
!   write(18) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )

!   close(18)


!end subroutine filer_bodyforce


