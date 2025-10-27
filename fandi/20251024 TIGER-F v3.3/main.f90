! 12 Jan 2025 - FDS
!    y=1 ______________                                                                                 
!       /             /|       |  Author  : Zi-Hsuan Wei                                                 
!      /             / |       |  Version : 3.3                                                          
!     /____________ /  |       |  Mod by  : Fandi D. Suprianto                              
!     |  |         |   |       |  Web     : https://mestaiwantech.com                                    
!     |  |         |   |                                          
!     |  | x=y=z=0 |   |                                           
!     |  |_________|___|x=1                                        
!     |  /         |  /                                         
!     | /          | /                                        
!     |/___________|/                                         
!    z=1                                                  
!                                                       
program main
   use variables
   use mpi
   use omp_lib
   use openacc ! openACC library
   !use nvtx    ! Nsight profiler library
   implicit none



   nthreads = 4 !(OpenMP or OpenACC)
   
   !---------------------------- OPENMP ----------------------------!
   call omp_set_num_threads(nthreads)
   ! Use only physical cores for best performance
   ! Disable hyperthreading (if any) using 'export OMP_PROC_BIND=TRUE')
   !----------------------------------------------------------------!


   !----------------------------- MPI ------------------------------!	*** MPI_division.f90 ***
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
   call Mpi_division()					! Equal load by default
   !call Mpi_division_unbalance2gpu()	! Unbalanced load for 2 different GPUs
   !----------------------------------------------------------------!
   
   !--------for nz--------!
   istart = gstart(myid)  !
   iend = gend0(myid)     !
   igcount = gcount(myid) !
   !----------------------!

   !----------------------------for sendrecv------------------------!
   l_nbr = myid - 1                                                 !
   r_nbr = myid + 1                                                 !         
   if(myid == 0) then; l_nbr=MPI_PROC_NULL; endif                   !
   if(myid == (nproc-1)) then; r_nbr=MPI_PROC_NULL; endif           !
   !----------------------------------------------------------------!


   !---------------------------- OPENACC ---------------------------!
   call acc_set_num_cores(nthreads)
   call acc_set_device_num(myid,acc_device_nvidia)
   !----------------------------------------------------------------!


   !------------------------ Toggle Selector -----------------------!
   
   dimensionality					= 0						! 0) Non-Dimensional, 1) Dimensional

   LES_mode                         = 'ON'                  ! the LES mode.

   wall_model						= 'OFF'					! Apply wall model.	(under development)
   
   Damping_F						= 'OFF'					! Van Driest damping function.	(under development)
   
   pressure_solver					= 1						! 1)RB_SOR, 2)SOR-CPU
   
   correction_stage					= 0						! number of correction stage in Prediction-correction process
		StartCorrection_time		= 0.d0					! initialize time of prediction-correction

   resume							= 'OFF'					! resume simulation. (Provide the backup file: 'STATIC.bak'/'DYNAMIC.bak')
   

   VOS_by							= 1						! Create VOS by 1)geometry function, 2)raycasting2D, 3)raycasting3D

   dat_file	                     	= '1_0.1_30_w.DAT'   	! input for raycasting2D (DAT file)
						z_length2D 	= 1.d0 					! scaling the geometry size with polygon's origin as the base point, Enter the targeted z_length

   stl_file                         = 'sphere_ratio_005.stl'! input for raycasting3D (STL ASCII file)
						z_length3D 	= 1.d0 					! scaling the geometry size with CAD's origin as the base point, Enter the targeted z_length
      
   select_ref_area					= 1						! Reference area in virtual force calculation
					user_defined	= L_ch*lx				! 1)Z_castingarea, 2)Y_castingarea, 3)X_castingarea, or 4)define manually

   solid_motion						= 0						! Select dynamic coupling: 0) Stationary, 1) oneway_coupling, 2) twoway_coupling
		StartDynamic_time           = 0.d0 					! initialize time of dynamic coupling (if solid_motion = 1 or 2)

!   DBD                              = 0                     ! DBD actuator on : 1 ; DBD actuator off : 0

!   unDBD                            = 0                     ! Strouhal number of brust DBD actuator

   !----------------------------------------------------------------!
   

   !---------------- Parameters for the simulation -----------------!
   
   zeta                             = 1.0d-4                ! zeta for solving pressure matrix. 
!															  SOR with a dynamic omega is used. Set omega bounds in SOR.f90 when necessary. (default = 1.1-1.9)
   itmax                            = 5000                  ! maximum for zeta in Gauss Seidel subroutines

   initial_time                     = 0.d0                  ! initialize time of simulation

   total_time                       = 300.d0	 			! total time of simulation

   time_stepping_by					= 'CFL'					! 'DT' = constant dt, 'CFL' = constant CFL
					constant_dt		= 3.d-3					! input dt if time_stepping_by a constant dt
					constant_cfl	= 0.2d0					! input cfl if time_stepping_by a constant cfl (use 0.2 if using 2nd order AB)

   !----------------------------------------------------------------!


   !---------------------Output Data Management---------------------!

   coeff_start						= 0.1d0					! initialize time of storing force & pressure coefficients (CD_time & Cp)
   
   filer_cp							= 'OFF'					! create surface pressure distribution plot. (set in filer.f90)
			filerCp_int				= 5.d0					! time interval to write cp data (start = coeff_start)
			
   filer3d							= 'OFF'					! create 3D output (3d_xxxx.q). 
			startfiler3d_time		= 0.d0					! initialize time of writing filer3d data
			filer3d_int				= 1.d0					! time interval to write filer3d data
   
   filer2d							= 'ON'					! create 2D output (2d_xxxx.q). 'AVG'/'ON' : Spanwise avg ; 'MID' : Spanwise Midplane slice
			startfiler2d_time		= 0.d0					! initialize time of writing filer2d data
			filer2d_int				= 1.d0					! time interval to write filer2d data

   filer_virtual_force				= 'OFF'					! create 3D output of virtual force field.
			startfilerVF_time		= 0.d0					! initialize time of writing virtual force files
			filerVF_int				= 1.d0   				! time interval to write virtual force files
   
   filer_running_avg				= 'OFF'					! Calculate running average. 'STANDARD'/'ON' : uvw+TKE ; 'FULL' : uvwp+Reynold stress
			startfilerAvg_time		= 0.d0					! initialize time of writing average data
			filerAvg_int			= 100.d0   				! time interval to write average data (the last timestep is enough in most cases)

   backup_int						= 1.d0					! time interval to write backup data

   
   !----------------------------------------------------------------!  
   

   totalstarttime = MPI_WTIME()
   totallasttime = totalstarttime
   lastprogresstime = totalstarttime
  
   !----------------------------------------------------------------!	*** read_data.f90 ***
   call read_data()              ! define dx, dy, dz
   !----------------------------------------------------------------!

  
   !--------------------- Gridder selection ------------------------!	*** gridder.f90 ***
   if(Gridder=='non-uniform-sin4')then
      call gridder_sin4()              ! define unequal mesh proposed by Kuyper
!   else if(Gridder=='non-uniform-sin4-3D')then
!      call gridder_sin4_3D()           ! define unequal mesh proposed by Kuyper in x,y,z
  else if(Gridder=='uniform')then
     call gridder_equal()             ! define equal mesh
!   else if(Gridder=='ground-3D')then
!      call gridder_ground()  		   ! define unequal mesh with two sine domain in x,y,z for object lies on the bottom wall
!   else if(Gridder=='non-uniform-sin5-3D')then
!      call gridder_sin5_3D()           ! define unequal mesh with two sine domain in x,y,z
!   else if(Gridder=='non-uniform-sin5')then
!      call gridder_sin5()              ! define unequal mesh with two sine domain in x,y
  else if(Gridder=='non-uniform-fall')then
     call gridder_fall()
  else if(Gridder=='non-uniform-fall-wall')then
     call gridder_fall_wall()
!   else if(Gridder=='gridder_flatplate')then
!	  call gridder_triuniform()
!   else if(Gridder=='non-uniform')then
!      call gridder_unequal()          ! define unequal mesh 
!   else if(Gridder=='non-uniform-sin3')then
!      call gridder_sin3()             ! define unequal mesh using simple SIN function from 0 to pi/2
!   else if(Gridder=='non-uniform-sin2')then
!      call gridder_sin2()             ! define unequal mesh
!   else if(Gridder=='non-uniform-sin' )then
!      call gridder_sin()              ! define unequal mesh
  else if(Gridder=='read_txt' )then
     call gridder_read_txt()          ! grid coordinate is supplied (grid.txt) in X(i),Y(j),Z(k) format
   end if
   !----------------------------------------------------------------!

   !----------------------------------------------------------------!	*** initial_conditions.f90 ***
    if (time_stepping_by=='DT') then
		dt = constant_dt
	else if(time_stepping_by=='CFL') then
		dt = constant_cfl * lzSml/ nzSml / U_inf	!initial dt. dt is updated in calcul_new_velocity.f90
	endif
   
   call initial_conditions()        !call initial conditions
   !----------------------------------------------------------------!
   
   !-------------------- Loading backup data -----------------------!	*** read_data.f90 ***
   if (resume=='OFF') then
	   time = initial_time
   else if (resume=='ON') then													
	  call read_backupfile()        				!input first .bak file			--> (used when resuming simulation)				
	  initial_time = lastsaved_time
	  time = lastsaved_time
   endif
   !----------------------------------------------------------------!


   !--------------------- Intervals handler ------------------------!
	istep = 0
    prog_int = 0.1d0*(total_time-initial_time)		! 10% time interval to calculate simulation progess
	next_progress= initial_time+prog_int

	next_filer3d  = MAX(INT(initial_time/filer3d_int)*filer3d_int + filer3d_int , startfiler3d_time)
	next_filer2d  = MAX(INT(initial_time/filer2d_int)*filer2d_int + filer2d_int, startfiler2d_time)
	next_filerVF  = MAX(INT(initial_time/filerVF_int)*filerVF_int + filerVF_int, startfilerVF_time)
	next_filerAvg = MAX(INT(initial_time/filerAvg_int)*filerAvg_int + filerAvg_int, startfilerAvg_time)
	next_filerWall= MAX(INT(initial_time/(10.*filer2d_int))*(10.*filer2d_int) + (10.*filer2d_int), initial_time)
	next_filerCp  = INT(initial_time/filerCp_int)*filerCp_int + filerCp_int	
	next_backup   = INT(initial_time/backup_int) *backup_int  + backup_int
	next_backupAvg= next_backup
   !----------------------------------------------------------------!


   !------------------ Loading avg backup data ---------------------!	*** filer.f90 ***
   if ((filer_running_avg=='STANDARD' .OR. filer_running_avg=='ON') .AND. &
		resume=='ON' .AND. time > startfilerAvg_time) then												
	  call read_backup_Avg()						!input first .bak file			--> (used when resuming simulation)				
   elseif (filer_running_avg=='FULL' .AND. resume=='ON' .AND. time > startfilerAvg_time) then
	  call read_backup_Avg_full1()
	  call read_backup_Avg_full2()
   endif
   !----------------------------------------------------------------!

!$acc data copy(X,Y,Z,Xs,Ys,Zs,iDx,iDy,iDz,ETA,p(:,:,istart-2:iend+2), &
!$acc	u0(:,:,istart-2:iend+2),v0(:,:,istart-2:iend+2),w0(:,:,istart-2:iend+2), &
!$acc 	 u(:,:,istart-2:iend+2), v(:,:,istart-2:iend+2), w(:,:,istart-2:iend+2))
   !----------------------------------------------------------------!	*** boundary_conditions.f90 ***
   call final_boundary_conditions()       !call boundary conditions
   !----------------------------------------------------------------!


   !------------------ Create initial VOS --------------------------!

   if (VOS_by == 2) then
		call read_dat()						  ! Read DAT file					*** vos_ray2d.f90 ***
vosstarttime = MPI_WTIME()
!call nvtxStartRange("VOS")
		call vos_ray2d()                      ! Create ETA from DAT file 		*** vos_ray2d.f90 ***
   elseif (VOS_by == 3) then
		call read_stl()						  ! Read STL file					*** vos_ray3d.f90 ***
vosstarttime = MPI_WTIME()
!call nvtxStartRange("VOS")
		call vos_ray3d()                      ! Create ETA from STL ASCII file 	*** vos_ray3d.f90 ***
   elseif (VOS_by == 1) then	
vosstarttime = MPI_WTIME()
!call nvtxStartRange("VOS")														*** vos_function.f90 ***
		call func_Cylinder()                  ! Create ETA of Cylinder  
!		call func_Cylplate_noblade()          ! Create ETA of Cylinder & flat plate
!		call func_Cylplate_1blade()           ! Create ETA of Cylinder & flat plate
!		call func_Cylplate_2blade()           ! Create ETA of Cylinder & flat plate    
!		call func_Cylplate_3blade()           ! Create ETA of Cylinder & flat plate 
!		call func_Sphere()                    ! Create ETA of Sphere
!		call airfoil_00xx()
!       call func_darrius_3blade
!		call square_cyl()
!		call flat_plate()
!		call func_torus3d
!		call func_rbc()
!		call func_64Spheres()
   endif
vosfinaltime = MPI_WTIME()
!call nvtxEndRange

!$acc end data


		!-----------------------calculating casting area----------------------!
		Z_castingarea=0.0
		do i=1,nx; do j=1,ny
		ETAs=0.0
		do k=1,nz 
			if ( ETA(i,j,k)>ETAs ) then
				ETAs=ETA(i,j,k)
			endif 
		enddo
		Z_castingarea=Z_castingarea+ETAs*iDx(i)*iDy(j)
		enddo; enddo

		Y_castingarea=0.0
		do i=1,nx; do k=1,nz
		ETAs=0.0
		do j=1,ny 
			if ( ETA(i,j,k)>ETAs ) then
				ETAs=ETA(i,j,k)
			endif 
		enddo
		Y_castingarea=Y_castingarea+ETAs*iDx(i)*iDz(k)
		enddo; enddo
	
		X_castingarea=0.0
		do j=1,ny; do k=1,nz
		ETAs=0.0
		do i=1,nx 
			if ( ETA(i,j,k)>ETAs ) then
				ETAs=ETA(i,j,k)
			endif 
		enddo
		X_castingarea=X_castingarea+ETAs*iDy(j)*iDz(k)
		enddo; enddo

		Solid_volume=0.0
		do k=1,nz; do j=1,ny; do i=1,nx
			solid_volume=solid_volume + ETA(i,j,k)*iDx(i)*iDy(j)*iDz(k)
		enddo; enddo; enddo


		if (select_ref_area == 1) then						! Reference area in virtual force calculation
			ref_area =  Z_castingarea
		elseif (select_ref_area == 2) then		
			ref_area =  Y_castingarea
		elseif (select_ref_area == 3) then		
			ref_area =  X_castingarea
		elseif (select_ref_area == 4) then		
			ref_area =  user_defined
		endif

	!vosfinaltime = MPI_WTIME()
		!---------------------------------------------------------------------!


	!------------calculating initial wall distance------------------!		*** wall_model.f90 *** (under development)
	! if (wall_model=='ON' .OR. Damping_F=='ON') then
		! !$acc data copy(Y,Ys,Zs,Dys,ETA,wdist,wsin,wcos)
			! call wall_distance()    ! calculate wall distance
		! !$acc end data
	! endif
   !----------------------------------------------------------------!   
   
   
   
   !------------------define plasma force field---------------------!
   !i=1
   !if(DBD == 1)then
   !   AOA = AOA1
   !   call plasma()
   !   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !do k=1,nz; do j=1,ny
      !   if(edelta(i,j,k)==1 .AND. myid==master)then
      !      write(*,*) F_tavex(i,j,k), F_tavey(i,j,k), edelta(i,j,k)
      !   end if
      !end do; end do

   !end if
   !----------------------------------------------------------------!


   !--------------------------Probe points--------------------------!		*** input points in variables_module.f90 ***
	! if(myid==master) then
		! call prober()												!		*** initial_conditions.f90 ***
	! endif
	! call MPI_BCAST(near_pt, 3*n_points, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   !----------------------------------------------------------------!
   
   
	!-----------------write initial output files--------------------!
	if(myid==master) then
      call acc_set_device_type(acc_device_host)
	  if((resume=='OFF') .OR. &
		 (filer3d=='ON' .AND. time > startfiler3d_time-dt)) then
		call filereachtime3d() 										!		*** filer.f90 ***
	  endif
	  
	  if((filer2d=='AVG' .OR. filer2d=='ON') .AND. &
		  time > startfiler2d_time-dt) then
		call filereachtime2d_Avg()				! write filer 2D			*** filer.f90 ***
	  elseif (filer2d=='MID' .AND. &
			  time > startfiler2d_time-dt) then
		call filereachtime2d_Mid()				! write filer 2D			*** filer.f90 ***
	  endif

	  if((filer_running_avg=='STANDARD' .OR. filer_running_avg=='ON') .AND. &
		  time > startfilerAvg_time-dt) then
		call filer3d_avg()						! write filer3d_avg			*** filer.f90 ***
	  elseif (filer_running_avg=='FULL' .AND. &
			  time > startfilerAvg_time-dt) then
		call filer3d_avg_full1()				! write filer3d_avg_full	*** filer.f90 ***
		call filer3d_avg_full2()
	  endif

      call filerInfo()												!		*** filer.f90 ***
          !---------------------
          !call filerProcess_cp() 
	  call acc_set_device_num(myid,acc_device_nvidia)
	endif
   !----------------------------------------------------------------!

   !call write_matrix_a()

    if(myid==master) then
	open (71,file='output.dat',position='append')
	write(71, '(/)')
	open (73,file='information.dat',position='append',form='formatted',access='stream')
    write(73,*) ' ' 
    write(73,*) '===============Simulation progress & Wall time==============='
	
	call date_and_time(date=date_str, time=time_str)
	write(73,'(1X,A4,"-",A2,"-",A2," ",A2,":",A2,":",A2,A)') &
		date_str(1:4), date_str(5:6), date_str(7:8), &
		time_str(1:2), time_str(3:4), time_str(5:6), ' Simulation started'
	flush(73)
	endif

	!open (72,file='vos_time.dat',position='append')

    do i = 0, nproc - 1
        if (myid == i) then
            open(71, file='output.dat', status='old', position='append')
            write(71, '(4(A,I5,3x))') 'Proc  =', myid, 'istart =', istart, 'iend =', iend, 'gcount =', igcount
			close(71)
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do

	open (71,file='output.dat',position='append')


    if(myid==master) then
		write(71,'(A,I5,3x)') ''
		write(71,'(A,I5,3x)') 'threads/Proc  =',nthreads
		write(71,800)
800 	FORMAT(63('-'), /)
	endif


	call MPI_BARRIER(MPI_COMM_WORLD, ierr) 



!$acc data copyin(X,Y,Z,Xs,Ys,Zs,iDx,iDy,iDz,Dxs,Dys,Dzs,ETA,p(:,:,istart-2:iend+2), &
!$acc	  P_Den(:,:,istart:iend,:),  Den_inv(:,:,istart:iend),     	  nut(:,:,istart:iend+1)    , &
!$acc	     u0(:,:,istart-2:iend+2),     v0(:,:,istart-2:iend+2),     w0(:,:,istart-2:iend+2), &
!$acc	      u(:,:,istart-2:iend+2),      v(:,:,istart-2:iend+2),      w(:,:,istart-2:iend+2)) &
!$acc create(u2(:,:,istart-2:iend+2),     v2(:,:,istart-2:iend+2),     w2(:,:,istart-2:iend+2), &
!$acc	u_star1(:,:,istart-2:iend+2),v_star1(:,:,istart-2:iend+2),w_star1(:,:,istart-2:iend+2), &
!$acc	 u_star(:,:,istart-2:iend+2), v_star(:,:,istart-2:iend+2), w_star(:,:,istart-2:iend+2), &
!$acc		 FX(:,:,istart-2:iend+2),     FY(:,:,istart-2:iend+2),     FZ(:,:,istart-2:iend+2))

!$acc enter data copyin(TKE(:,:,istart-2:iend+2), u_TimeAvg(:,:,istart-2:iend+2), &
!$acc	v_TimeAvg(:,:,istart-2:iend+2), w_TimeAvg(:,:,istart-2:iend+2), p_TimeAvg(:,:,istart-2:iend+2)) 	if(filer_running_avg=='STANDARD' .OR. filer_running_avg=='ON')

!$acc enter data copyin(u_TimeAvg(:,:,istart-2:iend+2), v_TimeAvg(:,:,istart-2:iend+2), &
!$acc	w_TimeAvg(:,:,istart-2:iend+2), p_TimeAvg(:,:,istart-2:iend+2), &
!$acc	RS_uu(:,:,istart-2:iend+2), RS_uv(:,:,istart-2:iend+2), RS_uw(:,:,istart-2:iend+2),  &
!$acc	RS_vv(:,:,istart-2:iend+2), RS_vw(:,:,istart-2:iend+2), RS_ww(:,:,istart-2:iend+2)) 				if(filer_running_avg=='FULL')

!$acc enter data copyin(wdist,wsin,wcos,fdamp(:,istart:iend),ypluss(:,istart:iend),upluss(:,istart:iend), &
!$acc	yplus_print(:,:,istart:iend),uplus_print(:,:,istart:iend),dist_print(:,:,istart:iend))	 			if(wall_model=='ON' .OR. Damping_F=='ON')

!$acc enter data copyin(p_no(:,:,istart-2:iend+2),p_pre(:,:,istart-2:iend+2,:)) &
!$acc create(FX1(:,:,istart-2:iend+2),FY1(:,:,istart-2:iend+2),FZ1(:,:,istart-2:iend+2)) 					if(correction_stage>0)

!-------------------------main loop on the timesteps--------------------------!
do while (time <= total_time)                                                                     

	if (time_stepping_by=='DT') then
		dt = constant_dt				!overwrititng dt from calcul_new_velocity.f90 
	endif
	
	inv_dt = 1.d0/dt	
	  
	  istep = istep+1
      time = time + dt
	  
   !time1start = MPI_WTIME()
  

    !*************************** Solid motion (dynamic coupling) ***********************************!      
      if( time .GE. StartDynamic_time .AND. solid_motion /= 0 )then

        !---------------------Select the coupling type-----------------------!		*** dynamic_coupling.f90 ***
		if (solid_motion == 1) then
			call oneway_coupling()  			! active solid motion
		elseif (solid_motion == 2) then
			call twoway_coupling()        		! passive solid motion
		endif

        !--------------------Create VOS at each timestep---------------------!

		if (VOS_by == 2) then
				call vos_ray2d()                    ! 								*** vos_ray2d.f90 ***
		elseif (VOS_by == 3) then
				call vos_ray3d()                    ! 								*** vos_ray3d.f90 ***
		elseif (VOS_by == 1) then
				call func_Cylinder() 				!								*** vos_function.f90 ***
		!		call func_Cylplate_noblade()
		!		call func_Cylplate_1blade()
		!		call func_Cylplate_2blade()  
		!		call func_Cylplate_3blade()
		!		call func_Sphere()
		!		call airfoil_00xx()
		!       call func_darrius_3blade
		!		call square_cyl()
		!		call flat_plate()
		!		call func_torus3d
		!		call func_rbc()
		!		call func_64Spheres()
		endif 
        !--------------------------------------------------------------------!

         !--------------------------------------------------------------------!
         !if(DBD == 1)then
         !   write(71,*)'DBD'
         !   call plasma()
         !   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         !end if
         !define plasma force field
         !--------------------------------------------------------------------!

      end if
    !***********************************************************************************************!  

   !time1end = MPI_WTIME()

      !if (unDBD >= 0.01d0 ) then
        ! write(71,*)'unDBD'
        ! call unsteady_plasma()
      !endif

      !if(myid==master)then 
      ! write(71,*) 'time = ',REAL(time) , '  , AOA = ', REAL(AOA), REAL(angular_vel)
	  !	write(71,*) 'Alpha = ', REAL(blade_alpha), '  , Blade angular vel = ', REAL(blade_omega)
      !endif

      !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
	  !$acc update self(u(:,:,istart-2:istart+1),v(:,:,istart-2:istart+1),w(:,:,istart-2:istart+1), &
	  !$acc				u(:,:,iend-1:iend+2),    v(:,:,iend-1:iend+2),    w(:,:,iend-1:iend+2)	 )  if(nproc>1)
      icount = 2*(nx+4)*(ny+4)
      itag = 100
      call MPI_SENDRECV( u(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         u(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 101
      call MPI_SENDRECV( v(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         v(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 102
      call MPI_SENDRECV( w(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         w(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )


      itag = 103
      call MPI_SENDRECV( u(-1,-1,iend-1),   icount, MPI_REAL8, r_nbr, itag, &
                         u(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 104
      call MPI_SENDRECV( v(-1,-1,iend-1),   icount, MPI_REAL8, r_nbr, itag, &
                         v(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 105
      call MPI_SENDRECV( w(-1,-1,iend-1),   icount, MPI_REAL8, r_nbr, itag, &
                         w(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )


      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	  !$acc update device(u(:,:,iend+1:iend+2),    v(:,:,iend+1:iend+2),    w(:,:,iend+1:iend+2),     &
	  !$acc				  u(:,:,istart-2:istart-1),v(:,:,istart-2:istart-1),w(:,:,istart-2:istart-1)) if(nproc>1)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   !time2start = MPI_WTIME()
      !------------------------------------------------------------------------------------------------------!

	  if (LES_mode=='ON') then
!call nvtxStartRange("Smagorinsky")
		 call CalculateSmagorinskyViscosity()    ! calculate Smagorinsky Viscosity		*** Smagorinsky.f90 ***
!call nvtxEndRange		 
	  endif
	  
	  if ((wall_model=='ON' .OR. Damping_F=='ON') ) then! .AND. time > 50.d0) then ! .AND. istep == INT((10.d0*dt)/dt)) then
!call nvtxStartRange("Wall model")
			! if( time .GE. StartDynamic_time )then
				! call wall_distance()    ! re-calculate wall distance
			! endif
		 call wall_mod_cyl()    ! calculate nut wall function and/or damping function		*** wall_model.f90 ***  (under development)
!call nvtxEndRange		 
	  endif	


      !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
	  !$acc update self(nut(:,:,istart), &
	  !$acc				nut(:,:,iend+1)) if(nproc>1)
      icount = (nx)*(ny)
      itag = 106
      call MPI_SENDRECV( nut(1,1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         nut(1,1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	  !$acc update device(nut(:,:,iend+1)) if(nproc>1)     
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


      !------------------------------------------------------------------------------------------------------!
!call nvtxStartRange("QUICK")
      !------------------------------------------------------------------------------------------------------!
      call discretisation_QUICK_centre()      ! calculate velocity field				*** QUICK.f90 ***
      !------------------------------------------------------------------------------------------------------!
!call nvtxEndRange
!call nvtxStartRange("AB")

      !------------------------------------------------------------------------------------------------------!
      call AdamsBashforth()                   ! Adams-Bashforth							*** AdamsBashforth.f90 ***
      !------------------------------------------------------------------------------------------------------!
!call nvtxEndRange
   !time2end = MPI_WTIME()

!call nvtxStartRange("u0_Boundary conditions")   
      call u0_boundary_conditions()			  !											*** boundary_conditions.f90 ***
!call nvtxEndRange

      !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
	  !$acc update self(w0(:,:,iend),		&
	  !$acc				w0(:,:,istart-1)  ) if(nproc>1)
      icount = (nx+4)*(ny+4)
      itag = 107
      call MPI_SENDRECV( w0(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
                         w0(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	  !$acc update device(w0(:,:,istart-1)) if(nproc>1)     
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


	  if (correction_stage .GT. 0 .AND. istep == INT(StartCorrection_time/dt)+1) then
		!$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
		!$acc parallel present(p_no(:,:,istart-2:iend+2),p_pre(:,:,istart-2:iend+2,:))
		!$acc loop independent collapse(3) gang vector 
		do k=istart-2,iend+2; do j=-1,ny+2; do i=-1,nx+2
			p_no(i,j,k)=p(i,j,k)
			p_pre(i,j,k,1)=p(i,j,k)
		end do; end do; enddo
		!$acc end parallel
		!$OMP END PARALLEL DO
	  endif

      if (correction_stage .GT. 0 .AND. time > StartCorrection_time) then
		!$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
		!$acc parallel present(p_no(:,:,istart-2:iend+2))
		!$acc loop independent collapse(3) gang vector 
		do k=istart-2,iend+2; do j=-1,ny+2; do i=-1,nx+2
			p(i,j,k)=p_no(i,j,k)
		end do; end do; enddo
		!$acc end parallel
		!$OMP END PARALLEL DO
	  endif
   !time3start = MPI_WTIME()	
!call nvtxStartRange("Pressure solver")    
      !------------------------------------------------------------------------------------------------------!	*** SOR.f90 ***
      ccc=0

      if (pressure_solver == 1) then
       call RB_SOR()			! calculate pressure field using red-black SOR --> nx MUST be an even number
      elseif (pressure_solver == 2) then
	   call SOR()			 	! calculate pressure field using traditional SOR	  
      endif
!      call gauss_seidel()        ! calculate pressure field using traditional SOR (Old)
!      call BICG_stab()           ! calculate pressure field
      !------------------------------------------------------------------------------------------------------!
   !time3end = MPI_WTIME()

!call nvtxEndRange
      if (correction_stage .GT. 0 .AND. time > StartCorrection_time) then
		!$OMP PARALLEL DO PRIVATE(i,j) collapse(nclps)
		!$acc parallel present(p_no(:,:,istart-2:iend+2))
		!$acc loop independent collapse(3) gang vector 
		do k=istart-2,iend+2; do j=-1,ny+2; do i=-1,nx+2
			p_no(i,j,k)=p(i,j,k)
		end do; end do; enddo
		!$acc end parallel
		!$OMP END PARALLEL DO
	  endif


      !--------------------Exit on diverged--------------------!			! Emergency stop if pressure solver diverges
		if (myid==master .AND. (pChangeMax_>huge(0.0) .OR. pChangeMax_/=pChangeMax_)) then
			call date_and_time(date=date_str, time=time_str)						
			write(73,*) ' '	
			write(73,'(1X,A4,"-",A2,"-",A2," ",A2,":",A2,":",A2,A,I7,A,F9.4,A)') &
				date_str(1:4), date_str(5:6), date_str(7:8), &
				time_str(1:2), time_str(3:4), time_str(5:6), &			
				' Pressure solver diverged at iteration ', istep, &
				' , time = ',time, ' sec'
			write(73,*) ' '	
			write(73,'(A,F9.4,A)') ' Last backup time = ',lastbackuptime,' sec'
			write(73, '(/)')
			stop
		end if		 
      !--------------------------------------------------------!
	  
	  

      !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
	  !$acc update self(p(:,:,istart),	&
	  !$acc				p(:,:,iend+1)) if(nproc>1)
      icount = (nx+4)*(ny+4)
      itag = 108
      call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	  !$acc update device(p(:,:,iend+1)) if(nproc>1)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


   !time4start = MPI_WTIME()
!call nvtxStartRange("calcul_new_vel")     
      !------------------------------------------------------------------------------------------------------!
      call calcul_new_velocity() ! update velocity field 							*** calcul_new_velocity.f90 ***
      !------------------------------------------------------------------------------------------------------!
!call nvtxEndRange

      if (correction_stage .GT. 0 .AND. time > StartCorrection_time) then
		call prediction_correction()	!											*** prediction_correction.f90 ***
      endif
   !time4end = MPI_WTIME() 

   !time5start = MPI_WTIME()
!call nvtxStartRange("virtual force integrator")
      !------------------------------------------------------------------------------------------------------!
	  ! Virtual force integrator													*** virtualForceIntegrator.f90 ***
      		call virtualForceIntegrator_nima() 	! CD and CL only
	  !	    call virtualForceTorqueIntegrator() ! CD, CL, CT, and CP
	  !	    call virtualForce_Magnus2blade() ! CD, CL, CT, and CP
	  !	    call virtualForce_Magnus3blade() ! CD, CL, CT, and CP
	  !		call virtualForceTorqueIntegrator3D() ! CD, CL, CT, and CP
	  !		call virtualForceTorque_frozen()	! CD, CL, CT, and CP for frozen rotor simulation
      !------------------------------------------------------------------------------------------------------!
!call nvtxEndRange


!call nvtxStartRange("updating velocity")
      !------------------------------------------------------------------------------------------------------!
      call Updating_velocity()   ! update velocity field 							*** calcul_new_velocity.f90 ***
      !------------------------------------------------------------------------------------------------------!
!call nvtxEndRange  


	  ! if ((wall_model=='ON' .OR. Damping_F=='ON') .AND. time > 50.d0) then ! .AND. istep == INT((10.d0*dt)/dt)) then
! !call nvtxStartRange("Wall model")
			! ! if( time .GE. StartDynamic_time )then
				! ! call wall_distance()    ! re-calculate wall distance
			! ! endif
		 ! call wall_mod_cyl()    ! calculate nut wall function and/or damping function		*** wall_model.f90 ***  (under development)
! !call nvtxEndRange		 
	  ! endif	

   
!call nvtxStartRange("final BC")
      !------------------------------------------------------------------------------------------------------!
      call final_boundary_conditions() ! recall boundary conditions to update them			*** boundary_conditions.f90 ***
      !------------------------------------------------------------------------------------------------------!
!call nvtxEndRange     
   !time5end = MPI_WTIME()     


        if (myid==master .AND. time > coeff_start) then
			call filerProcess()									! write CD_time				*** filer.f90 ***
        endif


      !---------- Calculate running average of u,v,w,p,TKE ----------!

		if (filer_running_avg/='OFF' .AND. time .GE. startfilerAvg_time) then

			!>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
			!$acc update self(w(:,:,iend),		&
			!$acc			  w(:,:,istart-1) ) if(nproc>1)
			icount = (nx+4)*(ny+4)
			itag = 113
			call MPI_SENDRECV( w(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   w(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
			call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			!$acc update device(w(:,:,istart-1)) if(nproc>1)     
			!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			if (filer_running_avg=='STANDARD' .OR. filer_running_avg=='ON') then
				call running_avg()
			elseif (filer_running_avg=='FULL') then
				call running_avg_full()
			endif
		endif


! ---------------------------------- Probe points -----------------------------------
        ! if (time > coeff_start) then
				! !>>>>>>>>>>>>>>>>> data transfer among MPI processes <<<<<<<<<<<<<<<<<<<
				! !$acc update self(u(:,:,istart-1:istart),v(:,:,istart-1:istart),w(:,:,istart-1:istart), &
				! !$acc			  u(:,:,iend:iend+1),    v(:,:,iend:iend+1),    w(:,:,iend:iend+1)   ) if(nproc>1)
				
				! icount = (nx+4)*(ny+4)
				! itag = 128
				! call MPI_SENDRECV( u(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
								   ! u(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
				! itag = 129
				! call MPI_SENDRECV( v(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
								   ! v(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
				! itag = 130
				! call MPI_SENDRECV( w(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
								   ! w(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
				! ! itag = xxx	******Already done (itag = 108)******
				! ! call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
								   ! ! p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )

				! itag = 131
				! call MPI_SENDRECV( u(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
								   ! u(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
				! itag = 132
				! call MPI_SENDRECV( v(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
								   ! v(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )      
				! itag = 133
				! call MPI_SENDRECV( w(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
								   ! w(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
				! call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
				! !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<            
			
			! do rb=1,n_points
				! if (near_pt(3*rb) .GE. istart .AND. near_pt(3*rb) .LE. iend) then
				! do k=1,2 ; do j=1,2 ; do i=1,2
					! !$acc update self (u(near_pt(3*rb-2)  +(i-1),near_pt(3*rb-1)+(j-1),near_pt(3*rb)+(k-1)), &
					! !$acc			   u(near_pt(3*rb-2)-1+(i-1),near_pt(3*rb-1)+(j-1),near_pt(3*rb)+(k-1)), &
					! !$acc			   v(near_pt(3*rb-2)+(i-1),near_pt(3*rb-1)  +(j-1),near_pt(3*rb)+(k-1)), &
					! !$acc			   v(near_pt(3*rb-2)+(i-1),near_pt(3*rb-1)-1+(j-1),near_pt(3*rb)+(k-1)), &
					! !$acc			   w(near_pt(3*rb-2)+(i-1),near_pt(3*rb-1)+(j-1),near_pt(3*rb)  +(k-1)), &
					! !$acc			   w(near_pt(3*rb-2)+(i-1),near_pt(3*rb-1)+(j-1),near_pt(3*rb)-1+(k-1)), &
					! !$acc			   p(near_pt(3*rb-2)+(i-1),near_pt(3*rb-1)+(j-1),near_pt(3*rb)+(k-1))    )
				! enddo; enddo; enddo
				! endif
				! call MPI_Barrier(MPI_COMM_WORLD, ierr)
			! enddo
			
			! call acc_set_device_type(acc_device_host)			
			! call filerProber()									!							*** filer.f90 ***
			! call acc_set_device_num(myid,acc_device_nvidia)
		! endif
! -----------------------------------------------------------------------------------



! ------------------------------- Create output files -------------------------------

!call nvtxStartRange("output files")
	if ((filer3d=='ON' .AND. time >= next_filer3d-1.0d-10) .OR. (filer2d/='OFF' .AND. time >= next_filer2d-1.0d-10) .OR. &
		(filer_cp=='ON' .AND. time >= next_filerCp-1.0d-10) .OR. time >= next_backup-1.0d-10) then        			!   modified 17_10_2024
	!$acc update self(u(:,:,istart:iend),v(:,:,istart:iend),w(:,:,istart:iend),p(:,:,istart:iend),ETA)
          
      !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
      if(myid>master)then
		 icount = igcount*(nx+4)*(ny+4)
         itag = 109
         call MPI_SEND( p(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 110
         call MPI_SEND( u(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 111
         call MPI_SEND( v(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 112
         call MPI_SEND( w(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if

      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount(i)*(nx+4)*(ny+4)
            itag = 109
            call MPI_RECV( p(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 110
            call MPI_RECV( u(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 111
            call MPI_RECV( v(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 112
            call MPI_RECV( w(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )     
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

           
      if(myid==master)then
	  call acc_set_device_type(acc_device_host)
   !time6start = MPI_WTIME()

        if (time > coeff_start .AND. time >= next_filerCp-1.0d-10) then
			if(filer_cp=='ON') then            
				call filerProcess_cp()				! write surface pressure distribution			*** filer.f90 ***
				next_filerCp = next_filerCp + filerCp_int
			endif
        endif


!call nvtxStartRange("write_result")         
        if (time > startfiler3d_time-dt .AND. time >= next_filer3d-1.0d-10) then
            !call filer_Reynoldstress()
			if(filer3d=='ON') then
				call filereachtime3d() 				! write filer 3D								*** filer.f90 ***
				next_filer3d = next_filer3d + filer3d_int
			endif
        end if 

        if (time > startfiler2d_time-dt .AND. time >= next_filer2d-1.0d-10) then
			if(filer2d=='AVG' .OR. filer2d=='ON') then
				call filereachtime2d_Avg()				! write filer 2D							*** filer.f90 ***
				next_filer2d = next_filer2d + filer2d_int
			elseif (filer2d=='MID') then
				call filereachtime2d_Mid()				! write filer 2D							*** filer.f90 ***
				next_filer2d = next_filer2d + filer2d_int
			endif
        end if 
!call nvtxEndRange 

!call nvtxStartRange("write_backup")		 
        if (time >= next_backup-1.0d-10) then
            call backupfile()						! write last data backup in case of power cut	*** filer.f90 ***
			lastbackuptime = time
			next_backup = next_backup + backup_int
        end if 
!call nvtxEndRange 
	  
	  call acc_set_device_num(myid,acc_device_nvidia)
      end if
	
	call MPI_BCAST(next_filer3d, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(next_filer2d, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(next_backup , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(next_filerCp, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	
   !time6end = MPI_WTIME()      
	end if
! -----------------------------------------------------------------------------------


! ------------------------------- Wall model output ----------------------------------

	if (wall_model=='ON' .AND. time > startfiler2d_time-dt .AND. time >= next_filerWall-1.0d-10) then        			!   modified 19_12_2024
	!$acc update self(yplus_print(:,:,istart:iend),uplus_print(:,:,istart:iend),dist_print(:,:,istart:iend),nut(:,:,istart:iend),ETA)
          
      !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
      if(myid>master)then
		 icount = igcount*(nx)*(ny)
         itag = 119
         call MPI_SEND(yplus_print(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 120
         call MPI_SEND(uplus_print(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 121
         call MPI_SEND(dist_print(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 122
         call MPI_SEND(nut(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if

      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount(i)*(nx)*(ny)
            itag = 119
            call MPI_RECV(yplus_print(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 120
            call MPI_RECV(uplus_print(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 121
            call MPI_RECV(dist_print(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 122
			call MPI_RECV(nut(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )			
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

           
      if(myid==master)then
	  call acc_set_device_type(acc_device_host)

            call filereachtime2d_wall()
			next_filerWall = next_filerWall + filer2d_int
	  call acc_set_device_num(myid,acc_device_nvidia)
      end if

	call MPI_BCAST(next_filerWall, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)    
	end if
! --------------------------------------------------------------------------------


! ------------------------------- Virtual force output ----------------------------------

	if (filer_virtual_force=='ON' .AND. time > startfilerVF_time-dt .AND. time >= next_filerVF-1.0d-10) then        			!   modified 06_02_2025
	!$acc update self(FX(:,:,istart:iend),FY(:,:,istart:iend),FZ(:,:,istart:iend),p(:,:,istart:iend),ETA)
          
      !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
      if(myid>master)then
		 icount = igcount*(nx+4)*(ny+4)
         itag = 114
         call MPI_SEND( p(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 115
         call MPI_SEND( FX(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 116
         call MPI_SEND( FY(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 117
         call MPI_SEND( FZ(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if

      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount(i)*(nx+4)*(ny+4)
            itag = 114
            call MPI_RECV( p(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 115
            call MPI_RECV( FX(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 116
            call MPI_RECV( FY(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 117
            call MPI_RECV( FZ(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )     
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

           
      if(myid==master)then
	  call acc_set_device_type(acc_device_host)

            call filereachtime3d_VF()
			next_filerVF = next_filerVF + filerVF_int
	  
	  call acc_set_device_num(myid,acc_device_nvidia)
      end if
	
	call MPI_BCAST(next_filerVF, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	end if
! --------------------------------------------------------------------------------



! -------------------------- STANDARD running Avg output files --------------------------

	if ((filer_running_avg=='STANDARD' .OR. filer_running_avg=='ON') .AND. time > startfilerAvg_time-dt .AND. &
		(time >= next_filerAvg-1.0d-10 .OR. time >= next_backupAvg-1.0d-10)) then        			!   modified 17_10_2024
	!$acc update self(u_TimeAvg(:,:,istart:iend),v_TimeAvg(:,:,istart:iend),w_TimeAvg(:,:,istart:iend), p_TimeAvg(:,:,istart:iend), TKE(:,:,istart:iend))
          
      !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
      if(myid>master)then
		 icount = igcount*(nx+4)*(ny+4)
         itag = 118
         call MPI_SEND(u_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 119
         call MPI_SEND(v_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 120
         call MPI_SEND(w_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 121
         call MPI_SEND(p_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 122
         call MPI_SEND(		 TKE(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if

      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount(i)*(nx+4)*(ny+4)
            itag = 118
            call MPI_RECV(u_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 119
            call MPI_RECV(v_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 120
            call MPI_RECV(w_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 121
			call MPI_RECV(p_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
			itag = 122
            call MPI_RECV(		TKE(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )			
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

           
      if(myid==master)then
	  call acc_set_device_type(acc_device_host)

        if (time >= next_filerAvg-1.0d-10) then 
			call filer3d_avg()						! calculate Kinetic Energy						*** filer.f90 ***
			next_filerAvg = next_filerAvg + filerAvg_int
        end if 
		
        if (time >= next_backupAvg-1.0d-10) then            
            call backupAvg()						! write last data backup in case of power cut	*** filer.f90 ***
			next_backupAvg = next_backupAvg + backup_int
		end if 
	  
	  call acc_set_device_num(myid,acc_device_nvidia)
      end if

	call MPI_BCAST(next_filerAvg, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr) 
	call MPI_BCAST(next_backupAvg, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)	
	end if
! --------------------------------------------------------------------------------


! -------------------------- FULL running Avg output files --------------------------

	if (filer_running_avg=='FULL' .AND. time > startfilerAvg_time-dt .AND. &
	   (time >= next_filerAvg-1.0d-10 .OR. time >= next_backupAvg-1.0d-10)) then        		
	!$acc update self(u_TimeAvg(:,:,istart:iend),v_TimeAvg(:,:,istart:iend),w_TimeAvg(:,:,istart:iend), &
	!$acc 	p_TimeAvg(:,:,istart:iend),RS_uu(:,:,istart:iend),RS_uv(:,:,istart:iend),RS_uw(:,:,istart:iend), &
	!$acc	RS_vv(:,:,istart:iend),RS_vw(:,:,istart:iend),RS_ww(:,:,istart:iend))
          
      !>>>>>>>>>>>>>>>> send results back to Master process <<<<<<<<<<<<<<<<<<
      if(myid>master)then
		 icount = igcount*(nx+4)*(ny+4)
         itag = 118
         call MPI_SEND(u_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 119
         call MPI_SEND(v_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 120
         call MPI_SEND(w_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 121
         call MPI_SEND(p_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 122
         call MPI_SEND(	   RS_uu(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 123
         call MPI_SEND(	   RS_uv(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 124
         call MPI_SEND(	   RS_uw(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 125
         call MPI_SEND(	   RS_vv(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 126
         call MPI_SEND(	   RS_vw(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 127
         call MPI_SEND(	   RS_ww(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if

      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount(i)*(nx+4)*(ny+4)
            itag = 118
            call MPI_RECV(u_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 119
            call MPI_RECV(v_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 120
            call MPI_RECV(w_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 121
			call MPI_RECV(p_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
			itag = 122
            call MPI_RECV(    RS_uu(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )			
			itag = 123
            call MPI_RECV(    RS_uv(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )			
			itag = 124
            call MPI_RECV(    RS_uw(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )			
			itag = 125
            call MPI_RECV(    RS_vv(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )			
			itag = 126
            call MPI_RECV(    RS_vw(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )			
			itag = 127
            call MPI_RECV(    RS_ww(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )			
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

           
      if(myid==master)then
	  call acc_set_device_type(acc_device_host)

        if (time >= next_filerAvg-1.0d-10) then 
			call filer3d_avg_full1()						
			call filer3d_avg_full2()
			next_filerAvg = next_filerAvg + filerAvg_int			
        end if 
		
        if (time >= next_backupAvg-1.0d-10) then            
            call backup_avg_full1()            
            call backup_avg_full2()
			next_backupAvg = next_backupAvg + backup_int			
        end if 
	  
	  call acc_set_device_num(myid,acc_device_nvidia)
      end if

	call MPI_BCAST(next_filerAvg, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr) 
	call MPI_BCAST(next_backupAvg, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)	    
	end if
! --------------------------------------------------------------------------------


      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   
      !----------calculate wall time----------!
      totalfinaltime = MPI_WTIME()
      totalcosttime = totalfinaltime-totalstarttime
      if(myid==master)then
         write(71,*) 'time = ', REAL(time), 'Cd = ', REAL(cDrag), ', Cl = ', REAL(cLift)
         write(71,*) 'step cost = ', REAL(totalfinaltime-totallasttime), 'total cost = ', REAL(totalcosttime/3600.) ,'(Hr)'
		 !write(71,*) 'dyn time  =', REAL(time1end-time1start), 'quickAB =', REAL(time2end-time2start), 'GS = ', REAL(time3end-time3start)
		 !write(71,*) 'vsol-PC =', REAL(time4end-time4start), 'cdcl-Vupdate = ', REAL(time5end-time5start)!, 'filer = ', REAL(time6end-time6start)
         write(71,*)
      end if
      totallasttime = totalfinaltime
      !----------calculate wall time----------!


      !---------------print simulation progress---------------!
		if (myid == master .AND. time-initial_time >= 2.d0 .AND. time-initial_time-dt < 2.d0) then
			timeapproxstart = MPI_WTIME()
		endif
		
		if (myid == master .AND. time-initial_time >= 10.d0 .AND. time-initial_time-dt < 10.d0) then
			timeapproxend = MPI_WTIME()
			write(73,'(1X,19("*"),A,F8.4,A)') ' Total cost time (Approx.)= ', &
				REAL((timeapproxend-timeapproxstart)*(total_time-initial_time)/8./3600.),' Hr'
			write(73,*) ' '
			write(73,*) '      Date     Time    Done        t*     Elapsed   Remaining'
			flush(73)		
		endif
		
		if (myid == master .AND. time >= next_progress-1.0d-10) then
			progress = INT(100*(time-initial_time)/(total_time-initial_time))
			call date_and_time(date=date_str, time=time_str)
			write(73,'(1X,A4,"-",A2,"-",A2," ",A2,":",A2,":",A2,A,I5,A,I6,A,2(F9.4,A))') &
				date_str(1:4), date_str(5:6), date_str(7:8), &
				time_str(1:2), time_str(3:4), time_str(5:6), ' ', &
				progress,' %',NINT(time),' sec',REAL(totalcosttime/3600.) ,' Hr', &
				REAL((totalfinaltime-lastprogresstime)*(100-progress)/10./3600.),' Hr'
			flush(73)
			lastprogresstime = totalfinaltime 
			next_progress = next_progress + prog_int
		end if
      !--------------------------------------------------------!

!call nvtxEndRange

   end do

!$acc exit data if(correction_stage>0)
!$acc exit data if(wall_model=='ON' .OR. Damping_F=='ON')
!$acc exit data if(filer_running_avg=='FULL')
!$acc exit data if(filer_running_avg=='STANDARD' .OR. filer_running_avg=='ON')
!$acc end data

!--------------------End of main loop on the timesteps----------------------!


   if(myid==master)then
      close(73)
	  call filer_final()
   end if

   call MPI_FINALIZE(ierr) 

end program main
