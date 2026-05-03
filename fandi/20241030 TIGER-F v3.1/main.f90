! 27 Oct 2024 - FDS
!    y=1 ______________                                                                                 
!       /             /|       |  Author  : Zi-Hsuan Wei                                                 
!      /             / |       |  Version : 3.1                                                          
!     /____________ /  |       |  Mod by  : Fandi D. Suprianto                              
!     |  |         |   |       |  Web     : http://smetana.me.ntust.edu.tw/                                    
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



   nthreads = 14 !(OpenMP or OpenACC)
   
   !---------------------------- OPENMP ----------------------------!
   call omp_set_num_threads(nthreads)
   ! Use only physical cores for best performance
   ! Disable hyperthreading (if any) using 'export OMP_PROC_BIND=TRUE')
   !----------------------------------------------------------------!


   !----------------------------- MPI ------------------------------!	*** MPI_division.f90 ***
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
   call Mpi_division()	! Equal load by default
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

   LES                              = 1                     ! the LES mode. 1 : on ; 0 : off

!   wall_model						= 0						! Apply wall model 1 : on ; 0 : off				(under development)
   
!   Damping_F						= 0						! Van Driest damping function 1 : on ; 0 : off 	(under development)
   
   pressure_solver					= 1						! 1)RB_SOR, 2)SOR-CPU
   
   correction_stage					= 0						! number of correction stage in Prediction-correction process
		StartCorrection_time		= 0.d0					! initialize time of prediction-correction

   steadiness                       = 2                     ! steady : 1 ; unsteady : 2

   resume							= 0						! resume simulation. 1:on ; 0:off (Provide the backup file: 'STATIC.bak'/'DYNAMIC.bak')
   

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
   
   omega                            = 1.7d0                 ! Set value for SOR method

   zeta                             = 1.0d-4                ! zeta for solving pressure matrix

   itmax                            = 5000                  ! maximum for zeta in Gauss Seidel subroutines

   initial_time                     = 0.d0                  ! initialize time of simulation

   total_time                       = 200.d0	 			! total time of simulation
    
!   zeta_vel                         = 1.0d-5                ! steady criteria (if steadiness = 1)

!   istea_int						 = 0.1d0				 ! time interval to check steady state

   !----------------------------------------------------------------!


   !---------------------Output Data Management---------------------!

   coeff_start						= 0.1d0					! initialize time of storing force & pressure coefficients (CD_time & Cp)

   filer3d							= 0						! create 3D output (3d_xxxx.q). 1 : on ; 0 : off
			startfiler3d_time		= 0.d0					! initialize time of storing filer3d data
			isto3d_int				= 1.d0					! time interval to storing filer3d data
   
   filer2d							= 1						! create 2D output (2d_xxxx.q). 1 : on ; 0 : off (Spanwise avg - Used only in the case of extruded solids)
			startfiler2d_time		= 0.d0					! initialize time of storing filer2d data
			isto2d_int				= 1.d0					! time interval to storing filer2d data
   
   filer_cp							= 0						! create surface pressure distribution plot. 1 : on ; 0 : off (set in filer.f90)
			istocp_int				= 5.d0					! time interval to storing cp data (start = coeff_start)
   
   backup_int						= 1.d0					! time interval to store backup data


   calculate_avg					= 0						! Calculate running average of u,v,w,p,TKE for all timesteps
			startfilerAvg_time		= 0.d0					! initialize time of storing average data
			istoAvg_int				= 5.d0   				! time interval to storing average data (the last timestep is enough in most cases)
   
   !----------------------------------------------------------------!  
   

   totalstarttime = MPI_WTIME()
   totallasttime = totalstarttime
   
  
   !----------------------------------------------------------------!	*** read_data.f90 ***
   call read_data()              ! define dx, dy, dz
   !----------------------------------------------------------------!

  
   !--------------------- Gridder selection ------------------------!	*** gridder.f90 ***
   if(Gridder=='non-uniform-sin4')then
      call gridder_sin4()              ! define unequal mesh proposed by Kuyper
   else if(Gridder=='non-uniform-sin4-3D')then
      call gridder_sin4_3D()           ! define unequal mesh proposed by Kuyper in x,y,z
   else if(Gridder=='uniform')then
      call gridder_equal()             ! define equal mesh
   else if(Gridder=='ground-3D')then
      call gridder_ground()  		   ! define unequal mesh with two sine domain in x,y,z for object lies on the bottom wall
   else if(Gridder=='non-uniform-sin5-3D')then
      call gridder_sin5_3D()           ! define unequal mesh with two sine domain in x,y,z
   else if(Gridder=='non-uniform-sin5')then
      call gridder_sin5()              ! define unequal mesh with two sine domain in x,y
   else if(Gridder=='non-uniform-fall')then
      call gridder_fall()
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
!   else if(Gridder=='read_txt' )then
!      call gridder_read_txt()          ! grid coordinate is supplied (grid.txt) in X(i),Y(j),Z(k) format
   end if
   !----------------------------------------------------------------!

   !----------------------------------------------------------------!	*** initial_conditions.f90 ***
   call initial_conditions()        !call initial conditions
   !----------------------------------------------------------------!
   
	istep = 0
	istep2= 0

   !-------------------- Loading backup data -----------------------!	*** read_data.f90 ***
   if (resume == 1) then													
	  call read_backupfile()        				!input first .bak file			--> (used when resuming simulation)				
	  nstep = nstep - NINT(lastsaved_time/dt)
	  initial_time = lastsaved_time
	  time = lastsaved_time
   endif
   !----------------------------------------------------------------!

   !------------------ Loading avg backup data ---------------------!	*** filer.f90 ***
   if (calculate_avg==1 .AND. resume == 1) then													
	  call read_backup_uvwp_Avg()					!input first .bak file			--> (used when resuming simulation)				
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
		do j=1,ny; do i=1,nx
		ETAs=0.0
		do k=1,nz 
			if ( ETA(i,j,k)>ETAs ) then
				ETAs=ETA(i,j,k)
			endif 
		enddo
		Z_castingarea=Z_castingarea+ETAs*iDx(i)*iDy(j)
		enddo; enddo

		Y_castingarea=0.0
		do k=1,nz; do i=1,nx
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
	! if (wall_model ==1 .OR. Damping_F ==1) then
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


   if(myid==master) then
	  call acc_set_device_type(acc_device_host)
		!call prober()

      if(filer3d==1) then
		call filereachtime3d() 										!		*** filer.f90 ***
	  endif
	  if(filer2d==1) then
		call filereachtime2d() 										!		*** filer.f90 ***
	  endif 
      !call filer_bodyforce()
      call filerInfo()												!		*** filer.f90 ***
          !---------------------
          !call filerProcess_cp() 
          !call filer_Reynoldstress()
	  call acc_set_device_num(myid,acc_device_nvidia)
   endif


   !call write_matrix_a()

	open (71,file='output.dat',position='append')

	!open (72,file='vos_time.dat',position='append')

	if(myid==nproc-1) then
		write(71,'(4(A,I5,3x))') ' Last Proc  =',myid+1,'istart =',istart,'iend =',iend,'gcount =',igcount	
	endif	 
   
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)   

    if(myid==master) then	
		write(71,'(A,I5,3x)') ' N threads  =',nthreads
		write(71,800)
800 	FORMAT(' ',63('-'))
	endif


	call MPI_BARRIER(MPI_COMM_WORLD, ierr) 


!$acc data copyin(X,Y,Z,Xs,Ys,Zs,iDx,iDy,iDz,Dxs,Dys,Dzs,ETA,p(:,:,istart-2:iend+2), &
!$acc	  P_Den(:,:,istart:iend,:),  Den_inv(:,:,istart:iend),     Viseff(:,:,istart:iend)    , &
!$acc	     u0(:,:,istart-2:iend+2),     v0(:,:,istart-2:iend+2),     w0(:,:,istart-2:iend+2), &
!$acc	      u(:,:,istart-2:iend+2),      v(:,:,istart-2:iend+2),      w(:,:,istart-2:iend+2)) &
!$acc create(u2(:,:,istart-2:iend+2),     v2(:,:,istart-2:iend+2),     w2(:,:,istart-2:iend+2), &
!$acc	u_star1(:,:,istart-2:iend+2),v_star1(:,:,istart-2:iend+2),w_star1(:,:,istart-2:iend+2), &
!$acc	 u_star(:,:,istart-2:iend+2), v_star(:,:,istart-2:iend+2), w_star(:,:,istart-2:iend+2), &
!$acc		 FX(:,:,istart-2:iend+2),     FY(:,:,istart-2:iend+2),     FZ(:,:,istart-2:iend+2))

!$acc enter data copyin(TKE(:,:,istart-2:iend+2), u_TimeAvg(:,:,istart-2:iend+2), &
!$acc	v_TimeAvg(:,:,istart-2:iend+2), w_TimeAvg(:,:,istart-2:iend+2), p_TimeAvg(:,:,istart-2:iend+2)) 	if(calculate_avg==1)

!!$acc enter data copyin(wdist,wsin,wcos,fdamp(:,istart:iend),ypluss(:,istart:iend),upluss(:,istart:iend)) 	if(wall_model==1 .OR. Damping_F==1)

!$acc enter data copyin(p_no(:,:,istart-2:iend+2),p_pre(:,:,istart-2:iend+2,:)) &
!$acc create(FX1(:,:,istart-2:iend+2),FY1(:,:,istart-2:iend+2),FZ1(:,:,istart-2:iend+2)) 					if(correction_stage>0)

!-------------------------main loop on the timesteps--------------------------!
   do istep=1,nstep                                                                      

      time = initial_time + dt*istep
      
	  
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

      !----------data transformation among nodes----------!
	  !$acc update self(u(:,:,istart-2:istart+1),v(:,:,istart-2:istart+1),w(:,:,istart-2:istart+1)) if(nproc>1)
	  !$acc update self(u(:,:,iend-1:iend+2),    v(:,:,iend-1:iend+2),    w(:,:,iend-1:iend+2)) if(nproc>1)
      icount = 2*(nx+4)*(ny+4)
      itag = 110
      call MPI_SENDRECV( u(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         u(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 120
      call MPI_SENDRECV( v(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         v(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 130
      call MPI_SENDRECV( w(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         w(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )


      itag = 140
      call MPI_SENDRECV( u(-1,-1,iend-1),   icount, MPI_REAL8, r_nbr, itag, &
                         u(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 150
      call MPI_SENDRECV( v(-1,-1,iend-1),   icount, MPI_REAL8, r_nbr, itag, &
                         v(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 160
      call MPI_SENDRECV( w(-1,-1,iend-1),   icount, MPI_REAL8, r_nbr, itag, &
                         w(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )


      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	  !$acc update device(u(:,:,iend+1:iend+2),    v(:,:,iend+1:iend+2),    w(:,:,iend+1:iend+2)) if(nproc>1)
	  !$acc update device(u(:,:,istart-2:istart-1),v(:,:,istart-2:istart-1),w(:,:,istart-2:istart-1)) if(nproc>1)
      !----------data transformation among nodes----------!

   !time2start = MPI_WTIME()
      !------------------------------------------------------------------------------------------------------!

	  if (LES == 1) then
!call nvtxStartRange("Smagorinsky")
		 call CalculateSmagorinskyViscosity()    ! calculate Smagorinsky Viscosity		*** Smagorinsky.f90 ***
!call nvtxEndRange		 
	  endif
	  

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

      !----------data transformation among nodes----------!
	  !$acc update self(w0(:,:,iend)) if(nproc>1)
	  !$acc update self(w0(:,:,istart-1)) if(nproc>1)
      icount = (nx+4)*(ny+4)
      itag = 240
      call MPI_SENDRECV( w0(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
                         w0(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	  !$acc update device(w0(:,:,istart-1)) if(nproc>1)     
	  !----------data transformation among nodes----------!

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
	  !		call virtualForceTorqueIntegrator3D() ! CD, CL, CT, and CP
	  !		call virtualForceTorque_frozen()	! CD, CL, CT, and CP for frozen rotor simulation
      !------------------------------------------------------------------------------------------------------!
!call nvtxEndRange


!call nvtxStartRange("updating velocity")
      !------------------------------------------------------------------------------------------------------!
      call Updating_velocity()   ! update velocity field 							*** calcul_new_velocity.f90 ***
      !------------------------------------------------------------------------------------------------------!
!call nvtxEndRange  


	  ! if ((wall_model ==1 .OR. Damping_F ==1) .AND. time > 50.d0) then ! .AND. istep == INT((10.d0*dt)/dt)) then
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

		if (calculate_avg==1 .AND. time .GE. startfilerAvg_time) then

			!----------data transformation among nodes----------!
			!$acc update self(w(:,:,iend)) if(nproc>1)
			!$acc update self(w(:,:,istart-1)) if(nproc>1)
			icount = (nx+4)*(ny+4)
			itag = 250
			call MPI_SENDRECV( w(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
							   w(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
			call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			!$acc update device(w(:,:,istart-1)) if(nproc>1)     
			!----------data transformation among nodes----------!

			call running_avg()
			istep2 = istep2 + 1	
		endif


! ------------------------------- Create output files -------------------------------

      !----------data collect among nodes for filer----------!

!call nvtxStartRange("output files")
	if ((filer3d==1 .AND. mod(istep,isto3d)==0) .OR. (filer2d==1 .AND. mod(istep,isto2d)==0) .OR. &
		(filer_cp==1 .AND. mod(istep,istocp)==0) .OR. mod(istep,ibackup)==0) then        			!   modified 17_10_2024
	!$acc update self(u(:,:,istart-2:iend+2),v(:,:,istart-2:iend+2),w(:,:,istart-2:iend+2),p(:,:,istart-2:iend+2),ETA)
          
      !Send my results back to the master
      if(myid>master)then
		 icount = igcount*(nx+4)*(ny+4)
         itag = 10
         call MPI_SEND( p(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 20
         call MPI_SEND( u(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 30
         call MPI_SEND( v(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 40
         call MPI_SEND( w(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if

      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount(i)*(nx+4)*(ny+4)
            itag = 10
            call MPI_RECV( p(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 20
            call MPI_RECV( u(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 30
            call MPI_RECV( v(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 40
            call MPI_RECV( w(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )     
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !Wait to receive results from each task
	  
      !----------data collect among nodes for filer----------!
      !end if

           
      if(myid==master)then
	  call acc_set_device_type(acc_device_host)
   !time6start = MPI_WTIME()

        if (time > coeff_start .AND. mod(istep,istocp)==0) then
			if(filer_cp==1) then            
				call filerProcess_cp()				! write surface pressure distribution			*** filer.f90 ***
			endif
        endif


!call nvtxStartRange("write_result")         
        if (time > startfiler3d_time-dt .AND. mod(istep,isto3d)==0) then
            !call filer_Reynoldstress()
			if(filer3d==1) then
				call filereachtime3d() 				! write filer 3D if isto3d = istep				*** filer.f90 ***
			endif
        end if 

        if (time > startfiler2d_time-dt .AND. mod(istep,isto2d)==0) then
			if(filer2d==1) then
				call filereachtime2d()				! write filer 2D if isto2d = istep				*** filer.f90 ***
			endif
        end if 
!call nvtxEndRange 

!call nvtxStartRange("write_backup")		 
        if (mod(istep,ibackup)==0) then
            call backupfile()						! write last data backup in case of power cut	*** filer.f90 ***
        end if 
!call nvtxEndRange 
        ! if (mod(istep,istea)==0) then ! write results if isto = istep
          !  call check_steady()
        ! end if
	  
	  call acc_set_device_num(myid,acc_device_nvidia)
      end if

   !time6end = MPI_WTIME()      
	end if
! -----------------------------------------------------------------------------------


! ------------------------------- TKE output files ----------------------------------

	if (calculate_avg==1 .AND. time > startfilerAvg_time-dt .AND. mod(istep,isto3davg)==0) then        			!   modified 17_10_2024
	!$acc update self(u_TimeAvg(:,:,istart-2:iend+2),v_TimeAvg(:,:,istart-2:iend+2),w_TimeAvg(:,:,istart-2:iend+2), p_TimeAvg(:,:,istart-2:iend+2), TKE(:,:,istart-2:iend+2))
          
      !Send my results back to the master
      if(myid>master)then
		 icount = igcount*(nx+4)*(ny+4)
         itag = 50
         call MPI_SEND(u_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 60
         call MPI_SEND(v_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 70
         call MPI_SEND(w_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 80
         call MPI_SEND(p_TimeAvg(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 90
         call MPI_SEND(		 TKE(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if

      if(myid==master)then
         do i = 1, (nproc-1)
            icount = gcount(i)*(nx+4)*(ny+4)
            itag = 50
            call MPI_RECV(u_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 60
            call MPI_RECV(v_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 70
            call MPI_RECV(w_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 80
			call MPI_RECV(p_TimeAvg(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
			itag = 90
            call MPI_RECV(		TKE(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )			
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !Wait to receive results from each task
	  
      !----------data collect among nodes for filer----------!
      !end if

           
      if(myid==master)then
	  call acc_set_device_type(acc_device_host)

			call filer3d_avg()						! calculate Kinetic Energy						*** filer.f90 ***

        if (mod(istep,ibackup)==0) then            
            call backupAvg()						! write last data backup in case of power cut	*** filer.f90 ***
        end if 
	  
	  call acc_set_device_num(myid,acc_device_nvidia)
      end if
    
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


      icount=1
      call MPI_BCAST ( VelocityDifference, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !---------------------exit time loop---------------------!
         if ( VelocityDifference < zeta_vel .AND. steadiness==1 ) then; exit; endif
      !--------------------------------------------------------!
!call nvtxEndRange

   end do

!$acc exit data if(correction_stage>0)
!!$acc exit data if(wall_model==1 .OR. Damping_F==1)
!$acc exit data if(calculate_avg==1)
!$acc end data

!--------------------End of main loop on the timesteps----------------------!


   if(myid==master)then
      call filer_final()
   end if

   call MPI_FINALIZE(ierr) 

end program main
