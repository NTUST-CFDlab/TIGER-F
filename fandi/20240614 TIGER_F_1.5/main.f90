! 22 Aug 2023 - FDS
!    y=1 ______________                                                                                 
!       /             /|       |  Author  : Zi-Hsuan Wei                                                 
!      /             / |       |  Version : 1.8                                                          
!     /____________ /  |       |  Web     : http://smetana.me.ntust.edu.tw/                              
!     |  |         |   |                                                        
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
   implicit none



   !---------------------------- OPENMP ----------------------------!
   nthreads = 20
   call omp_set_num_threads(nthreads)
   !---------------------------- OPENMP ----------------------------!

   !----------------------------- MPI ------------------------------!
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
   call Mpi_division()
   !----------------------------- MPI ------------------------------!
   
   !--------for nz--------!
   istart = gstart(myid)  !
   iend = gend0(myid)     !
   igcount = gcount(myid) !
   !--------for nz--------!

   !----------------------------for sendrecv----------------------------!
   l_nbr = myid - 1                                                     !
   r_nbr = myid + 1                                                     !         
   if(myid == 0) then; l_nbr=MPI_PROC_NULL; endif                       !
   if(myid == (nproc-1)) then; r_nbr=MPI_PROC_NULL; endif               !
   !----------------------------for sendrecv----------------------------!


   !------------------------ Toggle Selector -----------------------!
   
   dimensionality					= 0						! 0) Non-Dimensional, 1) Dimensional

   LES                              = 1                     ! the LES mode. 1 : on ; 0 : off
   
   pressure_solver					= 1						! 1)SOR, 2)RB_SOR
   
   correction_stage					= 0						! number of correction stage in Prediction-correction process

   steadiness                       = 2                     ! steady : 1 ; unsteady : 2

   resume							= 0						! resume simulation. 1:on ; 0:off (Provide the backup file: 'Static0.Q'/'OOXX.Q')
   

   VOS_by							= 1						! Create VOS by 1)geometry function, 2)raycasting2D, 3)raycasting3D

   dat_file	                     	= '1_0.1_30_w.DAT'   	! input for raycasting2D (DAT file)
							scale 	= 1/400.d0 				! geometry scale factor wrt. the input file --> for raycasting2d

   stl_file                         = 'sphere_ratio_005.stl'! input for raycasting3D (STL ASCII file)
   
   select_ref_area					= 4						! Reference area in virtual force calculation
					user_defined	= L_ch*lx				! 1)Z_castingarea, 2)Y_castingarea, 3)X_castingarea, or 4)define manually

!   DBD                              = 0                     ! DBD actuator on : 1 ; DBD actuator off : 0

!   unDBD                            = 0                     ! Strouhal number of brust DBD actuator

   !------------------------ Toggle Selector -----------------------!
   

   !---------------- Parameters for the simulation -----------------!
   
   omega                            = 1.7d0                 ! Set value for SOR method

   zeta                             = 1.0d-4                 ! zeta for solving pressure matrix

   itmax                            = 5000                  ! maximum for zeta in Gauss Seidel subroutines

   initial_time                     = 0.d0                  ! initialize time of simulation

   total_time                       = 200.d0	 			! total time of simulation
    
!   zeta_vel                         = 1.0d-5                 ! steady criteria (if steadiness = 1)

!   istea_int						 = 0.1d0				 ! time interval to check steady state
   
   StartDynamic_time                = 500.d0 				! initialize time of solid motion

   !---------------- Parameters for the simulation -----------------!


   !---------------------Output Data Management---------------------!

   coeff_start						= 0.1d0					! initialize time of storing force & pressure coefficients (CD_time & Cp)

   filer3d							= 0						! create 3D output (3d_xxxx.q). 1 : on ; 0 : off
			startfiler3d_time		= 600.d0				! initialize time of storing filer3d data
			isto3d_int				= 1.d0					! time interval to storing filer3d data
   
   filer2d							= 1						! create 2D output (2d_xxxx.q). 1 : on ; 0 : off (Spanwise avg - Used only in the case of extruded solids)
			startfiler2d_time		= 0.d0					! initialize time of storing filer2d data
			isto2d_int				= 1.d0					! time interval to storing filer2d data
   
   filer_cp							= 0						! create surface pressure distribution plot. 1 : on ; 0 : off (set in filer.f90)
			istocp_int				= 5.d0					! time interval to storing cp data (start = coeff_start)
   
   backup_int						= 1.d0					! time interval to store backup data
   
   !---------------------Output Data Management---------------------!   
   

   totalstarttime = MPI_WTIME()
   totallasttime = totalstarttime
   
  
   !----------------------------------------------------------------!
   call read_data()              ! define dx, dy, dz
   !----------------------------------------------------------------!

  
   !--------------------- Gridder selection ------------------------!
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

   !----------------------------------------------------------------!
   call initial_conditions()        !call initial conditions
   !----------------------------------------------------------------!
	istep = 0

   !-------------------- Loading backup data -----------------------!
   if (resume == 1) then													
	  call read_backupfile()        				!input first q file			--> (used when resuming simulation)					
	  num = NINT(lastsaved_time/(ibackup*dt))		!				
	  nstep = nstep - NINT(lastsaved_time/dt)
	  initial_time = lastsaved_time
	  time = lastsaved_time
   endif
   !----------------------------------------------------------------!

   !----------------------------------------------------------------!
   call final_boundary_conditions()       !call boundary conditions
   !----------------------------------------------------------------!

   !----------------------- Create VOS -----------------------------!

   if (VOS_by == 2) then
		call read_dat()						  ! Read DAT file
		call vos_ray2d()                      ! Create ETA from DAT file 
   elseif (VOS_by == 3) then
		call vos_ray3d()                      ! Create ETA from STL ASCII file 
   elseif (VOS_by == 1) then
		call func_Cylinder()                     ! Create ETA of Cylinder      
!		call func_Sphere()                       ! Create ETA of Sphere

   endif 

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
        
	vosfinaltime = MPI_WTIME()

   !----------------------------------------------------------------!

   
   !----------------------------------------------------------------!
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
   !define plasma force field
   !----------------------------------------------------------------!


   if(myid==master) then

		!call prober()

      if(filer3d==1) then
		call filereachtime3d() 
	  endif
	  if(filer2d==1) then
	    call filereachtime2d()
	  endif
		 num = num + 1	  
      !call filer_bodyforce()
      call filerInfo()
          !---------------------
          call filerProcess_cp() 
          !call filer_Reynoldstress()

   endif


   !call write_matrix_a()

	open (71,file='output.dat',position='append')


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

!-------------------------main loop on the timesteps--------------------------!
   do istep=1,nstep                                                                      

      
      time = initial_time + dt*istep

      
      if( time .GE. StartDynamic_time )then
		 
         !--------------------------------------------------------------------!
         !call twoway_dynamic() 
         !--------------------------------------------------------------------!
		 
         !--------------------------------------------------------------------!
         !call rotor_dynamic()                   ! Change AOA
         !--------------------------------------------------------------------!
		 
         !--------------------------------------------------------------------!
         !call blade_dynamic()                   ! blade rotation
         !--------------------------------------------------------------------!

         !--------------------------------------------------------------------!
		! if (VOS_by == 2) then
				! call vos_ray2d()                      ! Create ETA from DAT file 
		! elseif (VOS_by == 3) then
				! call vos_ray3d()                      ! Create ETA from STL ASCII file 
		! elseif (VOS_by == 1) then	
				! call func_Cylinder()                     ! Create ETA of Cylinder      
		! !		call func_Sphere()                       ! Create ETA of Sphere
		! endif 
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


      !if (unDBD >= 0.01d0 ) then
        ! write(71,*)'unDBD'
        ! call unsteady_plasma()
      !endif

      !if(myid==master)then 
      ! write(71,*) 'time = ',REAL(time) , '  , AOA = ', REAL(AOA), REAL(angular_vel)
	  !	write(71,*) 'Alpha = ', REAL(blade_alpha), '  , Blade angular vel = ', REAL(blade_omega)
      !endif

      !----------data transformation among nodes----------!
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



      ! icount = 2
      ! itag = 170
      ! call MPI_SENDRECV( iDz(istart), icount, MPI_REAL8, l_nbr, itag, &
                         ! iDz(iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      ! itag = 180
      ! call MPI_SENDRECV( Dzs(istart), icount, MPI_REAL8, l_nbr, itag, &
                         ! Dzs(iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )


      ! itag = 190
      ! call MPI_SENDRECV( iDz(iend-1),   icount, MPI_REAL8, r_nbr, itag, &
                         ! iDz(istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      ! itag = 200
      ! call MPI_SENDRECV( Dzs(iend-1),   icount, MPI_REAL8, r_nbr, itag, &
                         ! Dzs(istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                         
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      !------------------------------------------------------------------------------------------------------!
	  if (LES == 1) then
		 call CalculateSmagorinskyViscosity()    ! calculate Smagorinsky Viscosity
	  endif
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call discretisation_QUICK_centre()      ! calculate velocity field
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call AdamsBashforth()                   ! Adams-Bashforth
      !------------------------------------------------------------------------------------------------------!

   
      call u0_boundary_conditions()

      !----------data transformation among nodes----------!
      icount = (nx+4)*(ny+4)
      itag = 240
      call MPI_SENDRECV( w0(-1,-1,iend),     icount, MPI_REAL8, r_nbr, itag, &
                         w0(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      if (correction_stage .GT. 0 .AND. time > 2.d0) then
		do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2
			p(i,j,k)=p_no(i,j,k)
		end do; end do; enddo
	  endif
	  
      !------------------------------------------------------------------------------------------------------!
      ccc=0

      if (pressure_solver == 1) then
       call SOR()				 ! calculate pressure field using traditional SOR (faster)
      elseif (pressure_solver == 2) then
	   call RB_SOR()			 ! calculate pressure field using red-black SOR --> nx MUST be an even number	  
      endif
!      call gauss_seidel()        ! calculate pressure field using traditional SOR (Old)
!      call BICG_stab()           ! calculate pressure field
      !------------------------------------------------------------------------------------------------------!

      if (correction_stage .GT. 0 .AND. time > 2.d0) then
		do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2
			p_no(i,j,k)=p(i,j,k)
		end do; end do; enddo
	  endif



      ! !----------data transformation among nodes----------! --> Moved to GS.f90
      ! icount = (nx+4)*(ny+4)

      ! itag = 270
      ! call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         ! p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! !----------data transformation among nodes----------!
  
      !------------------------------------------------------------------------------------------------------!
      call calcul_new_velocity() ! update velocity field 
      !------------------------------------------------------------------------------------------------------!

      if (correction_stage .GT. 0 .AND. time > 2.d0) then
		call prediction_correction()
      endif


      !------------------------------------------------------------------------------------------------------!
	  ! Force integrator (comment this part if prediction-correction is used)
      		call virtualForceIntegrator_nima() 	! CD and CL only
	  !	    call virtualForceTorqueIntegrator() ! CD, CL, CT, and CP
	  !		call virtualForceTorque_frozen()	! CD, CL, CT, and CP for frozen rotor simulation
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call Updating_velocity()   ! update velocity field 
      !------------------------------------------------------------------------------------------------------!
     
      !------------------------------------------------------------------------------------------------------!
      call final_boundary_conditions() ! recall boundary conditions to update them
      !------------------------------------------------------------------------------------------------------!    



      !----------data collect among nodes for filer----------!

      if (mod(istep,isto3d)==0 .OR. mod(istep,isto2d)==0 .OR. mod(istep,istocp)==0) then        !   modified 29_03_2023
          
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
      end if

           
      if(myid==master)then


        if (time > coeff_start) then
            call filerProcess()									! write CD_time
			!call filerProber()									! point probes
			if(filer_cp==1 .AND. mod(istep,istocp)==0) then            
				call filerProcess_cp()							! write surface pressure distribution
			endif
        endif

         
        if (time > startfiler3d_time-dt .AND. mod(istep,isto3d)==0) then     ! write filer 3D if isto3d = istep
            !call filer_Reynoldstress()
			if(filer3d==1) then
				call filereachtime3d() 
			endif
        end if 

        if (time > startfiler2d_time-dt .AND. mod(istep,isto2d)==0) then     ! write filer 2D if isto2d = istep
			if(filer2d==1) then
				call filereachtime2d() 
			endif
        end if 
		 
        if (mod(istep,ibackup)==0) then            ! write last data backup in case of power cut
            call filerrerun()
			num = num + 1
        end if 

        ! if (mod(istep,istea)==0) then ! write results if isto = istep
          !  call check_steady()
        ! end if

      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     


      !----------calculate wall time----------!
      totalfinaltime = MPI_WTIME()
      totalcosttime = totalfinaltime-totalstarttime
      if(myid==master)then
         write(71,*) 'step cost = ', REAL(totalfinaltime-totallasttime), 'total cost = ', REAL(totalcosttime/3600.) ,'(Hr)'
         write(71,*) 'Cd = ', REAL(cDrag), ', Cl = ', REAL(cLift)
         write(71,*)
      end if
      totallasttime = totalfinaltime
      !----------calculate wall time----------!


      icount=1
      call MPI_BCAST ( VelocityDifference, icount, MPI_REAL8, master,MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !---------------------exit time loop---------------------!
         if ( VelocityDifference < zeta_vel .AND. steadiness==1 ) then; exit; endif
      !---------------------exit time loop---------------------!



   end do
   
!-------------------------main loop on the timesteps----------------------!


   if(myid==master)then
      call filer_final()
   end if

   call MPI_FINALIZE(ierr) 

end program main
