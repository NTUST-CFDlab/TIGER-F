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



   !------------------- OPENMP ------------------------!
   nthreads = 7
   call omp_set_num_threads(nthreads)
   !------------------- OPENMP ------------------------!

   !------------------------ MPI ------------------------!
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
   call Mpi_division()
   !------------------------ MPI ------------------------!
   
   !--------for nz--------!
   istart = gstart(myid)  !
   iend = gend0(myid)     !
   igcount = gcount(myid) !
   !--------for nz--------!



   !-----------------Parameters for the simulation------------------!
   
   omega                            = 1.7d0                ! Set value for SOR method

   zeta                             = 1.e-4                 ! zeta for solving pressure matrix

   itmax                            = 5000                  ! maximum for zeta in Gauss Seidel subroutines

   zeta_vel                         = 1.e-5                 ! convergence condition for velocity field

   initial_time                     = 0.d0                  ! initialize time of simulation

!   AOA1                             = 8.d0                  ! initialize AOA of simulation

!   AOA_amp                          = 10.D0                 ! initialize AOA of simulation

   nstep                            = INT(300.d0/dt)         ! number of timesteps for the simulation

   isto                             = INT(1.d0/dt)          ! data stored every 'isto' steps

   istea                            = INT(0.1d0/dt)         ! data stored every 'isto' steps

   LES                              = 1                     ! the LES mode. 1 : on ; 0 : off

   steadiness                       = 2                     ! steady : 1 ; unsteady : 2 

   totalcosttime                    = 0                     ! initialize wall time
   
   inputfile                        = 'AOA0.Q'              ! input first q file

!   StartDynamic                     = INT(500.d0/dt)        ! 

!   NACA_filename                    = 'NACA0012new.DAT'     ! 

!   DBD                              = 0                     ! DBD actuator on : 1 ; DBD actuator off : 0

!   unDBD                            = 0                     ! Strohul number of brust DBD actuator
   
   VelocityDifference               = 1024
   !-----------------Parameters for the simulation------------------!


   totalstarttime = MPI_WTIME()
   totallasttime = totalstarttime
   
   !--------------------------------------------------------------------!
   call reading_data()              ! define dx, dy, dz, nu = 1./Re
   !--------------------------------------------------------------------!

   
   !--------------------------------------------------------------------!
   if(Gridder=='non-uniform')then
      call gridder_unequal()           ! define unequal mesh 
   else if(Gridder=='non-uniform-sin5')then
      call gridder_sin5()             ! define unequal mesh with two non uniform domain 
   else if(Gridder=='non-uniform-sin4')then
      call gridder_sin4()             ! define unequal mesh proposed by Kuyper
   else if(Gridder=='non-uniform-sin3')then
      call gridder_sin3()             ! define unequal mesh using simple SIN function from 0 to pi/2
   else if(Gridder=='non-uniform-sin2')then
      call gridder_sin2()             ! define unequal mesh
   else if(Gridder=='non-uniform-sin')then
      call gridder_sin()              ! define unequal mesh
   else if(Gridder=='uniform')then
      call gridder_equal()             ! define equal mesh
   end if
   !--------------------------------------------------------------------!
   
   open (71,file='output.dat',position='append')

   !--------------------------------------------------------------------!
   call initial_conditions()        !call initial conditions
   !--------------------------------------------------------------------!

   !--------------------------------------------------------------------!
   !call reading_variables()        !input first q file
   !--------------------------------------------------------------------!

   !--------------------------------------------------------------------!
   call final_boundary_conditions()       !call boundary conditions
   !--------------------------------------------------------------------!

   !--------------------------------------------------------------------!
!   AOA = AOA1
!   call RayCasting()                     ! Create ETA of airfoil         
   !--------------------------------------------------------------------!

   call DFIB_Cylinder()                   ! Create ETA of Cylinder        
   !--------------------------------------------------------------------!
   !--------------------------------------------------------------------!
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
   !--------------------------------------------------------------------!


   if(myid==master)then

      call filereachtime()     
      call filer_bodyforce()
      call filerInfo()
	  !---------------------
	  call filerProcess_cp() 
	  !call filer_Reynoldstress()

   endif

   !call write_matrix_a()
    

   !----------------------------for sendrecv----------------------------!
   l_nbr = myid - 1                                                     !
   r_nbr = myid + 1                                                     !         
   if(myid == 0) then; l_nbr=MPI_PROC_NULL; endif                       !
   if(myid == (nproc-1)) then; r_nbr=MPI_PROC_NULL; endif               !
   !----------------------------for sendrecv----------------------------!

  
   
   write(71,*) myid, 'istart = ', istart, 'iend = ', iend, 'gcount = ', igcount


   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!-------------------------main loop on the timesteps----------------------!
   do istep=1,nstep                                                                      

      
   
      time = initial_time + dt*istep
      
      
      
      !if( istep > StartDynamic )then

         !--------------------------------------------------------------------!
         !call dynamic_AOA()                     ! Change AOA
         !--------------------------------------------------------------------!

         !--------------------------------------------------------------------!
         !call RayCasting()                      ! Creat ETA of airfoil
         !--------------------------------------------------------------------!

         !--------------------------------------------------------------------!
         !if(DBD == 1)then
         !   write(71,*)'DBD'
         !   call plasma()
         !   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         !end if
         !define plasma force field
         !--------------------------------------------------------------------!

      !end if

      !if (unDBD >= 0.01d0 ) then
        ! write(71,*)'unDBD'
        ! call unsteady_plasma()
      !endif

      !if(myid==master)then 
      !   write(71,*) 'time = ',REAL(time) , '  , AOA = ', REAL(AOA), REAL(angular_vel)
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
      call MPI_SENDRECV( u(-1,-1,iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         u(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 150
      call MPI_SENDRECV( v(-1,-1,iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         v(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 160
      call MPI_SENDRECV( w(-1,-1,iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         w(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )



      icount = 2
      itag = 170
      call MPI_SENDRECV( iDz(istart), icount, MPI_REAL8, l_nbr, itag, &
                         iDz(iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 180
      call MPI_SENDRECV( Dzs(istart), icount, MPI_REAL8, l_nbr, itag, &
                         Dzs(iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )


      itag = 190
      call MPI_SENDRECV( iDz(iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         iDz(istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 200
      call MPI_SENDRECV( Dzs(iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         Dzs(istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                         
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      !------------------------------------------------------------------------------------------------------!
      call CalculateSmagorinskyViscosity()    ! calculate Smagorinsky Viscosity
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
      call MPI_SENDRECV( w0(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
                         w0(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      !do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2
!
      !   p(i,j,k)=p_no(i,j,k)
!
      !end do; end do; enddo

      !------------------------------------------------------------------------------------------------------!
      ccc=0
      call gauss_seidel()       ! calculate pressure field
      !call BICG_stab()       ! calculate pressure field
      !------------------------------------------------------------------------------------------------------!


      !do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2
!
      !   p_no(i,j,k)=p(i,j,k)
!
      !end do; end do; enddo



      !----------data transformation among nodes----------!
      icount = (nx+4)*(ny+4)

      itag = 270
      call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      !------------------------------------------------------------------------------------------------------!
      call calcul_new_velocity() ! update velocity field 
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call virtualForceIntegrator_nima() ! COEFFICIENTS: CD and Center_CL,Distance,Velocity for solid  
										 !(comment this line if prediction-correction is used)
      !------------------------------------------------------------------------------------------------------!

      !if (time > 2.d0) then
         !call initial_Prediction()
         !call Prediction_Force()
      !endif

      !------------------------------------------------------------------------------------------------------!
      call Updating_velocity()   ! update velocity field 
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call final_boundary_conditions() ! recall boundary conditions to update them
      !------------------------------------------------------------------------------------------------------!
      



      !----------data collect among nodes for filer----------!

	  if (mod(istep,isto)==0) then        !   modified 17_09_2021
	  
      icount = igcount*(nx+4)*(ny+4)
      !Send my results back to the master
      if(myid>master)then
         itag = 10
         call MPI_SEND( p(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 20
         call MPI_SEND( u(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 30
         call MPI_SEND( v(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 40
         call MPI_SEND( w(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if
      !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !Wait to receive results from each task
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
      !----------data collect among nodes for filer----------!
	  end if

           
      if(myid==master)then

        ! if(mod(istep,iforce_sto)==0) then
           ! call filerProcess()
         !endif

          if (time > 1.0) then
            call filerProcess()
         endif


        if (time > 5.d0 .AND. mod(istep,isto)==0) then
            call filerProcess_cp()
        endif
         
         if (mod(istep,isto)==0) then           ! write results if isto = istep
			!call filer_Reynoldstress()
			call filereachtime()
         end if 

        ! if (mod(istep,istea)==0) then ! write results if isto = istep
          !  call check_steady()
        ! end if

         ! if (time > 0.1d0 .AND. cDrag > 2.0 ) then ! write results if isto = istep
         ! !if (cDrag > 0.2 ) then ! write results if isto = istep

         !      !$OMP PARALLEL DO PRIVATE(i,j)  
         !      do k=1,nz; do j=1,ny; do i=1,nx
         !         uc(i,j,k) = 0.5d0*(u(i,j,k)+u(i-1,j,k))
         !         vc(i,j,k) = 0.5d0*(v(i,j,k)+v(i,j-1,k))
         !         wc(i,j,k) = 0.5d0*(w(i,j,k)+w(i,j,k-1))
         !      enddo; enddo; enddo
         !      !$OMP END PARALLEL DO

         !      !$OMP PARALLEL DO PRIVATE(i,j)  
         !      do k=1,nz; do j=1,ny; do i=1,nx
         !         Qout(i,j,k,1)=p(i,j,k)
         !         Qout(i,j,k,2)=uc(i,j,k)
         !         Qout(i,j,k,3)=vc(i,j,k)
         !         Qout(i,j,k,4)=wc(i,j,k)
         !         Qout(i,j,k,5)=ETA(i,j,k)
         !      enddo; enddo; enddo
         !      !$OMP END PARALLEL DO

	      !      write(filename,'(I7.7)')INT(time*1.0E+5)
 	      !      fileformat = '.q'

         !      open (17,file=TRIM(filename)//fileformat,position='append',form='unformatted')
         !      write(17) nblocks
         !      write(17) nx, ny, nz
         !      write(17) temp, temp, temp, REAL(time)
         !      write(17) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )

         !      close(17)
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
   istep = istep - 1
!-------------------------main loop on the timesteps----------------------!


   if(myid==master)then
      call filer_final(); call filereachtime()
   end if

   call MPI_FINALIZE(ierr) 

end program main
