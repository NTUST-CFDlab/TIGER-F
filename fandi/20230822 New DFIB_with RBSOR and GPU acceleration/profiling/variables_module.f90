! 22 Aug 2023 - FDS
 MODULE variables
    use mpi
    implicit none 

    real*8, parameter            :: PI           = 4.d0*ATAN(1.d0)
	
    !This module gathers all structures and variable declarations. 
    !The different variables no need to be redefined and can this file can just be called in the subroutines.


    real*8, parameter            :: Re           = 3900.d0		  		  ! Calculate & insert Re manualy in dimensional simulation

    real*8, parameter            :: dt           = 2.5d0*1.0e-3 !8.d0*1.0e-4

    integer, parameter           :: ny           = 248 
			integer, parameter 	 :: nyLow 		 = 50					  !number of grid in lower side for non-uniform-sin4 (from y(0) to nySml)

    integer, parameter           :: nz           = 264
			integer, parameter 	 :: nzUps 		 = 47					  !number of grid in upstream for non-uniform-sin4, sin5 & ground grids (from z(0) to nzSml) 

    integer, parameter           :: nx           = 16 !45				  ! Must be an even number if RBSOR is used

    real*8, parameter            :: ly           = 30.d0

    real*8, parameter            :: lz           = 48.d0

    real*8, parameter            :: lx           = 2.d0 !PI

    character(len=20)            :: Gridder      = 'non-uniform-sin4'  !non-uniform, uniform, non-uniform-sin, non-uniform-sin2, non-uniform-sin3, non-uniform-sin4 (Kuyper)
                                                                          !non-uniform-sin5, non-uniform-sin5-3D, ground-3D

    !--------------------- Unequal grid ---------------------!

    real*8, parameter            :: GridderYc    = 15.d0            ! Center coordinate of small grid, GridderXc is set at lx/2. Set wrt. the geometry xc and yc.
    real*8, parameter            :: GridderZc    = 12.5d0

    real*8, parameter            :: lySml        = 2.22d0
    integer, parameter           :: nySml        = 148

    real*8, parameter            :: lzSml        = 2.22d0
    integer, parameter           :: nzSml        = 148

    real*8, parameter            :: lxSml        = 1.05               ! for Gridder: non-uniform-sin5-3D, ground-3D
    integer, parameter           :: nxSml        = 105

    real*8, parameter            :: lyMid        = 10.0                ! for Gridder: non-uniform, non-uniform-sin, non-uniform-sin2, non-uniform-sin5, non-uniform-sin5-3D, ground-3D
    integer, parameter           :: nyMid        = 10                 ! 'Only' nyMid

    real*8, parameter            :: lzMid        = 10.0                ! for Gridder: non-uniform, non-uniform-sin, non-uniform-sin2, non-uniform-sin5, non-uniform-sin5-3D, ground-3D
    integer, parameter           :: nzMid        = 10                 ! 'Only' nzMid

    real*8, parameter            :: lxMid        = 10.0                ! for Gridder: non-uniform-sin5-3D, ground-3D
    integer, parameter           :: nxMid        = 10                 ! 'Only' nyMid

    real*8                       :: dxSml, dxMid, dx

    real*8                       :: dySml, dyMid, dy

    real*8                       :: dzSml, dzMid, dz

    real*8,  parameter           :: thickness    = lySml/nySml*1.d0  ! For y+ calculation
		
	


    !------------------------------B.Cs-----------------------------!

    !    y=1 ______________                                                                                 
    !       /             /|                                                     
    !      /       N     / |                                                          
    !     /____________ /  |                                
    !     |  |         |   |                                                        
    !     |  | B       |   |                                          
    !   W |  | x=y=z=0 | E |                                           
    !     |  |_________|___|x=1                                        
    !     |  /         |  /                                         
    !     | /     S    | /                                        
    !     |/___________|/                                         
    !    z=1 F     
    !   Neumann     du/dn = 0
    !   Dirichlet   u = U_inf
    !   no-slip     u = 0
    !   Periodic    Un=U1
    !	*Symmetry bc: 'Neumann' for tangential planes + 'no-slip' for normal planes

    character(len=20)            :: WestWall_u         = 'Periodic'
    character(len=20)            :: WestWall_v         = 'Periodic'
    character(len=20)            :: WestWall_w         = 'Periodic'

    character(len=20)            :: EastWall_u         = 'Periodic'
    character(len=20)            :: EastWall_v         = 'Periodic'
    character(len=20)            :: EastWall_w         = 'Periodic'

    character(len=20)            :: SouthWall_u        = 'Neumann'
    character(len=20)            :: SouthWall_v        = 'no-slip'
    character(len=20)            :: SouthWall_w        = 'Neumann'
    
    character(len=20)            :: NorthWall_u        = 'Neumann'
    character(len=20)            :: NorthWall_v        = 'no-slip'
    character(len=20)            :: NorthWall_w        = 'Neumann'

    character(len=20)            :: BackWall_u         = 'no-slip'
    character(len=20)            :: BackWall_v         = 'no-slip'
    character(len=20)            :: BackWall_w         = 'Dirichlet'

    character(len=20)            :: FrontWall_u        = 'Neumann'
    character(len=20)            :: FrontWall_v        = 'Neumann'
    character(len=20)            :: FrontWall_w        = 'Neumann'


    !---------------------Solid reference position------------------!

    real*8, parameter            :: xc                  = lx/2.d0               ! Geometry's origin

    real*8, parameter            :: yc                  = ly/2.d0                	! Geometry's origin

    real*8, parameter            :: zc                  = 12.d0                 ! Geometry's origin

    real*8, parameter            :: offset_z            = 0.d0  ! offset distance from origin to center of rotation -->    o - - - - - - - - - - - - o
    ! 											                                                                     (zc,yc)<----- offset_z (+)----->center of rotation


	!-------------------- domain clipping for filer3d --------------!
	
	real*8, parameter            :: x_start          = 0.d0			! Customize the filer output range
	real*8, parameter            :: x_end            = lx
	
	real*8, parameter            :: y_start          = 0.5*ly-6.d0
	real*8, parameter            :: y_end            = 0.5*ly+6.d0
	
	real*8, parameter            :: z_start          = 10.d0		
	real*8, parameter            :: z_end            = 20.d0

	integer                      :: k_start,k_end,j_start,j_end,i_start,i_end



    !--------------------- Physical variable -----------------------!

    real*8, parameter            :: U_inf              = 1.d0  ! Free stream velocity = inlet velocity (1.d0 for dimensionless simulation)

    real*8, parameter            :: den_flu            = 1.d0  ! Fluid density (1.d0 for dimensionless simulation)
	
    real*8                       :: L_ch			   = 1.d0  ! Characteristic length used in Re calculation (1.d0 for dimensionless simulation)

    real*8                       :: nu, inv_rho, inv_dt 
	

    !---------------- Rotor dynamic model function -----------------!

	real*8, parameter			 :: rotor_tsr			= 0.d0					! + = CCW
	
    real*8, parameter            :: rotor_r             = 0.d0					! rotor radius

    real*8						 :: AOA1                = 0.d0                  ! initial AOA (0 = windward, CCW)

    real*8						 :: AOA_amp             = 0.D0                 ! AOA amplitude

	real*8						 :: rotor_omega

    real*8, parameter            :: reduce_frequency_k  = 0.d0

    real*8                       :: AOA, dis_Y=0.d0, ratio, angular_vel=0.d0


    !--------------- Blade dynamic model function ------------------!

    real*8, parameter            :: blade_alpha			= 0.d0					! Cylinder spin ratio, + = CCW
	
    real*8, parameter            :: blade_r             = 0.5d0					! blade radius
	
	real*8						 :: blade_omega


    !------------------ VOS by Geometry function -------------------!
	
    real*8, parameter            :: r                   = 0.5d0	

    integer, parameter           :: nSubGrids_f 		= 100 !   --> for function
	


    !--------------------- VOS by RayCasting -----------------------!

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: ETA

    real*8                                           :: min_dist = 0.03 ! minimum distance between parts --> for raycasting2d
	
    real*8                                           :: scale
	
	real*8, dimension(:), allocatable 				 :: iaz, iay

	real*8											 :: yc_t, zc_t !transformation of geometry's center
	
	integer											 :: poly

    character(len=50)                                :: dat_file, stl_file ! 							 

    integer, parameter                               :: nSubGrids_2d = 100 !  --> for raycasting2d
	
    integer, parameter                               :: nSubGrids_3d = 10 !   --> for raycasting3d

    integer                                          :: A, B, C, E

    integer, dimension(1:nz,1:ny)                    :: points

    integer, dimension(1:nz,1:ny)                    :: intersection, sub_intersection

    real*8                                           :: m_pa, m_ab

    real*8, dimension(1:nz,1:ny)                     :: ETA_1

    real*8                                           :: Y_castingarea, Z_castingarea, X_castingarea

    integer                                          :: position



    !--------------------------- Plasma ----------------------------!

    real*8 ,dimension(1:nx,1:ny,1:nz)                :: F_tavex, F_tavey

    real*8 ,dimension(1:nx,1:ny,1:nz)                :: edelta, EE

    real*8                                           :: PlasmaZc, PlasmaYc

    real*8                                           :: unDBD_cycle



    !-------------- Dimensional Plasma Parameter (SI) --------------!
    !real*8, parameter :: theta                          = 3000.d0          
    !real*8, parameter :: roc                            = 1.0e17
    !real*8, parameter :: poto                           = 5656.85d0          
    !real*8, parameter :: Eb                             = 3.0E6
    !real*8, parameter :: delta_t                        = 6.7E-5  
    !real*8, parameter :: alfa                           = 1.0
    !real*8, parameter :: ec                             = 1.6e-19
    
    !real*8, parameter :: reference_L                    = 0.127d0                                           !!! Chord length 127mm
    !real*8, parameter :: Freestream                     = Re / reference_L * 1.47d0 * 1.0E-5               !!! Chord length 127mm


    !------------------ Dimensionaless Parameter -------------------!
    !real*8, parameter :: len_d                          = 0.0075                                            ! Dimensionless Plasma parameter
    !real*8, parameter :: len_a                          = 0.02
    !real*8, parameter :: len_b                          = 0.0118  


    !real*8 :: PlasmaVy1, PlasmaVy2, PlasmaVy3, PlasmaVz1, PlasmaVz2, PlasmaVz3
    !real*8 :: PlasmaZr, PlasmaYr, PlasmaZu, PlasmaYu
    !real*8 :: c1, c2, c3, E0, k1, k2
	

    !-------------------- iteration variable -----------------------!

    integer                                          :: ik, k, i, j, rb, isto, istea, istep, nstep, ccc, ibackup
	
    integer                                          :: isto3d, isto2d, istocp
	
    real*8                                           :: isto_int, istea_int, backup_int, coeff_start
	
	real*8                                           :: startfiler3d_time, startfiler2d_time

	real*8                                           :: isto3d_int, isto2d_int, istocp_int

    real*8                                           :: time, initial_time, total_time, StartDynamic_time

    real*8                                           :: VelocityDifference
	
	real*8                                           :: p_spansum, uc_spansum, vc_spansum, wc_spansum
   
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: p, p_old, p_new, p_no

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2,1:2)   :: p_pre

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: u, v, w
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: u1, v1, w1

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: u0, v0, w0
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: u_star, v_star, w_star
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: last_velocity

    real*8, dimension(1:nx,1:ny,1:nz)                :: div,uc,vc,wc,pre
	
    real*8, dimension(1,1:ny,1:nz)                	 :: p_spanavg, uc_spanavg, vc_spanavg, wc_spanavg


  
   
    !------------------------ Gauss Seidel -------------------------!

    real*8                                  :: pNew, pChange, omega, pChangeMax, itmax

    real*8                                  :: pChangeMax_

    real*8, dimension(1:nx,1:ny,1:nz)       :: mChange, Den, Den_inv

    real*8, dimension(1:nx,1:ny,1:nz,1:6)   :: P_Den


 
    !-------------------------- QUICK ------------------------------!

    real*8                              :: u_tilde_x1, u_tilde_x2, u_tilde_y1, u_tilde_y2, u_tilde_z1, u_tilde_z2 
    
    real*8                              :: v_tilde_x1, v_tilde_x2, v_tilde_y1, v_tilde_y2, v_tilde_z1, v_tilde_z2
    
    real*8                              :: w_tilde_x1, w_tilde_x2, w_tilde_y1, w_tilde_y2, w_tilde_z1, w_tilde_z2
    
    real*8                              :: ue, uw, un, us, uf, ub, vnu, vsu, wfu, wbu
    
    real*8                              :: ve, vw, vn, vs, vf, vb, uev, uwv, wfv, wbv
    
    real*8                              :: we, ww, wn, ws, wf, wb, uew, uww, vnw, vsw



    !---------------------------- BICG -----------------------------!

    !integer, parameter                   :: s=nx*ny*nz, nxy=nx*ny
!
    !real*8, dimension(1:3,1:7*s)         :: aa
!
    !real*8, dimension(1:s)               :: r0, bm, rm, p0, uu, xm, xmp, ap, ss, as
!
    !real*8                               :: alpha, om, beta, gamap, sub, gama, apr, asas, ass



    !---------------------------- LES ------------------------------!

    real*8, parameter                   :: Cs = 0.1d0

    real*8                              :: nut, delta, mutsgs

    real*8, dimension(nx,ny,nz,3)       :: dudx, dvdx, dwdx

    real*8, dimension(nx,ny,nz)         :: Viseff

       !-------------Damping Wall Function-------------------------------!
                               
		real*8                            ::    Cf, cf_d
  
		real*8                            ::    tau, tau_d
 
		real*8                            ::    ufric, ufric_d
 
		real*8                            ::    Yplus, yplus_d
 
		real*8                            ::    fwall, fwall_d 


    !---------------------------- Grid -----------------------------!

    !Initial grid coordinates for evaluating grid lengths
    real*8, dimension (-1:nx+3)         :: X

    real*8, dimension (-1:ny+3)         :: Y

    real*8, dimension (-1:nz+3)         :: Z 

    !Grid lengths
    real*8, dimension (-1:nx+2)         :: iDx

    real*8, dimension (-1:nx+2)         :: Dxs

    real*8, dimension (-1:ny+2)         :: iDy

    real*8, dimension (-1:ny+2)         :: Dys

    real*8, dimension (-1:nz+2)         :: iDz

    real*8, dimension (-1:nz+2)         :: Dzs

    !Midpoints of grid coordinates
    real*8, dimension (1:nx)            :: Xs

    real*8, dimension (1:ny)            :: Ys

    real*8, dimension (1:nz)            :: Zs



    !------------------- virtualForceIntegrator --------------------!

    real*8                                         :: u_solid = 0
    
    real*8                                         :: v_solid = 0
    
    real*8                                         :: w_solid = 0
	
	real*8                              		   :: ref_area, area_ref
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)     :: u2, v2, w2 
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)     :: FX, FY, FZ, FX1, FY1, FZ1
    
    real*8                                         :: totalFX, totalFY, totalFZ
	
    real*8                                         :: totalTX, totalTY, totalTorq
    
    real*8                                         :: cDrag ,cLift, cTorq, cPower

    real*8, dimension(1:nx,1:ny)                   :: FXz, FYz, FZz
    
    real*8, dimension(1:nx)                        :: FXy, FYy, FZy




    !--------------------------- output ----------------------------!
    
    character(len=20)                     :: filename, fileformat

    integer, parameter                    :: nblocks = 1
    
    real, dimension(1:nx,1:ny,1:nz)       :: Xout, Yout, Zout
    
    real, dimension(1:nx,1:ny,1:nz,5)     :: Qout
    
    real                                  :: temp = 1.0    ! mach, alpha, reyn, time 
    
    integer                               :: h,num
    


    !--------------------------- input -----------------------------!
    
    integer                               :: inblocks
    
    integer                               :: inx
    
    integer                               :: iny
    
    integer                               :: inz 
	
	real*8								  :: lastsaved_time
    
    character(len=20)                     :: inputfile
    


  
    !--------------------------- OPENMP ----------------------------!

    integer                      :: nthreads
	
	integer, parameter			 :: nclps = 1 ! use 2 collapsed loops if nz/nproc/nthread < 5


    !----------------------------- MPI -----------------------------!

    integer                      :: nproc, myid, ierr, dest

    integer                      :: status(MPI_STATUS_SIZE)

    integer, parameter           :: master=0

    integer                      :: Zdv, Zr

    integer, dimension(1:1024)   :: gstart ,gend, gend0, gcount

    integer                      :: l_nbr, r_nbr, icount, iend, istart, itag, igcount
	
	real						 :: igcount_min




    !--------------------- calculate wall time ---------------------!
    
    real*8                   :: totalstarttime

    real*8                   :: totalfinaltime

    real*8                   :: totalcosttime

    real*8                   :: totallasttime
	
	real*8					 :: time1start, time2start, time3start, time4start, time5start, time6start

	real*8					 :: time1end, time2end, time3end, time4end, time5end, time6end



    !--------------------------- toggle ----------------------------!

    integer                  :: steadiness

	integer					 :: LES

    integer                  :: resume
	
    integer                  :: filer3d, filer2d, filer_cp

    real*8                   :: zeta_vel

    real*8                   :: zeta

    integer                  :: VOS_by
	
    integer                  :: select_ref_area

    real*8                   :: user_defined

    integer                  :: DBD

    real*8                   :: unDBD

    integer                  :: unDBD_onoff

    integer                  :: Prediction



    !--------------------- time Step Schemer -----------------------!

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star1, v_star1, w_star1

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star2, v_star2, w_star2

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star3, v_star3, w_star3

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star4, v_star4, w_star4

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_starm, v_starm, w_starm


end module
