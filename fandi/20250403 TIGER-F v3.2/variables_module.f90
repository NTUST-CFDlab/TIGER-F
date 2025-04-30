! 12 Jan 2025 - FDS
	
!This module gathers all structures and variable declarations. 

MODULE variables
    use mpi
    implicit none 

    real*8, parameter            :: PI           = 4.d0*ATAN(1.d0)


    !--------------------- Physical parameters ---------------------!

    real*8, parameter            :: U_inf              = 1.d0   ! Free stream velocity = inlet velocity (1.d0 for dimensionless simulation)

    real*8, parameter            :: den_flu            = 1.d0  ! Fluid density (1.d0 for dimensionless simulation)
	
	real*8, parameter            :: den_sol            = 100.d0  ! Solid density (only needed in 2way)
	
	real*8, parameter            :: i_sol              = 500.d0  ! moment of inertia of solid (only needed in 2way)
	
    real*8, parameter            :: L_ch			   = 1.d0  ! Characteristic length used in Re calculation (1.d0 for dimensionless simulation)


    real*8, parameter            :: Re				   = 3900.d0	! Reynolds number (for a dimensionless simulation)
	real*8						 :: nu= (U_inf*L_ch)/Re				! comment if dimensional
	
    ! real*8, parameter            :: nu           	   = 1.d-4	! Kinematic viscosity (for a dimensional simulation)
	! real*8					     :: Re= (U_inf*L_ch)/nu					! comment if dimensionless

    real*8, parameter            :: dt           	   = 3.d-3
	
	real*8, parameter 			 :: g                  = 981.d0

    !----------------------- Gridder settings ----------------------!
	
    integer, parameter           :: ny           = 124
			integer, parameter 	 :: nyLow 		 = 37					  !number of grid in lower side for non-uniform-sin4 (from y(0) to nySml)

    integer, parameter           :: nz           = 183
			integer, parameter 	 :: nzUps 		 = 37					  !number of grid in upstream for non-uniform-sin4, sin5 & ground grids (from z(0) to nzSml) 

    integer, parameter           :: nx           = 12    				  !Must be an even number if RBSOR is used
			integer, parameter 	 :: nxLow 		 = 2					  !number of grid in lower side for non-uniform-fall

    real*8, parameter            :: ly           = 24.d0

    real*8, parameter            :: lz           = 36.d0

    real*8, parameter            :: lx           = 3.d0

    character(len=30)            :: Gridder      = 'non-uniform-sin4'  	!non-uniform, uniform, non-uniform-sin, non-uniform-sin2, non-uniform-sin3, non-uniform-sin4 (Kuyper)
                                                                        !non-uniform-sin5, non-uniform-sin5-3D, ground-3D, non-uniform-fall, non-uniform-fall-wall,read_txt (format: nx+1 ny+1 nz+1)

    !--------------------- Unequal grid ---------------------!

    real*8, parameter            :: GridderYc    = 12.d0            ! Center coordinate of fine grid, GridderXc is set at lx/2. Set wrt. the geometry x0 and y0.
    real*8, parameter            :: GridderZc    = 12.9d0

    real*8, parameter            :: lySml        = 1.5d0
    integer, parameter           :: nySml        = 50 ! 1423

    real*8, parameter            :: lzSml        = 3.d0
    integer, parameter           :: nzSml        = 100 ! 1423

    real*8, parameter            :: lxSml        = 1.05d0             ! for Gridder: non-uniform-sin5-3D, ground-3D
    integer, parameter           :: nxSml        = 105

    real*8, parameter            :: lyMid        = 2.d0               ! for Gridder: non-uniform, non-uniform-sin, non-uniform-sin2, non-uniform-sin5, non-uniform-sin5-3D, ground-3D
    integer, parameter           :: nyMid        = 4                 ! 'Only' nyMid

    real*8, parameter            :: lzMid        = 2.d0              ! for Gridder: non-uniform, non-uniform-sin, non-uniform-sin2, non-uniform-sin5, non-uniform-sin5-3D, ground-3D
    integer, parameter           :: nzMid        = 4                 ! 'Only' nzMid

    real*8, parameter            :: lxMid        = 2.d0              ! for Gridder: non-uniform-sin5-3D, ground-3D
    integer, parameter           :: nxMid        = 4                  ! 'Only' nxMid

    real*8                       :: dxSml, dxMid, dx

    real*8                       :: dySml, dyMid, dy

    real*8                       :: dzSml, dzMid, dz
	
	


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
	!	moving_ref	u = Usolid
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

	real*8, parameter            :: p_ref			   = 0.d0	! Set pressure reference (at the domain's exit flow p(nx,ny,nz+1))
	

    !---------------------Solid reference position------------------!

    real*8, parameter            :: x0                  = lx/2.d0               ! Geometry's origin

    real*8, parameter            :: y0                  = ly/2.d0              	! Geometry's origin

    real*8, parameter            :: z0                  = 12.d0                 ! Geometry's origin

    real*8, parameter            :: offset_z            = 0.d0  ! offset distance from origin to center of rotation -->    o - - - - - - - - - - - - o
    ! 											                                                                     (z0,y0)<----- offset_z (+)----->center of rotation


    !--------------- Blade dynamic model function ------------------!

    real*8, parameter            :: blade_alpha			= 0.d0					! Cylinder spin ratio, + = CCW
	
    real*8, parameter            :: blade_r             = 0.d0					! blade radius

	real*8						 :: blade_omega


    !---------------- Rotor dynamic model function -----------------!

	real*8, parameter			 :: rotor_tsr			= 0.d0					! + = CCW
	
    real*8, parameter            :: rotor_r             = 0.d0					! rotor radius

    real*8						 :: AOA1                = 0.d0                  ! initial AOA (0 = z+, + = CCW)

    real*8						 :: AOA_amp             = 0.d0                  ! AOA amplitude

	real*8						 :: rotor_omega

    real*8, parameter            :: reduce_frequency_k  = 0.d0

    real*8                       :: AOA, dis_Y=0.d0, ratio, angular_vel=0.d0

	real*8            			 :: T_gen	= 0.54d0*(0.5d0*den_flu*2.d0*(rotor_r+blade_r)*lx*U_inf*U_inf*(rotor_r+blade_r)) ! Generator torque (set constant for self start analysis)
	
	real*8						 :: P_gen	= 0.4d0 *(0.5d0*den_flu*2.d0*(rotor_r+blade_r)*lx*U_inf**3)	! Generator power (set constant for steady analysis)



    !--------------- VIV of a cylinder ------------------!

    ! real*8 :: Ur_star 			= (Re - 900.d0)*16.d0/(15000.d0 - 900.d0) + 1.d0
	
	! real*8 :: k_spring 			= (2.d0 * PI / ((Re - 900.d0)*16.d0/(15000.d0 - 900.d0) + 1.d0)) ** 2
	
	! real*8 :: c_damper 			= 4.d0 * PI * 0.0054d0 / ((Re - 900.d0)*16.d0/(15000.d0 - 900.d0) + 1.d0)



    !---------------------Solid motion------------------!

	real*8						 :: x0_t, y0_t, z0_t !transformation of geometry's center due to 2way

	real*8						 :: x0_t1, y0_t1, z0_t1 !transformation of geometry's center
	
	real*8						 :: x0_t2, y0_t2, z0_t2 !transformation of geometry's center

	real*8						 :: x0_t3, y0_t3, z0_t3 !transformation of geometry's center
	
    real*8                       :: u_solid, u_solid1, du_solid
    
    real*8                       :: v_solid, v_solid1, dv_solid
    
    real*8                       :: w_solid, w_solid1, dw_solid
	
	real*8						 :: sol_speed



	!-------------------- domain clipping for filer3d --------------!
	
	real*8, parameter            :: x_start          = 0.d0			! Customize the filer output range
	real*8, parameter            :: x_end            = lx
	
	real*8, parameter            :: y_start          = 0.d0		!y0-6.d0
	real*8, parameter            :: y_end            = ly		!y0+6.d0
	
	real*8, parameter            :: z_start          = 0.d0		!z0-5.d0		
	real*8, parameter            :: z_end            = lz		!z0+15.d0

	integer                      :: k_start,k_end,j_start,j_end,i_start,i_end


    !-------------------- Common VOS parameters -----------------------!

    real*8, dimension(:,:,:), allocatable, pinned   :: ETA
	
	integer          			 :: iBgnVOS, iEndVOS , jBgnVOS, jEndVOS , kBgnVOS, kEndVOS	
	
    real*8                       :: Y_castingarea, Z_castingarea, X_castingarea, solid_volume
	

    !------------------ VOS by Geometry function ----------------------!

    real*8, parameter            :: r                   = 0.5d0	

    integer, parameter           :: nSubGrids_f 		= 100 !   --> for function

	integer, dimension(1:nz,1:ny)      	:: ETA_2D
	
	integer, dimension(1:nz,1:ny,1:nx)  :: ETA_3D
	

    !--------------------- VOS by RayCasting 2D -----------------------!

    integer, parameter                :: nSubGrids_2d = 100 !  --> for raycasting2d

    real*8                            :: min_dist = 2.d0*lySml/nySml ! minimum distance between parts --> for raycasting2d
	
    real*8                            :: z_length2D
	
	real*8, dimension(:), allocatable :: iaz, iay
	
	integer							  :: poly

    character(len=50)                 :: dat_file						 

    integer                           :: A, B, C, E

    integer, dimension(1:nz,1:ny)     :: points

    integer, dimension(1:nz,1:ny)     :: intersection, sub_intersection

    real*8                            :: m_pa, m_ab

    real*8, dimension(1:nz,1:ny)      :: ETA_1

    integer                           :: position


    !--------------------- VOS by RayCasting 3D -----------------------!

    integer, parameter                :: nSubGrids_3d = 10 !   --> for raycasting3d

	integer*4                         :: nf=99999999		! max limit number of facets

	real*8                            :: z_length3D
	
	real*8                            :: ETAs
	
    character(len=50)                 :: stl_file
	
	real*4  ,allocatable              :: iFN(:,:)
	
	real*4  ,allocatable              :: iPx(:), iPy(:), iPz(:)


    !--------------------------- Wall function ---------------------!
	
	real*8, parameter	:: Kappa = 0.41d0 			! Von karman's constant

	real*8, parameter	:: E_wm = 9.8d0				! An empirical wall coefficient (for smooth walls)
	
	real*8, parameter	:: invE = 1.d0/9.8d0		! An empirical wall coefficient (for smooth walls)
	
	real*8, parameter	:: ROOTVSMALL = 1.d-20 		! to avoid a division by zero
	
	integer, parameter	:: NRmaxIter = 10 			! Max number of Newton-Raphson interation
	
	real*8, parameter	:: NRtolerance = 0.01d0 	! tolerance for Newton-Raphson interation

	integer 										 :: wp_k, void, wp_j_up, wp_j_low, k_wp, j_wp, j_wpp
	
	real*8											 :: wdistance,wdistance_,max_dist,wsin_,wcos_,hyp
	
	!real*8, dimension(:), allocatable   			 :: wp_up,wp_low,z_wp

	!real*8, dimension(:,:), allocatable 			 :: y_wp_up, y_wp_low

    real*8, parameter   :: c_ref = 3.d0
	
    real*8, parameter   :: yplusLam = 10.d0

	real*8, dimension(1:ny,1:nz)					 :: wdist,wsin,wcos,ypluss,upluss,fdamp
	
	integer 										 :: NRit,countYplus,countYplus_
	
	real*8											 :: NRerr,Ucp,Uwp,Urel,uTau,uTau_	! Urel is the relative velocity wrt wall
	
	real*8											 :: kUu,fkUu,NR_f,NR_df
	
	real*8											 :: yplusss,yplusMAX,yplusMIN
	
	real*8											 :: yplusss_,yplusMAX_,yplusMIN_

    real*8, dimension(:,:,:), allocatable, pinned 	 :: yplus_print,uplus_print,dist_print
	

    ! !--------------------------- Plasma ----------------------------!

    ! real*8 ,dimension(1:nx,1:ny,1:nz)                :: F_tavex, F_tavey

    ! real*8 ,dimension(1:nx,1:ny,1:nz)                :: edelta, EE

    ! real*8                                           :: PlasmaZc, PlasmaYc

    ! real*8                                           :: unDBD_cycle



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

    integer                                          :: ik, k, i, j, rb, isto, istea, istep, nstep, ccc, ibackup, istep2
	
    real*8                                           :: time, initial_time, total_time, StartDynamic_time
	
    real*8                       	    			 :: inv_dt = 1.d0/dt

    real*8                                           :: VelocityDifference, Re_t
	
	real*8                                           :: p_spansum, uc_spansum, vc_spansum, wc_spansum
   
    real*8, dimension(:,:,:), allocatable, pinned    :: p   
    !real*4, dimension(:,:,:), allocatable, pinned    :: p  					!fp32
	
	real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: p_old, p_new, p_no
	!real*4, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: p_old, p_new, p_no 	!fp32

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2,1:2)   :: p_pre
    !real*4, dimension(-1:nx+2,-1:ny+2,-1:nz+2,1:2)   :: p_pre					!fp32

    real*8, dimension(:,:,:), allocatable, pinned    :: u, v, w
    
    real*8, dimension(:,:,:), allocatable, pinned    :: u0, v0, w0
	
    real*8, dimension(:,:,:), allocatable, pinned    :: u_TimeAvg, v_TimeAvg, w_TimeAvg, TKE, p_TimeAvg

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: u1, v1, w1
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: u2, v2, w2 

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: u_star, v_star, w_star
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: last_velocity

    real*8, dimension(1:nx,1:ny,1:nz)                :: div,uc,vc,wc,pre
	
    real*8, dimension(1,1:ny,1:nz)                	 :: p_spanavg, uc_spanavg, vc_spanavg, wc_spanavg

	real*8, dimension(1,1:ny,1:nz)     				 :: KE_spanavg, KE_timeavg, KE_timestd  
   
    !------------------------ SOR -------------------------!

    integer                                 :: itmax

    real*8                                  :: pNew, pChange, omega, pChangeMax, pChangeMax_
    !real*4                                  :: pNew, pChange, omega, pChangeMax, pChangeMax_	!fp32

    real*8, dimension(1:nx,1:ny,1:nz)       :: mChange, Den, Den_inv
    !real*4, dimension(1:nx,1:ny,1:nz)       :: mChange, Den, Den_inv							!fp32

    real*8, dimension(1:nx,1:ny,1:nz,1:6)   :: P_Den
    !real*4, dimension(1:nx,1:ny,1:nz,1:6)   :: P_Den											!fp32

	real*8 									:: omega_max, omega_min, d_omega, log_zeta, k_zeta_zeta, x0_sigmoid
	
	real*8, parameter 						:: k_zeta = 5.d0		! zeta multiplier to determine the zeta upper bound in dynamic omega
! 																	  set k_zeta to 0.5 for static omega = omega_max

	real*8 									:: log_k_zeta = LOG(k_zeta)

	real*8 									:: k_sigmoid = 5.d0		! Sigmoid steepness parameter

 
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

    real*8                              :: strainRateTensor,delta, Viseff, idds

    real*8, dimension(nx,ny,nz,3)       :: dudx, dvdx, dwdx

    real*8, dimension(:,:,:), allocatable, pinned    :: nut

       !-------------Wall Damping Function-------------------------------!
                               
		real*8                            ::    Cf
  
		real*8                            ::    tau
 
		real*8                            ::    ufric
 
		real*8                            ::    Yplus
 
		real*8                            ::    fwall 


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
	
	real*8                              		   :: ref_area, area_ref
    
    real*8, dimension(:,:,:), allocatable, pinned  :: FX, FY, FZ

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)     :: FX1, FY1, FZ1
    
    real*8                                         :: totalFX, totalFY, totalFZ

    real*8                                         :: totalFX_, totalFY_, totalFZ_
	
    real*8                                         :: totalTX, totalTY, totalTZ, totalTorq

    real*8                                         :: totalTX_, totalTY_, totalTZ_, totalTorq_
    
	real*8                                         :: totalTorqx, totalTorqy, totalTorqz
	
	
	real*8                                         :: totalT_XZ, totalT_XY
	
	real*8                                         :: totalT_YZ, totalT_YX
	
	real*8                                         :: totalT_ZY, totalT_ZX
	
	real*8                                         :: totalT_XZ_, totalT_XY_
	
	real*8                                         :: totalT_YZ_, totalT_YX_
	
	real*8                                         :: totalT_ZY_, totalT_ZX_


    real*8                                         :: cDrag ,cLift, cTorq, cPower

    real*8, dimension(1:nx,1:ny)                   :: FXz, FYz, FZz
    
    real*8, dimension(1:nx)                        :: FXy, FYy, FZy

	real*8                                         :: rotate_sx, rotate_sy, rotate_sz


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
	
	integer, parameter			 :: nclps = 1 ! use 2 collapsed loops if nz/nproc/nthread < 5 (Only for OpenMP)


    !-------------------------- MPI (MAIN) -------------------------!

    integer                      :: nproc, myid, ierr, dest

    integer                      :: status(MPI_STATUS_SIZE)

    integer, parameter           :: master=0

    integer, dimension(:), allocatable :: gstart ,gend, gend0, gcount

    integer                      :: Zdv, Zr

    integer                      :: l_nbr, r_nbr, icount, iend, istart, itag, igcount
	
	real						 :: igcount_min


    !-------------------------- MPI (VoS) -------------------------!
	
    integer, dimension(:), allocatable :: gstart_vos ,gend_vos, gend0_vos, gcount_vos

    integer                      :: Zdv_vos, Zr_vos

    integer                      :: iend_vos, istart_vos, igcount_vos	




    !--------------------- calculate wall time ---------------------!
    
    real*8                   :: totalstarttime, vosstarttime

    real*8                   :: totalfinaltime, vosfinaltime

    real*8                   :: totalcosttime

    real*8                   :: totallasttime
	
	real*8					 :: time1start, time2start, time3start, time4start, time5start, time6start

	real*8					 :: time1end, time2end, time3end, time4end, time5end, time6end



    !--------------------------- toggle ----------------------------!

    integer                  :: steadiness

	integer					 :: LES

    integer                  :: wall_model
	
    integer                  :: Damping_F

    integer                  :: resume

    real*8                   :: zeta_vel

    real*8                   :: zeta

    integer                  :: VOS_by
	
    integer                  :: select_ref_area

    real*8                   :: user_defined

    integer                  :: solid_motion

    integer                  :: DBD

    real*8                   :: unDBD

    integer                  :: unDBD_onoff
	
    integer                  :: pressure_solver	

    integer                  :: correction_stage
	
    real*8                   ::	StartCorrection_time
	
	integer					 :: dimensionality

    integer                  :: filer3d, filer2d, filer_cp

    integer                  :: filer_virtual_force,filer_running_avg

    integer                  :: isto3d, isto2d, istocp, istoVF, istoAvg
	
    real*8					 :: isto_int, istea_int, backup_int, coeff_start
	
	real*8                   :: startfiler3d_time, startfiler2d_time, startfiler_VF_time, startfilerAvg_time

	real*8                   :: isto3d_int, isto2d_int, istocp_int, istoVF_int, istoAvg_int


    !--------------------- time Step Schemer -----------------------!

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star1, v_star1, w_star1

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star2, v_star2, w_star2

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star3, v_star3, w_star3

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star4, v_star4, w_star4

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_starm, v_starm, w_starm



end module
