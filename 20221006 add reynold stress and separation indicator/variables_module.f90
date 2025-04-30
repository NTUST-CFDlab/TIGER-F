 MODULE variables
    use mpi
    implicit none 

    !This module gathers all structures and variable declarations. 
    !The different variables no need to be redefined and can this file can just be called in the subroutines.


    real*8, parameter            :: Re           = 3900.d0

    real*8, parameter            :: dt           = 3.d0*1.0e-3

    integer, parameter           :: nx           = 36

    integer, parameter           :: ny           = 134

    integer, parameter           :: nz           = 168

    real*8, parameter            :: lx           = 4.d0*DATAN(1.d0)

    real*8, parameter            :: ly           = 20.d0

    real*8, parameter            :: lz           = 20.d0

    character(len=20)            :: Gridder      = 'non-uniform-sin4'  !non-uniform, uniform, non-uniform-sin, non-uniform-sin2, non-uniform-sin3, non-uniform-sin4 (Kuyper)

    !--------------------- Unequal grid ---------------------!

    real*8, parameter            :: GridderYc        = ly/2.d0		! Center coordinate of small grid
    real*8, parameter            :: GridderZc        = 5.6d0

    real*8, parameter            :: lySml            = 4.0
    integer, parameter           :: nySml            = 80

    real*8, parameter            :: lyMid            = 5.0			! for Gridder: non-uniform, non-uniform-sin, non-uniform-sin2
    integer, parameter           :: nyMid            = 30

    real*8, parameter            :: lzSml            = 4.0
    integer, parameter           :: nzSml            = 80 

    real*8, parameter            :: lzMid            = 5.0			! for Gridder: non-uniform, non-uniform-sin, non-uniform-sin2
    integer, parameter           :: nzMid            = 30

    real*8                       :: dySml, dyMid, dy

    real*8                       :: dxSml, dx

    real*8                       :: dzSml, dzMid, dz
	
	real*8,  parameter           :: thickness          = lySml/nySml*1.d0  ! For y+ calculation

    !--------------------- Unequal grid ---------------------!



    
    !---------------------------Dimensional Plasma Parameter (SI)---------------------------!
    !real*8, parameter :: theta                          = 3000.d0          
    !real*8, parameter :: roc                            = 1.0e17
    !real*8, parameter :: poto                           = 5656.85d0          
    !real*8, parameter :: Eb                             = 3.0E6
    !real*8, parameter :: delta_t                        = 6.7E-5  
    !real*8, parameter :: alfa                           = 1.0
    !real*8, parameter :: ec                             = 1.6e-19
    
    !real*8, parameter :: reference_L                    = 0.127d0                                           !!! Chord length 127mm
    !real*8, parameter :: Freestream                     = Re / reference_L * 1.47d0 * 1.0E-5               !!! Chord length 127mm


    !---------------------------Dimensionaless Parameter---------------------------!
    !real*8, parameter :: len_d                          = 0.0075                                            ! Dimensionless Plasma parameter
    !real*8, parameter :: len_a                          = 0.02
    !real*8, parameter :: len_b                          = 0.0118  


    !real*8 :: PlasmaVy1, PlasmaVy2, PlasmaVy3, PlasmaVz1, PlasmaVz2, PlasmaVz3
    !real*8 :: PlasmaZr, PlasmaYr, PlasmaZu, PlasmaYu
    !real*8 :: c1, c2, c3, E0, k1, k2




    !----------------Dynamic airfoil model function---------------------!

    real*8, parameter            :: xc                  = lx/2.d0

    real*8, parameter            :: yc                  = ly/2.d0

    real*8, parameter            :: zc                  = 5.d0

    real*8, parameter            :: reduce_frequency_k  = 0.4d0

    real*8, parameter            :: r                   = 0.5d0

    real*8, parameter            :: PI                  = 4.d0*ATAN(1.d0)

    real*8                       :: AOA, AOA1, AOA_amp, StartDynamic, dis_Y=0.d0, ratio, angular_vel=0.d0


    !----------------Dynamic airfoil model funtion---------------------!



    !----------------------------B.Cs---------------------------!

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
    !   Dirichlet   u = 1
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



    !----------------------------B.Cs---------------------------!



    !-----------------------RayCasting-------------------------!

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)       :: ETA

    integer, parameter                               :: poly = 5040  ! equal POINT

    character(len=50)                                :: NACA_filename

    integer, parameter                               :: nSubGrids = 60

    integer                                          :: A, B, C, E

    real*8, dimension(1:poly)                        :: az, ay, iaz, iay

    integer, dimension(1:nz,1:ny)                    :: points

    integer, dimension(1:nz,1:ny)                    :: intersection, sub_intersection

    real*8                                           :: m_pa, m_ab

    real*8, dimension(1:nz,1:ny)                     :: ETA_1

    integer                                          :: position


    !-----------------------RayCasting-------------------------!



    !-------------------------Plasma---------------------------!

    real*8 ,dimension(1:nx,1:ny,1:nz)                :: F_tavex, F_tavey

    real*8 ,dimension(1:nx,1:ny,1:nz)                :: edelta, EE

    real*8                                           :: PlasmaZc, PlasmaYc

    real*8                                           :: unDBD_cycle

    !-------------------------Plasma---------------------------!



    !-------------------------iteration variable---------------------------!

    integer                                          :: ik, k, i, j, isto, istea, istep, nstep, ccc

    real*8                                           :: time, initial_time

    real*8                                           :: VelocityDifference

    !-------------------------iteration variable---------------------------!
   
    

    !---------------------------Physical variable-----------------------!
	
    real*8, parameter 			 :: U_inf		   		= 1.d0  ! Free stream velocity
	
    real*8, parameter            :: den_flu             = 1.d0
	
    real*8			             						:: nu 
	
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: p, p_old, p_no

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2,1:2)      :: p_pre

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u, v, w
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u1, v1, w1

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u0, v0, w0
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star, v_star, w_star
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: last_velocity

    real*8, dimension(1:nx,1:ny,1:nz)                   :: div,uc,vc,wc,pre

    !---------------------------Physical variable-----------------------!

  
    
    !---------------------------Gauss Seidel-----------------------!

    real*8                                  :: pNew, pChange, omega, pChangeMax, itmax

    real*8                                  :: pChangeMax_

    real*8, dimension(1:nx,1:ny,1:nz)       :: mChange, Den

    real*8, dimension(1:nx,1:ny,1:nz,1:6)   :: P_Den
    !---------------------------Gauss Seidel-----------------------!



    !---------------------------QUICK---------------------------!

    real*8                              :: u_tilde_x1, u_tilde_x2, u_tilde_y1, u_tilde_y2, u_tilde_z1, u_tilde_z2 
    
    real*8                              :: v_tilde_x1, v_tilde_x2, v_tilde_y1, v_tilde_y2, v_tilde_z1, v_tilde_z2
    
    real*8                              :: w_tilde_x1, w_tilde_x2, w_tilde_y1, w_tilde_y2, w_tilde_z1, w_tilde_z2
    
    real*8                              :: ue, uw, un, us, uf, ub, vnu, vsu, wfu, wbu
    
    real*8                              :: ve, vw, vn, vs, vf, vb, uev, uwv, wfv, wbv
    
    real*8                              :: we, ww, wn, ws, wf, wb, uew, uww, vnw, vsw

    !---------------------------QUICK---------------------------!



    !---------------------------BICG----------------------------!

    !integer, parameter                   :: s=nx*ny*nz, nxy=nx*ny
!
    !real*8, dimension(1:3,1:7*s)         :: aa
!
    !real*8, dimension(1:s)               :: r0, bm, rm, p0, uu, xm, xmp, ap, ss, as
!
    !real*8                               :: alpha, om, beta, gamap, sub, gama, apr, asas, ass

    !---------------------------BICG----------------------------!



    !---------------------------LES-----------------------------!

    real*8, parameter                   :: Cs = 0.1d0

    real*8                              :: nut, delta, mutsgs

    real*8, dimension(nx,ny,nz,3)       :: dudx, dvdx, dwdx

    real*8, dimension(nx,ny,nz)         :: Viseff

       !-------------Damping Wall Function-------------------------------!
	                                    
	
	  real*8                              ::    Cf
	  
	  real*8                              ::    tau 
	 
	  real*8                              ::    ufric 
	 
	  real*8                              ::    Yplus 
	 
      real*8                              ::    fwall 

    !---------------------------LES-----------------------------!




    !---------------------------Grid----------------------------!

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

    !---------------------------Grid----------------------------!



    !------------------virtualForceIntegrator------------------!

    real*8                                         :: u_solid = 0
    
    real*8                                         :: v_solid = 0
    
    real*8                                         :: w_solid = 0
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)     :: u2, v2, w2 
    
    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)     :: FX, FY, FZ, FX1, FY1, FZ1
    
    real*8                                         :: totalFX, totalFY, totalFZ
    
    real*8                                         :: cDrag ,cLift

    real*8, dimension(1:nx,1:ny)                   :: FXz, FYz, FZz
    
    real*8, dimension(1:nx)                        :: FXy, FYy, FZy

    !------------------virtualForceIntegrator------------------!



    !--------------------- output ---------------------!
    
    character(len=20)                     :: filename, fileformat

    integer, parameter                    :: nblocks = 1
    
    real, dimension(1:nx,1:ny,1:nz)       :: Xout, Yout, Zout
    
    real, dimension(1:nx,1:ny,1:nz,5)     :: Qout
    
    real                                  :: temp = 1.0    ! mach, alpha, reyn, time 
    
    integer                               :: h,num
    
    !--------------------- output ---------------------!



    !---------------------- input ---------------------!
    
    integer                               :: inblocks
    
    integer                               :: inx
    
    integer                               :: iny
    
    integer                               :: inz 
    
    character(len=20)                     :: inputfile
    
    !---------------------- input ---------------------!


    
    !------------------- OPENMP ------------------------!

    integer                      :: nthreads

    !------------------- OPENMP ------------------------!




    !----------------------------MPI----------------------------!

    integer                      :: nproc, myid, ierr, dest

    integer                      :: status(MPI_STATUS_SIZE)

    integer, parameter           :: master=0

    integer                      :: Zdv, Zr

    integer, dimension(1:1024)   :: gstart ,gend, gend0, gcount

    integer                      :: l_nbr, r_nbr, icount, iend, istart, itag, igcount

    !----------------------------MPI----------------------------!




    !--------------------calculate wall time--------------------!
    
    real*8                   :: totalstarttime

    real*8                   :: totalfinaltime

    real*8                   :: totalcosttime

    real*8                   :: totallasttime

    !--------------------calculate wall time--------------------!



    !-------------------------chooser---------------------------!

    integer                  :: steadiness
    
    integer                  :: LES

    real*8                   :: zeta_vel

    real*8                   :: zeta

    integer                  :: DBD

    real*8                   :: unDBD

    integer                  :: unDBD_onoff

    integer                  :: Prediction

    !-------------------------chooser---------------------------!


    !-------------------time Step Schemer-----------------------!

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star1, v_star1, w_star1

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star2, v_star2, w_star2

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star3, v_star3, w_star3

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_star4, v_star4, w_star4

    real*8, dimension(-1:nx+2,-1:ny+2,-1:nz+2)          :: u_starm, v_starm, w_starm


    !-------------------time Step Schemer-----------------------!

end module
