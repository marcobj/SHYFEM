!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_para

  integer, save :: rmode = 13 ! Ensemble Kalman Filter with SVD pseudo inversion of SS'+ EE'
  !22 Square root algorithm with SVD pseudo inversion of SS'+(N-1)R
  !23 Square root algorithm with SVD pseudo inversion of SS'+ EE'
  !10 exact update scheme for diagonal obs-err-cov-mat

  logical, parameter :: verbose = .true. ! Prints diagnostic output

  integer, parameter :: inflate = 2  ! Inflation
                                     ! 0 = off
				     ! 1 = multiplicative
				     ! 2 = adaptive according to Evensen 2009
  real, parameter :: infmult = 1.    ! If inflate=1 -> infmult in the range 1.01-1.1 may be reasonable inflation factors 
		                                       ! dependent on the ensemble size.
						       ! If inflate=2 -> you can set it to 1.0 since infmult only serves as an adjustment of
						       ! the adaptive multiplication factor
                                   

  ! Set these parameters for local analysis. Important.				       
  !
  integer, save :: is_local = 0 !Local analysis. 0 disable, 1 local analysis. Specify radii in the info files.

  ! set this to 1 to include the model errors in the analysis (to test) or equal 2 to
  ! include model parameters (todo)
  !
  integer, parameter :: mode_an = 0

  ! This coefficient limits the maximum for the innovations that are too large
  real, parameter :: inn_alpha = 5.

!------------
! Settings for the creation of the initial ensemble
  integer,parameter :: fmult_init = 10             !mult factor to determine the supersampling
  real, parameter :: theta_init = 0.              !rotation of the random fields (0 East, anticlockwise)
  real, parameter :: sigma_init_z = 0.03               !standard deviation of zeta
  real, parameter :: sigma_init_t = 1.                 !standard deviation of temperature
  real, parameter :: sigma_init_s = 1.                 !standard deviation of salinity
  logical, parameter ::  sample_fix_init = .true.

!
!------------
! Settings for the observation perturbations (mod_enkf)
!
  ! decay time for the red noise of the observations (sec). 
  ! Set lower than 0 to disable (white noise)
  !
  double precision, parameter :: TTAU_0D = -1
  double precision, parameter :: TTAU_2D = 3*3600.

!------------
! Settings to manage the observations (mod_manage_obs)
!
  ! multiplication factor for the ens std, to resize
  ! the obs std (default = 2). Set <= 0 to disable it.
  ! See Sakov et al. 2012.
  !
  real, parameter :: KSTD = 2

  ! parameters for the analysis
  real, parameter :: truncation = 0.99 ! truncation of the SVD eigenvalues
  logical, parameter :: lrandrot = .false. ! True if additional random rotation is used
  logical, parameter :: lupdate_randrot = .true. 
  !logical, parameter :: lupdate_randrot = .false. ! False to avoid randomness
  logical, parameter :: lsymsqrt = .true. ! true if symmetrical square root of Sakov is used (should be used)

  ! time interval (sec) to select an observation
  ! with respect to the analysis time
  !
  double precision, parameter :: TEPS = 300.

  ! standard flag for a missing observation
  !
  real, parameter :: OFLAG = -999.

  ! min-max values for the observation and model check
  !
  real, parameter :: TEM_MIN = -20.0d0
  real, parameter :: TEM_MAX = 60.0d0
  real, parameter :: SAL_MIN = -0.5d0
  real, parameter :: SAL_MAX = 60.0d0
  real, parameter :: SSH_MIN = -8.0d0
  real, parameter :: SSH_MAX = 8.0d0
  real, parameter :: VEL_MIN = -90000.0d0
  real, parameter :: VEL_MAX = 90000.0d0	!used also for u and v, min with -

end module mod_para
