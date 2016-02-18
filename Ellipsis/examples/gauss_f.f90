program gauss
  use gauss_post
  implicit none
  integer, parameter :: ndim=gauss_numDims
  double precision , dimension(ndim) :: st
  double precision , dimension(ndim) :: stp_sz
  double precision scl_fct
  character(len=128), parameter :: fl_pfx="gauss_f"
  integer seed
  integer fb_int
  integer max_stp
  integer resume
  logical diag
  ! intrinsic function trim
  character trim
  integer nburn,nsamp
  integer doMaxLike

  print*, " "
  print*, "================================================="
  print*, "|    Example program in F95: 2                  |"
  print*, "|    Correlated N-D Gaussian                    |"
  print*, "================================================="


  ! diagonal covariance matrix?
  diag=.true.
  doMaxLike = 0

  ! initialise the multivariate Gaussian posterior
  call initialiseGaussData(diag)

  ! start point and step sizes for each parameter
  st=gauss_startPoint
  stp_sz=gauss_stepSz


  scl_fct=1.      ! dimensionality scaling factor (a value between 0 and 1)
                  ! ideal choise is the one which gives you ~68% acceptance rate

  fb_int=1000     ! feed back to console interval

  max_stp=10      ! maximum number of steps in the leapforg (10 is fine)

  resume=0        ! resume from previous run? 0(no) or 1(yes)

  seed=1234       ! random number generator seed

  nburn=0         ! number of samples to be burned

  nsamp=1000      ! number of samples to be drawn
                  ! note that if nsamp>0 the sampler will
                  ! stop after nsamp samples. This doesn't mean that the algorithm
                  ! has converged. Check the hanson values to be sure.

  !
  ! open file for mcmc extract
  ! The accepted samples are written to this file.

  ! openf90 gives me an error if recl is not huge
  ! gfortran,intel,sunf90 are fine
  if(resume .eq. 1) then 
     open(unit=gauss_lunExt,file= (trim(fl_pfx)//".extract.dat"),POSITION='APPEND',&
          recl=1500000)
  else
     open(unit=gauss_lunExt,file= (trim(fl_pfx)//".extract.dat"),recl=1500000)
  endif

  ! call GHS
  call run_guided_hmc(ndim,st,scl_fct,max_stp,stp_sz,fl_pfx,seed,resume,&
       fb_int,negLogPostAndGrad,writeExtract,nburn,nsamp, doMaxLike);

  ! close file
  close(gauss_lunExt)

end program gauss
