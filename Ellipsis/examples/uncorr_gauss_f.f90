program uncorr_gauss

  implicit none

  integer, parameter :: ndim=1000
  double precision , dimension(ndim) :: st
  double precision , dimension(ndim) :: stp_sz
  double precision scl_fct
  character(len=128), parameter :: fl_pfx="uncorr_gauss_f"
  integer seed
  integer fb_int
  integer max_stp
  integer resume
  integer i
  integer nburn,nsamp
  external nd_uncorr_gauss_neg_log_post
  external write_gauss_ghs_extract
  integer doMaxLike


  character trim

  print*, " "
  print*, "================================================="
  print*, "|    Example program in F95: 1                  |"
  print*, "|    Uncorrelated N-D Gaussian                  |"
  print*, "================================================="

  ! start point and step sizes for each parameter
  ! We need to start from a point as close to the 
  ! peak as possible. Here it is (0,0,...,0).
  ! Step-sizes should be approximately the width of
  ! the posterior distributiion of the parameter.
  ! Here it is (1,1,...,1).

  do i=1,ndim
     st(i)=0.
     stp_sz(i)=1.
  enddo


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

  doMaxLike = 0

  !
  ! open file for mcmc extract
  ! The accepted samples are written to this file.
  ! Note that the unit=20 should same as the unit
  ! opened in write_gauss_ghs_extract() subroutine
  ! Otherwise undefined behavior!
  !
 
  ! open the extract file in the right way
  ! openf90 gives me an error if recl is not huge
  ! gfortran,intel,sunf90 are fine
  if(resume .eq. 1) then
     open(unit=20,file= (trim(fl_pfx)//".extract.dat"),POSITION='APPEND',&
          recl=1500000)
  else
     open(unit=20,file= (trim(fl_pfx)//".extract.dat"),recl=1500000)
  endif

  ! call GHS
  call run_guided_hmc(ndim,st,scl_fct,max_stp,stp_sz,fl_pfx,seed,resume,&
       fb_int,nd_uncorr_gauss_neg_log_post,write_gauss_ghs_extract,nburn,nsamp, doMaxLike);

  ! close file
  close(20)

end program uncorr_gauss


!
! A subroutine for calculating n-d gaussian posteior
! The sampler requires -log(posterior)=Psi
! ndim -> number of dimensions
! x    -> input vector at which the posterior is be calculated
! g    -> gradient of the -log(posterior) = dPsi/dx
!
subroutine nd_uncorr_gauss_neg_log_post(ndim,x,v,g)
  implicit none
  integer ndim
  double precision, dimension(ndim) :: x
  double precision, dimension(ndim) :: g
  double precision v
  integer i

  v=0
  do i=1,ndim
     v=v+x(i)*x(i)
     g(i)=x(i)
  enddo
  v=0.5*v
end subroutine nd_uncorr_gauss_neg_log_post

!
! A subroutie to write the extract file
! I have assumed that the unit=20 is opened for
! writing (append) earlier. In general only write
! those parametes which are estimated. The files 
! can be really big depending on the dimensionality
! 
subroutine write_gauss_ghs_extract(ndim,x,v,g)
  implicit none
  integer ndim
  double precision, dimension(ndim) :: x
  double precision, dimension(ndim) :: g
  double precision v
  integer i

  do i=1,ndim-1
     write(20,'(E18.10,a)',advance='no') x(i), ','
  enddo
  write(20,'(E18.10)') x(ndim)

end subroutine write_gauss_ghs_extract
