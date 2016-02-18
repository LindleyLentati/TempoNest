module gauss_post

  use linearAlgebraUtils

  implicit none

  integer, parameter :: gauss_numDims=1000
  integer, parameter :: gauss_lunExt=20
  double precision gauss_mean(gauss_numDims),gauss_eigenValues(gauss_numDims)
  double precision gauss_startPoint(gauss_numDims),gauss_stepSz(gauss_numDims)
  double precision gauss_featureMatrix(gauss_numDims,gauss_numDims)
  double precision gauss_inverCovMatrix(gauss_numDims,gauss_numDims)
  contains

    !------------------------------------------------------------------------!
    ! initialise all the data needed for the gauss                           !
    !------------------------------------------------------------------------!
    subroutine initialiseGaussData(diagonal)
      implicit none
      logical diagonal
      double precision, dimension(gauss_numDims,gauss_numDims) :: Mat,MatTrans,posDefCovMat
      integer i,j
      double precision, dimension(gauss_numDims*gauss_numDims) ::  uniReal
      
      ! get uniform random numbers
      call random_number(uniReal);

      ! create a random postitive definite matrix
      do i=1,gauss_numDims
      
         gauss_mean(i)=0.
         
         ! start point for the sampler
         ! For a Gaussian case with mean=0 the peak of the distribution is 0
         gauss_startPoint(i)=0.
         
         do j=1,gauss_numDims
            !uniReal=rand_uni(rng_handle)
            Mat(i,j)=uniReal((i-1)*gauss_numDims+j)! row major?
            MatTrans(j,i)=Mat(i,j)
         enddo
      enddo
      
      ! create the positive definite covariane matrix C=M*M^T
      posDefCovMat=matmul(Mat,MatTrans)
      
      ! if we need a diagonal matrix
      if(diagonal) then
         print*, " "
         print*, "Making the covariance matrix diagonal"
         do i=1,gauss_numDims
            do j=1,gauss_numDims
               if(i .eq. j) then
                  posDefCovMat(i,j)=1.
               else
                  posDefCovMat(i,j)=0.
               endif
            enddo
         enddo
      endif
	
      ! find the inverse of the covariance matrix
      ! This is the Hessian of the -log(posterior)
      ! This gives a gaussian approximation at the peak of 
      ! the distribution. In a Gaussian case it is just the 
      ! inverse of the covariance matrix
      gauss_inverCovMatrix=posDefCovMat
      call SymPosDefMatInv(gauss_inverCovMatrix,gauss_numDims)
      
      ! find its eigen values (PCA)
      ! the feature matrix is a matrix of eigen values
      ! allows us to sample in principal coordinates
      gauss_featureMatrix=gauss_inverCovMatrix
      call SymPosDefEigenDecomp(gauss_featureMatrix,gauss_eigenValues,gauss_numDims)
      
      ! assign step size as the inverse square root of eigen values
      ! inverse square root of the eigen values gives the width of the 
      ! distribution in the principal coordinates
      do i=1,gauss_numDims
         gauss_stepSz(i)=1/sqrt(gauss_eigenValues(i))
      enddo
      
      ! rotate the start point to the principal coordinates
      ! we need to start from the peak in the principal coordinates
      call rotate2Principal(gauss_numDims,gauss_startPoint)
		
    end subroutine initialiseGaussData
    
    !------------------------------------------------------------------------!
    ! get the negative log posterior  and its gradients in one go            !
    !------------------------------------------------------------------------!
    subroutine negLogPostAndGrad(ndim,x,nlpostval,g)
      implicit none
      integer :: i,ndim
      double precision,dimension(ndim) ::x,g,xmu,cinvxmu,xphy,gphy
      double precision nlpostval
      
      if(ndim .ne. gauss_numDims) then 
         write(*,*) "Error @ initialiseGaussData() ! dimensions does not match."
         write(*,*) " "
         stop
      endif
      
      ! rotate x from pricipal coordinates to physical coordinates
      xphy=x
      call rotate2Physical(ndim,xphy)
      
      !(mu-x)
      xmu=(gauss_mean-xphy)
      
      !C^{-1}(mu-x)
      cinvxmu=matmul(gauss_inverCovMatrix,xmu)
      
      ! find the negative log post
      nlpostval=0.
      do i=1,ndim
         nlpostval=nlpostval+cinvxmu(i)*xmu(i)
      enddo
      nlpostval=0.5*nlpostval
      
      !print*, "postval=",nlpostval
      
      !find the gradeint
      gphy=-cinvxmu
      
      !rotate to pricipal
      g=gphy
      call rotate2Principal(ndim,g)
    end subroutine negLogPostAndGrad
    
    !------------------------------------------------------------------------!
    ! rotate a vector x to physical coordinates from principal coordinates   !
    !------------------------------------------------------------------------!
    subroutine rotate2Physical(ndim,x)
      implicit none
      integer ndim
      double precision,dimension(ndim) ::x
      
      if(ndim .ne. gauss_numDims) then 
         write(*,*) "Error @ rotate2Physical() ! dimensions does not match."
         write(*,*) " "
         stop
      endif
      
      x=matmul(gauss_featureMatrix,x)
    end subroutine rotate2Physical
    
    !------------------------------------------------------------------------!
    ! rotate a vector x to principal coordinates from physical coordinates   !
    !------------------------------------------------------------------------!
    subroutine rotate2Principal(ndim,x)
      implicit none
      integer ndim
      double precision,dimension(ndim) ::x
      
      if(ndim .ne. gauss_numDims) then 
         write(*,*) "Error @ rotate2Principal() ! dimensions does not match."
         write(*,*) " "
         stop
      endif
      
      x=matmul(x,gauss_featureMatrix)
    end subroutine rotate2Principal
    
    !------------------------------------------------------------------------!
    ! write the mcmc extract after converting to physical coordinates        !
    !------------------------------------------------------------------------!        
    subroutine writeExtract(ndim,x,neg_log_post,g)
      implicit none
      integer ndim
      double precision neg_log_post
      double precision x(ndim),xphy(ndim),g(ndim)
      integer i
      
      if(ndim .ne. gauss_numDims) then 
         write(*,*) "Error @ writeExtract() ! dimensions does not match."
         write(*,*) " "
         stop
      endif
      
      xphy=x

      ! move to the physical coordinates
      call rotate2Physical(ndim,xphy)
      
      ! write extact
      do i=1,ndim-1
         write(gauss_lunExt,'(E18.10,a)',advance='no') xphy(i), ','
      enddo
      write(gauss_lunExt,'(E18.10)') xphy(ndim)
      
    end subroutine writeExtract
    
  end module gauss_post
  
