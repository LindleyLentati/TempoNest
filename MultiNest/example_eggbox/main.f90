program main

	use params
      	use nestwrapper
        use RandomNS
      
	implicit none
      
      
	integer i
        call InitRandomNS(1,100)

	pi = 4d0*atan(1d0)
        allocate(datvec(100))
        do i = 1, 100
                datvec(i) = 0*sin(2*pi*i/100.0)+Gaussian1NS(0)*100
        end do	
               call killRandomNS      
	!setting priors
	spriorran(1:sdim,1)=0d0
      	spriorran(1:sdim,2)=10d0*pi
      
      	!no parameters to wrap around
      	nest_pWrap=0

      	call nest_Sample
      	stop
end
