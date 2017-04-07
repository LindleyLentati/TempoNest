module like
	
use params
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
	implicit none
      
	double precision Cube(nest_nPar),slhood, Datavector(21), Amp, prior
	integer i, minE, maxE
	
	slhood = 0d0
	
        !Cube(1) = -10*Cube(1)+5
        !Amp = 10d0**Cube(1)	
        

        Cube(1) = 1000*Cube(1)
        Amp = Cube(1)	
        
        prior = 0.5*(9.5-Amp)**2/5**2
        do i = 1, 100
                         slhood = slhood + 0.5d0*(datvec(i) - Amp*sin(2*pi*i/100.0))**2/(100.0*100.0)
	enddo

        slhood=-slhood-prior
	

end subroutine slikelihood
      
!=======================================================================

end module like
