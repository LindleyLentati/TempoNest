module like
	
use params
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
	implicit none
      
	double precision Cube(nest_nPar),slhood, Datavector(21)
	integer i, minE, maxE
	
	slhood = 0d0

        minE = 3
        maxE = 20
	
        Cube(1) = floor((maxE - minE)*Cube(1)+minE)
        Cube(2) = 10d0*Cube(2)

        Datavector = 0
        Datavector(11) = 8

	do i = 1, 20
                if(i .ne. Cube(1))then
		        slhood = slhood + 0.5d0*(Datavector(i))**2
                else if (i .eq. Cube(1))then
                         slhood = slhood + 0.5d0*(Datavector(i) - Cube(2))**2
                end if
	enddo

        slhood=-slhood
	

end subroutine slikelihood
      
!=======================================================================

end module like
