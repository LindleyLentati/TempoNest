import numpy as np
import matplotlib.pyplot as plt
coeffs=(np.loadtxt("Coeffs.dat").T)[2:]
Ncoeff=50

RefCoeffs=np.zeros(Ncoeff)
RFreq=138
refloc=2
Beta= 0.0360142951246998633
polyorder = 2

doplot = False


for i in range(Ncoeff):
	n=i
	coeffs[:][2][n::Ncoeff][refloc] = 10.0**-10
	line=np.zeros(len(coeffs[:][0][n::Ncoeff]))
	c=np.polynomial.polynomial.polyfit((coeffs[:][0][n::Ncoeff]-RFreq)/1000, coeffs[:][1][n::Ncoeff], polyorder, rcond=None, full=False, w=1.0/coeffs[:][2][n::Ncoeff])
	print i, c

	if doplot:
		for i in range(len(c)):
			line+=c[i]*((coeffs[:][0][n::Ncoeff]-RFreq)/1000)**i
		plt.errorbar((coeffs[:][0][n::Ncoeff]-RFreq)/1000, coeffs[:][1][n::Ncoeff], yerr=coeffs[:][2][n::Ncoeff],linestyle='')
		plt.errorbar((coeffs[:][0][n::Ncoeff]-RFreq)/1000, line)
		plt.show()

f=open("ProfileInfo.Initial.40C.dat", "w")
f.write(str(Ncoeff)+" "+str(Ncoeff)+"\n")
for i in range(Ncoeff):
	RefCoeffs[i] = coeffs[:][1][i::Ncoeff][refloc]
	print RefCoeffs[i]
	f.write(str(RefCoeffs[i])+"\n")

print Beta
f.write(str(Beta)+"\n")
for p in range(polyorder):
	for i in range(Ncoeff):
		n=i
		c=np.polynomial.polynomial.polyfit((coeffs[:][0][n::Ncoeff]-RFreq)/1000, coeffs[:][1][n::Ncoeff], polyorder, rcond=None, full=False, w=1.0/coeffs[:][2][n::Ncoeff])
		if(i==0):
			print 0
			f.write(str(0)+"\n")
		else:
			print c[p+1]
			f.write(str(c[p+1])+"\n")

f.close()		
