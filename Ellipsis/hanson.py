import numpy as np
import acor

root="Sim1-2B-Time-"
diag=np.loadtxt(root+'.diag.dat', delimiter=',')
samps = np.float128((np.loadtxt(root+'.extract.dat').T)[:-1])
grads = np.float128((np.loadtxt(root+'.gextract.dat').T)[:-1])




NSamp = len(samps[0])

E_x=diag[1]/NSamp
E_x2=diag[2]/NSamp
E_x3dpsi=diag[3]/NSamp
E_x2dpsi=diag[4]/NSamp
E_xdpsi=diag[5]/NSamp
E_dpsi=diag[6]/NSamp

var2=(E_x3dpsi-3.*E_x*E_x2dpsi+3.*E_x*E_x*E_xdpsi - E_x*E_x*E_x*E_dpsi)/3.
var1=E_x2-E_x*E_x


SampE_x=np.float128(np.zeros(len(samps)))
SampE_x2=np.float128(np.zeros(len(samps)))
SampE_dpsi = np.float128(np.zeros(len(samps)))
SampE_xdpsi = np.float128(np.zeros(len(samps)))
SampE_x2dpsi = np.float128(np.zeros(len(samps)))
SampE_x3dpsi = np.float128(np.zeros(len(samps)))

for i in range(len(samps)):
	SampE_x[i] = np.float128(np.mean(samps[i]))
	SampE_x2[i] = np.float128(np.mean(samps[i]**2))
	
	SampE_dpsi[i] = np.float128(np.mean(grads[i]))
	SampE_xdpsi[i] =  np.float128(np.mean(samps[i]*grads[i]))
	SampE_x2dpsi[i] = np.float128(np.mean((samps[i]**2)*grads[i]))
	SampE_x3dpsi[i] = np.float128(np.mean((samps[i]**3)*grads[i]))


Sampvar1=SampE_x2-SampE_x*SampE_x
Sampvar2=(SampE_x3dpsi-3.*SampE_x*SampE_x2dpsi+3.*SampE_x*SampE_x*SampE_xdpsi - SampE_x*SampE_x*SampE_x*SampE_dpsi)/3.

print "np: ", np.mean(samps[0]), np.std(samps[0])**2
print "Svar: ", SampE_x[0], Sampvar1[0]
print "diag: ", E_x[0], var1[0]

print "vars:", var1[0], var2[0], Sampvar1[0], Sampvar2[0]
H=Sampvar1/Sampvar2

acors=np.zeros(len(samps))
for i in range(len(samps)):
	acors[i] = acor.acor(np.float64(samps[i]))[0]


for p in range(len(samps)):
	index=p
	sum=0
	for i in range(len(samps[index])-1):
    		sum=sum+np.abs(samps[index][i+1]-samps[index][i])
	sum=sum/(NSamp-1)
	print p, sum, np.std(samps[index]), sum/np.std(samps[index])


