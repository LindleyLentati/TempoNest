import numpy as np
import acor

root="102050cm-NOEQ-FitTime-FitProfile-4p-NewDMXRef-J1713+0747-"

samps = open(root+".extract.dat").readlines()[:-1]

for i in range(len(samps)):
	samps[i] = samps[i].strip('\n').split()

samps=np.array(samps).T
Likes=samps[-1]
samps=np.float64(samps[:-1])

TimingParams = open(root+"T2scaling.txt").readlines()
for i in range(len(TimingParams)):
	TimingParams[i] = TimingParams[i].strip('\n').split()


NSamp = len(samps[0])
acors=np.zeros(len(samps))
means=np.zeros(len(samps))
stds=np.zeros(len(samps))

for i in range(len(samps)):
	acors[i] = acor.acor(np.float64(samps[i]))[0]
	means[i]=np.mean(samps[i])
	stds[i] = np.std(samps[i])

for i in range(len(TimingParams)):
	print i, TimingParams[i][0], means[i+1]*np.float64(TimingParams[i][-1])+np.float64(TimingParams[i][-2]), stds[i+1]*np.float64(TimingParams[i][-1])

for p in range(len(samps)):
	index=p
	sum=0
	for i in range(len(samps[index])-1):
    		sum=sum+np.abs(samps[index][i+1]-samps[index][i])
	sum=sum/(NSamp-1)
	print p, sum, np.std(samps[index]), sum/np.std(samps[index])


