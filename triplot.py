#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, LinearLocator, NullFormatter, NullLocator
import matplotlib.ticker
from optparse import OptionParser

"""
Given a 2D matrix of (marginalised) likelihood levels, this function returns
the 1, 2, 3- sigma levels. The 2D matrix is usually either a 2D histogram or a
likelihood scan

"""
def getsigmalevels(hist2d):
  # We will draw contours with these levels
  sigma1 = 0.68268949
  level1 = 0
  sigma2 = 0.95449974
  level2 = 0
  sigma3 = 0.99730024
  level3 = 0

  #
  lik = hist2d.reshape(hist2d.size)
  sortlik = np.sort(lik)

  # Figure out the 1sigma level
  dTotal = np.sum(sortlik)
  nIndex = sortlik.size
  dSum = 0
  while (dSum < dTotal * sigma1):
    nIndex -= 1
    dSum += sortlik[nIndex]
  level1 = sortlik[nIndex]

  # 2 sigma level
  nIndex = sortlik.size
  dSum = 0
  while (dSum < dTotal * sigma2):
    nIndex -= 1
    dSum += sortlik[nIndex]
  level2 = sortlik[nIndex]

  # 3 sigma level
  nIndex = sortlik.size
  dSum = 0
  while (dSum < dTotal * sigma3):
    nIndex -= 1
    dSum += sortlik[nIndex]
  level3 = sortlik[nIndex]

  return level1, level2, level3



def makesubplot2d(ax, samples1, samples2, weights=None):
    xmin = np.min(samples1)
    xmax = np.max(samples1)
    ymin = np.min(samples2)
    ymax = np.max(samples2)

    hist2d,xedges,yedges = np.histogram2d(samples1, samples2, weights=weights, \
            bins=40,range=[[xmin,xmax],[ymin,ymax]])
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
    
    xedges = np.delete(xedges, -1) + 0.5*(xedges[1] - xedges[0])
    yedges = np.delete(yedges, -1) + 0.5*(yedges[1] - yedges[0])
    
    level1, level2, level3 = getsigmalevels(hist2d)
    
    contourlevels = (level1, level2, level3)
    
    #contourcolors = ('darkblue', 'darkblue', 'darkblue')
    contourcolors = ('black', 'black', 'black')
    contourlinestyles = ('-', '--', ':')
    contourlinewidths = (2.0, 2.0, 2.0)
    contourlabels = [r'1 $\sigma$', r'2 $\sigma$',r'3 $\sigma$']
    
    line1 = plt.Line2D(range(10), range(10), linewidth=contourlinewidths[0], \
            linestyle=contourlinestyles[0], color=contourcolors[0])
    line2 = plt.Line2D(range(10), range(10), linewidth=contourlinewidths[1], \
            linestyle=contourlinestyles[1], color=contourcolors[1])
    line3 = plt.Line2D(range(10), range(10), linewidth=contourlinewidths[2], \
            linestyle=contourlinestyles[2], color=contourcolors[2])
    
    contall = (line1, line2, line3)
    contlabels = (contourlabels[0], contourlabels[1], contourlabels[2])

    c1 = ax.contour(xedges,yedges,hist2d.T,contourlevels, \
            colors=contourcolors, linestyles=contourlinestyles, \
            linewidths=contourlinewidths, zorder=2)
    
    
def makesubplot1d(ax, samples, weights=None):
    ax.hist(samples, 100, color='k', histtype='bar', linewidth=2.0)


# The mcmc chain (ASCII file, with columns the values of the parameters each step)
# Note that the first two columns are walker index, and loglikelihood value
# So the indices are used '+2' in the chain below
# Layout ASCII-file
#     Col 1        Col 2        Col 3           Col 4
#   walker id  loglikelihood  parameter 1   parameter 2
#   walker id  loglikelihood  parameter 1   parameter 2
#  .... etc.

parser = OptionParser()
parser.add_option("-f", "--infile",dest="root",metavar='INFILE')
(options,args)=parser.parse_args()


shortname=options.root
chainfilename = shortname+'-post_equal_weights.dat'
figurefilename = shortname+'-triplot.png'
parfilename=shortname+'-.paramnames'
chain = np.loadtxt(chainfilename)

parfile = open(parfilename)
lines=[line.strip() for line in parfile]
parlabels=[]
for i in range(len(lines)):
	lines[i]=lines[i].split(" ")
	parlabels.append(lines[i][1])
	


parplotlabels = []
parplotnums= []

for i in range(len(chain[0,:])-1):
	mean=np.sum(chain[:,i])
	stdev=np.std(chain[:,i])
	if stdev != 0 :
		parplotlabels.append(parlabels[i])
		parplotnums.append(i)
# The labels, and the indices, of the parameters

parameters = np.array(parplotnums)
parlabels=parplotlabels
# Create the plot array
f, axarr = plt.subplots(nrows=len(parameters), ncols=len(parameters))

for i in range(len(parameters)):
    # for j in len(parameters[np.where(i <= parameters)]:
    for j in range(len(parameters)):
        ii = i
        jj = len(parameters) - j - 1

        xmajorLocator = matplotlib.ticker.MaxNLocator(nbins=4,prune='both')#LinearLocator(3)
        ymajorLocator = matplotlib.ticker.MaxNLocator(nbins=4,prune='both')#LinearLocator(3)

        if j <= len(parameters)-i-1:
            axarr[jj][ii].xaxis.set_minor_locator(NullLocator())
            axarr[jj][ii].yaxis.set_minor_locator(NullLocator())
            axarr[jj][ii].xaxis.set_major_locator(NullLocator())
            axarr[jj][ii].yaxis.set_major_locator(NullLocator())

            axarr[jj][ii].xaxis.set_minor_formatter(NullFormatter())
            axarr[jj][ii].yaxis.set_minor_formatter(NullFormatter())
            axarr[jj][ii].xaxis.set_major_formatter(NullFormatter())
            axarr[jj][ii].yaxis.set_major_formatter(NullFormatter())
            xmajorFormatter = FormatStrFormatter('%g')
            ymajorFormatter = FormatStrFormatter('%g')

            if ii == jj:
                # Make a 1D plot
                makesubplot1d(axarr[ii][ii], chain[:,parameters[ii]])
	    else:
                # Make a 2D plot
                makesubplot2d(axarr[jj][ii], chain[:,parameters[ii]], \
                        chain[:,parameters[jj]])

            axarr[jj][ii].xaxis.set_major_locator(xmajorLocator)
            axarr[jj][ii].yaxis.set_major_locator(ymajorLocator)
        else:
            axarr[jj][ii].set_visible(False)
            #axarr[jj][ii].axis('off')

        if jj == len(parameters)-1:
            axarr[jj][ii].xaxis.set_major_formatter(xmajorFormatter)
            axarr[jj][ii].set_xlabel(parlabels[ii])

        if ii == 0:
            if jj == 0:
                axarr[jj][ii].yaxis.set_major_locator(NullLocator())
                axarr[jj][ii].set_ylabel('Post.')
            else:
                axarr[jj][ii].yaxis.set_major_formatter(ymajorFormatter)
                axarr[jj][ii].set_ylabel(parlabels[jj])


#f.subplots_adjust(hspace=0)
#plt.setp([a.get_xticklabels() for a in f.axes[:-0-2]], visible=False)
#plt.tight_layout() # Or equivalently,  "plt.tight_layout()"

#plt.savefig('pulsar-' + str(psr) + '.png')
#plt.savefig(figurefilename)
plt.show()

"""
# Fine-tune figure: make subplots close to each other and hide x ticks for all
# but the bottom plot
# Also add some space for the legend below the plots
f.subplots_adjust(bottom=0.22)
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-0-2]], visible=False)
"""
