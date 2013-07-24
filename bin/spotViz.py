#!/usr/bin/env python
import h5py, sys, numpy, random
from pbsuite.honey.Honey import *
import matplotlib.pyplot as plt


#transposon
#start = 1976517 - buffer
#end   = 1977304 + buffer
#tandem duplication
#start = 1096192 - buffer
#end   = 1096825 + buffer
#tandem deletion
#start = 4294209 - buffer
#end   = 4294314 + buffer    
#start, end = 304417-buffer, 304460+buffer #tails
#start, end = 1456604 - buffer, 1456650 + buffer
#start, end = 3255896 - buffer, 3255932 + buffer

def makeTransformPlots(key, data, start, end, buffer, binsize, fignum, normalize=True):
    orig = numpy.convolve(data, avgWindow, "same") #smooth
    if normalize:
        norm = orig / cov
    else:
        norm = orig
    slopeT = numpy.convolve(norm, slopWindow, "same") #slope transform
    abs = numpy.abs(slopeT) # absolute value
    smo = numpy.convolve(abs, avgWindow, "same") #smooth
    space = numpy.convolve(smo * orig, avgWindow, "same") # Take back to absolute space and smooth again
    win = range(start-buffer, end+buffer)
    plt.figure(fignum)
    
    plt.subplot(711)
    plt.plot(win, data, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(712)
    plt.plot(win, orig, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(713)
    plt.plot(win, norm, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(714)
    plt.plot(win, slopeT, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(715)
    plt.plot(win, abs, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(716)
    plt.plot(win, smo, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(717)
    plt.plot(win, space,  'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.show()
    plt.savefig("transform_%s.png" % key)

def makeLinePlots(data, start, end, buffer, binsize):
    
    avgWindow = numpy.ones(int(binsize))/float(binsize)
    slopWindow = numpy.zeros(binsize)
    slopWindow[:binsize/3] = -1
    slopWindow[-binsize/3:] = 1
    sumWindow = numpy.ones(binsize)
    s = start - buffer
    e = end   + buffer
    cov  = numpy.convolve( data[COV], avgWindow, "same")
    #mat  = signalTransform(data[MAT], cov, slopWindow, avgWindow)
    mis  = signalTransform(data[MIS], cov, slopWindow, avgWindow)
    ins  = signalTransform(data[INS], cov, slopWindow, avgWindow)
    insz = signalTransform(data[INSZ]*ins, cov*binsize, slopWindow, avgWindow)
    dele = signalTransform(data[DEL], cov, slopWindow, avgWindow)
    tai  = signalTransform(data[TAI], cov, slopWindow, avgWindow)
    #maq  = signalTransform(data[MAQ], cov, slopWindow, avgWindow)
    
    win = range(s, e)
    
    #matP  = plt.plot(win, mat, "g-", linewidth=1)
    misP  = plt.plot(win, mis, "r-", linewidth=1)
    insP  = plt.plot(win, ins, "c-", linewidth=1)
    inszP = plt.plot(win, insz, "c-", linewidth=2)
    deleP = plt.plot(win, dele, "b-", linewidth=1)
    tailP = plt.plot(win, tai, "m-", linewidth=1)
    #maqP  = plt.plot(win, maq, "k-", linewidth=1)
    
    ticks = range(s, e, (e-s)/5)[:-1]
    labels = range(s, e, (e-s)/5)[:-1]
    plt.xticks(ticks, labels, horizontalalignment="left", rotation=17)
    plt.xlabel("position")
    plt.ylabel("metric")
    plt.legend([misP, insP, inszP, deleP, tailP], ["MIS", "INS", "INSZ", "DEL", "TAI"])
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    plt.suptitle("%d bp sv (%d - %d)" % (end - start - (buffer*2), start+buffer, end-buffer))
    plt.show()
    plt.savefig("metrics.png")
    

if __name__ == '__main__':
    buffer = 500
    binsize = 50
    avgWindow = numpy.ones(int(binsize))/float(binsize)
    slopWindow = numpy.zeros(binsize)
    slopWindow[:binsize/3] = -1
    slopWindow[-binsize/3:] = 1
    
    h5    = h5py.File(sys.argv[1], 'r')
    cols  = h5.attrs["columns"]
    key   = sys.argv[2]
    start = int(sys.argv[3])
    end   = int(sys.argv[4])
    
    data = h5[key]["data"][: , start-buffer:end+buffer]
    h5.close()
    
    cov = numpy.convolve(data[COV], avgWindow, "same")
    
    makeLinePlots(data, start, end, buffer, binsize)
    for i,k in enumerate(cols):
        norm = not k == 'coverage'
        makeTransformPlots(k, data[i], start, end, buffer, binsize, i, norm)
    
