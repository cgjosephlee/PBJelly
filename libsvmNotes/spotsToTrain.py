import h5py, sys, random

from pbsuite.honey.Honey import *

#first iter, I'm only doing errors, numbasins, and matches - err and mat is as a function of coverage
#second iter, I need to average the data - this will give me fewer spots, but will be more accurate?

def create(fg, spot, data, sample = None):
    """
    creates a +1 training spot
    """
    if sample is not None:
        i = 0
        while i < sample:
            d = random.randint(spot[0], spot[1])
            print fg+" %s" % getData(data[:,d])
            i += 1
    else:
        for i in xrange(spot[0],spot[1]):
            print fg+" %s" % getData(data[:,i])
    

def getData(point):
    """
    make the data or the spot
    #only doing three right now
    todo - all individually - super classification of the three types
    (4 if you include concordant)
    """
    return ("1:%f " \
           "2:%f "  \
           "3:%f" )\
             % (point[MAT]/point[COV], point[TER]/point[COV], point[DIFF]/point[COV] )
        
if __name__ == '__main__':
    spots = h5py.File(sys.argv[1])
    data = spots[spots.keys()[0]]["data"]
    
    null1 = 0, 1096408
    spot1 = 1096408, 1096820
    null2 = 1096821, 1976507
    spot2 = 1976507, 1977314
    null3 = 1977315, 4294292
    spot3 = 4294292, 4294443
    null4 = 4294444, len(data[0])
    
    create("+1", spot1, data)
    create("+1", spot2, data)
    create("+1", spot3, data)
    create("-1", null1, data, 2500)
    create("-1", null2, data, 2500)
    create("-1", null3, data, 2500)
    create("-1", null4, data, 2500)
    
