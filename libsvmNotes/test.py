import sys, h5py
from svmutil import *
from pbsuite.honey.Honey import *

print "reading model"
#y, x = svm_read_problem("../heart_scale")
model = svm_load_model(sys.argv[1])
print "reading data"
file = h5py.File(sys.argv[2])
data = file[file.keys()[0]]["data"]
"""
>>> y, x = svm_read_problem('../heart_scale')
>>> m = svm_train(y[:200], x[:200], '-c 4')
>>> p_label, p_acc, p_val = svm_predict(y[200:], x[200:], m)

# Construct problem in python format
# Dense data
>>> y, x = [1,-1], [[1,0,1], [-1,0,-1]]
I'm going to need a method to create this data...
len(y) = numDataPts
len(x) = numDataPts
len(x[i]) = numObs per dataPt
"""
import copy

lower = -1.
upper = 1.

def normalize_instances(instances, ranges = None) :
    """
    This works for now, but I'll need it to work on numpy.arrays
    --- If I can get it out of the dictionary, I can supply the sparse matrix
        to svm_predict(y, x, m)
    """
    normalized_instances = copy.deepcopy(instances)
    if ranges == None :
        ranges_dict = dict()
    max_attribute = max(normalized_instances[0].keys())
    for attribute in range(1, max_attribute + 1) :  # we iterate on the attributes
        column = []
        for instance in normalized_instances:
            if attribute not in instance.keys():
                instance[attribute] = 0.
            column.append(instance[attribute])
        if ranges != None :
            minimum = ranges[attribute][0]
            maximum = ranges[attribute][1]
        else :
            minimum = min(column)
            maximum = max(column)
            ranges_dict[attribute] = [minimum, maximum]
        for i in range(len(column)) :
            if column[i] == minimum :
                column[i] = lower
            elif column[i] == maximum :
                column[i] = upper
            else :
                column[i] = lower + (upper-lower) * (column[i] - minimum) / (maximum - minimum)
        # Copying normalized values in memory
        for elem, instance in zip(column, normalized_instances):
            instance[attribute] = elem
    if ranges == None :
        return normalized_instances, ranges_dict
    else :
        return normalized_instances

#insertion
start, end =  1096192, 1096820 
#big deletion
start, end = 1976507, 1977314
#whole thing
#start, end = 0, 4500000
start -= 500
end += 500
def honH5ToSvmProblem(data):
    """Turn to the data format that libsvm is expecting"""
    x = []
    j = 0
    for i in data.value.T[start:end]:
        j += 1
        x.append({1:i[MAT]/i[COV], 2:i[TER]/i[COV], 3:i[DIFF]/i[COV]})
    print "parsed %d" % j
    return x, [0]*len(x)

print "hon2svm"
x, y = honH5ToSvmProblem(data)
#This is manually entered.... but I need something to read it in automatically
myscale = {1: [0, 1], 2: [0, 1.611111], 3:[-1, 1.611111]}
#prob = svm_problem(y, x)
normx = normalize_instances(x, myscale)
print normx
sparseNorm = map(lambda x: x.values(), normx)
print "predict"
p_label, p_acc, p_val = svm_predict(y, normx, model)
for pos,i in enumerate(range(len(p_label))):
    print pos+start, p_label[i],   p_val[i]

p_label, p_acc, p_val = svm_predict(y, sparseNorm, model)
#print p_label
#print p_acc
#print p_val
print "label accuracy meansquareerr squaredcorrelationcoeff probestimates"
for pos,i in enumerate(range(len(p_label))):
    print pos+start, p_label[i],   p_val[i]
