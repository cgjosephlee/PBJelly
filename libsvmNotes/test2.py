from svmutil import *
"""
I need to figure out the normalize stuff...
"""
import copy

lower = -1.
upper = 1.

def normalize_instances(instances, ranges = None) :
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

y, x = [1,-1], [{1:1, 3:1}, {1:-1,3:-1}]
data = svm_read_problem("mydata.txt")
print data
exit(0)
prob  = svm_problem(y, x)
param = svm_parameter('-c 4 -b 1')
m = svm_train(prob, param)
p_label, p_acc, p_val = svm_predict(y, x, m)
z = [{1:4}, {1:12}, {1:19}]

