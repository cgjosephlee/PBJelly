#!/hgsc_software/python/python-2.6/bin/python
import sys, json, math

def getStats(seqLengths):
    data = {}

    seqLengths.sort(reverse=True)
    
    data["numSeqs"] = len(seqLengths)
    data["totalLength"] = sum(seqLengths)
    tl = data["totalLength"]
    n50_mark = data["totalLength"] * .5
    n90_mark = data["totalLength"] * .90
    n95_mark = data["totalLength"] * .95
    
    data["n50"] = None
    data["n50_gt_count"] = None
    data["n90"] = None
    data["n90_gt_count"] = None
    data["n95"] = None
    data["n95_gt_count"] = None
    basesSeen = 0
    
    for pos,n in enumerate(seqLengths):
        basesSeen += n
        if data["n50"] is None and basesSeen > n50_mark:
            data["n50"] = n
            data["n50_gt_count"] = pos
        if data["n90"] is None and basesSeen > n90_mark:
            data["n90"] = n
            data["n90_gt_count"] = pos
        if data["n95"] is None and basesSeen > n95_mark:
            data["n95"] = n
            data["n95_gt_count"] = pos
            break
    #may not have gaps
    if data["numSeqs"] == 0:
        return data
    data["min"] = seqLengths[-1]
    data["FstQu"] = seqLengths[ int(math.floor(data["numSeqs"]*.75)) ]
    median = data["numSeqs"]*.50
    data["median"] = int( (seqLengths[ int(math.floor(median)) ] + \
                           seqLengths[ int(math.floor(median)) ]) / 2)
    data["mean"] = data["totalLength"]/data["numSeqs"]
    data["TrdQu"] = seqLengths[ int(math.floor(data["numSeqs"]*.25)) ] 
    data["max"] = seqLengths[0]

    return data
 
if __name__ == '__main__':
    data = map(float, sys.stdin.read().strip().split('\n'))
    ret = getStats(data)
    
    print "ReadStats: "
    print json.dumps(ret, indent=4)
