import sys

if __name__ == '__main__':
    fh = open(sys.argv[1],'r')
    h = fh.readline().strip()[1:]
    header = {}
    for pos, item in enumerate(h.split('\t')):
        header[item] = pos

    for line in fh.readlines():
        chrom, os, s, ins, ine, e, oe, info = line.strip().split('\t')
        data = {}
        for i in info.split(';'):
            key,val = i.split('=')
            data[key] = val
        #d = filter(lambda x: x.startswith('label'), info.split(';'))[0].split('=')[1]
        d = data['label']
        if os == '.':
            print "\t".join([chrom, ine, oe, d+'_e'])
        elif ine == '.':
            print "\t".join([chrom, os, ins, d+'_s'])
        else:
            print "\t".join([chrom, os, oe, d, '0', '+', ins, ine])
    fh.close()
            
            

