import sys

forg = sys.argv[-1]
fname = "test.out"

f = open(fname, 'r').readlines()
out_all = 0

for line in f:
    tmp = line.split()
    if int(tmp[2]) == 1: out_all += 1


print("{0:35}  {1:4}".format(forg, int(out_all))) #,int(brk)))
#print("{0:4}  {1:4}".format(int(out_all)) #,int(brk)))
