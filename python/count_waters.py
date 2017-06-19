import sys

nlines = int(sys.argv[1])
forg = sys.argv[-1]
fname = "wat.lines"

f = open(fname, 'r').readlines()
out, old, crr = -1, -1, 0

for line in f:
    tmp = line.split()
    ct = tmp[0].split(":")
    if out == -1:  # This is the first line
        out = int(ct[0])
        old = out - 1

    crr = int(ct[0])

    if abs(old-crr) > 330: 
        brk = crr - old
        old = crr - 1

    old = crr - 1

lst = crr

brk = float(brk+1)/3.0
out_all = float(out+(nlines-lst))/3.0

#print("{0:35}  {1:4}  {2:4}".format(forg, int(out_all),int(brk)))
print("{0:4}  {1:4}".format(int(out_all),int(brk)))
print("{0:4}  {1:4}".format(int(out_all),int(brk)))
