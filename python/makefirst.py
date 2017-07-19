import sys

lmp = sys.argv[1]
xyz = "first.xyz"

f = open(lmp, "r").readlines()
at_flag = 0; bd_flag = 0; xyzs = []

for line in f:
    if "Bond" in line: break

    if at_flag > 0:
        if at_flag > 1 and bd_flag == 0 and line != "\n":
            tmp = line.split()
            xyzs.append([int(tmp[2]),float(tmp[-3]),float(tmp[-2]),float(tmp[-1])])
        at_flag += 1
    if "Atoms" in line: at_flag = 1

f = open(xyz, "w")
f.write("{0}\nAtoms. Timestep: 0\n".format(len(xyzs)))
for l in range(len(xyzs)):
    f.write("{0} {1:.2f} {2:.2f} {3:.2f}\n".format(*xyzs[l]))

f.close()
