

def xyz_reader(filename, times = [], dims = []):
    ''' Method to create array of atom coordinates in xyz file'''
    time, atoms, atom, types = [], [], [], []
    grab_snap, t_ct = 0, 0
    print(filename, len(times))
    f=open(filename, "r")
    f.readline()
    line = f.readline()
    if times == [] or times[0][0] == int(line.split()[-1]):
        grab_snap = 1
        time.append(int(line.split()[-1]))
                                                                                   
    '''First, we read in the coordinates from the trajectory file and
       save them into atom=[[coords],[coords]...,[coords]]'''
    for line in f:
        if times != [] and t_ct >= len(times): break

        #Gets past # and atoms lines
        if line[0]=="A":
            if grab_snap == 1:
                atom+=[atoms] #; typ += [types]
                atoms, types = [], []

            if times == [] or times[t_ct][0] == int(line.split()[-1]):
                t_ct += 1
                print(line.split()[-1])
                time.append(int(line.split()[-1]))
            else:
                grab_snap = 0
        # Only grab snapshots that have volume data
        elif len(line.split()) > 1 and grab_snap == 1:
            #crds is a list of the coordinates in string formats
            crds=line.split()
            if t_ct == 0: types.append(int(crds[0]))
            if dims == []:
                atoms.append([float(crds[1]),float(crds[2]),float(crds[3])])
            else: # Move coords so that box is from 0 -> xmin
                atoms.append([float(crds[1])-dims[t_ct][0],
                              float(crds[2])-dims[t_ct][2],
                              float(crds[3])-dims[t_ct][4]])
    #For last timestep, since it doesn't have a line about Atoms after it.
    if grab_snap == 1:
        atom=atom+[atoms] #; typ = typ+[types]
    print(len(time), len(times))
    f.close()
    return time, atom, types
