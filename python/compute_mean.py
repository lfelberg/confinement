import sys
import numpy as np

def open_file(fname, datloc):
    '''Given filename and data column location, return np array of data'''
    f = open(fname, 'r').readlines()
    data = np.zeros(len(f))
    for lin in range(len(f)):
       #tmp = f[lin].split(",")
        tmp = f[lin].split(",")
        data[lin] = float(tmp[datloc])
    return data

def compute_stats(data):
    '''Given an array of data, print mean and std'''
    print("Mean: {0}, stdev: {1}".format(np.mean(data),np.std(data)))

def main():
    filename=sys.argv[1]
    loc = int(sys.argv[2])
    compute_stats(open_file(filename, loc))

if __name__=="__main__":
    main()

