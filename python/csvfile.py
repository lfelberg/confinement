import sys
import numpy as np
import csv

class CSVFile:
    '''A class for csv files'''
    csvfname = ''

    def __init__(self, fname):
         self.csvfname = fname
         self.get_dic(fname)

    def get_dic(self, fname):
        '''Method to get dictionary from csv file'''
        with open(fname) as csvfile:                                              
            reader, rw_ct, rw, dt = csv.DictReader(csvfile), 0, 0, []              
            rw = 0                                                                 
            for row in reader:                                                     
                dat, rk = [], row.keys()                                           
                if rw == 0: 
                    self.key = sorted(rk)
                for key in sorted(rk):                                             
                    if "_" in row[key]:
                        tmp = row[key].split('_')
                        for i in range(len(tmp)): dat.append(float(tmp[i]))
                    elif key != "":
                        dat.append(float(row[key]))                                    
                dt.append(dat)                                                     
                rw += 1 
        self.dat = np.array((dt)).T
        print("Keys: ", self.key, " and dat shape ", self.dat.shape)

    def find_keyword(self, word):
        '''Find locations of key dim, where given word is found'''
        locs = [ x for x in range(len(self.key)) if word in self.key[x]]
        return locs

    def find_not_keyword(self, word):
        '''Find locations of key dim, where given word is not found'''
        locs = [ x for x in range(len(self.key)) if word not in self.key[x]]
        return locs
