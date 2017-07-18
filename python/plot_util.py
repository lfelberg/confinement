import numpy as np

def get_gr_lst(csv_nm):
    ''' Given a filename for g(r) csv, find sep size and use that to get
        columns to plot'''
    sep_siz = float(csv_nm.split("_")[3]); dim = csv_nm.split("_")[2]; nfct=1.0

    if sep_siz == 6.: 
        nfct = 1.3; cls = [2]
    elif sep_siz >= 8.25 and sep_siz <= 8.75:  cls = [4]
    # flexible
    elif sep_siz == 91.:  
         cls = [2,30]; nfct = 1.4
    elif sep_siz == 92.:  
         cls = [4,28,32]; nfct = 1.5

    # benzene
    elif sep_siz == 9.:
         cls = [3,5,6]; nfct = 1.1

    # rigid
   #elif sep_siz == 9.:  cls = [5,15,18,27]

    elif sep_siz == 10. or sep_siz == 11.: cls = [3,6,]
    elif sep_siz == 12.: 
         cls = [7]; nfct = 1.1
    elif sep_siz == 131 or sep_siz == 141: cls =[5] # inner layer of 3L
    elif sep_siz == 132 or sep_siz == 142: cls =[3,7] # layer close to graph
    elif sep_siz == 143: cls =[3,8] # 14 sep with benzene
    elif sep_siz == 161: cls =[5,7] #inner layer of 3/4
    elif sep_siz == 162: cls =[3,9] #layer near graph
    elif sep_siz == 201: cls =[3,9] #inner layer of 3/4
    elif sep_siz == 202: cls =[5,7] #layer near graph

    else: cls = [2] # dont have this sep size saved
    return cls, nfct
