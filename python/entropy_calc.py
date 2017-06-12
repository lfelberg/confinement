import sys

from csvfile import CSVFile
from volfile import VolFile
from entropy_calc_util import STrans, SOrien

def main():
    ''' Given a csv file with 5 angles, oxygen pair distance and 2D cell vol
        calculate translational and rotational entropy
    '''
    angname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    dim = int(sys.argv[5]); ent_type = sys.argv[6]; 
    angC = CSVFile(angname); dr = []
    dis_loc = angC.find_keyword("dis"); vol_loc = angC.find_keyword("vol")

    if ent_type == "trans" or ent_type == "both":
        if "gr" in angname:
            for i in range(7,len(sys.argv)): dr.append(sys.argv[i])
            s_t = STrans(1, 1, 0.03, dim); ent_t = s_t.trans_gr(angC.dat, dr)
        else:
            dis_dt = angC.dat[dis_loc]
            s_t=STrans(dis_dt.shape[0],dis_dt.shape[1], 0.03, dim)
           #ent_t = s_t.trans_entropy(dis_dt, angC.dat[vol_loc,0]*0.925)
            ent_t = s_t.trans_entropy(dis_dt, angC.dat[vol_loc,0])
        print("Translational entropy (cal/mol/K): {0:.7f}".format(ent_t))
    if ent_type == "orien" or ent_type == "both":
        nord = int(sys.argv[7]); other_loc = angC.find_not_keyword("dis")
        oth_key = [angC.key[i] for i in other_loc]

        if "trans_gr" in angname:
            for i in range(7,len(sys.argv)): dr.append(sys.argv[i])
            sorC = SOrien(1, 1, nord, 0.10, dim)
            e_or = sorC.orien_gang(dr, angC.dat, dr)
        else:
            dis_dt = angC.dat[dis_loc]
            sorC=SOrien(dis_dt.shape[0],dis_dt.shape[1], nord, 0.10, dim)
            e_or=sorC.orien_order_entropy(dis_dt,angC.dat[other_loc],
                                          angC.dat[vol_loc,0])

if __name__=="__main__":
    main()

