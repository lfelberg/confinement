import sys

from csvfile import CSVFile
from volfile import VolFile
from entropy_calc_util import orien_order_entropy
from entropy_calc_util import STrans

def main():
    ''' Given a csv file with 5 angles, oxygen pair distance and 2D cell vol
        calculate translational and rotational entropy
    '''
    angname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    ent_type = sys.argv[5]; nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    angC = CSVFile(angname)
    dis_loc = angC.find_keyword("dis"); vol_loc = angC.find_keyword("vol")

    if ent_type == "trans" or ent_type == "both":
        if "gr" in angname:
            s_t = STrans(1,1, 0.03, 2); s_t.trans_gr(angC.dat)
            ent_t = s_t.trans_gr(angC.dat)
        else:
            dis_dt = angC.dat[dis_loc]
            s_t=STrans(dis_dt.shape[0],dis_dt.shape[1], 0.03, 2)
            ent_t = s_t.trans_entropy(dis_dt, angC.dat[vol_loc,0]*0.925)
        print("Translational entropy (cal/mol/K): {0:.7f}".format(ent_t))
    if ent_type == "orien" or ent_type == "both":
        nord = int(sys.argv[6]); other_loc = angC.find_not_keyword("dis")
        oth_key = [angC.key[i] for i in other_loc]
        ent_or = orien_order_entropy(nord,oth_key,angC.dat[dis_loc],
                                angC.dat[other_loc],angC.dat[vol_loc,0],2)

if __name__=="__main__":
    main()
