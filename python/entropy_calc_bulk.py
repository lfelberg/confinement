import sys

from csvfile import CSVFile
from volfile import VolFile
from entropy_calc_util import trans_entropy,orien_order_entropy

def main():
    ''' Given a csv file with 5 angles, oxygen pair distance and dist from wall
        calculate translational and rotational entropy
    '''
    angname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    ent_type = sys.argv[5]; nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    angC = CSVFile(angname)
    dis_loc = angC.find_keyword("dis"); vol_loc = angC.find_keyword("vol")

    if ent_type == "trans" or ent_type == "both":
        if "gr" in angname:    ent_t = trans_gr(angC.dat) 
        else:  ent_t = trans_entropy(angC.dat[dis_loc], angC.dat[vol_loc,0])
        print("Translational entropy (cal/mol/K): {0:.7f}".format(ent_t))
    if ent_type == "orien" or ent_type == "both":
        nord = int(sys.argv[6]); other_loc = angC.find_not_keyword("dis")
        oth_key = [angC.key[i] for i in other_loc]
        if nord == 1:
            ent_or1 = orien_order_entropy(1,oth_key,angC.dat[dis_loc],
                                    angC.dat[other_loc],angC.dat[vol_loc,0])
        if nord == 2:
            ent_or2 = orien_order_entropy(2,oth_key,angC.dat[dis_loc],
                                    angC.dat[other_loc],angC.dat[vol_loc,0])
        if nord == 3:
            ent_or3 = orien_order_entropy(3,oth_key,angC.dat[dis_loc],
                                    angC.dat[other_loc],angC.dat[vol_loc,0])
        if nord == 4:
            ent_or4 = orien_order_entropy(4,oth_key,angC.dat[dis_loc],
                                    angC.dat[other_loc],angC.dat[vol_loc,0])
        if nord == 5:
            ent_or5 = orien_order_entropy(4,oth_key,angC.dat[dis_loc],
                                    angC.dat[other_loc],angC.dat[vol_loc,0])

if __name__=="__main__":
    main()
