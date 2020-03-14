#Arg1 - Path of juicertool
#Arg2 - .hic file 
#Arg3 - TADs bed file
#Arg4 - normalisation (refer dump)
#Arg5 - normalisation method (refer dump)
#Arg6 - resolution

import subprocess
import sys


def construct_command_line(juicer_path, hic_file, coordinate, normalisation, norm_method, res, out_file):
    comnd = ["java","-jar", juicer_path, "dump", normalisation, norm_method,hic_file,coordinate,coordinate,"BP",res,"-d",out_file]
    return(comnd)


juice_path = sys.argv[1]
hic_file = sys.argv[2]
tads = sys.argv[3]
norm = sys.argv[4]
norm_meth = sys.argv[5]
res = sys.argv[6]

with open(tads) as f:
    for l in f.readlines():
        l = l.rstrip()
        l = l.split()
        coordinate = l[0].replace("chr","") + ":" + l[1] + ":" + l[2]
        outfile = "TAD" + "_" + l[0] + "_" + l[1] + "_" + l[2] + "_cliq" + l[3] + ".mat"
        command=construct_command_line(juice_path,hic_file,coordinate,norm,norm_meth,res,outfile)
        subprocess.call(command)
