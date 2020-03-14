import sys
import numpy as np

orient_file = sys.argv[1]

tad_dict = {}
comp_score_dict = {}
with open(orient_file) as cf:
    for line in cf:
        line = line.rstrip()
        line = line.split()
        key = line[0] + "\t" + line[1] + "\t" + line[2]+ "\t" + line[3] + "\t" + line[4]
        total = line[6]
        comp_score = line[5]
        comp_score = comp_score.split(",")
        for i in comp_score:
            i = i.split(":")
            comp_score_dict[i[0]] = float(i[1])
        tad_dict[key] = [comp_score_dict,float(total)]
        comp_score_dict = {}

for k in tad_dict.keys():
    vals = []
    if "." in tad_dict[k][0]:
        tad_dict[k].append(tad_dict[k][0]["."]/tad_dict[k][1])
        vals.append(tad_dict[k][0]["."]/tad_dict[k][1])
    else:
        tad_dict[k].append(0)
        vals.append(0)
    if "+" in tad_dict[k][0]:
        tad_dict[k].append(tad_dict[k][0]["+"]/tad_dict[k][1])
        vals.append(tad_dict[k][0]["+"]/tad_dict[k][1])
    else:
        tad_dict[k].append(0)
        vals.append(0)
    if "-" in tad_dict[k][0]:
        tad_dict[k].append(tad_dict[k][0]["-"]/tad_dict[k][1])
        vals.append(tad_dict[k][0]["-"]/tad_dict[k][1])
    else:
        tad_dict[k].append(0)
        vals.append(0)

    vals = np.asarray(vals)
    max_ind = np.where(vals == vals.max())[0]
    if len(max_ind) > 1:
        tad_dict[k].append("multi")
    else:
        tad_dict[k].append(str(int(max_ind[0])))
    #tad_dict[k].append(max_ind)

for k in tad_dict.keys():
    #val = k + "\t" + str(tad_dict[k][2]) + "\t" + str(tad_dict[k][3]) + "\t" + str(tad_dict[k][4]) + "\t" + str(tad_dict[k][5])
    if tad_dict[k][5] == "0":
        orient = "NoCTCF"
    elif tad_dict[k][5] == "1":
        orient = "F"
    elif tad_dict[k][5] == "2":
        orient = "R"
    elif tad_dict[k][5] == "multi":
        orient = "NA"
    val = k + "\t" + orient
    print(val)
