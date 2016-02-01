from os import path, listdir,  sep
from collections import OrderedDict
import numpy as np
import sys
import matplotlib.pyplot as plt
    
def get_info(path):
    # Read file
    f = open(path, "r")
    text = f.readlines()
    f.close()

    # Get info
    info = {}
    for s in text:
        if s != "\n":
            key, value = s.split(": ")
            try:
                v = eval(value[:-1])
            except:
                v = value[:-1]
            if isinstance(v, list):
                v = v[0]
            info[key] = v

    return info

def find_SD(data):
    data = np.asarray(data)
    mean = np.mean(data)
    var = np.sum([(x-mean)**2 for x in data]) / (float(len(data))-1)
    SD = np.sqrt(var)
    return mean, SD


def plot_(info):
    fig = plt.figure()
    b= {"bif": "bif", 1: 1, 2: 2}
    for m in ["mean", "max"]:
        for bif in ["bif", 1, 2]:
            normal = False
            anu = False
            for key, value in info.items():
                color = "b" if key.startswith("P") else "r"
                for key, v in value.items():
                    if key == "phi1_bif_half":
                        phi1 = float(v)
                    elif key == "phi2_bif_half":
                        phi2 = float(v)
                    #elif key == "curvature_mean_bif":
                    #    cur = float(v)
                    elif key == "curvature_%s_new_%s" % (m, bif):
                        cur = float(v)
                    #elif key == "sex":
                    #    color = "b" if v == "M" else "r"
                b_ = b[bif]
                if b_ == "bif":
                    phi = phi2+phi1
                    phi_t = r"$\varphi_1 + \varphi_2$"
                else:
                    phi = phi1 if b_ == 1 else phi2
                    phi_t = r"$\varphi_1$" if b_ == 1 else r"$\varphi_2$"
                
                if not normal:
                    plt.plot([phi], [cur], marker="o", color=color,
                            label="Normal")
                    #plt.plot([phi], [cur_new], marker='^', color=color)
                    normal = True
                elif not anu:
                    plt.plot([phi], [cur], marker="o", color=color, label="Aneurysmal")
                    #plt.plot([phi], [cur_new], marker='^', color=color)
                    anu = True
                else:
                    plt.plot([phi], [cur], marker="o", color=color)
                    #plt.plot([phi], [cur_new], marker='^', color=color)
        
            plt.xlabel(phi_t)
            plt.ylabel(r"%s curvature" % m)
            plt.legend()
            plt.hold("on")
    
            plt.title(r"%s curvature in bifurcation %s" % (m, phi_t))
            plt.savefig("plots/phi_%s_%s.png" % (m, bif))
            plt.show()
    
    
def getParameter(data, key):
    val = []
    for value in data.itervalues():
        val.append(value[key])

    return val


if __name__ == "__main__":
    info = OrderedDict()
    folder_path = "/home/aslak/master/src/aneurysms"
    folders = []
    # Anu cases LAT
    #folders += [path.join(folder_path, c) for c in listdir(folder_path) if c.startswith("C0")]

    # Healty cases
    #folders += [path.join(folder_path, c) for c in listdir(folder_path) if
    #        c.startswith("P0")]

    # Anu cases TER ICA/MCA
    #folders += [path.join(folder_path, "new_cases", c) for c in \
    #                listdir(path.join(folder_path, "new_cases")) if c.startswith("C0")]

    # Manipulates area cases
    folders += [path.join(folder_path, c) for c in listdir(folder_path) if c.startswith("A0")]
    folders.sort(key=lambda x: float(x.split("/")[-1][1:]))

    for folder in folders:
        if "C0093" in folder or "C0087" in folder: continue
        info[folder.split("/")[-1]] = get_info(path.join(folder, "manifest.txt"))

    unruptured = []
    ruptured = []
    healty = []

    """
    for key, val in info.iteritems():
        try:
            tmp = val["ruptureStatus"]
            if tmp == "U":
                unruptured.append(val["max_min_ratio_area"])
            else:
                ruptured.append(val["max_min_ratio_area"])
        except KeyError:
            print "Healty", key
            healty.append(val["max_min_ratio_area"])


    color = ["gx", "kv", "c+"]
    label = ["Unruptured", "Ruptured", "Healty"]
    fig, ax = plt.subplots(figsize=(15, 10))
    #ax.set_xtics
    for i, data in enumerate([unruptured, ruptured, healty]):
        mean, SD = find_SD(data)
        first = True
        for p in data:
            if first:
                plt.plot([i], [p], color[i], label=label[i])
                plt.hold("on")
                first = False
            else:
                plt.plot([i], [p], color[i])

        plt.errorbar([i], [mean], yerr=[[-SD], [SD]], fmt="o", color=color[i][0]) 
    plt.xlim([-0.1, 2.1])
    plt.savefig("test_diff.png")

    sys.exit(0)
    """

    R = getParameter(info, "max_min_ratio_area")
    origin = getParameter(info, "origin")
    origin = [o.split(sep)[-3] for o in origin]
    base = "/home/aslak/master/src/aneurysms"
    for i in range(len(R) // 3):
        print "Case", origin[i*3]
        r_original =  get_info(path.join(base, origin[i*3], "manifest.txt"))["max_min_ratio_area"]
        R_tmp = R[i*3:(i+1)*3]
        index_min = R_tmp.index(min(R_tmp))
        index_max = R_tmp.index(max(R_tmp))
        b = range(3)
        b.pop(b.index(index_max))
        b.pop(b.index(index_min))
        index_mean = b[0]
        print R_tmp, r_original
        #print "Index mean:", index_mean, "Index min", index_min, "Index max", index_max
        print "- SD diff", abs(r_original - R_tmp[index_min])
        print "+ SD diff", abs(r_original - R_tmp[index_mean])
        print "+2SD diff", abs(r_original - R_tmp[index_max])
        print ""
    
    C = getParameter(info, "min_circleness")
    E = getParameter(info, "local_max_stepest_disent")
    r_max = getParameter(info, "max_mean_ratio_area")
    r_min = getParameter(info, "min_mean_ratio_area")
    
    beta = []
    for value in info.itervalues():
        folder_path = value["origin"].split("/")
        beta.append(float(folder_path[-1].split("_")[2]))

    status = [0 for i in range(len(beta))]
    for i in range(len(beta) // 3):
        beta_tmp = beta[i*3:(i+1)*3]

        index = beta_tmp.index(min(beta_tmp))
        status[i*3 + index] = 2
        index = beta_tmp.index(max(beta_tmp))
        status[i*3 + index] = 1
    

    """
    status = []
    for value in info.itervalues():
        if value["aneurysmType"] is None:
            status.append(0)
        else:
            stat = 1 if value["ruptureStatus"] == "U" else 2
            status.append(stat)
    """

    print "Total number of cases:", len(R)
    
    num_2sd = sum([s == 1 for s in status])
    num_msd = sum([s == 2 for s in status])
    num_psd = sum([s == 0 for s in status])
    leg = ["R", "Circleness", "A'(s)", "R_max", "R_min"]
    #top_leg = ["Healty (10)", "Unruptured (%s)" % num_unruptured, "Ruptured (%s)" % num_ruptured]
    top_leg = ["  + SD (%s)" % num_2sd, "   +2SD (%s)" % num_msd, "   - SD (%s)" % num_psd]
    
    for i in range(3):
        print "=" *15
        print top_leg[i]
        print "="*15
        for j, data in enumerate([R, C, E, r_max, r_min]):
            tmp = []
            for k, stat in enumerate(status):
                if i == stat:
                    tmp.append(data[k])
    
            mean, SD = find_SD(tmp)
            print "%s %01.04f %01.05f" % (leg[j], mean, SD)

        print ""
