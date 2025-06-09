#!/home/iaw/soft/conda/2024.06.1/envs/pytorch3.9/bin/python
##################################################################################
##python runRMS.py A.csv
##T
##Zhai Jihang               HENU
##Linux python3.7
##################################################################################

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys

def draw(df  ,c_lis,center_lis) ->  None:

    # 设置字体格式
    plt.rcParams["axes.labelweight"] ="bold"
    plt.rcParams["font.family"]="Times New Roman"
    plt.rcParams["font.weight"]="bold"
    plt.rcParams["font.size"]=10

    # 创建你画布
    fig,ax = plt.subplots()
    # 设置边框
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    # 设置标题
    ax.set_title("Cluster By RMSD")
    x,y = df["fra"],df["rms"]
    y_max = max(y)
    y_min = min(y)
    x_max = max(x)
    x_min = min(x)

    colors = ["#2878B5","#F8AC8C","#C82423","#FF8884","#9AC9D8","#54B345","#32B897","#05B9E2","#8983BF","#C76DA2"]
    markers = ["o","v","^","s","*","D","p","<",">","P","X"]


    plt.scatter(x,y,c=[colors[i] for i in c_lis], linestyle='-',marker = 'o',s=2.5,alpha=0.4,linewidth=0.5)
    for i in center_lis:
        plt.scatter(x[i-1],y[i-1],c=colors[c_lis[i-1]], linestyle='-',marker = '*',s=80,label = "Center[{}]-Frames[{}]".format(center_lis.index(i)+1,i),alpha=1,linewidth=1,edgecolors="black")
        print("Drawing: Center[{}]-Frames[{}]".format(center_lis.index(i)+1,i))
    ax.set_xlabel("FRAME [Index]")
    ax.set_ylabel("RMSD [Å]")
    ax.set_ylim([0,y_max*2])
    ax.set_xlim([x_min-2,x_max+2])
    plt.legend()
    plt.savefig("./cluster.jpg",dpi=300)
    plt.savefig("./cluster.tiff",dpi=300)

def readcolordat(p):
    count = 1
    out = []
    with open(p,mode="r") as F:

        for line in F.readlines():
            if count:
                count = 0
            else:
                out.append(eval(line.rstrip("\n").split()[1]))
    return out
def main():

    p = sys.argv[1]
    p2 = sys.argv[2]
    center = [eval(i) for i in sys.argv[3].split(",")]
 

    
    df=pd.read_csv(p,header=None,names=["fra","rms"],sep=",")

#  #Frame     Cnum_00002
#      1            1
#      2            1
#      3            1
#      4            1
#      5            1
#      6            1
    color_lis = readcolordat(p2)

    draw(df,color_lis,center)

if __name__ == "__main__":
    main()


