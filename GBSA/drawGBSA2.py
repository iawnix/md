#!/home/iaw/soft/conda/2024.06.1/envs/pytorch3.9/bin/python
from importlib.resources import path
from turtle import color
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import sys

def drawDEC(path,res):
    data = pd.read_csv(path)
    if eval(res) in data["num"].to_list():
        data =  data.drop(index = data[(data["num"] == eval(res))].index.tolist())

    plt.rcParams["axes.labelweight"] ="bold"
    plt.rcParams["font.family"]="Times New Roman"
    plt.rcParams["font.weight"]="bold"
    plt.rcParams["font.size"]=10

    color = ["#0F3D71","#3B88A6","#94BA9E","#C7C391"]
    data = data.sort_values(by = "num", ascending=True)
    data["vdW + non"] = data["vdW"] + data["nPS"]
    data["ele + pol"] = data["Ele"] + data["PS"]
    plt.figure(figsize=(10,7),dpi=300)

    plt.xticks([i+0.15 for i in range(len(data["num"]))],data["num"])

    ax=plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_position(('data',0))
    ax.spines['left'].set_position(('axes',0))
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.tick_params("x",pad = -15)
    #ax.set_xlabel("Residue [index]")
    ax.set_ylabel("Binding Free Energy [kcal/mol]")
    plt.bar(range(len(data["num"])),data["vdW + non"],width=0.1,label='vdW + non',color=color[0])
    plt.bar([i +0.1 for i in range(len(data["num"]))],data["ele + pol"],width=0.1,label='ele + pol',color=color[1])
    plt.bar([i +0.1*2 for i in range(len(data["num"]))],data["Side"],width=0.1,label='Side',color=color[2])
    plt.bar([i +0.1*3 for i in range(len(data["num"]))],data["Bone"],width=0.1,label='Bone',color=color[3])

    plt.legend(loc = 3)
    plt.savefig("./dec.jpg")
    plt.savefig("./dec.tiff",dpi=300)


def drawTol(path,res):
    data = pd.read_csv(path,names=["name","num","e"])
    if eval(res) in data["num"].to_list():
        data =  data.drop(index = data[(data["num"] == eval(res))].index.tolist())

    data = data.sort_values(by = "num", ascending=True)
    
    plt.figure(figsize=(7,5),dpi=300)

    num1,num2 = data["num"].max(),data["num"].min()
    ni = (num1- num2) // 10

    plt.xticks([i for i in range(num2,num1+ni,ni)],[i for i in range(num2,num1+ni,ni)],color='blue')

    ax=plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.spines['bottom'].set_position(('data',0))
    ax.spines['left'].set_position(('axes',0))

    ax.tick_params("x",pad = -15)

    plt.bar(range(len(data["num"])),data["e"],width=1,label='Total')
    plt.plot(range(len(data["num"])),[-1 for i in range(len(data["num"]))],alpha = 0.3,color="black", linestyle = "--")
    plt.text(-10,-1+0.2,s="-1",ha = "left")
    plt.legend()
    
    plt.savefig("./tol.jpg",dpi=300)
    plt.savefig("./tol.tiff",dpi=300)

def main():
    path = sys.argv[1]
    res = sys.argv[2]
    if "KEY.csv" in path:
        drawDEC(path,res)
    
    elif "Tol.csv" in path:
        drawTol(path,res)

if __name__ == "__main__":
    main()
