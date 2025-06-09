#!/home/iaw/soft/conda/2024.06.1/envs/pytorch3.9/bin/python
#################################################################################################################################
##python -In PDBNAMEn.pdb -Gs GaussianVersion -OBJ Lig -Charge NUM -Key 1 OR 0 -Where whereisATOMnanme -mode 1 OR 2 OR 3 OR 4
##IAW  Zhaijihang HENU 
##Linux python3.7
#################################################################################################################################

import subprocess
import os
import sys
import argparse


def Pipe(cmd1 : str, cmd2 : str) -> str:
    ret = subprocess.Popen(cmd1,bufsize=-1,shell=True,encoding="utf-8",stdout = subprocess.PIPE,stderr=subprocess.PIPE)
    out = ret.communicate(input=None)
    out1,error1 = out[0],out[1]
    code1 = ret.returncode
    if error1 != "":
        if not code1:								# 解决UBUNTU上浮点数的Note导致程序跳出
            print("Sucessful: {},But: {}".format(cmd1,error1))
            return out1
        else:
            print("Error: {}".format(error1))
            sys.exit(1)
    else:
        ret_ = subprocess.Popen(cmd2,bufsize=-1,shell=True,encoding="utf-8",stdout = subprocess.PIPE,stdin = subprocess.PIPE,stderr=subprocess.PIPE)
        out_ = ret_.communicate(input=out1)
        out2,error2 = out_[0],out_[1]
        code2 = ret_.returncode
        if error1 != "":
            if not code2:								# 解决UBUNTU上浮点数的Note导致程序跳出
                print("Sucessful: {},But: {}".format(cmd2,error2))
                return out2
            else:
                print("Error: {}".format(error2))
                sys.exit(1)
        else:
            print("Sucessful: {}".format(cmd2))
            return out2
        
def runCMD(cmd:str) -> str:
    ret1 = subprocess.Popen(cmd,bufsize=-1,shell=True,encoding="utf-8",stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    ret1_ = ret1.communicate(input=None)
    out1,error1 = ret1_[0],ret1_[1]
    code = ret1.returncode
    if error1 != "":
        if not code:								# 解决UBUNTU上浮点数的Note导致程序跳出
            print("Sucessful: {},But: {}".format(cmd,error1))
            return out1
        else:
            print("Error: {}".format(error1))
            sys.exit(1)
    else:
        print("Sucessful: {}".format(cmd))
        return out1

def cmdout2list(string : str) -> list:
    var = []
    for i in string.split("\n"):
        if i != "":
            var.append(i)
    return var


def GAUSS(OBJ:str,gauss:str):
    
    if not (os.system("{} {}.com {}.out".format(gauss,OBJ,OBJ))):
        print("Gaussian Finish")
    else:
        print("Error: Please check Gaussian Input")
        sys.exit(1)

def GetOBJfromPDB(file:str,OBJ:str)->None:
    #"grep 'ATOM\|HETATM' ./file > out.pdb"
    
    out = "{}.pdb".format(OBJ)
    grep1 = "grep 'ATOM\|HETATM' ./{}".format(file)
    grep2 = "grep '{}'  > {}".format(OBJ,out)
    Pipe(grep1,grep2)

def PdbToCom(OBJ:str,Charge:int,KEY:str) -> None:
    #"antechamber -i BTB.pdb -fi pdb -o BTB.com -nc Charge -fo gcrt"
    if KEY != "":
        cmd_antechamber = 'antechamber -i {}.pdb -fi pdb -o {}.com -nc {} -fo gcrt -rn {} -gk "{}"'.format(OBJ,OBJ,Charge,OBJ,KEY)
    else:
        cmd_antechamber = 'antechamber -i {}.pdb -fi pdb -o {}.com -nc {} -fo gcrt -rn {}'.format(OBJ,OBJ,Charge,OBJ)
    sedCOM = lambda obj:"sed -i '{}s/{}/%chk={}/g' ./{}.com".format(2,"%chk=molecule",obj,obj)
    runCMD(cmd_antechamber)
    runCMD(sedCOM(OBJ))

def OutToMolFrcmod(OBJ:str) -> None:

    cmd_antechamber = "antechamber -i {}.out -fi gout -o {}.mol2 -fo mol2 -c resp -rn {} -pf y".format(OBJ,OBJ,OBJ)
    cmd_antechamber_AC = "antechamber -i {}.out -fi gout -o {}.ac -fo ac -at amber -c resp -pf y".format(OBJ,OBJ)       # 增加AC文件
    cmd_parmchk2 = "parmchk2 -i {}.mol2 -f mol2 -o {}.frcmod".format(OBJ,OBJ)
    
    runCMD(cmd_antechamber)
    runCMD(cmd_parmchk2)
    runCMD(cmd_antechamber_AC)

def Alter(file2 : str, OBJ : str, whereisATOMname : int):
    file1 = "{}.mol2".format(OBJ)

    sign1 = "@<TRIPOS>ATOM"
    sign2 = "@<TRIPOS>BOND"
    
    mol2name = []
    mol2start = 0
    mol2end = 0

    pdbname = []
    pdbnum= []

    grep = lambda a,b:"grep -n '{}' ./{}".format(a,b)
    cut = "cut -d : -f 1"
    sed = lambda a,b: "sed -n '{}p' {} ".format(a,b)
    Ised = lambda num,old,new,f:"sed -i '{}s/{}/{}/g' ./{}".format(num,old,new,f) 
    awk = lambda a: "awk '{}print ${}{}' ".format("{",a,"}")

    mol2start = cmdout2list(Pipe(grep(sign1,file1),cut))[0]
    mol2end = cmdout2list(Pipe(grep(sign2,file1),cut))[0]
    for i in range(eval(mol2start)+1,eval(mol2end)):
        mol2name.append(cmdout2list(Pipe(sed(i,file1),awk(2)))[0])

    pdbnum = cmdout2list(Pipe(grep(OBJ,file2),cut))
    pdbname = cmdout2list(Pipe(grep(OBJ,file2),awk(whereisATOMname)))

    if len(pdbnum) == len(pdbname) and len(pdbname) == len(mol2name):
        for i in range(len(pdbnum)):
            
            ## 20220407 更正替换字符串位数不一致引起的pdb文件异常
            if len(pdbname[i]) >= len(mol2name[i]):
                while len(mol2name[i]) != len(pdbname[i]):
                    mol2name[i] += " "
            else:
                while len(mol2name[i]) != len(pdbname[i]):
                    pdbname[i] += " "
            runCMD(Ised(pdbnum[i],pdbname[i],mol2name[i],file2))
    else:
        print(len(pdbnum),len(pdbname),len(mol2name))
        print("ERROR pleace check alter")
        sys.exit(1)

def Parm():

    parser = argparse.ArgumentParser(description='Calculating the charge of small molecules in PDB')
    parser.add_argument('-In', type=str, nargs=1,help='PDB')
    #20220506 新增切换不同的高斯版本
    parser.add_argument('-Gs',type=str, nargs=1, help='Gaussian version')
    parser.add_argument('-OBJ', type=str, nargs=1,help='反应物')
    parser.add_argument('-Charge', type=int, nargs=1,help='反应物的小分子')
    parser.add_argument("-Key",type=int, nargs=1,help='是否使用新的计算resp方法')
    parser.add_argument("-Where",type=int, nargs=1,help='pdb中原子名称所在的列')
    parser.add_argument('-mode', type=int, nargs=1,help='模式选择: 1为准备计算GAUSS; 2为进行高斯计算 3为高斯结果转化 4为修改PDB')
    
    return parser.parse_args()

def main():
    parm = Parm()
    KEY = " opt b3lyp/6-31g* scrf=(smd,solvent=water) pop=mk iop(6/33=2,6/42=6)"
    if parm.mode[0] == 1:
        # creat OBJ.pdb
        GetOBJfromPDB(parm.In[0],parm.OBJ[0])
        # OBJ.pdb To OBJ.com
        if parm.Key[0] == 1:
            PdbToCom(parm.OBJ[0],parm.Charge[0],KEY)
        elif parm.Key[0] == 0:
            PdbToCom(parm.OBJ[0],parm.Charge[0],"")
        else:
            print("ERROR: Key argument, only 1 or 2")
            sys.exit(1)
    
    elif  parm.mode[0] == 2:
        # run GAUSS
        GAUSS(parm.OBJ[0],parm.Gs[0])
    elif  parm.mode[0] == 3:
        # creat OBJ.mol2 OBJ.fremod OBJ.ac
        OutToMolFrcmod(parm.OBJ[0])
    elif  parm.mode[0] == 4:
        # alter in.PDB 
        Alter(parm.In[0],parm.OBJ[0],parm.Where[0])
    else:
        print("ERROR: please check mode")
        sys.exit(1)

if __name__ == "__main__":
    main()
