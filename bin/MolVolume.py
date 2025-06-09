#!/home/iaw/soft/conda/2024.06.1/envs/AmberTools23/bin/python
from rich import print as rprint
from rich.progress import track
import math
from typing import List, Set, Dict
import numpy as np
from numpy import ndarray as arr
import parmed
from pathlib import Path
import pathlib
from biopandas.pdb import PandasPdb
import argparse
from multiprocessing import Pool
import sys
from rich.progress import Progress
from tqdm import tqdm

def NumWat(v: float, _p: float) -> int:
    Na: float = 6.022*1e23
    M_mol: float = 18
    v: float = v * (1e-8) *(1e-8) *(1e-8)                           # cm³
    M: float = _p * v
    num: float = M/M_mol * Na
    rprint("The number of Water is {:.2f}".format(num))

    return math.floor(num)

class initAtoms1():
    """
        只需要读取pos 以及原子类型
    """
    def __init__(self, top: str, pdb: str) -> None:
        top_fp: pathlib.PosixPath = Path(top).absolute()
        pdb_fp: pathlib.PosixPath = Path(pdb).absolute()

        self.parm = parmed.amber.LoadParm(top_fp.as_posix())
        self.pdb = PandasPdb().read_pdb(pdb_fp.as_posix())


    def __vdW__(self):
        out: List[float] = [] 
        for atm in self.parm.atoms:
            out.append(atm.rmin)
        return np.array(out)
    
    def __isinBox__(self, box: List[float], pos: arr) -> bool:
        x_min, x_max = min(box[0], box[3]), max(box[0], box[3])
        y_min, y_max = min(box[1], box[4]), max(box[1], box[4])
        z_min, z_max = min(box[2], box[5]), max(box[2], box[5])
        inside: bool = (x_min <= pos[0] <= x_max) and (y_min <= pos[1] <= y_max) and (z_min <= pos[2] <= z_max)
        
        return inside

    def __call__(self, box: List[float]) -> Set:

        pos = self.pdb.df["ATOM"][["x_coord","y_coord","z_coord"]]
        vdw = self.__vdW__()
        pos = np.vstack([pos["x_coord"].to_numpy().ravel(), pos["y_coord"].to_numpy().ravel(), pos["z_coord"].to_numpy().ravel()]).T
      
        # 获取box内的原子

        _del_i: List[int] = []

        for i in range(pos.shape[0]):
            i_pos = pos[i]
            if not self.__isinBox__(box, i_pos):
                _del_i.append(i)
        
        _del_i = list(set(_del_i))

        pos = np.delete(pos, list(set(_del_i)), axis=0)
        vdw = np.delete(vdw, list(set(_del_i)), axis=0)   
        
        return (pos, vdw)

class initAtoms2(initAtoms1):
    def __init__(self, top: str, pdb: str) -> None:
        super(initAtoms2, self).__init__(top, pdb)

    def __call__(self, space: float) -> Set:
        pos: arr[arr] = self.pdb.df["ATOM"][["x_coord","y_coord","z_coord"]]
        vdw: arr = self.__vdW__()
        pos = np.vstack([pos["x_coord"].to_numpy().ravel(), pos["y_coord"].to_numpy().ravel(), pos["z_coord"].to_numpy().ravel()]).T
      
        # 获取这些原子的xmin, ymin, zmin以及xmax, ymax, zmax
        xmin, ymin, zmin = pos.min(axis=0)
        xmax, ymax, zmax = pos.max(axis=0)
        # 需要外扩
        box: List[float] = [  math.floor(xmin - space*2)
                            , math.floor(ymin - space*2)
                            , math.floor(zmin - space*2)
                            , math.ceil(xmax + space*2)
                            , math.ceil(ymax + space*2)
                            , math.ceil(zmax + space*2)]
        
        rprint("The sum of Sphere Volume is {:.2f}".format(np.sum((4/3)*math.pi*np.power(vdw, 3))))
        return (pos, vdw, box)


class Grid():
    """
        box: List[float] -> x, y, z, x1, y1, z1             # Å
        space: float -> 格点长度                               # Å
    """
    def __init__(self, box: List[float], space: float, atoms: Set) -> None:
        
        self.pos_s = atoms[0]
        self.vdw_s = atoms[1]

        self.space: float = space
        
        self._x: float = np.round(abs(box[0] - box[3]), 1)
        self._y: float = np.round(abs(box[1] - box[4]), 1)
        self._z: float = np.round(abs(box[2] - box[5]), 1)
        
        self.Nx: int = int(self._x / space)
        self.Ny: int = int(self._y / space)
        self.Nz: int = int(self._z / space)

        #rprint(self._z, space, self.Nz)
        
        _axis_x: arr = np.linspace(
                                              np.round(box[0], 1)
                                            , np.round(box[3], 1)
                                            , self.Nx + 1)
        _axis_y: arr = np.linspace(
                                              np.round(box[1], 1)
                                            , np.round(box[4], 1)
                                            , self.Ny + 1)
        _axis_z: arr = np.linspace(
                                              np.round(box[2], 1)
                                            , np.round(box[5], 1)
                                            , self.Nz + 1)
        
        _axis_x = [(_axis_x[i] + _axis_x[i+1])/2 for i in range(_axis_x.shape[0]-1)]
        _axis_y = [(_axis_y[i] + _axis_y[i+1])/2 for i in range(_axis_y.shape[0]-1)]
        _axis_z = [(_axis_z[i] + _axis_z[i+1])/2 for i in range(_axis_z.shape[0]-1)]

        #rprint( np.round(box[2], 1), np.round(box[5], 1),self._z , self.Ny + 1, _axis_z)
        #sys.exit(-1)
        X, Y, Z = np.meshgrid(_axis_x, _axis_y, _axis_z, indexing='ij')

        self.points: arr = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

        self.n_points = self.points.shape[0]
        self.GridVolume = self._x * self._y * self._z
        
        rprint("Number of points of grid: {}".format(self.n_points))
        with open("./2.pml", "w+") as F:
            for i in range(self.points.shape[0]):
                F.writelines("cmd.pseudoatom('{}', pos=[{:.4f}, {:.4f}, {:.4f}])\n".format(i, self.points[i][0], self.points[i][1], self.points[i][2]))
        
        #rprint(self.exclude_point(np.array([55, 45, 77]), [19.469, 102.288, 92.732], 1.359))
        #sys.exit(0)


    def __worker__(self, i_pos: int):
        del_point: List[int] = []
        
        # 计算斜对角的一半, 这里有个前提是立方体的边与各个坐标轴是平行关系 
        h = np.sqrt((self.space/2)**3)
        
        v = np.abs(self.points - self.pos_s[i_pos])
        u = v - h
        
        # 归零
        u[u<=0] = 0
        
        # 获取符合要求的索引
        b = np.sum(np.multiply(u, u), axis = 1)
        return np.where(b <= (self.vdw_s[i_pos])**2)[0].tolist()

    def __call__(self, n_cpu) -> Dict:
        
        del_point: List[int] = []
    
        # 加入多线程加速
        pool = Pool(n_cpu)
        rprint("[green]Please wait...")
        result = pool.map(self.__worker__, list(range(self.pos_s.shape[0])))
        pool.close()
        pool.join()
        
        for i in result:
            del_point.extend(i)

        if del_point:
            self.points = np.delete(self.points, list(set(del_point)), axis=0)

        return {"n_reserve_points":len(self.points), "n_points":self.n_points, "Vgrid":self.GridVolume}

def pos2Pymol(fp: str, pos: arr):
     with open(fp, "w+") as F:
        for i in range(pos.shape[0]):
            F.writelines("cmd.pseudoatom('p{}', pos=[{:.4f}, {:.4f}, {:.4f}])\n".format(i, pos[i][0], pos[i][1], pos[i][2]))
        
        # group

def Parm():
    parser = argparse.ArgumentParser(description=
                                     "The author is very lazy and doesn't want to write anything\n"
                                     "Author: IAW [ECNU]"
                                     )
    parser.add_argument("-cpu",type=str, nargs=1, help="The number of cpus")
    parser.add_argument("-top",type=str, nargs=1, help="SmilesFilePath")
    parser.add_argument("-pdb",type=str, nargs=1, help="PqbqtFilePath")
    parser.add_argument("-box",type=str, nargs=1, help="\"x1,y1,z1,x2,y2,z2\"")
    parser.add_argument("-space",type=str, nargs=1, help="haha")
    return parser.parse_args()


def main():
    myP = Parm()
    n_cpu: int = eval(myP.cpu[0])
    fp_top: pathlib.PosixPath = Path(myP.top[0]).absolute()
    fp_pdb: pathlib.PosixPath = Path(myP.pdb[0]).absolute()

    space: float = eval(myP.space[0])
    
    if myP.box[0] != "":
        box: List[float, float, float, float, float, float] = [eval(i) for i in myP.box[0].split(",")]
        atoms_in_box = initAtoms1(fp_top, fp_pdb)(box)
    else:
        atoms_in_box = initAtoms2(fp_top, fp_pdb)(space)
        box: List[float, float, float, float, float, float] = atoms_in_box[2]

    pos2Pymol(fp = "./1.pml", pos = atoms_in_box[0])
    rprint("Number of atoms to be counted: {}".format(len(atoms_in_box[0])))
    result: Dict = Grid(box = box, space=space, atoms = atoms_in_box)(n_cpu)
    rprint(result)
    
    # 计算分子的体积
    rprint("The volume of mol is : {:.2f}".format(( (result["n_points"] - result["n_reserve_points"]) / result["n_points"] ) * result["Vgrid"]))

    # 计算可以填充水分子的数目
    v: float = ( result["n_reserve_points"] / result["n_points"] ) * result["Vgrid"]
    NumWat(v, 1) 

if __name__ == "__main__":
    main()
