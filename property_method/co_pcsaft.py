# -*- coding: utf-8 -*-
"""
共聚物的状态方程, 待完善
Copolymer PC-SAFT EOS
These functions implement Copolymer PC-SAFT equation of state. 
"""


import math
import numpy as np
import pcsaft


kb = 1.380648465952442093e-23  # Boltzmann constant, J K^-1
PI = 3.141592653589793
N_AV = 6.022140857e23          #  Avagadro's number
E_CHRG = 1.6021766208e-19      #  elementary charge, units of coulomb
perm_vac = 8.854187817e-22     # permittivity in vacuum, C V^-1 Angstrom^-1
R = kb*N_AV                    # Gas constant
π = np.pi

# Table 1. Universal Model Constants for Equations 18 and 19
a0 = np.array([0.910563145, 0.636128145, 2.686134789, -26.54736249, 97.75920878, -159.5915409, 91.29777408   ])
a1 = np.array([-0.308401692, 0.186053116, -2.503004726, 21.41979363, -65.25588533, 83.31868048, -33.74692293 ])
a2 = np.array([-0.090614835, 0.452784281, 0.596270073, -1.724182913, -4.130211253, 13.77663187, -8.672847037 ])
b0 = np.array([0.724094694, 2.238279186, -4.002584949, -21.00357682, 26.85564136, 206.5513384, -355.6023561  ])
b1 = np.array([-0.575549808, 0.699509552, 3.892567339, -17.21547165, 192.6722645, -161.8264617, -165.2076935 ])
b2 = np.array([0.097688312, -0.255757498, -9.155856153, 20.64207597, -38.80443005, 93.62677408, -29.66690559 ])


# 自定义异常类
class InputError(Exception):
    # Exception raised for errors in the input.
    def __init__(self, message):
        self.message = message

# 自定义异常类
class SolutionError(Exception):
    # Exception raised when a solver does not return a value.
    def __init__(self, message):
        self.message = message

# 检查输入是否合理
def  check_input(x, name1, var1, name2, var2):
    ''' Perform a few basic checks to make sure the input is reasonable. '''
    if abs(np.sum(x) - 1) > 1e-7:
        raise InputError('The mole fractions do not sum to 1. x = {}'.format(x))
    if var1 <= 0:
        raise InputError('The {} must be a positive number. {} = {}'.format(name1, name1, var1))
    if var2 <= 0:
        raise InputError('The {} must be a positive number. {} = {}'.format(name2, name2, var2))


# 二元交互参数(公式来自aspen文档"PC-SAFT Binary Parameters" )
def pcsaft_kij(t, a=0, b=0, c=0, d=0, e=0, t_ref=298.15):
    t_r = t/t_ref
    k_ij = a + b/t_r + c*math.log(t_r) + d*t_r + e*t_r**2
    return


def pcsaft_ares(t, rho, x, mole_args, seg_args=None, k_ij=None):
    """
    Calculate the residual Helmholtz energy.

    Parameters
    ----------
    t : float
        Temperature (K)
    rho : float
        Molar density (mol/m^3)
    x : ndarray, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    mole_args : dict
        A dictionary containing PC-SAFT parameters of conventional molecule and polymer that can be passed for
        use in PC-SAFT:

        components : ndarray, shape (n,)
            Name of each component (excluding segments)
        type : ndarray, shape (n,)
            type of each component ("conventional" or "polymer", excluding segments)
        m : ndarray, shape (n,)
            Segment number for each component.
            The parameter of polymer need to be calculated according to the segment composition.
        s : ndarray, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
            The parameter of polymer need to be calculated according to the segment composition.
        e : ndarray, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
            The parameter of polymer need to be calculated according to the segment composition.
    seg_args : dict
        A dictionary containing PC-SAFT parameters of segments that can be passed for use in PC-SAFT:

        segments : ndarray, shape (ns,)
            Name of segments. It has a length of ns, where ns is the number of segments in the system.
        type : ndarray, shape (ns,)
            type of each component ("repeat" or "end" or "branch3" or "branch4")
        r : ndarray, shape (ns,)
            Segment ratio parameter for each component.
        s : ndarray, shape (ns,)
            Segment diameter for each segment. Units of Angstrom.
        e : ndarray, shape (ns,)
            Dispersion energy of each segment. Units of K.  
        mw : ndarray, shape (ns,)
            Relative molecular mass of each segment
    



    Returns
    -------
    ares : float
        Residual Helmholtz energy (J/mol)
    """

    components = mole_args["components"]    # 组分名称列表(不包含链段)
    type = mole_args["type"]                # 组分类型列表(不包含链段)
    nc = len(components)    # 组分数(不包含链段)
    # 如果没有链段, 则按照pcsaft计算
    if seg_args == None:
        if k_ij != None:
            
        return pcsaft.pcsaft_ares(t, rho, x, mole_args)

    segments = seg_args["segments"]         # 链段名称列表
    
    ns = len(segments)    # 链段数

    r = np.zeros((nc, ns))     # 链段比率参数r_iα
    σ = np.zeros((nc, ns))     # 链段直径σ_iα
    ε = np.zeros((nc, ns))     # 链段能量参数ε_iα
    SFRAC = np.zeros((nc, ns)) # 在共聚物i上链段α的摩尔分数
    DPN = np.zeros(nc)         # 共聚物的数均聚合度
    mw = np.zeros((nc, ns))    # 链段的相对分子质量
    M = np.zeros((nc, ns))     # 链段的数均分子量M_iα = DPN[i]*SFRAC[i][α]*mw[i][α]
    m_iα = np.zeros((nc, ns))  # 链段数m_iα = r_iα*M_iα
    m = np.zeros(nc)           # 分子的链段数m_i = sum(m_iα)
    z = np.zeros((nc, ns))     # 链段分数z_iα = m_iα / m_i

    # 给分子和共聚物的基本参数m, σ, ε赋值
    for i in range(nc):
        if type[i] == "conventional":
            m_iα[i][0] = mole_args["m"][i]
            σ[i][0] = mole_args["s"][i]
            ε[i][0] = mole_args["e"][i]
        elif type[i] == "polymer":
            DPN[i] = polymer[components[i]]["DPN"]       # 共聚物i的数均聚合度
            x[i] = x[i]/DPN[i]                           # 将链段的摩尔分数换算为聚合物的摩尔分数
            SFRAC[i] = polymer[components[i]]["SFRAC"]   # 共聚物i的链段组成(一维数组)
            mw[i] = seg_args["mw"]                       # 共聚物i的各链段分子质量(一维数组)
            r[i] = seg_args["r"]
            σ[i] = seg_args["s"]
            ε[i] = seg_args["e"]
            # for α in range(ns):
            #     M[i][α] = DPN[i]*SFRAC[i][α]*mw[i][α]    # 链段的数均分子量M_iα
            #     m_iα[i][α] = M[i][α]*r[i][α]
            M[i] = DPN[i]*SFRAC[i]*mw[i]    # 链段的数均分子量M_iα = DPN[i]*SFRAC[i][α]*mw[i][α]
            m_iα[i] = M[i]*r[i]             # 链段数m_iα = r_iα*M_iα
        else:
            raise InputError("组分type输入错误, 只能为conventional和polymer")
            
        m[i] = np.sum(m_iα[i])     # m_i = sum(m_iα)
        z[i] = m_iα[i]/m[i]        # z_iα = m_iα / m_i

    x = x/sum(x)          # 归一化
    m_avg = np.sum(x*m)   # 平均链段数 
    d = σ*(1-0.12*np.exp(-3*ε/t))  # 温度相关的链段直径
    ζ = np.zeros(4)        # 核心的中间变量, eq A.10
    ρ = rho*N_AV/1.0e30    # 分子数密度, 单位1/A^(3)

    zd = np.zeros(nc)
    for n in range(4):
        for i in range(nc):
            zd[i] = np.sum(z[i]*d[i]**n)    # eq A.10
        ζ[n] = π/6*ρ*np.sum(x*m*zd)         # eq A.10
    ares_hs = 1/ζ[0]*(3*ζ[1]*ζ[2]/(1-ζ[3]) + ζ[2]**3/(ζ[3]*(1-ζ[3])**2) \
                + (ζ[2]**3/ζ[3]**2 - ζ[0])*np.log(1-ζ[3]))

    
    ghs = np.zeros((nc, ns, nc, ns))  # 径向分布函数ghs
    B = np.zeros((nc, ns, nc, ns))    # 键分数B
    ares_chain = 0       # 成键项剩余自由能

    # 参考Aspen文档"Copolymer PC-SAFT Dispersion Term"
    k_iαjβ = np.zeros((nc, ns, nc, ns))   # 链/链段的二元交互参数k_iαjβ
    σ_iαjβ = np.zeros((nc, ns, nc, ns))   # 交叉链段直径σ_iαjβ
    ε_iαjβ = np.zeros((nc, ns, nc, ns))   # 交叉链段能量参数ε_iαjβ

    m2εσ3 = 0    # 计算色散项的中间变量，Aspen用X表示该变量
    m2ε2σ3 = 0   # 计算色散项的中间变量，Aspen用Y表示该变量

    for i in range(nc):
        for α in range(ns):
            for j in range(nc):
                for β in range(ns):
                    # 计算径向分布函数ghs
                    if d[i][α] == 0 or d[j][β] == 0:
                        ghs[i][α][j][β] = 0     # 实际上没有这个索引，设为0
                    else: 
                        ghs[i][α][j][β] = 1/(1-ζ[3]) + (d[i][α]*d[j][β]/(d[i][α]+d[j][β]))*3*ζ[2]/(1-ζ[3])**2 + \
                            (d[i][α]*d[j][β]/(d[i][α]+d[j][β]))**2*2*ζ[2]**2/(1-ζ[3])**3
                    
                    # 计算成键分数B_iαiβ
                    # 对于常规组分, 第一项为1, 其余为0
                    if i == j and type[i] == "conventional":
                        B[i][0][j][0] = 1
                    # 对于共聚物, 参考Aspen文档"Hard-chain Fluids and Chain Connectivity"
                    if i == j and type[i] == "polymer":
                        if polymer[components[i]]["type"] == "random":
                            B[i][α][j][β] = z[i][α]*z[j][β]
                        elif polymer[components[i]]["type"] == "block": 
                            raise InputError("block类型还未实现")
                        elif polymer[components[i]]["type"] == "alternating":
                            raise InputError("alternating类型还未实现")
                        else:
                            raise InputError("聚合物类型输入错误, 只能为random或block或alternating")
                    
                    # 计算剩余自由能的成键项
                    if i == j:
                        ares_chain = ares_chain -x[i]*(m[i]-1)*B[i][α][j][β]*np.log(ghs[i][α][j][β])

                    # 计算色散项
                    # 给二元交互参数k_iαjβ赋值
                    if type[i] == "conventional":
                        if type[j] == "conventional":
                            k_iαjβ[i][0][j][0] = k_ij.get((components[i], components[j]),0)
                        elif type[j] == "polymer":  
                            k_iαjβ[i][0][j][β] = k_ij.get((components[i], segments[β]),0)
                    elif type[i] == "polymer":
                        if type[j] == "conventional":
                            k_iαjβ[i][α][j][0] = k_ij.get((segments[α], components[j]),0)
                        elif type[j] == "polymer":  
                            k_iαjβ[i][α][j][β] = k_ij.get((segments[α], segments[β]),0)

                    ε_iαjβ[i][α][j][β] = (1-k_iαjβ[i][α][j][β])*np.sqrt(ε[i][α]*ε[j][β])
                    σ_iαjβ[i][α][j][β] = (σ[i][α] + σ[j][β])/2
                    m2εσ3 = m2εσ3 + x[i]*x[j]*m[i]*m[j]*z[i][α]*z[j][β]*ε_iαjβ[i][α][j][β]/t*σ_iαjβ[i][α][j][β]**3
                    m2ε2σ3 = m2ε2σ3 + x[i]*x[j]*m[i]*m[j]*z[i][α]*z[j][β]*(ε_iαjβ[i][α][j][β]/t)**2*σ_iαjβ[i][α][j][β]**3   

    ares_hc = m_avg*ares_hs + ares_chain

    a = a0 + (m_avg-1)/m_avg*a1 + (m_avg-1)/m_avg*(m_avg-2)/m_avg*a2  # eq A.16
    b = b0 + (m_avg-1)/m_avg*b1 + (m_avg-1)/m_avg*(m_avg-2)/m_avg*b2  # eq A.17
    η = ζ[3]
    idx = np.arange(7)
    I1 = np.sum(a*η**idx)   # eq A.14
    I2 = np.sum(b*η**idx)   # eq A.15
    # eq A.13
    C1 = 1/(1 + m_avg*(8*η-2*η**2)/(1-η)**4 + (1-m_avg)*(20*η-27*η**2+12*η**3-2*η**4)/((1-η)*(2-η))**2)

    ares_disp = -2*π*ρ*I1*m2εσ3 - π*ρ*m_avg*C1*I2*m2ε2σ3    # eq A.12

    # Association term -------------------------------------------------------
    # To implement
    ares_assoc = 0

    # Dipole term (Gross and Vrabec term) --------------------------------------
    # To implement
    ares_polar = 0
    
    # Ion term ---------------------------------------------------------------
    # To implement
    ares_ion = 0

    # 汇总各项结果，得到体系的剩余Helmholtz自由能
    ares = ares_hc + ares_disp + ares_polar + ares_assoc + ares_ion

    return ares




if __name__ == "__main__":
    # 形式一: 列表
    # 分子参数： "TICL4", "TEA", "C2H4", "C2H6" ,"C3H6", "C3H8", "H2", "N2", "PP"
    mole_args = {
        "components": np.array(['TiCl4', 'TEA', 'C2H4', 'C2H6', 'C3H6', 'C3H8', 'H2', 'N2', 'PP']),
        "type": np.array(["conventional", "conventional", "conventional", "conventional", "conventional", "conventional", "conventional", "conventional", "polymer"]),
        "m": np.array([20, 20, 1.593, 1.6069, 1.9597, 2.002, 0.8285, 1.2053, None]),
        "s": np.array([3.54, 3.54,	3.445, 3.5206, 3.5356, 3.6184, 2.9729,	3.313, None]),
        "e": np.array([210, 210, 176.47, 191.42, 207.19, 208.11, 12.53,	90.96, None]),
        "w": np.array([0.283732, 0.841783, 0.0862484, 0.099493,	0.137588, 0.152291,	-0.215993,	0.0377215,	0]), 
        "Tc": np.array([638, 678.15, 282.34, 305.32, 364.85, 369.83, 33.19,	126.2, 2000]), 
        "Pc": np.array([4660950, 8930000, 5041000, 4872000,	4600000, 4248000, 1313000, 3400000,	5000000]), 
        "mw": np.array([189.6908, 114.16664, 28.05376, 30.06964, 42.08064, 44.09652, 2.01588, 28.01348,	42.08064])}

    # 链段参数：C2H4-R, C3H6-R 
    seg_args = {"segments": np.array(["C2H4-R", "C3H6-R"]),
                "type":np.array(["repeat", "repeat"]),
                "r": np.array([0.04132, 0.02528]),
                "s": np.array([3.4751, 4.1473]),
                "e": np.array([267.1854, 298.6392]),
                "mw": np.array([28.05376, 42.08064])}

    # 二元交互参数
    k_ij = {
        ("C2H4", "C3H6"): -0.001,
        ("C3H6", "C2H4"): -0.001,
        ("C2H4", "N2"): 0.075,
        ("N2", "C2H4"): 0.075,
        ("C3H6", "C2H4-R"): 0.029,
        ("C2H4-R", "C3H6"): 0.029,
        ("C3H6", "C3H6-R"): 0.029,
        ("C3H6-R", "C3H6"): 0.029,
        ("C3H8", "C2H4-R"): 0.0206,
        ("C2H4-R", "C3H8"): 0.0206,
        ("C2H4-R", "C3H6-R"): -0.009,
        ("C3H6-R", "C2H4-R"): -0.009,
    }

    # 组分： "TICL4", "TEA", "C2H4", "C2H6" ,"C3H6", "C3H8", "H2", "N2", "PP"
    x = np.array([2.35514E-06, 0.00016953, 0.000311436, 5.84603E-06, 0.00104029, 8.75865E-05, 2.55071E-09, 0.000467398, 0.9979156])
    print(sum(x))
    # 指明聚合物的信息：
    # 聚合物物种：数均聚合度, 链段分数, 链段的连接类型(随机random/嵌段block/交替alternating)
    polymer = {"PP": {"DPN": 2934.836, "SFRAC": np.array([0.1112102, 0.8887898]), "type":"random"}}
    t = 273.15
    ρ = 0.01/2934.84    # 液相摩尔密度(换算为聚合物的)
