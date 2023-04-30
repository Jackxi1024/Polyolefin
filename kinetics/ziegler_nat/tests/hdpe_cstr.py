""" CSTR反应器生产HDPE的模拟 """

# 导包
import numpy as np
from polymer_model.database.pure_components import args   # 导入物性数据
from polymer_model.kinetics.ziegler_nat.temp.cstr3 import *     # 导入CSTR反应器


# 反应器1的计算
# 设置组分
components = {
    "C2H4": {"type": "conventional"},
    "H2": {"type": "conventional"},
    "Hexane": {"type": "conventional"},
    "HDPE": {"type": "polymer"},
    "TiCl4": {"type": "conventional"},
    "TEA": {"type": "conventional"},
    "C2H4-R": {"type": "segment"},
}

# 设置聚合物链段
segments = {"C2H4-R":{"type": "repeat"}}

# 设置催化剂
catalyst = {"TiCl4": {"type": "Z-N", 
                    "site_types_num":4,   # 位点数
                    "site_conc":0.1,      # 位点浓度, mol/kgcat
            }}

polymers = {
    "segments": segments,
    "catalyst": catalyst,
}

# 设置进料参数: 温度、压力、摩尔流量、摩尔分数、催化剂各位点流量、入口聚合物各阶矩流量
feed = {
    "T": 70+273.15,
    "P": 20265000,
    "Mole Flow": 233.4925,  
    "Mole Fraction": {
        "C2H4": 0.2544396,
        "H2": 0.000177044,
        "Hexane": 0.7452519,
        "HDPE": 0,
        "TiCl4": 3.76296E-05,
        "TEA": 9.37838E-05,
    },
    "TiCl4": {
        "CPSFLOW":0.000166667, 
        "CDSFLOW":0, 
        "CISFLOW":np.array([0, 0, 0, 0]), 
        "CVSFLOW":np.array([0, 0, 0, 0])
    },
    "HDPE": {
        "LSZMOM": np.array([0, 0, 0, 0]),
        "DSZMOM": np.array([0, 0, 0, 0]),
        "LSFMOM": np.array([0, 0, 0, 0]),
        "DSFMOM": np.array([0, 0, 0, 0]),
        "LSSMOM": np.array([0, 0, 0, 0]),
        "DSSMOM": np.array([0, 0, 0, 0]),
    }
}

# 设置Z-N反应物种
species = { "polymer": "HDPE",
            "tdb segment": None,
            "monomers": "C2H4",
            "segments": "C2H4-R",
            "precatalyst": None,
            "catalyst": "TiCl4",
            "cocatalysts": "TEA",
            "solvents": "Hexane",
            "transfer agent": None,
            "hydrogens": "H2",
            "poisons": None,
            "elec don": None,
            "byproduct": None,
        }


# 设置Ziegler-Natta反应动力学参数
# 反应类型, 催化位点, 组分1, 组分2, 前指因子, 活化能, 反应级数, 终端双键分数, 参考温度
r1 = [["Act-Spon", 1, "TiCl4", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 2, "TiCl4", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 3, "TiCl4", None, 0, 0, 1, None, 1e35],
    ["Act-Spon", 4, "TiCl4", None, 0, 0, 1, None, 1e35],
    ["Act-Cocat", 1, "TiCl4", "TEA", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 2, "TiCl4", "TEA", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 3, "TiCl4", "TEA", 0, 0, 1, None, 1e35],
    ["Act-Cocat", 4, "TiCl4", "TEA", 0, 0, 1, None, 1e35],
    ["Chain-Ini", 1, "C2H4", None, 255/1000, 0, 1, None, 1e35],
    ["Chain-Ini", 2, "C2H4", None, 90/1000,   0, 1, None, 1e35],
    ["Chain-Ini", 3, "C2H4", None, 0,  0, 1, None, 1e35],
    ["Chain-Ini", 4, "C2H4", None, 0,  0, 1, None, 1e35],
    ["Propagation", 1, "C2H4-R", "C2H4", 255/1000,  0, 1, None, 1e35],
    ["Propagation", 2, "C2H4-R", "C2H4", 90/1000,  0, 1, None, 1e35],
    ["Propagation", 3, "C2H4-R", "C2H4", 0,  0, 1, None, 1e35],
    ["Propagation", 4, "C2H4-R", "C2H4", 0,  0, 1, None, 1e35],
    ["Chat-Mon", 1, "C2H4-R", "C2H4", 0.09/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 2, "C2H4-R", "C2H4", 0.24/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 3, "C2H4-R", "C2H4", 0,  0, 1, None, 1e35],
    ["Chat-Mon", 4, "C2H4-R", "C2H4", 0,  0, 1, None, 1e35],
    ["Chat-H2", 1, "C2H4-R", "H2", 5.55/1000,  0, 1, None, 1e35],
    ["Chat-H2", 2, "C2H4-R", "H2", 18.5/1000,  0, 1, None, 1e35],
    ["Chat-H2", 3, "C2H4-R", "H2", 0,  0, 1, None, 1e35],
    ["Chat-H2", 4, "C2H4-R", "H2", 0,  0, 1, None, 1e35],
    ["Chat-Spon", 1, "C2H4-R", None, 0.004,  0, 1, None, 1e35],
    ["Chat-Spon", 2, "C2H4-R", None, 0.012,  0, 1, None, 1e35],
    ["Chat-Spon", 3, "C2H4-R", None, 0,  0, 1, None, 1e35],
    ["Chat-Spon", 4, "C2H4-R", None, 0,  0, 1, None, 1e35],
    ["Deact-Spon", 1, None, None, 0.0001,  0, 1, None, 1e35],
    ["Deact-Spon", 2, None, None, 0.0006,  0, 1, None, 1e35],
    ["Deact-Spon", 3, None, None, 0,  0, 1, None, 1e35],
    ["Deact-Spon", 4, None, None, 0,  0, 1, None, 1e35]]


# 设置反应器参数
specs = {
    "P": 20265000,
    "T": 433.15,
    "Valid phases": "liquid",
    "Specification type": "Reactor volume",
    "Reactor volume": 60,
    "Reactions": {"Type":"Z-N", 
                    "Species":species, 
                    "Reactions":r1,
                    "Reacting phase": "liquid"
                }
}


# 运行计算
print("\nCSTR-1 Product")
product1, duty1 = cstr(components, args, polymers, [feed], specs)
print("product: ",product1[0])
print(duty1)

# 计算链长分布
NC, Ns, sol = cstr_run(components,args, polymers, [feed], specs)
n = get_points(Np=100, upper=100000)
eqs, τ1, SFMOM1, mw = model(sol.x, components, args, polymers, [feed], specs) 
w1 = cld(τ1, SFMOM1, n, GPC=True)
print("局部链长分布: \n", w1)
# 输出每个位点前8个数据, 其中第一行为各位点总和
print(w1[:, 0:8])   

# 绘制链长分布图
cld_plot(n, w1)

# 计算分子量分布
Mn = n*mw
w1 = mwd(τ1, SFMOM1, mw, Mn)
print("局部分子量分布: ")
# 输出每个位点前8个数据, 其中第一行为各位点总和
print(w1[:, 0:8])   

# 绘制分子量分布图
mwd_plot(Mn,w1)


# 反应器2的计算
# 补充催化剂进料
feed2 = {
    "T": 70+273.15,
    "P": 20265000,
    "Mole Flow": 45.85611,  
    "Mole Fraction": {
        "C2H4": 0.4318562,
        "H2": 0.00600988,
        "Hexane": 0.5620224,
        "HDPE": 0,
        "TiCl4": 3.19341E-05,
        "TEA": 7.95889E-05,
    },
    "TiCl4": {
        "CPSFLOW":2.77778E-05, 
        "CDSFLOW":0, 
        "CISFLOW":np.array([0, 0, 0, 0]), 
        "CVSFLOW":np.array([0, 0, 0, 0])
    },
}


# 设置Ziegler-Natta反应动力学参数
# 反应类型, 催化位点, 组分1, 组分2, 前指因子, 活化能, 反应级数, 终端双键分数, 参考温度
r1 = [["Act-Spon", 1, "TiCl4", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 2, "TiCl4", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 3, "TiCl4", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 4, "TiCl4", None, 0.08, 0, 1, None, 1e35],
    ["Act-Cocat", 1, "TiCl4", "TEA", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 2, "TiCl4", "TEA", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 3, "TiCl4", "TEA", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 4, "TiCl4", "TEA", 0.15/1000, 0, 1, None, 1e35],
    ["Chain-Ini", 1, "C2H4", None, 255/1000, 0, 1, None, 1e35],
    ["Chain-Ini", 2, "C2H4", None, 90/1000,   0, 1, None, 1e35],
    ["Chain-Ini", 3, "C2H4", None, 255/1000,  0, 1, None, 1e35],
    ["Chain-Ini", 4, "C2H4", None, 90/1000,  0, 1, None, 1e35],
    ["Propagation", 1, "C2H4-R", "C2H4", 255/1000,  0, 1, None, 1e35],
    ["Propagation", 2, "C2H4-R", "C2H4", 90/1000,  0, 1, None, 1e35],
    ["Propagation", 3, "C2H4-R", "C2H4", 255/1000,  0, 1, None, 1e35],
    ["Propagation", 4, "C2H4-R", "C2H4", 90/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 1, "C2H4-R", "C2H4", 0.09/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 2, "C2H4-R", "C2H4", 0.24/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 3, "C2H4-R", "C2H4", 0.09/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 4, "C2H4-R", "C2H4", 0.24/1000,  0, 1, None, 1e35],
    ["Chat-H2", 1, "C2H4-R", "H2", 5.55/1000,  0, 1, None, 1e35],
    ["Chat-H2", 2, "C2H4-R", "H2", 18.5/1000,  0, 1, None, 1e35],
    ["Chat-H2", 3, "C2H4-R", "H2", 5.55/1000,  0, 1, None, 1e35],
    ["Chat-H2", 4, "C2H4-R", "H2", 18.5/1000,  0, 1, None, 1e35],
    ["Chat-Spon", 1, "C2H4-R", None, 0.004,  0, 1, None, 1e35],
    ["Chat-Spon", 2, "C2H4-R", None, 0.012,  0, 1, None, 1e35],
    ["Chat-Spon", 3, "C2H4-R", None, 0.004,  0, 1, None, 1e35],
    ["Chat-Spon", 4, "C2H4-R", None, 0.012,  0, 1, None, 1e35],
    ["Deact-Spon", 1, None, None, 0.0001,  0, 1, None, 1e35],
    ["Deact-Spon", 2, None, None, 0.0006,  0, 1, None, 1e35],
    ["Deact-Spon", 3, None, None, 0.0001,  0, 1, None, 1e35],
    ["Deact-Spon", 4, None, None, 0.0006,  0, 1, None, 1e35]]


# 设置反应器参数
specs = {
    "P": 20265000,
    "T": 433.15,
    "Valid phases": "liquid",
    "Specification type": "Reactor volume",
    "Reactor volume": 60,
    "Reactions": {"Type":"Z-N", 
                    "Species":species, 
                    "Reactions":r1,
                    "Reacting phase": "liquid"
                }
}


# 运行计算
print("\nCSTR-2 Product")
product2, duty2 = cstr(components, args, polymers, [product1[0], feed2], specs)
print("product: ",product2[0])
print(duty2)

# 计算链长分布
NC, Ns, sol = cstr_run(components,args, polymers, [product1[0], feed2], specs)
n = get_points(Np=100, upper=100000)
eqs, τ2, SFMOM2, mw = model(sol.x, components, args, polymers, [feed], specs) 
w2 = cld(τ2, SFMOM2, n, GPC=True)
print("局部链长分布: \n", w2)
# 输出每个位点前8个数据, 其中第一行为各位点总和
print(w2[:, 0:8])   

# 绘制链长分布图
cld_plot(n, w2)

# 计算分子量分布
Mn = n*mw
w2 = mwd(τ2, SFMOM2, mw, Mn)
print("局部分子量分布: ")
# 输出每个位点前8个数据, 其中第一行为各位点总和
print(w2[:, 0:8])   

# 绘制分子量分布图
mwd_plot(Mn,w2)


# 计算反应器2的累计链长分布
w = cld_cum(τ1, SFMOM1, τ2, SFMOM2, n, GPC=True)
print("累计链长分布: ")
print(w[:, 0:8])

# 绘制累计链长分布图
cld_plot(n, w)

# 计算反应器2的累计分子量分布
w = mwd_cum(τ1, SFMOM1, τ2, SFMOM2, mw, Mn)
print("累计分子量分布: ")
print(w[:, 0:8])

# 绘制累计分子量分布图
mwd_plot(Mn, w)