""" CSTR反应器生产聚丙烯的模拟 """

# 导包
import numpy as np
from polymer_model.database.pure_components import args   # 导入物性数据
from polymer_model.kinetics.ziegler_nat.cstr import *     # 导入CSTR反应器


# 设置组分参数：C3H6、C3H8、H2、N2、PP、TiCl4、TEA、C3H6-R
components = {"C3H6": {"type": "conventional"},
    "C3H8": {"type": "conventional"},
    "H2": {"type": "conventional"},
    "N2": {"type": "conventional"},
    "PP": {"type": "polymer"},
    "TiCl4": {"type": "conventional"},
    "TEA": {"type": "conventional"},
    "C3H6-R": {"type": "segment"},
    "H2O": {"type": "conventional"},
}

# 设置聚合物链段
segments = {"C3H6-R":{"type": "repeat"}}

# 设置催化剂
catalyst = {"TiCl4": {"type": "Z-N",      # 催化剂类型
                    "site_types_num":4,   # 位点数
                    "site_conc":0.45,      # 位点浓度, mol/kgcat
            }}

# 设置聚合物
polymers = {
    "segments": segments,
    "catalyst": catalyst,
}

# 设置反应器进料
feed = {'T': 330.7119, 'P': 3000000, 'Mole Flow': 7012.333, 
        'Mole Fraction': {'C3H6': 0.7695529, 'C3H8': 0.1446299, 
        'H2': 0.0786014, 'N2': 0.00721171, 'PP': 0.0, 
        'TiCl4': 6.26484E-07, 'TEA': 3.46973E-06, 'H2O': 0.0}, 
        'TiCl4': {'CPSFLOW': 0.000375, 
                'CDSFLOW': 0, 
                'CISFLOW': np.array([0, 0, 0, 0]), 
                'CVSFLOW': np.array([0, 0, 0, 0])}}

# 设置Ziegler-Natta反应物种
species = { "polymer": "PP",          # 聚合物
            "tdb segment": None,      # 终端双键(用于计算支链)
            "monomers": "C3H6",       # 单体
            "segments": "C3H6-R",     # 链段
            "precatalyst": None,      # 预催化剂
            "catalyst": "TiCl4",      # 催化剂
            "cocatalysts": "TEA",     # 助催化剂
            "solvents": None,         # 溶剂
            "transfer agent": None,   # 链转移剂
            "hydrogens": "H2",        # 氢
            "poisons": None,          # 毒物
            "elec don": None,         # 电子供体
            "byproduct": None,        # 副产物
        }

# 设置Z-N反应：反应类型, 催化位点, 组分1, 组分2, 前指因子, 活化能, 反应级数, 终端双键分数, 参考温度
r1 = [["Act-Spon", 1, "TiCl4", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Act-Spon", 2, "TiCl4", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Act-Spon", 3, "TiCl4", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Act-Spon", 4, "TiCl4", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Chain-Ini", 1, "C3H6", None, 108.85/1000, 3.0145e4, 1, None, 343.15],
    ["Chain-Ini", 2, "C3H6", None, 24.5/1000,   3.0145e4, 1, None, 343.15],
    ["Chain-Ini", 3, "C3H6", None, 170.8/1000,  3.0145e4, 1, None, 343.15],
    ["Chain-Ini", 4, "C3H6", None, 60.55/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 1, "C3H6-R", "C3H6", 108.85/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 2, "C3H6-R", "C3H6", 24.5/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 3, "C3H6-R", "C3H6", 170.8/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 4, "C3H6-R", "C3H6", 60.55/1000,  3.0145e4, 1, None, 343.15],
    ["Chat-Mon", 1, "C3H6-R", "C3H6", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Mon", 2, "C3H6-R", "C3H6", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Mon", 3, "C3H6-R", "C3H6", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Mon", 4, "C3H6-R", "C3H6", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Cocat", 1, "C3H6-R", "TEA", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-Cocat", 2, "C3H6-R", "TEA", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-Cocat", 3, "C3H6-R", "TEA", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-Cocat", 4, "C3H6-R", "TEA", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-H2", 1, "C3H6-R", "H2", 4.8/1000,  4.47988e4, 1, None, 343.15],
    ["Chat-H2", 2, "C3H6-R", "H2", 8.88/1000,  4.47988e4, 1, None, 343.15],
    ["Chat-H2", 3, "C3H6-R", "H2", 2.64/1000,  4.47988e4, 1, None, 343.15],
    ["Chat-H2", 4, "C3H6-R", "H2", 6.6/1000,  4.47988e4, 1, None, 343.15],
    ["Deact-Spon", 1, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Deact-Spon", 2, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Deact-Spon", 3, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Deact-Spon", 4, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Atact-Prop", 1, "C3H6-R", "C3H6", 8/1000,  3.0145e4, 1, None, 343.15],
    ["Atact-Prop", 2, "C3H6-R", "C3H6", 1/1000,  3.0145e4, 1, None, 343.15],
    ["Atact-Prop", 3, "C3H6-R", "C3H6", 0.1/1000,  3.0145e4, 1, None, 343.15],
    ["Atact-Prop", 4, "C3H6-R", "C3H6", 0.1/1000,  3.0145e4, 1, None, 343.15]]

# 设置反应器参数
specs = {
    "P": 3000000,
    "T": 333.15,
    "Valid phases": "vapor-liquid",
    "Specification type": "Reactor volume & Phase volume",
    "Reactor volume": 90,
    "Condensed phase volume": 60, 
    "Reactions": {"Type":"Z-N", 
                    "Species":species, 
                    "Reactions":r1, 
                    "Reacting phase":"liquid"},
    "Streams": {"vapor", "liquid"},
}


# 运行计算
print("\nCSTR Product")
product, duty = cstr(components, args, polymers, [feed], specs)
print("liquid: ",product[0])
print("vapor: ",product[1])
print(duty)

# 计算链长分布
NC, Ns, sol = cstr_run(components,args, polymers, [feed], specs)
n = get_points(Np=100, upper=100000)
eqs, τ, SFMOM, mw = model(sol.x, components, args, polymers, [feed], specs) 
w = cld(τ, SFMOM, n, GPC=True)
print("局部链长分布: \n", w)

# 绘制链长分布图
cld_plot(n, w)

# 计算分子量分布
Mn = n*mw
w = mwd(τ, SFMOM, mw, Mn)
print("局部分子量分布: \n", w)

# 绘制分子量分布图
mwd_plot(Mn,w)