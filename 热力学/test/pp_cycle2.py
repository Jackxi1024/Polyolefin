"""
聚丙烯生产全流程的模拟, 目的是测试循环的求解
Wegstein法
"""

# 导包
from polymer_model.thermodynamics.thermo import *      # 导入混合、加热、压缩模块
from polymer_model.kinetics.ziegler_nat.cstr import *  # 导入反应模块
from polymer_model.thermodynamics.args import args 
from polymer_model.thermodynamics.reaction import species, r1   # 导入反应包
from polymer_model.solver import *   # 导入求解器
from polymer_model.utility import SolutionError
from scipy import optimize             # 用于求解CSTR模型(非线性方程组)

# 全局属性设置
# 设置组分参数：C3H6、C3H8、H2、N2、PP、TiCl4、TEA、C3H6-R 
components = {
    "C3H6": {"type": "conventional"},
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

NC = 8  # 真实组分数

# 以下物流是混合器FMIX1的进出口物流
# 设置流股参数: 温度、压力、摩尔流量、摩尔分数、
# 如果流股中有聚合物，则需要设置聚合物的数均聚合度DPN和链段分数SFRAC
N2FEED = {"T": 303.15, "P": 3000000, "Mole Flow": 2.37587,  
    "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2": 0,
        "N2": 1, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
}

H2FEED = {"T": 303.15, "P": 3000000, "Mole Flow": 11.68827,  
    "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2": 1,
        "N2": 0, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
}

CAT = {"T": 303.15, "P": 3000000, "Mole Flow": 1.964734,  
    "Mole Fraction": { "C3H6": 0.9958403, "C3H8": 0.00192372, "H2": 0,
        "N2": 0, "PP": 0, "TiCl4": 0.00223599, "TEA": 0, "H2O":0},
    "TiCl4": {
            "CPSFLOW":0.000375, 
            "CDSFLOW":0., 
            "CISFLOW":np.array([0., 0., 0., 0.]), 
            "CVSFLOW":np.array([0., 0., 0., 0.])
        },
    "Valid phases": "liquid"
}

COCAT = {"T": 303.15, "P": 3000000, "Mole Flow": 0.0243309,  
    "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2": 0,
        "N2": 0, "PP": 0, "TiCl4": 0, "TEA": 1, "H2O":0},
    "Valid phases": "liquid"
}

C3FEED = {"T": 303.15, "P": 3000000, "Mole Flow": 2004.329,  
    "Mole Fraction": { "C3H6": 0.8331688, "C3H8": 0.1668312, "H2": 0,
        "N2": 0, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
}


cycle_in = dict()
cycle_in["Mole Fraction"] = dict()

# 初始化循环物流
# cycle_in = {"T": 334.8567, "P": 3200000, "Mole Flow": 4991.95,  
#     "Mole Fraction": { "C3H6": 0.7460931, "C3H8": 0.1361802, "H2": 0.1080722,
#         "N2": 0.00965455, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
# }


def fun(z):

    cycle_in["T"] = z[0]
    cycle_in["P"] = z[1]
    cycle_in["Mole Flow"] = np.sum(z[2:])
    cycle_in["Mole Fraction"]["C3H6"] = z[2]/cycle_in["Mole Flow"]
    cycle_in["Mole Fraction"]["C3H8"] = z[3]/cycle_in["Mole Flow"]
    cycle_in["Mole Fraction"]["H2"] = z[4]/cycle_in["Mole Flow"]
    cycle_in["Mole Fraction"]["N2"] = z[5]/cycle_in["Mole Flow"]
    z[6] = 0
    cycle_in["Mole Fraction"]["PP"] = z[6]/cycle_in["Mole Flow"]
    cycle_in["Mole Fraction"]["TiCl4"] = z[7]/cycle_in["Mole Flow"]
    cycle_in["Mole Fraction"]["TEA"] = z[8]/cycle_in["Mole Flow"]
    cycle_in["Mole Fraction"]["H2O"] = z[9]/cycle_in["Mole Flow"]

    print("cycle in")
    print(cycle_in)

    # 混合器的入口物流
    streams = [N2FEED, H2FEED, CAT, COCAT, C3FEED, cycle_in]

    # 混合器出口物流
    RFeed1 = mix(components, catalyst, args, streams)

    # CSTR反应器
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
    product, duty = cstr(components, args, polymers, [RFeed1], specs)
    POWDER1 = product[0]
    VAP_A = product[1]

    # 压缩机
    pressure_spec = {"Pressure increase": 200000}
    VAP_B, T_isentropic, isentropic_compression_work, true_compression_work, net_work = compress_isentropic(components, catalyst, args, [VAP_A], pressure_spec)

    # 多流股换热器
    # 设置进出口
    hot_streams = [VAP_B]
    H2O_IN = {"T": 20+273.15, "P": 3000000, "Mole Flow": 13067.61,  
        "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2": 0,
            "N2": 0, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":1}
    }
    cold_streams = [H2O_IN]
    specs = {
        "hot":[{"spec":None, "P":0, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01},],
        "cold": [{"spec":{"T":30+273.15}, "P":0, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01},]
    }
    # 运行计算
    hot_out, hot_duty, cold_out, cold_duty, duty = multi_stream_heat_exchang(components, catalyst, args, hot_streams, cold_streams, specs)

    cycle_out = hot_out[0]
    print("cycle_out: \n", cycle_out)

    res = np.zeros(10)
    res[0] = cycle_out["T"] - z[0]
    res[1] = cycle_out["P"] - z[1]
    mole_flow = cycle_out["Mole Flow"]
    res[2] = cycle_out["Mole Fraction"]["C3H6"]*mole_flow - z[2]
    res[3] = cycle_out["Mole Fraction"]["C3H8"]*mole_flow - z[3]
    res[4] = cycle_out["Mole Fraction"]["H2"]*mole_flow - z[4]
    res[5] = cycle_out["Mole Fraction"]["N2"]*mole_flow - z[5]
    res[6] = cycle_out["Mole Fraction"]["PP"]*mole_flow - z[6]
    res[7] = cycle_out["Mole Fraction"]["TiCl4"]*mole_flow - z[7]
    res[8] = cycle_out["Mole Fraction"]["TEA"]*mole_flow - z[8]
    res[9] = cycle_out["Mole Fraction"]["H2O"]*mole_flow - z[9]
    print("res: ", res)
    return res


# x0 = [334.8567, 3200000, 3724.461, 679.8047, 539.4908, 48.19503, 0, 0, 0, 0]
x0 = np.array([334, 3200000, 3700, 700, 500, 50, 0, 0, 0, 0]) 

# sol = solve(fun, x0, args=(), method='wegstein', max_iter=100, tol=1e-3, w_min=0, w_max=5)
# if sol["success"]:
#     z = sol["x"]
#     iter = sol["iter"]
#     print("z: ", z)
#     print("iter: ", iter)
# else:
#     raise SolutionError("求解失败")

sol = optimize.root(fun, x0, method='broyden1', tol=1e-5)
if sol.success :
    print("Number of iterations: ",sol.nfev)
    print("sol.x: ", sol.x)
else : 
    raise SolutionError(sol.message)

print("cycle in")
print(cycle_in)

# 混合器的入口物流
streams = [N2FEED, H2FEED, CAT, COCAT, C3FEED, cycle_in]

# 混合器出口物流
RFeed1 = mix(components, catalyst, args, streams)

# CSTR反应器
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
product, duty = cstr(components, args, polymers, [RFeed1], specs)
POWDER1 = product[0]
VAP_A = product[1]

# 压缩机
pressure_spec = {"Pressure increase": 200000}
VAP_B, T_isentropic, isentropic_compression_work, true_compression_work, net_work = compress_isentropic(components, catalyst, args, [VAP_A], pressure_spec)

# 多流股换热器
# 设置进出口
hot_streams = [VAP_B]
H2O_IN = {"T": 20+273.15, "P": 3000000, "Mole Flow": 13067.61,  
    "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2": 0,
        "N2": 0, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":1}
}
cold_streams = [H2O_IN]
specs = {
    "hot":[{"spec":None, "P":0, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01},],
    "cold": [{"spec":{"T":30+273.15}, "P":0, "Valid phases": "vapor-liquid", "Max iter":100, "Tol":0.01},]
}
# 运行计算
hot_out, hot_duty, cold_out, cold_duty, duty = multi_stream_heat_exchang(components, catalyst, args, hot_streams, cold_streams, specs)

cycle_out = hot_out[0]
print("cycle_out: \n", cycle_out)

# 混合器出口物流
print("\nRFeed1: ")
print(RFeed1)
print("Mole Entropy: ")
print(molar_enthalpy_of_mixture(components, args, RFeed1)[2])


# 反应器出口物流
print("\nCSTR Product")
print("气相出口: \n", VAP_A)
print("液相出口: \n", POWDER1)


# 压缩机
print("\nVAP-B: ")
print(VAP_B)
print("T_isentropic: ", T_isentropic)
print("等熵压缩功: ", isentropic_compression_work)
print("真实压缩功: ", true_compression_work)
print("压缩机净功率: ", net_work)

# 多流股换热器
print("\n多流股换热器: ")
print("热物流出口: \n", hot_out)
print("冷物流出口: \n", cold_out)
print("热物流放热: ", hot_duty)
print("冷物流吸热: ", cold_duty)


print("\n闪蒸罐 STRIP1")
flash_type = {"T": 65+273.15, "P":500000}
POWDER2, GAS1, beta, duty = flash(components, catalyst, args, [POWDER1], flash_type, max_iter=300)
print("气化分数: ", beta)
print("热负荷: ", duty)
print("气相: \n", GAS1)
print("液相: \n", POWDER2)


print("\n闪蒸罐 STRIP2")
# N2抽提
STRIP_N2 = {"T": 303.15, "P": 300000, "Mole Flow": 4.957931,  
    "Mole Fraction": { "C3H6": 0, "C3H8": 0, "H2": 0,
        "N2": 1, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
}
flash_type = {"T": 60+273.15, "P":100000}
product, GAS2, beta, duty = flash(components, catalyst, args, [POWDER2, STRIP_N2], flash_type, max_iter=300)
print("气化分数: ", beta)
print("热负荷: ", duty)
print("气相: \n", GAS2)
print("液相: \n", product)



