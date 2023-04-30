"""
聚丙烯生产全流程的计算, 目的是测试每个单元块的相互连接
"""

# 导包
from polymer_model.thermodynamics.thermo import *      # 导入混合、加热、压缩模块
from polymer_model.kinetics.ziegler_nat.cstr import *  # 导入反应模块
from polymer_model.database.pure_components import args  # 导入纯组分数据
from polymer_model.thermodynamics.reaction import species, r1   # 导入反应包
import time

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
            "CDSFLOW":0, 
            "CISFLOW":np.array([0, 0, 0, 0]), 
            "CVSFLOW":np.array([0, 0, 0, 0])
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

# 初始化循环物流
CYCGAGB = {"T": 334.8567, "P": 3200000, "Mole Flow": 4991.95,  
    "Mole Fraction": { "C3H6": 0.7460931, "C3H8": 0.1361802, "H2": 0.1080722,
        "N2": 0.00965455, "PP": 0, "TiCl4": 0, "TEA": 0, "H2O":0}
}

# tol = 0.001
# error = 2*tol
# while error > tol:

# 混合器的入口物流
streams = [N2FEED, H2FEED, CAT, COCAT, C3FEED, CYCGAGB]

# 混合器出口物流
print("\nRFeed1: ")
RFeed1 = mix(components, catalyst, args, streams)
print(RFeed1)

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
print("\nCSTR Product")
product, duty = cstr(components, args, polymers, [RFeed1], specs)
POWDER1 = product[0]
VAP_A = product[1]
print("气相出口: \n", VAP_A)
print("液相出口: \n", POWDER1)
print("负荷: ", duty)

# 压缩机
pressure_spec = {"Pressure increase": 200000}

print("\nVAP-B: ")
VAP_B, T_isentropic, isentropic_compression_work, true_compression_work, net_work = compress_isentropic(components, catalyst, args, [VAP_A], pressure_spec)
print(VAP_B)
print("T_isentropic: ", T_isentropic)
print("等熵压缩功: ", isentropic_compression_work)
print("真实压缩功: ", true_compression_work)
print("压缩机净功率: ", net_work)

# 多流股换热器
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
print("\n多流股换热器: ")
hot_out, hot_duty, cold_out, cold_duty, duty = multi_stream_heat_exchang(components, catalyst, args, hot_streams, cold_streams, specs)
print("热物流出口: \n", hot_out)
print("冷物流出口: \n", cold_out)
print("热物流放热: ", hot_duty)
print("冷物流吸热: ", cold_duty)
print("换热器负荷: ", duty)


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



