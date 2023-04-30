""" 
对聚丙烯生产过程进行模拟, 测试本程序各部件的正确性
闭环
"""

# 导包
import numpy as np
from polymer_model.flowsheet.flowsheet import Flowsheet
from polymer_model.unit_model.mixer import Mixer
from polymer_model.unit_model.cstr import CSTR
from polymer_model.unit_model.compressor import Compressor
from polymer_model.unit_model.mheatx import MHeatX
from polymer_model.unit_model.flash import Flash
from polymer_model.unit_model.stream import MaterialStream


# 流程的组分：C3H6、C3H8、H2、N2、PP、TiCl4、TEA、C3H6-R, H2O
components = {
    "Propylene": {"type": "conventional"},
    "Propane": {"type": "conventional"},
    "Hydrogen": {"type": "conventional"},
    "Nitrogen": {"type": "conventional"},
    "Polypropylene": {"type": "polymer"},
    "Titanium Tetrachloride": {"type": "conventional"},
    "Triethyl Aluminum": {"type": "conventional"},
    "Propylene-R": {"type": "segment"},
    "Water": {"type": "conventional"},
}


# 流程的聚合物
polymers = {
    # 链段
    "segments": {"Propylene-R":{"type": "repeat"},
                }, 
    # 催化剂
    "catalyst": {"Titanium Tetrachloride": {"type": "Z-N",       # 催化剂类型
                        "site_types_num":4,  # 位点数
                        "site_conc":0.45,    # 位点浓度, mol/kgcat
                        }
                },
    # 聚合物分子量分布
    "Distribution": {"Number of points":100, "Upper limit":100000, "GPC":True}
}

# 物性方法
property_method = "PC-SAFT"

# 创建流程
fs = Flowsheet(components, polymers, property_method)

# 添加单元模块
# 添加多流股换热器
exc1 = MHeatX()
exc1.set_name("Exc1")
fs.add_block(exc1)
# 添加混合器
mix1 = Mixer()
mix1.set_name("Mix1")
fs.add_block(mix1)
# 添加反应器
cstr1 = CSTR()
cstr1.set_name("CSTR1")
fs.add_block(cstr1)
# 添加压缩机
comp1 = Compressor()
comp1.set_name("Comp1")
fs.add_block(comp1)
# 添加闪蒸罐1
flash1 = Flash()
flash1.set_name("Flash1")
fs.add_block(flash1)
# 添加闪蒸罐2
flash2 = Flash()
flash2.set_name("Flash2")
fs.add_block(flash2)

# 添加流股
# 混合器进料流股
N2_feed = MaterialStream(source=None, destination=mix1.inlet)
N2_feed.set_name("N2 Feed")
fs.add_stream(N2_feed)

H2_feed = MaterialStream(source=None, destination=mix1.inlet)
H2_feed.set_name("H2 Feed")
fs.add_stream(H2_feed)

Cat_feed = MaterialStream(source=None, destination=mix1.inlet)
Cat_feed.set_name("Cat Feed")
fs.add_stream(Cat_feed)

Cocat_feed = MaterialStream(source=None, destination=mix1.inlet)
Cocat_feed.set_name("Cocat Feed")
fs.add_stream(Cocat_feed)

C3_feed = MaterialStream(source=None, destination=mix1.inlet)
C3_feed.set_name("C3 Feed")
fs.add_stream(C3_feed)

cycle_gas = MaterialStream(source=exc1.hot_outlet, destination=mix1.inlet)
cycle_gas.set_name("Cycle Gas")
fs.add_stream(cycle_gas)

# 反应器进料流股
rfeed1 = MaterialStream(source=mix1.outlet, destination=cstr1.inlet)
rfeed1.set_name("Reactor Feed1")
fs.add_stream(rfeed1)

# 压缩机进料流股
vap_a = MaterialStream(source=cstr1.vapor_outlet, destination=comp1.inlet)
vap_a.set_name("Vap-A")
fs.add_stream(vap_a)

# 多流股反应器进出口流股
vap_b = MaterialStream(source=comp1.outlet, destination=exc1.hot_inlet)
vap_b.set_name("Vap-B")
fs.add_stream(vap_b)

# cycle_gas2 = MaterialStream(source=exc1.hot_outlet, destination=None)
# cycle_gas2.set_name("Cycle Gas2")
# fs.add_stream(cycle_gas2)

H2O_in = MaterialStream(source=None, destination=exc1.cold_inlet)
H2O_in.set_name("H2O-In")
fs.add_stream(H2O_in)

H2O_out = MaterialStream(source=exc1.cold_outlet, destination=None)
H2O_out.set_name("H2O-Out")
fs.add_stream(H2O_out)

# 闪蒸罐1进出料流股
powder1 = MaterialStream(source=cstr1.liquid_outlet, destination=flash1.inlet)
powder1.set_name("Powder1")
fs.add_stream(powder1)

gas1 = MaterialStream(source=flash1.vapor_outlet, destination=None)
powder1.set_name("Gas1")
fs.add_stream(gas1)

powder2 = MaterialStream(source=flash1.liquid_outlet, destination=flash2.inlet)
powder1.set_name("Powder2")
fs.add_stream(powder2)

# 闪蒸罐2进出料流股
N2_strip = MaterialStream(source=None, destination=flash2.inlet)
powder1.set_name("Strip N2")
fs.add_stream(N2_strip)

gas2 = MaterialStream(source=flash2.vapor_outlet, destination=None)
gas2.set_name("Gas2")
fs.add_stream(gas2)

product = MaterialStream(source=flash2.liquid_outlet, destination=None)
product.set_name("PP product")
fs.add_stream(product)

# 设置各进料流股参数
# 设置混合器进料流股参数
# C3H6、C3H8、H2、N2、PP、TiCl4、TEA、H2O
z1 = np.array([0, 0, 0, 1, 0, 0, 0, 0])
N2_feed.input(303.15, 3000000, None, 0.0665563, z1, "Mass & Mass Frac")

z2 = np.array([0, 0, 1, 0, 0, 0, 0, 0])
H2_feed.input(303.15, 3000000, None, 0.0235621, z2, "Mass & Mass Frac")

z3 = np.array([0.988, 0.002, 0, 0, 0, 0.01, 0, 0])
component_attribute = {
    "Titanium Tetrachloride": {
        "CPSFRAC": 1,
        "CDSFRAC": 0.,
        "CVSFRAC": np.array([0., 0., 0., 0.,]),
        "CISFRAC": np.array([0., 0., 0., 0.,]),

    },
}
Cat_feed.input(303.15, 3000000, None, 0.0833333, z3, "Mass & Mass Frac", component_attribute, valid_phases="liquid")
z4 = np.array([0, 0, 0, 0, 0, 0, 1, 0])
Cocat_feed.input(303.15, 3000000, None, 0.00277778, z4, "Mass & Mass Frac", valid_phases="liquid")
z5 = np.array([0.8265628, 0.1734372, 0, 0, 0, 0, 0, 0])
C3_feed.input(303.15, 3000000, None, 85.01753, z5, "Mass & Mass Frac")
# z6 = np.array([0.8286228, 0.1584892, 0.0057499, 0.00713806, 0, 0, 0, 0])
# cycle_gas.input(334.8567, 3200000, None, 189.1424, z6, "Mass & Mass Frac")

# 设置冷却水参数
z2 = np.array([0, 0, 0, 0, 0, 0, 0, 1])
H2O_in.input(293.15, 3000000, None, 235.4167, z2, "Mass & Mass Frac")

# 设置脱挥气参数
z2 = np.array([0, 0, 0, 0.1388889, 0, 0, 0, 0])
N2_strip.input(303.15, 300000, None, None, z2, "Mass & Mass Flow", valid_phases="vapor")

# 设置各单元模块操作参数
mix1.input()

# 设置反应器
# 设置Ziegler-Natta反应物种
species = { "polymer": "Polypropylene",          # 聚合物
            "tdb segment": None,      # 终端双键(用于计算支链)
            "monomers": ["Propylene"],       # 单体
            "segments": {"Propylene":"Propylene-R"},     # 链段
            "precatalyst": None,      # 预催化剂
            "catalyst": ["Titanium Tetrachloride"],      # 催化剂
            "cocatalysts": ["Triethyl Aluminum"],     # 助催化剂
            "solvents": None,         # 溶剂
            "transfer agent": None,   # 链转移剂
            "hydrogens": ["Hydrogen"],  # 氢
            "poisons": None,          # 毒物
            "elec don": None,         # 电子供体
            "byproduct": None,        # 副产物
        }

# 设置Z-N反应：反应类型, 催化位点, 组分1, 组分2, 前指因子, 活化能, 反应级数, 终端双键分数, 参考温度
r1 = [["Act-Spon", 1, "Titanium Tetrachloride", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Act-Spon", 2, "Titanium Tetrachloride", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Act-Spon", 3, "Titanium Tetrachloride", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Act-Spon", 4, "Titanium Tetrachloride", None, 0.0013, 3.19872e4, 1, None, 343.15],
    ["Chain-Ini", 1, "Propylene", None, 108.85/1000, 3.0145e4, 1, None, 343.15],
    ["Chain-Ini", 2, "Propylene", None, 24.5/1000,   3.0145e4, 1, None, 343.15],
    ["Chain-Ini", 3, "Propylene", None, 170.8/1000,  3.0145e4, 1, None, 343.15],
    ["Chain-Ini", 4, "Propylene", None, 60.55/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 1, "Propylene-R", "Propylene", 108.85/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 2, "Propylene-R", "Propylene", 24.5/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 3, "Propylene-R", "Propylene", 170.8/1000,  3.0145e4, 1, None, 343.15],
    ["Propagation", 4, "Propylene-R", "Propylene", 60.55/1000,  3.0145e4, 1, None, 343.15],
    ["Chat-Mon", 1, "Propylene-R", "Propylene", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Mon", 2, "Propylene-R", "Propylene", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Mon", 3, "Propylene-R", "Propylene", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Mon", 4, "Propylene-R", "Propylene", 0.012/1000,  5.2e4, 1, None, 343.15],
    ["Chat-Cocat", 1, "Propylene-R", "Triethyl Aluminum", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-Cocat", 2, "Propylene-R", "Triethyl Aluminum", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-Cocat", 3, "Propylene-R", "Triethyl Aluminum", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-Cocat", 4, "Propylene-R", "Triethyl Aluminum", 0.12/1000,  5.02416e4, 1, None, 343.15],
    ["Chat-H2", 1, "Propylene-R", "Hydrogen", 4.8/1000,  4.47988e4, 1, None, 343.15],
    ["Chat-H2", 2, "Propylene-R", "Hydrogen", 8.88/1000,  4.47988e4, 1, None, 343.15],
    ["Chat-H2", 3, "Propylene-R", "Hydrogen", 2.64/1000,  4.47988e4, 1, None, 343.15],
    ["Chat-H2", 4, "Propylene-R", "Hydrogen", 6.6/1000,  4.47988e4, 1, None, 343.15],
    ["Deact-Spon", 1, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Deact-Spon", 2, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Deact-Spon", 3, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Deact-Spon", 4, None, None, 0.001,  4.1868e3, 1, None, 343.15],
    ["Atact-Prop", 1, "Propylene-R", "Propylene", 8/1000,  3.0145e4, 1, None, 343.15],
    ["Atact-Prop", 2, "Propylene-R", "Propylene", 1/1000,  3.0145e4, 1, None, 343.15],
    ["Atact-Prop", 3, "Propylene-R", "Propylene", 0.1/1000,  3.0145e4, 1, None, 343.15],
    ["Atact-Prop", 4, "Propylene-R", "Propylene", 0.1/1000,  3.0145e4, 1, None, 343.15]]

# 设置反应
reaction={"Type":"Ziegler-Nat", 
            "Species": species, 
            "Reactions": r1, 
            "Reacting phase": "liquid"}

# 设置反应器体积和各相体积(或者反应器停留时间、各相停留时间)
holdup = {"Valid phases": "vapor-liquid",
        "Specification type": "Reactor volume & Phase volume", 
        "Reactor volume": 90, 
        "Liquid volume": 60}


# 设置反应器
cstr1.input(spec={"Pressure":3e6 ,"Temperature":333.15}, holdup=holdup, reactions=reaction,
            flash_params={"Max iter":100, "Tol":1e-7})

comp1.input(specs={"Pressure increase": 200000})

# 设置多流股换热器的操作参数
specs = {
    vap_b:{"Outlet": cycle_gas, "Valid phases":"vapor-liquid", "spec":None, 
            "Pressure":0, "Duty estimate":None, "Max iter":100, "Tol":1e-6},

    H2O_in:{"Outlet": H2O_out, "Valid phases":"vapor-liquid", "spec":{"Temperature":303.15}, 
            "Pressure":0, "Duty estimate":None, "Max iter":100, "Tol":1e-6},
}

exc1.input(specs)

# 设置闪蒸罐
flash_specs = {"Temperature": 338.15, "Pressure": 500000}
flash1.input(flash_specs, max_iter=300)

flash_specs = {"Temperature": 333.15, "Pressure": 100000}
flash2.input(flash_specs, max_iter=500)


# 运行流程
# 人为提供撕裂物流和初值, 如果不提供, 则通过内置算法自动确定
fs.run(estimate={cycle_gas: np.array([300, 3000000, 3500, 700, 500, 50, 0, 0, 0, 0])})

# 查看结果
# fs.streams_results()
cycle_gas.print_result()

product.print_result()



