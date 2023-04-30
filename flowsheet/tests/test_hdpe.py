"""
高压液相法生产HDPE, 2个CSTR反应器g
Aspen文件: ./Aspen/hdpe/hdpe.apw
"""

from polymer_model.flowsheet.flowsheet import Flowsheet
from polymer_model.unit_model.stream import MaterialStream
from polymer_model.unit_model.cstr import CSTR
from polymer_model.unit_model.flash import Flash
import time
import numpy as np

# 流程的组分：
components = {
    "Titanium Tetrachloride": {"type": "conventional"},
    "Triethyl Aluminum": {"type": "conventional"},
    "Ethylene": {"type": "conventional"},
    "Hydrogen": {"type": "conventional"},
    "N-hexane": {"type": "conventional"},
    "HDPE": {"type": "polymer"},
    "Ethylene-R": {"type": "segment"},
}

# 流程的聚合物
polymers = {
    # 链段
    "segments": {"Ethylene-R":{"type": "repeat"}}, 
    # 催化剂
    "catalyst": {"Titanium Tetrachloride": {"type": "Z-N",       # 催化剂类型
                        "site_types_num":4,  # 位点数
                        "site_conc":0.1,    # 位点浓度, mol/kgcat
                        }},
    # 聚合物分子量分布
    "Distribution": {"Number of points":100, "Upper limit":100000, "GPC":True}
}

# 物性方法
property_method = "PC-SAFT"

# 创建流程
fs = Flowsheet(components, polymers, property_method)

# 添加单元模块
# 添加反应器1
cstr1 = CSTR()
cstr1.set_name("CSTR1")
fs.add_block(cstr1)

# 添加反应器2
cstr2 = CSTR()
cstr2.set_name("CSTR2")
fs.add_block(cstr2)

# 添加闪蒸罐
flash = Flash()
flash.set_name("Flash")
fs.add_block(flash)

# 添加流股
feed1 = MaterialStream(source=None, destination=cstr1.inlet)
feed1.set_name("Feed1")
fs.add_stream(feed1)

product1 = MaterialStream(source=cstr1.outlet, destination=cstr2.inlet)
product1.set_name("Product1")
fs.add_stream(product1)

feed2 = MaterialStream(source=None, destination=cstr2.inlet)
feed2.set_name("Feed2")
fs.add_stream(feed2)

product2 = MaterialStream(source=cstr2.outlet, destination=flash.inlet)
product2.set_name("Product2")
fs.add_stream(product2)

gas = MaterialStream(source=flash.vapor_outlet, destination=None)
gas.set_name("Gas")
fs.add_stream(gas)

product = MaterialStream(source=flash.liquid_outlet, destination=None)
product.set_name("Polymer")
fs.add_stream(product)


# 设置进料流股参数
# 设置反应器1进料流股
z1 = np.array([0.0001, 0.00015, 0.1, 5E-06, 0.899745, 0])
component_attribute = {
    "Titanium Tetrachloride": {
        "CPSFRAC": 1.,
        "CDSFRAC": 0.,
        "CVSFRAC": np.array([0., 0., 0., 0.,]),
        "CISFRAC": np.array([0., 0., 0., 0.,]),
    },
}
feed1.input(343.15, 20265000, None, 100/6, z1, "Mass & Mass Frac", component_attribute)

# 设置反应器2进料流股
z2 = np.array([0.0001, 0.00015, 0.2, 0.0002, 0.79955, 0])
component_attribute = {
    "Titanium Tetrachloride": {
        "CPSFRAC": 1.,
        "CDSFRAC": 0.,
        "CVSFRAC": np.array([0., 0., 0., 0.,]),
        "CISFRAC": np.array([0., 0., 0., 0.,]),
    },
}
feed2.input(343.15, 20265000, None, 100/36, z2, "Mass & Mass Frac", component_attribute)

# 设置反应
# 设置Ziegler-Natta反应物种
species = { "polymer": "HDPE",
            "tdb segment": None,
            "monomers": ["Ethylene"],
            "segments": {"Ethylene": "Ethylene-R"},
            "precatalyst": None,
            "catalyst": ["Titanium Tetrachloride"],      # 催化剂
            "cocatalysts": ["Triethyl Aluminum"],     # 助催化剂
            "solvents": ["N-hexane"],
            "transfer agent": None,
            "hydrogens": ["Hydrogen"],
            "poisons": None,
            "elec don": None,
            "byproduct": None,
        }


# 设置Z-N反应：反应类型, 催化位点, 组分1, 组分2, 前指因子, 活化能, 反应级数, 终端双键分数, 参考温度
r1 = [["Act-Spon", 1, "Titanium Tetrachloride", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 2, "Titanium Tetrachloride", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 3, "Titanium Tetrachloride", None, 0, 0, 1, None, 1e35],
    ["Act-Spon", 4, "Titanium Tetrachloride", None, 0, 0, 1, None, 1e35],
    ["Act-Cocat", 1, "Titanium Tetrachloride", "Triethyl Aluminum", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 2, "Titanium Tetrachloride", "Triethyl Aluminum", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 3, "Titanium Tetrachloride", "Triethyl Aluminum", 0, 0, 1, None, 1e35],
    ["Act-Cocat", 4, "Titanium Tetrachloride", "Triethyl Aluminum", 0, 0, 1, None, 1e35],
    ["Chain-Ini", 1, "Ethylene", None, 255/1000, 0, 1, None, 1e35],
    ["Chain-Ini", 2, "Ethylene", None, 90/1000,   0, 1, None, 1e35],
    ["Chain-Ini", 3, "Ethylene", None, 0,  0, 1, None, 1e35],
    ["Chain-Ini", 4, "Ethylene", None, 0,  0, 1, None, 1e35],
    ["Propagation", 1, "Ethylene-R", "Ethylene", 255/1000,  0, 1, None, 1e35],
    ["Propagation", 2, "Ethylene-R", "Ethylene", 90/1000,  0, 1, None, 1e35],
    ["Propagation", 3, "Ethylene-R", "Ethylene", 0,  0, 1, None, 1e35],
    ["Propagation", 4, "Ethylene-R", "Ethylene", 0,  0, 1, None, 1e35],
    ["Chat-Mon", 1, "Ethylene-R", "Ethylene", 0.09/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 2, "Ethylene-R", "Ethylene", 0.24/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 3, "Ethylene-R", "Ethylene", 0,  0, 1, None, 1e35],
    ["Chat-Mon", 4, "Ethylene-R", "Ethylene", 0,  0, 1, None, 1e35],
    ["Chat-H2", 1, "Ethylene-R", "Hydrogen", 5.55/1000,  0, 1, None, 1e35],
    ["Chat-H2", 2, "Ethylene-R", "Hydrogen", 18.5/1000,  0, 1, None, 1e35],
    ["Chat-H2", 3, "Ethylene-R", "Hydrogen", 0,  0, 1, None, 1e35],
    ["Chat-H2", 4, "Ethylene-R", "Hydrogen", 0,  0, 1, None, 1e35],
    ["Chat-Spon", 1, "Ethylene-R", None, 0.004,  0, 1, None, 1e35],
    ["Chat-Spon", 2, "Ethylene-R", None, 0.012,  0, 1, None, 1e35],
    ["Chat-Spon", 3, "Ethylene-R", None, 0,  0, 1, None, 1e35],
    ["Chat-Spon", 4, "Ethylene-R", None, 0,  0, 1, None, 1e35],
    ["Deact-Spon", 1, None, None, 0.0001,  0, 1, None, 1e35],
    ["Deact-Spon", 2, None, None, 0.0006,  0, 1, None, 1e35],
    ["Deact-Spon", 3, None, None, 0,  0, 1, None, 1e35],
    ["Deact-Spon", 4, None, None, 0,  0, 1, None, 1e35]]

r2 = [["Act-Spon", 1, "Titanium Tetrachloride", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 2, "Titanium Tetrachloride", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 3, "Titanium Tetrachloride", None, 0.08, 0, 1, None, 1e35],
    ["Act-Spon", 4, "Titanium Tetrachloride", None, 0.08, 0, 1, None, 1e35],
    ["Act-Cocat", 1, "Titanium Tetrachloride", "Triethyl Aluminum", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 2, "Titanium Tetrachloride", "Triethyl Aluminum", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 3, "Titanium Tetrachloride", "Triethyl Aluminum", 0.15/1000, 0, 1, None, 1e35],
    ["Act-Cocat", 4, "Titanium Tetrachloride", "Triethyl Aluminum", 0.15/1000, 0, 1, None, 1e35],
    ["Chain-Ini", 1, "Ethylene", None, 255/1000, 0, 1, None, 1e35],
    ["Chain-Ini", 2, "Ethylene", None, 90/1000,   0, 1, None, 1e35],
    ["Chain-Ini", 3, "Ethylene", None, 255/1000,  0, 1, None, 1e35],
    ["Chain-Ini", 4, "Ethylene", None, 90/1000,  0, 1, None, 1e35],
    ["Propagation", 1, "Ethylene-R", "Ethylene", 255/1000,  0, 1, None, 1e35],
    ["Propagation", 2, "Ethylene-R", "Ethylene", 90/1000,  0, 1, None, 1e35],
    ["Propagation", 3, "Ethylene-R", "Ethylene", 255/1000,  0, 1, None, 1e35],
    ["Propagation", 4, "Ethylene-R", "Ethylene", 90/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 1, "Ethylene-R", "Ethylene", 0.09/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 2, "Ethylene-R", "Ethylene", 0.24/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 3, "Ethylene-R", "Ethylene", 0.09/1000,  0, 1, None, 1e35],
    ["Chat-Mon", 4, "Ethylene-R", "Ethylene", 0.24/1000,  0, 1, None, 1e35],
    ["Chat-H2", 1, "Ethylene-R", "Hydrogen", 5.55/1000,  0, 1, None, 1e35],
    ["Chat-H2", 2, "Ethylene-R", "Hydrogen", 18.5/1000,  0, 1, None, 1e35],
    ["Chat-H2", 3, "Ethylene-R", "Hydrogen", 5.55/1000,  0, 1, None, 1e35],
    ["Chat-H2", 4, "Ethylene-R", "Hydrogen", 18.5/1000,  0, 1, None, 1e35],
    ["Chat-Spon", 1, "Ethylene-R", None, 0.004,  0, 1, None, 1e35],
    ["Chat-Spon", 2, "Ethylene-R", None, 0.012,  0, 1, None, 1e35],
    ["Chat-Spon", 3, "Ethylene-R", None, 0.004,  0, 1, None, 1e35],
    ["Chat-Spon", 4, "Ethylene-R", None, 0.012,  0, 1, None, 1e35],
    ["Deact-Spon", 1, None, None, 0.0001,  0, 1, None, 1e35],
    ["Deact-Spon", 2, None, None, 0.0006,  0, 1, None, 1e35],
    ["Deact-Spon", 3, None, None, 0.0001,  0, 1, None, 1e35],
    ["Deact-Spon", 4, None, None, 0.0006,  0, 1, None, 1e35]]

# 设置反应器1
# 设置反应
reaction={"Type":"Ziegler-Nat", 
            "Species": species, 
            "Reactions": r1, 
            "Reacting phase": "liquid"}

# 设置反应器体积和各相体积(或者反应器停留时间、各相停留时间)
holdup = {"Valid phases": "liquid",
        "Specification type": "Reactor volume", 
        "Reactor volume": 60}

cstr1.input(spec={"Pressure":20265000 ,"Temperature":433.15}, holdup=holdup, reactions=reaction,
            flash_params={"Max iter":100, "Tol":1e-7})

# 设置反应器2
# 设置反应
reaction={"Type":"Ziegler-Nat", 
            "Species": species, 
            "Reactions": r2, 
            "Reacting phase": "liquid"}

cstr2.input(spec={"Pressure":20265000 ,"Temperature":433.15}, holdup=holdup, reactions=reaction,
            flash_params={"Max iter":100, "Tol":1e-7})

# 设置闪蒸罐
flash_specs = {"Temperature": 433.15, "Pressure": 1013250}
flash.input(flash_specs, max_iter=500)


# 运行
start = time.time()
fs.run()
end = time.time()

# 输出各反应器结果
cstr1.print_results()
cstr2.print_results()
flash.print_results()

# 输出各流股结果
cstr1.print_stream_results()
cstr2.print_stream_results()
flash.print_stream_results()

# 输出各反应器生成的聚合物
cstr1.print_local_cld()
cstr2.print_local_cld()

# 输出流股的聚合物的分布信息
product1.print_cld()
product1.plot_cld()

product2.print_cld()
product2.plot_cld()

product.print_cld()
product.plot_cld()


print("运行时间："+str(end - start)+"秒")