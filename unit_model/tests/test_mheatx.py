# 导包
from polymer_model.flowsheet.flowsheet import Flowsheet
from polymer_model.unit_model.stream import MaterialStream
from polymer_model.unit_model.mheatx import MHeatX
import time
import numpy as np

# 流程的组分：C3H6、C3H8、H2O
components = {
    "Propylene": {"type": "conventional"},
    "Propane": {"type": "conventional"},
    "Water": {"type": "conventional"},
}

# 物性方法
property_method = "PC-SAFT"

# 创建流程
fs = Flowsheet(components, None, property_method)

# 添加单元模块
# 添加混合器
exc = MHeatX()
fs.add_block(exc)

# 添加流股
s1 = MaterialStream(source=None, destination=exc.hot_inlet)
fs.add_stream(s1)
s2 = MaterialStream(source=None, destination=exc.hot_inlet)
fs.add_stream(s2)
s3 = MaterialStream(source=exc.hot_outlet, destination=None)
fs.add_stream(s3)
s4 = MaterialStream(source=exc.hot_outlet, destination=None)
fs.add_stream(s4)
s5 = MaterialStream(source=None, destination=exc.cold_inlet)
fs.add_stream(s5)
s6 = MaterialStream(source=None, destination=exc.cold_inlet)
fs.add_stream(s6)
s7 = MaterialStream(source=exc.cold_outlet, destination=None)
fs.add_stream(s7)
s8 = MaterialStream(source=exc.cold_outlet, destination=None)
fs.add_stream(s8)

# 设置混合器进料流股参数
z1 = np.array([0.8, 0.2, 0])
s1.input(373.15, 202650, None, 1/3.6, z1, "Mass & Mass Frac")

z2 = np.array([0.5, 0.5, 0])
s2.input(363.15, 506625, None, 2/3.6, z2, "Mass & Mass Frac")

z5 = np.array([0, 0, 1])
s5.input(293.15, 101325, None, 1/3.6, z5, "Mass & Mass Frac")

z6 = np.array([0, 0, 1])
s6.input(298.15, 101325, None, 2/3.6, z6, "Mass & Mass Frac")

# 设置多流股换热器的操作参数
specs = {
    s1:{"Outlet": s3, "Valid phases":"vapor-liquid", "spec":{"Temperature":333.15}, 
            "Pressure":0, "Duty estimate":None, "Max iter":100, "Tol":1e-6},

    s2:{"Outlet": s4, "Valid phases":"vapor-liquid", "spec":{"Temperature":333.15}, 
            "Pressure":0, "Duty estimate":None, "Max iter":100, "Tol":1e-6},

    s5:{"Outlet": s7, "Valid phases":"vapor-liquid", "spec":None, 
            "Pressure":0, "Duty estimate":None, "Max iter":100, "Tol":1e-6},

    s6:{"Outlet": s8, "Valid phases":"vapor-liquid", "spec":None, 
            "Pressure":0, "Duty estimate":None, "Max iter":100, "Tol":1e-6},
}

exc.input(specs)

# 运行
start = time.time()
exc.run()
end = time.time()

# 输出结果
exc.print_results()

# 输出结果
exc.print_stream_results()

print("运行时间："+str(end - start)+"秒")