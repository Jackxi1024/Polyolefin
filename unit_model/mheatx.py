"""
多流股换热器: 
    有1个冷流股入口、1个冷流股出口、1个热流股入口、1个热流股出口
    每个入口和出口均可有多条流股

多流股换热器的run方法: 
(1) 指定某些流股的入口和出口参数, 通过能量平衡, 计算未指定的流股出料温度

目前还未实现冷热流股的温度匹配(即Zone Anysis)
"""

from polymer_model.unit_model.block import Block, Port
from polymer_model.utility import SolutionError      # 异常类
import prettytable as pt  # 表格工具, 在终端以表格形式输出结果

class MHeatX(Block):
    """
    Multi stream heat exchanger

    Attributes
    ----------

    """

    def __init__(self):
        super().__init__()
        self.hot_inlet = Port(self)    # 多流股换热器的热流股入口对象
        self.hot_outlet = Port(self)   # 多流股换热器的热流股出口对象
        self.cold_inlet = Port(self)   # 多流股换热器的冷流股入口对象
        self.cold_outlet = Port(self)  # 多流股换热器的冷流股出口对象


    def input(self, specs={}, convergence_params={"Max iter":100, "Tol":1e-7}):
        """
        输入


        specs : dict.
            Specifications of outlet cold and hot stream.
            For example:

            exc = MHeatX()
            gas_in = MaterialStream(source=None, destination=exc.hot_inlet)
            gas_out = MaterialStream(source=None, destination=exc.hot_outlet)
            H2O_in = MaterialStream(source=None, destination=exc.cold_inlet)
            H2O_out = MaterialStream(source=None, destination=exc.cold_inlet)

            specs = {
                gas_in:{"Outlet": gas_out, "Valid phases":"vapor-liquid", "spec":{"Temperature":30}, 
                        "Pressure":0, "Duty estimate":None, "Max iter":100, "Tol":1e-6},
                H2O_in:{"Outlet": H2O_out, "Valid phases":"vapor-liquid", "spec":None, "Pressure":0, 
                        "Duty estimate":None, "Max iter":100, "Tol":1e-6},
            }

        """
        self.specs = specs    # 多流股换热器的操作条件
        self.convergence_params = convergence_params   # 计算出口温度的收敛参数

    
    def run(self):
        """
        运行
        """

        h_tol = self.convergence_params["Tol"]   # 计算焓平衡的容差

        # 每条冷热流股的出口和传热量 
        hot_duty = {}
        cold_duty = {}
        undetermined = []  # 未指定的流股
        total_duty = 0  # 换热器总负荷

        # 判断是否所有的热流股都指定了出口
        flag = True
        for stream in self.hot_inlet.streams:
            if self.specs[stream].get("spec") == None:
                flag = False

        # 如果所有热流股都指定了出口
        if flag:
            T_max = 0   # 入口热流股温度的最大值，也是出口温度的上界
            # 计算每条热流股的换热量
            for stream in self.hot_inlet.streams:
                # 获取每条热流股的入口压力、焓值
                if stream.status == False:
                    stream.run()
                T_in = stream.temperature
                P_in = stream.pressure
                H_in = stream.enthalpy_flow
                if T_in > T_max:
                    T_max = T_in
                
                # 获取每条热流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                outlet_stream = self.specs[stream]["Outlet"]   # 对应的出口流股
                P_out = self.specs[stream]["Pressure"]      # 出口压力
                if P_out <= 0:
                    P_out = P_in - P_out
                valid_phases = self.specs[stream]["Valid phases"]  # 出口有效相
                max_iter = self.specs[stream]["Max iter"]
                tol = self.specs[stream]["Tol"]
                if self.specs[stream]["spec"].get("Temperature") != None:
                    T_out = self.specs[stream]["spec"].get("Temperature")
                    outlet_stream.input(T_out, P_out, None, None, stream.component_mole_flow, "Mole & Mole Flow", 
                                    stream.component_attributes, valid_phases, max_iter, tol)
                    outlet_stream.run()
                    hot_duty[stream] = outlet_stream.enthalpy_flow - H_in  # 热流股换热量
            
            T_min = 0  # 未指定的冷流股入口温度最大值，也是其出口温度的下界
            # 获取冷流股的出口规范, 并计算冷流股的吸热量
            for stream in self.cold_inlet.streams:
                # 获取每条冷流股的入口压力、焓值
                if stream.status == False:
                    stream.run()
                T_in = stream.temperature
                P_in = stream.pressure
                H_in = stream.enthalpy_flow
                # 指定了流股的出口规范
                if self.specs[stream].get("spec") != None:
                    # 获取每条冷流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                    outlet_stream = self.specs[stream]["Outlet"]   # 对应的出口流股
                    P_out = self.specs[stream]["Pressure"]      # 出口压力
                    if P_out <= 0:
                        P_out = P_in - P_out
                    valid_phases = self.specs[stream]["Valid phases"]  # 出口有效相
                    max_iter = self.specs[stream]["Max iter"]
                    tol = self.specs[stream]["Tol"]
                    # 如果指定了出口温度
                    if self.specs[stream]["spec"].get("Temperature") != None:
                        T_out = self.specs[stream]["spec"].get("Temperature")
                        outlet_stream.input(T_out, P_out, None, None, stream.component_mole_flow, "Mole & Mole Flow", 
                                        stream.component_attributes, valid_phases, max_iter, tol)
                        outlet_stream.run()
                        cold_duty[stream] = outlet_stream.enthalpy_flow - H_in  # 热流股换热量
                # 未指定了流股的出口规范
                else:
                    if T_in > T_min:
                        T_min = T_in
                    # 第i条流股未指定出口
                    undetermined.append(stream)  
            # 总负荷
            total_duty = -sum(hot_duty.values())
            # 剩余负荷
            remaining_duty = total_duty - sum(cold_duty.values())

            # 二分法迭代, 计算未指定流股的出口温度
            while True:
                T_out = (T_min+T_max)/2
                calc_duty = 0
                for stream in undetermined:
                    # 获取每条冷流股的入口压力、焓值
                    if stream.status == False:
                        stream.run()
                    P_in = stream.pressure
                    H_in = stream.enthalpy_flow
                    
                    # 获取每条冷流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                    outlet_stream = self.specs[stream]["Outlet"]   # 对应的出口流股
                    P_out = self.specs[stream]["Pressure"]      # 出口压力
                    if P_out <= 0:
                        P_out = P_in - P_out
                    valid_phases = self.specs[stream]["Valid phases"]  # 出口有效相
                    max_iter = self.specs[stream]["Max iter"]
                    tol = self.specs[stream]["Tol"]

                    outlet_stream.input(T_out, P_out, None, None, stream.component_mole_flow, "Mole & Mole Flow", 
                                    stream.component_attributes, valid_phases, max_iter, tol)
                    outlet_stream.run()
                    cold_duty[stream] = outlet_stream.enthalpy_flow - H_in  # 冷流股换热量
                    calc_duty += cold_duty[stream]  
                    
                error = calc_duty - remaining_duty
                if abs(error/remaining_duty) <= h_tol:
                    break
                if error > 0:
                    T_max = T_out
                else:
                    T_min = T_out
        
        # 如果有热流股未指定出口，那么所有的冷流股应该指定了出口
        else:
            T_min = 2000    # 入口冷流股温度的最小值，也是出口温度的下界
            # 计算每条冷流股的换热量
            for stream in self.cold_inlet.streams:
                # 获取每条冷流股的入口压力、焓值
                if stream.status == False:
                    stream.run()
                T_in = stream.temperature
                P_in = stream.pressure
                H_in = stream.enthalpy_flow
                if T_in < T_min:
                    T_min = T_in
                
                # 获取每条冷流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                outlet_stream = self.specs[stream]["Outlet"]   # 对应的出口流股
                P_out = self.specs[stream]["Pressure"]      # 出口压力
                if P_out <= 0:
                    P_out = P_in - P_out
                valid_phases = self.specs[stream]["Valid phases"]  # 出口有效相
                max_iter = self.specs[stream]["Max iter"]
                tol = self.specs[stream]["Tol"]
                if self.specs[stream]["spec"].get("Temperature") != None:
                    T_out = self.specs[stream]["spec"].get("Temperature")
                    outlet_stream.input(T_out, P_out, None, None, stream.component_mole_flow, "Mole & Mole Flow", 
                                    stream.component_attributes, valid_phases, max_iter, tol)
                    outlet_stream.run()
                    cold_duty[stream] = outlet_stream.enthalpy_flow - H_in  # 热流股换热量
          
            T_max = 2000  # 未指定的热流股入口温度最小值，也是其出口温度的上界
            # 获取热流股的出口规范, 并计算冷流股的吸热量
            for stream in self.hot_inlet.streams:
                # 获取每条热流股的入口压力、焓值
                if stream.status == False:
                    stream.run()
                T_in = stream.temperature
                P_in = stream.pressure
                H_in = stream.enthalpy_flow
                # 指定了流股的出口规范
                if self.specs[stream].get("spec") != None:
                    # 获取每条热流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                    outlet_stream = self.specs[stream]["Outlet"]   # 对应的出口流股
                    P_out = self.specs[stream]["Pressure"]      # 出口压力
                    if P_out <= 0:
                        P_out = P_in - P_out
                    valid_phases = self.specs[stream]["Valid phases"]  # 出口有效相
                    max_iter = self.specs[stream]["Max iter"]
                    tol = self.specs[stream]["Tol"]
                    # 如果指定了出口温度
                    if self.specs[stream]["spec"].get("Temperature") != None:
                        T_out = self.specs[stream]["spec"].get("Temperature")
                        outlet_stream.input(T_out, P_out, None, None, stream.component_mole_flow, "Mole & Mole Flow", 
                                        stream.component_attributes, valid_phases, max_iter, tol)
                        outlet_stream.run()
                        hot_duty[stream] = outlet_stream.enthalpy_flow - H_in  # 热流股换热量
                # 未指定了流股的出口规范
                else:
                    if T_in < T_max:
                        T_max = T_in
                    # 第i条流股未指定出口
                    undetermined.append(stream)  
            # 总负荷
            total_duty = -sum(cold_duty.values())
            # 剩余负荷
            remaining_duty = total_duty - sum(hot_duty.values())

            # 二分法迭代, 计算未指定流股的出口温度
            while True:
                T_out = (T_min+T_max)/2
                calc_duty = 0
                for stream in undetermined:
                    # 获取每条流股的入口压力、焓值
                    if stream.status == False:
                        stream.run()
                    P_in = stream.pressure
                    H_in = stream.enthalpy_flow
                    
                    # 获取每条流股的出口规范、压力、进料条件、有效相、最大迭代数、容差
                    outlet_stream = self.specs[stream]["Outlet"]   # 对应的出口流股
                    P_out = self.specs[stream]["Pressure"]      # 出口压力
                    if P_out <= 0:
                        P_out = P_in - P_out
                    valid_phases = self.specs[stream]["Valid phases"]  # 出口有效相
                    max_iter = self.specs[stream]["Max iter"]
                    tol = self.specs[stream]["Tol"]

                    outlet_stream.input(T_out, P_out, None, None, stream.component_mole_flow, "Mole & Mole Flow", 
                                    stream.component_attributes, valid_phases, max_iter, tol)
                    outlet_stream.run()
                    hot_duty[stream] = outlet_stream.enthalpy_flow - H_in  # 冷流股换热量
                    calc_duty += hot_duty[stream]  
                    
                error = calc_duty - remaining_duty
                if abs(error/remaining_duty) <= h_tol:
                    break
                if error > 0:
                    T_max = T_out
                else:
                    T_min = T_out

        self.hot_duty = hot_duty
        self.cold_duty = cold_duty
        self.duty = total_duty

        self.status = True


    def print_results(self):
        """ 
        Print results of the multi stream heat exchanger.
        """

        if self.status == False:
            raise SolutionError("计算失败！")
        
        print("Duty: ", self.duty)   
        
        # 以表格的形式输出
        tb = pt.PrettyTable()

        # 创建列标题(参考Aspen)
        column_header = ["Exchanger side", "Outlet stream", 
                        "Inlet temperature", "Inlet pressure", "Inlet vapor fraction", 
                        "Outlet temperature", "Outlet pressure", "Outlet vapor fraction", "duty"]
        # 设置表格的列标题
        tb.add_column("Inlet stream", column_header)
        
        # 添加流股的数据
        for stream in self.cold_inlet.streams:
            outlet_stream = self.specs[stream]["Outlet"]
            column_data = ["cold", outlet_stream.name,
                    stream.temperature, stream.pressure, stream.vapor_fraction,
                    outlet_stream.temperature, outlet_stream.pressure, outlet_stream.vapor_fraction,
                    self.cold_duty[stream]]

            # 往表格中添加数据
            tb.add_column(stream.name, column_data)

        # 添加流股的数据
        for stream in self.hot_inlet.streams:
            outlet_stream = self.specs[stream]["Outlet"]
            column_data = ["hot", outlet_stream.name,
                    stream.temperature, stream.pressure, stream.vapor_fraction,
                    outlet_stream.temperature, outlet_stream.pressure, outlet_stream.vapor_fraction,
                    self.hot_duty[stream]]

            # 往表格中添加数据
            tb.add_column(stream.name, column_data)

        # 输出表格
        print("Streams Results: ")
        print(tb)





if __name__ == "__main__":
    from polymer_model.flowsheet.flowsheet import Flowsheet
    from polymer_model.unit_model.stream import MaterialStream
    import time
    import numpy as np

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

    # 物性方法
    property_method = "PC-SAFT"

    # 创建流程
    fs = Flowsheet(components, None, property_method)

    # 添加单元模块
    # 添加混合器
    exc = MHeatX()
    exc.set_name("Exc")
    fs.add_block(exc)

    # 添加流股
    gas_in = MaterialStream(source=None, destination=exc.hot_inlet)
    gas_in.set_name("Gas-In")
    fs.add_stream(gas_in)

    gas_out = MaterialStream(source=exc.hot_outlet, destination=None)
    gas_out.set_name("Gas-Out")
    fs.add_stream(gas_out)

    H2O_in = MaterialStream(source=None, destination=exc.cold_inlet)
    H2O_in.set_name("H2O_in")
    fs.add_stream(H2O_in)

    H2O_out = MaterialStream(source=exc.cold_outlet, destination=None)
    H2O_out.set_name("H2O_out")
    fs.add_stream(H2O_out)

    # 设置混合器进料流股参数
    # C3H6、C3H8、H2、N2、PP、TiCl4、TEA、H2O
    z1 = np.array([0.7460931, 0.1361802, 0.1080722, 0.00965455, 0, 0, 0, 0])
    gas_in.input(337.593, 3200000, None, 4991.95, z1, "Mole & Mole Frac", valid_phases="vapor")

    z2 = np.array([0, 0, 0, 0, 0, 0, 0, 1])
    H2O_in.input(293.15, 3000000, None, 235.4167, z2, "Mass & Mass Frac")

    # 设置多流股换热器的操作参数
    specs = {
        gas_in:{"Outlet": gas_out, "Valid phases":"vapor-liquid", "spec":None, 
                "Pressure":0, "Duty estimate":None, "Max iter":100, "Tol":1e-6},

        H2O_in:{"Outlet": H2O_out, "Valid phases":"vapor-liquid", "spec":{"Temperature":303.15}, 
                "Pressure":0, "Duty estimate":None, "Max iter":100, "Tol":1e-6},
    }

    exc.input(specs)

    # 运行
    start = time.time()
    exc.run()
    end = time.time()

    # 输出结果
    exc.print_results()
    print("\n")
    exc.print_stream_results()

    print("运行时间："+str(end - start)+"秒")


