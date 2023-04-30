"""
压缩机: 有1个入口和一个出口, 入口可以有多条进料, 出口只有一条出料
混合器的run方法: 
(1) 通过物料衡算, 计算出料流量和组成
(2) 出口压力为设定值
(3) 通过压缩规定(等熵压缩/多变压缩), 计算出口流股温度、压缩机做功等
"""


from polymer_model.unit_model.block import Block, Port

class Compressor(Block):
    """
    Compressor

    Parameters
    ----------
    type : str, default="Isentropic" 
            The compression type of compressor:  "Isentropic" or "Polytropic"
    specs: dict, 
            Operating conditions of compressor.
            Optional keywords : "Discharge pressure"、"Pressure increase"、"Pressures ratio"
    efficiencies : dict, default={"Isentropic": 0.72, "Polytropic":0.72, "Mechanical":1}
            Efficiencies of compressor
            When Type is "Isentropic",  isentropic efficiency and mechanical efficiency are required.
            When Type is "Polytropic",  polytropic efficiency and mechanical efficiency are required.
    flash_params : dict, default={"Valid phases": "vapor", "Max iter": 100, "Tol":1e-6}
            Parameters required for flash calculation.
            "Valid phases"、"Max iter" and "Tol" are required.
    entropy_balance_params : dict, default={"Max iter": 100, "Tol":1e-6}
            Parameters required for entropy balance calculation.
            "Max iter" and "Tol" are required.

    Attributes
    ----------
    type : str, default="Isentropic" 
            The compression type of compressor:  "Isentropic" or "Polytropic"
    """


    def __init__(self):
        super().__init__()
        self.inlet = Port(self)    # 压缩机的入口对象
        self.outlet = Port(self)   # 压缩机的出口对象


    def input(self, type="Isentropic", specs={}, efficiencies={"Isentropic": 0.72, "Polytropic":0.72, "Mechanical":1},
            flash_params={"Valid phases": "vapor", "Max iter": 100, "Tol":1e-6},
            entropy_balance_params={"Max iter": 100, "Tol":1e-6}):
        """ 
        输入 
        """

        self.type = type
        self.specs = specs
        self.efficiencies = efficiencies
        self.flash_params = flash_params
        self.valid_phases = flash_params["Valid phases"]
        self.max_iter = flash_params["Max iter"]
        self.tol = flash_params["Tol"]
        self.entropy_balance_params = entropy_balance_params
    

    def isentropic_compress(self):
        """  等熵压缩 """
        # 计算入口物流的摩尔焓
        feed = self.mix(self.inlet.streams)
        T_in = feed.temperature    
        H_in = feed.enthalpy_flow
        S_in = feed.entropy_flow

        P_in = feed.pressure   # 入口压力
        if self.specs.get("Discharge pressure") != None:
            P_out = self.specs.get("Discharge pressure")
        elif self.specs.get("Pressure increase") != None:
            P_out = P_in + self.specs.get("Pressure increase")
        elif self.specs.get("Pressure ratio") != None:
            P_out = P_in*self.specs.get("Pressure ratio")
        else:
            raise NotImplementedError("The specification has not been implemented!")

        # 初始化出口物流
        product = self.outlet.streams[0]  # 出料流股
        product.distributions = feed.distributions

        # 迭代, 更新等熵出口温度的上下界
        T_min = T_in
        T_max = T_min*1.1
        while True:
            product.input(T_max, P_out, None, None, feed.component_mole_flow, "Mole & Mole Flow", 
                        feed.component_attributes, self.valid_phases, self.max_iter, self.tol)
            product.run()
            if product.entropy_flow < S_in:
                T_min = T_max
                T_max = 1.1*T_max
            else:
                break
        
        # 二分法迭代, 计算等熵出口温度
        s_tol = self.entropy_balance_params["Tol"]
        error = 2*s_tol
        while True:
            T_isen = (T_min + T_max)/2   # 等熵出口温度
            product.input(T_isen, P_out, None, None, feed.component_mole_flow, "Mole & Mole Flow", 
                        feed.component_attributes, self.valid_phases, self.max_iter, self.tol)
            product.run()            
            
            S_isen = product.entropy_flow
            error = S_isen - S_in

            if abs(error/S_in) <= s_tol:
                self.isentropic_outlet_temperature = T_isen
                break
            if error > 0:
                T_max = T_isen
            else:
                T_min = T_isen
        

        # 计算等熵压缩功(Aspen用Isentropic power requirement表示)
        self.isentropic_compress_work = product.enthalpy_flow - H_in
        # 计算实际压缩功(Aspen用Indicated horsepower表示) 
        self.isentropic_efficiency = self.efficiencies["Isentropic"]
        self.efficiency = self.isentropic_efficiency
        self.actual_compress_work = self.isentropic_compress_work/self.isentropic_efficiency
        # 计算压缩机所需机械功, aspen用Brake horsepower表示
        self.mechanical_efficiency = self.efficiencies["Mechanical"]
        self.mechanical_work = self.actual_compress_work/self.mechanical_efficiency
        # 计算压缩机所需净功率, 当没有其他功流输入时, 其等同于Brake horsepower
        self.net_work = self.mechanical_work
        
        # 计算出口焓值
        H_out = feed.enthalpy_flow + self.actual_compress_work

        # 根据出口焓值, 迭代计算真实的出口温度
        # 迭代, 更新出口温度的上下界
        T_min = T_isen
        T_max = T_min+1
        while True:
            product.input(T_max, P_out, None, None, feed.component_mole_flow, "Mole & Mole Flow", 
                        feed.component_attributes, self.valid_phases, self.max_iter, self.tol)
            product.run()     
            if product.enthalpy_flow < H_out:
                T_min = T_max
                T_max = T_min+1
            else:
                break

        # 二分法迭代, 计算出口温度
        h_tol = s_tol
        while True:
            T_out = (T_min + T_max)/2   # 等熵出口温度
            product.input(T_out, P_out, None, None, feed.component_mole_flow, "Mole & Mole Flow", 
                        feed.component_attributes, self.valid_phases, self.max_iter, self.tol)
            product.run() 
            error = product.enthalpy_flow - H_out
            if abs(error/H_out) <= h_tol:
                break
            if error > 0:
                T_max = T_out
            else:
                T_min = T_out

    
    def run(self):
        """ 
        运行 
        """

        if self.type == "Isentropic":
            self.isentropic_compress()
        else:
            raise NotImplementedError("The type has not been implemented!")
        
        self.status = True


    def print_results(self):
        """ 
        Print results of the mixer.
        """

        print("Compressor model:\t", self.type)
        print("Isentropic power:\t", self.isentropic_compress_work)
        print("Indicated horsepower:\t", self.actual_compress_work)
        print("Brake horsepower:\t", self.mechanical_work)
        print("Net work required:\t", self.net_work)
        print("Efficiency:\t", self.efficiency)
        print("Mechanical efficiency:\t", self.mechanical_efficiency)   
        print("Isentropic outlet temperature:\t", self.isentropic_outlet_temperature)
        
        product = self.outlet.streams[0]  # 出料流股
        print("Outlet temperature: ", product.temperature)
        print("Outlet pressure: ", product.pressure)
        print("Vapor fraction: ", product.vapor_fraction)



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
    comp1 = Compressor()
    comp1.set_name("Comp1")
    fs.add_block(comp1)

    # 添加流股
    vap_a = MaterialStream(source=None, destination=comp1.inlet)
    vap_a.set_name("Vap-A")
    fs.add_stream(vap_a)

    vap_b = MaterialStream(source=comp1.outlet, destination=None)
    vap_b.set_name("Vap-B")
    fs.add_stream(vap_b)

    # 设置混合器进料流股参数
    # C3H6、C3H8、H2、N2、PP、TiCl4、TEA、H2O
    z1 = np.array([0.7460931, 0.1361802, 0.1080722, 0.00965455, 0, 0, 0, 0])
    vap_a.input(333.15, 3000000, None, 4991.95, z1, "Mole & Mole Frac", valid_phases="vapor")

    # 设置各单元模块操作参数
    comp1.input(specs={"Pressure increase": 200000})

    # 运行
    start = time.time()
    comp1.run()

    end = time.time()


    # 输出结果
    comp1.print_results()
    print("\n")
    comp1.print_stream_results()

    print("运行时间："+str(end - start)+"秒")


