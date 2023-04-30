"""
PFR反应器: 

PFR反应器的run方法: 
    
"""

import numpy as np
from polymer_model.unit_model.block import Block, Port


class PFR(Block):
    """
    PFR

    Attributes
    ----------

    """

    def __init__(self):
        super().__init__()
        self.inlet = Port(self)      # 混合器的入口对象
        self.outlet = Port(self)     # 混合器的出口对象
        self.vapor_outlet = Port(self)    # 混合器的气相出口对象
        self.liquid_outlet = Port(self)   # 混合器的液相出口对象


    def input(self, pressure=0, valid_phases="vapor-liquid", T_guess=None, max_iter=30, tol=0.0001):
        """ 
        输入
        """
        self.pressure = pressure
        self.valid_phases = valid_phases
        self.T_guess = T_guess
        self.max_iter = max_iter
        self.tol = tol
    
    def run(self):
        """ 
        运行
        """
        pass

    
    def print_results(self):
        """ 
        Print results of the mixer.
        """

        print("Outlet temperature: ", self.outlet_temperature)
        print("Outlet pressure: ", self.outlet_pressure)
        print("Vapor fraction: ", self.vapor_fraction)




if __name__ == "__main__":
    from polymer_model.flowsheet.flowsheet import Flowsheet
    from polymer_model.unit_model.stream import MaterialStream