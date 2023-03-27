"""
物性方法封装工具: 将各个状态方程/活度方程封装在一起
"""

import pcsaft  

class PropertyMethod:
    """ 
    物性方法封装类, 其他热力学函数在调用此类方法时, 会根据实际的物性方法返回计算结果 
    """

    def __init__(self, property_method, property_args) -> None:
        self.property_method = property_method
        self.property_args = property_args
    
    def molar_density(self, temperature, pressure, mole_frac, phase):
        if self.property_method == "PC-SAFT":
            if phase in ("l", "liquid", "Liquid"):
                phase = "liq"
            elif phase in ("v", "vapor", "Vapor"):
                phase = "vap"
            return pcsaft.pcsaft_den(temperature, pressure, mole_frac, 
                                    self.property_args, phase)
        # To be completed
        elif self.property_method == "IDEAL":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        # To be completed
        elif self.property_method == "PENG-ROB":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        # To be completed
        elif self.property_method == "SRK":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        else:
            raise NotImplementedError(self.property_method+" has not been implemented yet!")

    def fugacity_coefficients(self, temperature, pressure, mole_frac, phase):
        if self.property_method == "PC-SAFT":
            if phase in ("l", "liquid", "Liquid"):
                phase = "liq"
            elif phase in ("v", "vapor", "Vapor"):
                phase = "vap"
            rho = pcsaft.pcsaft_den(temperature, pressure, mole_frac, 
                                    self.property_args, phase)
            fugcoef = pcsaft.pcsaft_fugcoef(temperature, rho, mole_frac, self.property_args)
            return fugcoef
        # To be completed
        elif self.property_method == "IDEAL":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        # To be completed
        elif self.property_method == "PENG-ROB":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        # To be completed
        elif self.property_method == "SRK":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        else:
            raise NotImplementedError(self.property_method+" has not been implemented yet!")        

    def molar_residual_enthalpy(self, temperature, pressure, mole_frac, phase):
        if self.property_method == "PC-SAFT":
            if phase in ("l", "liquid", "Liquid"):
                phase = "liq"
            elif phase in ("v", "vapor", "Vapor"):
                phase = "vap"
            rho = pcsaft.pcsaft_den(temperature, pressure, mole_frac, 
                                    self.property_args, phase)
            hres = pcsaft.pcsaft_hres(temperature, rho, mole_frac, self.property_args)
            return hres
        # To be completed
        elif self.property_method == "IDEAL":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        # To be completed
        elif self.property_method == "PENG-ROB":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        # To be completed
        elif self.property_method == "SRK":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        else:
            raise NotImplementedError(self.property_method+" has not been implemented yet!")

    def molar_residual_entropy(self, temperature, pressure, mole_frac, phase):
        if self.property_method == "PC-SAFT":
            if phase in ("l", "liquid", "Liquid"):
                phase = "liq"
            elif phase in ("v", "vapor", "Vapor"):
                phase = "vap"
            rho = pcsaft.pcsaft_den(temperature, pressure, mole_frac, 
                                    self.property_args, phase)
            sres = pcsaft.pcsaft_sres(temperature, rho, mole_frac, self.property_args)
            return sres
        # To be completed
        elif self.property_method == "IDEAL":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        # To be completed
        elif self.property_method == "PENG-ROB":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        # To be completed
        elif self.property_method == "SRK":
            raise NotImplementedError(self.property_method+" has not been implemented yet!")
        else:
            raise NotImplementedError(self.property_method+" has not been implemented yet!")


# test
if __name__ == "__main__":
    import numpy as np

    property_method = "IDEAL"   # 物性方法
    m = np.array([0.8285,])
    s = np.array([2.9729, ])
    e = np.array([12.53, ])
    property_args = {"m":m, "s":s, "e":e}  # 物性参数

    eos = PropertyMethod(property_method, property_args)

    T = 303.15
    P = 3000000
    z = np.array([1])
    
    print("Molar Densidty: ", )
    print(eos.molar_density(T, P, z, "liquid"))
    print("fugacity_coefficients: ", )
    print(eos.fugacity_coefficients(T, P, z, "liquid"))
    print("Molar Enthalpy: ", )
    print(eos.molar_residual_enthalpy(T, P, z, "liquid"))
    print("Molar Entropy: ", )
    print(eos.molar_residual_entropy(T, P, z, "liquid"))
    
