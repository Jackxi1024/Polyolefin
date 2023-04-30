"""
计算焓/熵/热容的公式
1.纯物质焓、熵计算
2.单相混合物焓熵计算
3.多相混合物焓熵计算
4.相态判断(vap, liq, mix)
"""

import numpy as np
import pcsaft
from math import sinh, cosh, log, exp
from scipy import integrate
import json
from polymer_model.thermodynamics.flash_algorithm import *
from polymer_model.utility import InputError

# 将json格式的数据库转化为字典格式
# database_json = open('modified_chemsep_database.json').read()
# database = json.loads(database_json)
# no = {"eqno":{"@value":"no"}}
# for i in database: 
#     type = database[i].get("RPPHeatCapacityCp",no)["eqno"]["@value"]

R = 8.314

# 计算纯组分的饱和蒸气压
def vapor_pressure(T, eqno, params):
    # 蒸气压方程
    if eqno == 101:
        k1, k2, k3, k4, k5, Tmin, Tmax = params
        if T >= Tmin and T <= Tmax:
            return exp(k1 + k2/T + k3*log(T) + k4*T**k5)
    # 安托因方程
    elif eqno == 10:
        k1, k2, k3, Tmin, Tmax = params
        if T >= Tmin and T <= Tmax:
            return exp(k1 - k2/(T+k3))
    


# 汽化焓, 公式106
def enth_vap(T, Tc, params):
    Tr = T/Tc
    if T >= params[5] and T <= params[6]:
        return params[0]*(1-Tr)**(params[1] + params[2]*Tr + params[3]*Tr**2 + params[4]*Tr**3)

# 组分的理想气体摩尔比热容(J/mol/K)
def cp_ig_comp(T, eqno, k):
    # 链段/聚合物的理想气体摩尔比热容
    if eqno == 0:
        k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11 = k
        if T < k7:
            return k9 + k10*pow(T,k11)
        elif T <= k8:
            return k1*pow(T,0) + k2*pow(T,1) + k3*pow(T,2) + k4*pow(T,3) + k5*pow(T,4) + k6*pow(T,5)
        else:
            dCp_dT = k2 + 2*k3*k8 + 3*k4*pow(k8,2) + 4*k5*pow(k8,3) + 5*k6*pow(k8,4)
            return cp_ig_comp(k8,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11) + dCp_dT*(T-k8)
    # 小分子的理想气体摩尔比热容
    elif eqno == 4:
        k1,k2,k3,k4 = k
        return k1 + k2*T + k3*pow(T,2) + k4*pow(T,3)
    elif eqno == 5:
        k1,k2,k3,k4,k5 = k
        return k1 + k2*T + k3*pow(T,2) + k4*pow(T,3) + k5*pow(T,4)
    elif eqno == 107:
        k1,k2,k3,k4,k5,k6,k7 = k
        # Tangent extrapolation
        if T >= k6 and T <= k7:
            return k1 + k2*(k3/T/sinh(k3/T))**2 + k4*(k5/T/cosh(k5/T))**2

# 比热容/温度，用于计算熵
def cp_ig_comp_S(T, eqno, k):
    return cp_ig_comp(T, eqno, k)/T

# 组分的理想气体焓
def enth_ig_comp(T, enth_ig_form, eqno, k, Tref=298.15):
    if T >= Tref:
        delt_H = integrate.quad(cp_ig_comp, Tref, T, args=(eqno, k))[0]
    else:
        delt_H = -integrate.quad(cp_ig_comp, T, Tref, args=(eqno, k))[0]
    return enth_ig_form + delt_H

# 组分的理想熵
def entr_ig_comp(T, P, entr_ig_form, eqno, k, Tref=298.15, Pref=101325):
    if T >= Tref:
        delt_S = integrate.quad(cp_ig_comp_S, Tref, T, args=(eqno, k))[0]
    else:
        delt_S = -integrate.quad(cp_ig_comp_S, T, Tref, args=(eqno, k))[0]
    return entr_ig_form + delt_S

# 混合物的理想焓
def enth_ig_mixture(T, x, enth_ig_form, eqno, k, Tref):
    # 理想焓
    H_ig = np.zeros(len(x))
    for i in range(len(x)):
        H_ig[i] = enth_ig_comp(T, enth_ig_form[i], eqno[i], k[i], Tref[i]) 
    return H_ig*x


# 混合物的理想熵
def entr_ig_mixture(T, P, x, entr_ig_form, eqno, k, Tref, Pref):
    # 理想熵
    S_ig = np.zeros(len(x))
    for i in range(len(x)):
        S_ig[i] = entr_ig_comp(T, P, entr_ig_form[i], eqno[i], k[i], Tref[i], Pref[i]) 
    return S_ig*x


# 单相混合物的焓
def enth_phase(T, P , phase, x, pyargs, enth_ig_form, eqno, k, Tref):
    """ 
    Calculate molar enthalpy of single phase mixture
    
    Parameters
    ----------
    T : float
        Temperature of mixture
    P : float
        pressure of mixture
    phase : str,  "vap" or "liq"
    x : ndarray,  mole fraction of mixture
    pyargs : dict,  physical parameters required by pcsaft package
    enth_ig_form : ndarray,  molar enthalpy of formation of ideal gas of each component
    eqno : ndarray,  Formula number for calculating molar specific heat capacity of ideal gas 
                    for each component
    k : ndarray,  parameters for calculating the molar specific heat capacity of ideal gas 
                    for each component
    Tref : ndarray,  reference temperature default 298.15 K
 
    Returns
    ----------       
    float,  molar enthalpy of single phase mixture
    """
    # 理想焓
    H_ig = np.zeros(len(x))
    for i in range(len(x)):
        H_ig[i] = enth_ig_comp(T, enth_ig_form[i], eqno[i], k[i], Tref[i]) 

    # 剩余焓
    pcsaft_den = pcsaft.pcsaft_den(T, P, x, pyargs, phase)
    H_res = pcsaft.pcsaft_hres(T, pcsaft_den, x, pyargs)
    return H_ig*x + H_res

# 单相混合物的熵
def entr_phase(T, P , phase, x, pyargs, entr_ig_form, eqno, k, Tref, Pref):
    """ 
    Calculate molar entropy of single phase mixture
    
    Parameters
    ----------
    T : float,  temperature of mixture
    P : float,  pressure of mixture
    phase : str,  "vap" or "liq"
    x : ndarray,  mole fraction of mixture
    pyargs : dict,  physical parameters required by pcsaft package
    enth_ig_form : ndarray,  molar enthalpy of formation of ideal gas of each component
    entr_ig_form : ndarray,  molar entropy of formation of ideal gas of each component
    eqno : ndarray,  Formula number for calculating molar specific heat capacity of ideal gas 
                    for each component
    k : ndarray,  parameters for calculating the molar specific heat capacity of ideal gas 
                    for each component
    Tref : float,  reference temperature default 298.15 K
    Pref : float,  reference pressure default 101325 Pa
 
    Returns
    ----------       
    float,  molar entropy of single phase mixture
    """
    # 理想熵
    S_ig = np.zeros(len(x))
    for i in range(len(x)):
        S_ig[i] = entr_ig_comp(T, P, entr_ig_form[i], eqno[i], k[i], Tref[i], Pref[i]) 

    # 剩余熵
    pcsaft_den = pcsaft.pcsaft_den(T, P, x, pyargs, phase)
    S_res = pcsaft.pcsaft_sres(T, pcsaft_den, x, pyargs)
    return S_ig*x + S_res

# 计算单相混合物的吉布斯自由能
def gibbs_phase(T, P, phase, x, pyargs, enth_ig_form, entr_ig_form, eqno, k, Tref, Pref):
    """ 
    Calculate Gibbs free energy of single phase mixture
    
    Parameters
    ----------
    T : float,  temperature of mixture
    P : float,  pressure of mixture
    phase : str,  "vap" or "liq"
    x : ndarray,  mole fraction of mixture
    pyargs : dict,  physical parameters required by pcsaft package
    enth_ig_form : ndarray,  molar enthalpy of formation of ideal gas of each component
    entr_ig_form : ndarray,  molar entropy of formation of ideal gas of each component
    eqno : ndarray,  Formula number for calculating molar specific heat capacity of ideal gas 
                    for each component
    k : ndarray,  parameters for calculating the molar specific heat capacity of ideal gas 
                    for each component
    Tref : float,  reference temperature default 298.15 K
    Pref : float,  reference pressure default 101325 Pa
 
    Returns
    ----------       
    float, molar gibbs free energy of single phase mixture
    """

    H = enth_phase(T, P , phase, x, pyargs, enth_ig_form, eqno, k, Tref)
    S = entr_phase(T, P , phase, x, pyargs, entr_ig_form, eqno, k, Tref, Pref)
    return H - T*S

# 判断流股的相态
def phase(T, P, beta, x, *args):
    if beta == None:
        Pb = bubble_pressure(T, x, args)
        Pd = dew_pressure(T, x, args)
        if P > Pd:
            return "liq"
        elif P < Pd and P > Pb:
            return "mix"
        elif P < Pb:
            return "vap"
        elif P == Pb or P == Pd:
            return "liquid"
    else:
        if beta == 0:
            return "liq"
        elif beta > 0 and beta < 1:
            return "mix"
        elif beta == 1:
            return "vap"
        else:
            raise InputError("Value is out of range;it must be between 0 and 1")

if __name__ == "__main__":
    # 以丙烯为例
    T = 333.1399   # K
    P = 202650     # K
    enth_ig_form = 20230   # J/mol
    entr_ig_form = -142.2438369948013  # 单位: J/mol/K
    eqno = 107
    k = np.array([43.852, 150.600, 1398.8, 74.754, 616.46, 130, 1500])
    H_ig = enth_ig_comp(T, enth_ig_form, eqno, k)
    S_ig = entr_ig_comp(T, P, entr_ig_form, eqno, k)
    print(H_ig)  
    print(S_ig)    

    print(50*"*")
    # 以乙丙共聚为例
    # 乙烯链段
    T = 350  # K
    P = 202650     # K
    enth_ig_form = -44000   # J/mol
    entr_ig_form = 204  # 单位: J/mol/K
    eqno = 0
    k = np.array([-39.748, 0.4, -0.0004998, 0.0000002298, 0, 0, 280, 1000, 36.0292, 1.0708E-57, 23.41042])
    H_ig1 = enth_ig_comp(T, enth_ig_form, eqno, k)
    S_ig1 = entr_ig_comp(T, P, entr_ig_form, eqno, k)-8.314*T*log(0.2)

    # 丙烯链段
    enth_ig_form = -70700   # J/mol
    entr_ig_form = 317  # 单位: J/mol/K
    eqno = 0
    k = np.array([-42.339, 0.50092, -0.0005574, 0.0000002412, 0, 0, 280, 1000, 36.0292, 0.00000161263, 2.927166])
    H_ig2 = enth_ig_comp(T, enth_ig_form, eqno, k)
    S_ig2 = entr_ig_comp(T, P, entr_ig_form, eqno, k)-8.314*T*log(0.8)
    H_ig = 0.2*H_ig1 + 0.8*H_ig2
    S_ig = 0.2*S_ig1 + 0.8*S_ig2
    print(H_ig*1.071416)
    print(S_ig*1.071416)

    T = 298.15
    P = 1013250
    entr_ig_form = -317
    eqno = 0
    k = np.array([-42.339, 0.50092, -0.0005574, 0.0000002412, 0, 0, 280, 1000, 36.0292, 0.00000161263, 2.927166])
    S_ig = entr_ig_comp(T, P, entr_ig_form, eqno, k, Tref=298.15, Pref=101325)
    print("S_ig: ",S_ig)
    import pcsaft
    x = np.array([1.0])
    PP_MWN = 1000*42.0806
    pyargs= {"m":np.array([0.02528*PP_MWN]), "s":np.array([4.1473]), "e":np.array([298.6392])}
    rho_liq = pcsaft.pcsaft_den(T, P, x, pyargs, "liq")
    H_res_liq = pcsaft.pcsaft_hres(T, rho_liq, x, pyargs)
    S_res_liq = pcsaft.pcsaft_sres(T, rho_liq, x, pyargs)
    print("S_res", S_res_liq)

    # TEA
    print(50*"*")
    T = 303.15
    P = 3e6
    enth_ig_form = -163600 
    eqno = 107
    k = np.array([112.060, 384.330, 1533, 242.740, 707.5, 298.15, 1500])
    h_ig = enth_ig_comp(T, enth_ig_form, eqno, k, Tref=298.15)
    print("h_ig: ",h_ig)
    pyargs= {"m":np.array([20]), "s":np.array([3.5356]), "e":np.array([207.19])}
    x = np.array([1.0])
    rho_liq = pcsaft.pcsaft_den(T, P, x, pyargs, "liq")
    print(rho_liq)
    rho_liq = 1607.867
    H_res_liq = pcsaft.pcsaft_hres(T, rho_liq, x, pyargs)
    print("h_res: ",H_res_liq)
    print(h_ig+H_res_liq)

    print("\n氮气饱和蒸气压: ")
    T = 126.2
    eqno = 101
    params = np.array([42.32946, -965.9771, -4.321774, 7.97271e-05, 2.0, 60.81, 126.26])
    print(vapor_pressure(T, eqno, params))

    print("\nH2O饱和蒸气压: ")
    T = 373.15
    eqno = 101
    params = np.array([74.55502, -7295.586, -7.442448, 4.2881e-06, 2.0, 263.15, 647.29])
    print(vapor_pressure(T, eqno, params))