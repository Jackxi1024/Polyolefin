"""
The flash algorithm includes: 
    Calculation of bubble point and dew point、stability analysis、
    TP flash、

Reference: 
Michelsen & Mollerup - Thermodynamic models:Fundamentals and Computational Aspects

"""

import numpy as np
import math
import pcsaft
import sys
sys.path.append('E:\\AZFT\\ploymer_model')
from ploymer_eos.pcsaft_elec import *
from scipy import optimize             # 求解VLE模型(非线性方程组)
from scipy.optimize import differential_evolution, Bounds, minimize, NonlinearConstraint
from utility import SolutionError


def bubble_pressure(T, x, args, max_iter=100):
    """
    Given the temperature and liquid phase composition, calculate the bubble point pressure
    
    Parameters
    ----------
    T ：float, Given temperature(K)
    x : ndarray, shape (n,), liquid phase composition in molar fractions
    args : A dictionary includes critical temperature "Tc", critical pressure "Pc" and eccentricity factor "w"
        If pc-saft equation is used, segment number "m", segment diameter "s" and 
        dispersion energy "e" are also included. 
        Take the mixed components of propylene, propane and hydrogen as an example:
        args = {"Tc": np.array([364.85, 369.83, 33.19]), 
                "Pc": np.array([4600000, 4248000, 1313000]), 
                "w": np.array([0.137588, 0.152291, -0.215993]),
                "m": np.array([1.9597, 2.002, 0.9846]),
                "s": np.array([3.5356, 3.6184, 2.8263]),
                "e": np.array([207.19, 208.11, 20.893])}
    max_iter : Maximum number of iterations

    Returns
    -------
    P : float, bubble point pressure(Pa)
    y :  ndarray,shape (n,), vapour phase composition
    iter_inner : Number of inner iterations until convergence
    iter  : Number of outer iterations until convergence
    """

    # Get input args
    Tc = args["Tc"]
    Pc = args["Pc"]
    omega = args["w"]
    pcsaft_args = {"m":args["m"], "s":args["s"], "e":args["e"]}
    
    # Initialize pressure and K-values (through Wilson equation)
    P = 0
    for i in range(len(x)):
        P = P + x[i]*Pc[i]*math.exp(5.373*(1+omega[i])*(1-Tc[i]/T))

    K = np.zeros(len(x))
    for i in range(len(x)): 
        K[i] = Pc[i]/P*math.exp(5.373*(1+omega[i])*(1-Tc[i]/T))

    # Initial gas phase composition yi
    y_old = K*x
    y = y_old

    # Initialize variable to converge (sum of vapour compositions=1) and iteration counter
    sum_y = 2
    iter = 0
    tempy = np.zeros(max_iter+1)   # Record "sum_y" of the iterative process
    tempp = np.zeros(max_iter+1)   # Record pressure of the iterative process

    # Outer loop: use Newton iterative method to modify the pressure make sum_y == 1
    while (abs(sum_y-1)>1e-7) & (iter<max_iter) :
        if iter == 0:
            pass
        elif iter == 1:
            if sum_y > 1 :
                P = P*1.01
            else :
                P = P*0.99
        else:
            Pnew = tempp[iter]+(1-tempy[iter])*(tempp[iter]-tempp[iter-1])/(tempy[iter]-tempy[iter-1])
            if Pnew < 0 :
                P = 0.9*P
            else :
                P = Pnew

        # Inner loop: use the equation of state to calculate K-value to make yi converge
        iter_inner = 0
        dy = np.ones(len(x))

        rho_l = pcsaft.pcsaft_den(T, P, x, pcsaft_args, "liq")
        fugcoef_l = pcsaft.pcsaft_fugcoef(T, rho_l, x, pcsaft_args)
        while (np.linalg.norm(dy)>1e-7) & (iter_inner<max_iter) :
            iter_inner = iter_inner+1
            y = y/sum(y)
            rho_v = pcsaft.pcsaft_den(T, P, y, pcsaft_args, "vap")
            fugcoef_v = pcsaft.pcsaft_fugcoef(T, rho_v, y, pcsaft_args)
            K = fugcoef_l/fugcoef_v   # 聚合物气相逸度系数可能为0, 相平衡常数为Nan
            np.nan_to_num(K,False,0)  # 因此将K中的Nan替换为0
            y = K*x
            dy = y-y_old
            y_old = y
        
        sum_y = sum(y_old)
        iter = iter+1
        tempy[iter] = sum_y
        tempp[iter] = P
    return P, y, iter_inner, iter

def bubble_temperature(P, x, args, max_iter=100):
    """
    Given the pressure and liquid phase composition, calculate the bubble point temperature
    
    Parameters
    ----------
    P ：float, Given pressure(Pa)
    x : ndarray, shape (n,), liquid phase composition in molar fractions
    args : A dictionary includes critical temperature "Tc", critical pressure "Pc" and eccentricity factor "w"
        If pc-saft equation is used, segment number "m", segment diameter "s" and 
        dispersion energy "e" are also included. 
        Take the mixed components of propylene, propane and hydrogen as an example:
        args = {"Tc": np.array([364.85, 369.83, 33.19]), 
                "Pc": np.array([4600000, 4248000, 1313000]), 
                "w": np.array([0.137588, 0.152291, -0.215993]),
                "m": np.array([1.9597, 2.002, 0.9846]),
                "s": np.array([3.5356, 3.6184, 2.8263]),
                "e": np.array([207.19, 208.11, 20.893])}

    Returns
    -------
    T : float, bubble point temperature(K)
    y :  ndarray,shape (n,), vapour phase composition
    iter_inner : Number of inner iterations until convergence
    iter  : Number of outer iterations until convergence
    """

    # Get input args
    Tc = args["Tc"]
    Pc = args["Pc"]
    omega = args["w"]
    pcsaft_args = {"m":args["m"], "s":args["s"], "e":args["e"]}

    # Initialize Temperature and K-values
    T = 350
    dT = 1
    iter = 0
    K = np.zeros(len(x))
    while (abs(dT)>1e-7) & (iter<max_iter):
        for i in range(len(x)):
            K[i]=Pc[i]/P*math.exp(5.373*(1+omega[i])*(1-Tc[i]/T) )
        f = sum(K*x)-1
        dfdT = 1/T**2*sum(5.373*(1+omega)*Tc*K*x)
        dT = -f/dfdT
        Tnew = T+dT
        if Tnew < 0:
            T = 0.9*T
        else:
            T = Tnew
        iter = iter+1
    
    # Initial gas phase composition yi
    y_old = K*x
    y = y_old

    # Initialize variable to converge (sum of vapour compositions=1) and iteration counter
    sum_y = 2
    iter = 0
    tempy = np.zeros(max_iter+1)  # Record "sum_y" of the iterative process
    tempt = np.zeros(max_iter+1)  # Record temperature of the iterative process

    # Outer loop: use Newton iterative method to modify the temperature make sum_y == 1
    while (abs(sum_y-1)>1e-7) & (iter<max_iter) :
        if iter == 0:
            pass
        elif iter == 1:
            if sum_y < 1 :
                T = T*1.01
            else :
                T = T*0.99
        else:
            Tnew = tempt[iter]+(1-tempy[iter])*(tempt[iter]-tempt[iter-1])/(tempy[iter]-tempy[iter-1])
            if Tnew < 0 :
                T = 0.9*T
            else :
                T = Tnew
        
        # Inner loop: use the equation of state to calculate K-value to make yi converge
        iter_inner = 0
        dy = np.ones(len(x))
        
        rho_l = pcsaft.pcsaft_den(T, P, x, pcsaft_args, "liq")
        fugcoef_l = pcsaft.pcsaft_fugcoef(T, rho_l, x, pcsaft_args)
        while (np.linalg.norm(dy)>1e-7) & (iter_inner<max_iter) :
            iter_inner = iter_inner+1
            y = y/sum(y)
            rho_v = pcsaft.pcsaft_den(T, P, y, pcsaft_args, "vap")
            fugcoef_v = pcsaft.pcsaft_fugcoef(T, rho_v, y, pcsaft_args)
            K = fugcoef_l/fugcoef_v   # 聚合物气相逸度系数可能为0, 相平衡常数为Nan
            np.nan_to_num(K,False,0)  # 因此将K中的Nan替换为0
            y = K*x
            dy = y-y_old
            y_old = y
            
        sum_y = sum(y_old)
        iter = iter+1
        tempy[iter] = sum_y
        tempt[iter] = T
    return T, y, iter_inner, iter

def dew_pressure(T, y, args, max_iter=100):
    """"
    Given the temperature and vapor phase composition, calculate the dew point pressure
    
    Parameters
    ----------
    T ：float, Given temperature(K)
    y : ndarray, shape (n,), vapor phase composition in molar fractions
    args : A dictionary includes critical temperature "Tc", critical pressure "Pc" and eccentricity factor "w"
        If pc-saft equation is used, segment number "m", segment diameter "s" and 
        dispersion energy "e" are also included. 
        Take the mixed components of propylene, propane and hydrogen as an example:
        args = {"Tc": np.array([364.85, 369.83, 33.19]), 
                "Pc": np.array([4600000, 4248000, 1313000]), 
                "w": np.array([0.137588, 0.152291, -0.215993]),
                "m": np.array([1.9597, 2.002, 0.9846]),
                "s": np.array([3.5356, 3.6184, 2.8263]),
                "e": np.array([207.19, 208.11, 20.893])}

    Returns
    -------
    P : float, dew point pressure(Pa)
    x :  ndarray,shape (n,), liquid phase composition
    iter_inner : Number of inner iterations until convergence
    iter  : Number of outer iterations until convergence
    """

    # Get input args
    Tc = args["Tc"]
    Pc = args["Pc"]
    omega = args["w"]
    pcsaft_args = {"m":args["m"], "s":args["s"], "e":args["e"]}
    
    # Initialize pressure and K-values (through Wilson equation)
    P = 101325
    dP = 1
    iter = 0
    K = np.zeros(len(y))
    while (abs(dP)>1e-7) & (iter<max_iter):
        for i in range(len(y)):
            K[i]=Pc[i]/P*math.exp(5.373*(1+omega[i])*(1-Tc[i]/T))
        f = sum(y/K)-1
        dfdP = sum(y/K/P) 
        dP = -f/dfdP
        Pnew = P + dP
        if Pnew < 0:
            P = 0.9*P
        else:
            P = Pnew
        iter = iter+1
 
    # Initial liquid phase composition yi
    x_old = y/K
    x = x_old 

    # Initialize variable to converge (sum of vapour compositions=1) and iteration counter
    sum_x = 2
    iter = 0
    tempx = np.zeros(max_iter+1)   # Record "sum_x" of the iterative process
    tempp = np.zeros(max_iter+1)   # Record pressure of the iterative process

    # Outer loop: use Newton iterative method to modify the pressure make sum_x == 1
    while (abs(sum_x-1)>1e-7) & (iter<max_iter) :
        if iter == 0:
            pass
        elif iter == 1:
            if sum_x > 1 :
                P = P*0.99
            else :
                P = P*1.01
        else:
            Pnew = tempp[iter]+(1-tempx[iter])*(tempp[iter]-tempp[iter-1])/(tempx[iter]-tempx[iter-1])
            if Pnew < 0 :
                P = 0.9*P
            else :
                P = Pnew

        # Inner loop: use the equation of state to calculate K-value to make xi converge
        iter_inner = 0
        dx = np.ones(len(y))

        rho_v = pcsaft.pcsaft_den(T, P, y, pcsaft_args, "vap")
        fugcoef_v = pcsaft.pcsaft_fugcoef(T, rho_v, y, pcsaft_args)
        while (np.linalg.norm(dx)>1e-7) & (iter_inner<max_iter) :
            iter_inner = iter_inner+1
            x = x/sum(x)
            rho_l = pcsaft.pcsaft_den(T, P, x, pcsaft_args, "liq")
            fugcoef_l = pcsaft.pcsaft_fugcoef(T, rho_l, x, pcsaft_args)
            K = fugcoef_l/fugcoef_v   # 聚合物气相逸度系数可能为0, 相平衡常数为Nan
            np.nan_to_num(K,False,0)  # 因此将K中的Nan替换为0
            x = y/K
            dx = x-x_old
            x_old = x
        
        sum_x = sum(x_old)
        iter = iter+1
        tempx[iter] = sum_x
        tempp[iter] = P
    return P, x, iter_inner, iter

def dew_temperature(P, y, args, max_iter=100):
    """
    Given the pressure and vapor phase composition, calculate the dew point temperature.
    
    Parameters
    ----------
    P ：float, Given pressure(Pa)
    y : ndarray, shape (n,), vapor phase composition in molar fractions
    args : A dictionary includes critical temperature "Tc", critical pressure "Pc" and eccentricity factor "w"
        If pc-saft equation is used, segment number "m", segment diameter "s" and 
        dispersion energy "e" are also included. 
        Take the mixed components of propylene, propane and hydrogen as an example:
        args = {"Tc": np.array([364.85, 369.83, 33.19]), 
                "Pc": np.array([4600000, 4248000, 1313000]), 
                "w": np.array([0.137588, 0.152291, -0.215993]),
                "m": np.array([1.9597, 2.002, 0.9846]),
                "s": np.array([3.5356, 3.6184, 2.8263]),
                "e": np.array([207.19, 208.11, 20.893])}

    Returns
    -------
    T : float, dew point temperature(K)
    x :  ndarray,shape (n,), liquid phase composition
    iter_inner : Number of inner iterations until convergence
    iter  : Number of outer iterations until convergence
    """

    # Get input args
    Tc = args["Tc"]
    Pc = args["Pc"]
    omega = args["w"]
    pcsaft_args = {"m":args["m"], "s":args["s"], "e":args["e"]}

    # Initialize Temperature and K-values
    T = 350
    dT = 1
    iter = 0
    K = np.zeros(len(y))
    while (abs(dT)>1e-7) & (iter<max_iter):
        for i in range(len(y)):
            K[i]=Pc[i]/P*math.exp(5.373*(1+omega[i])*(1-Tc[i]/T))
        f = sum(y/K)-1
        dfdT = -1/T**2*sum(5.373*(1+omega)*y*Tc/K)
        dT = -f/dfdT
        Tnew = T+dT
        if Tnew < 0:
            T = 0.9*T
        else:
            T = Tnew
        iter = iter+1
    
    # Initial gas phase composition yi
    x_old = y/K
    x = x_old

    # Initialize variable to converge (sum of vapour compositions=1) and iteration counter
    sum_x = 2
    iter = 0
    tempx = np.zeros(max_iter+1)  # Record "sum_x" of the iterative process
    tempt = np.zeros(max_iter+1)  # Record temperature of the iterative process

    # Outer loop: use Newton iterative method to modify the temperature make sum_x == 1
    while (abs(sum_x-1)>1e-7) & (iter<max_iter) :
        if iter == 0:
            pass
        elif iter == 1:
            if sum_x < 1 :
                T = T*0.99
            else :
                T = T*1.01
        else:
            Tnew = tempt[iter]+(1-tempx[iter])*(tempt[iter]-tempt[iter-1])/(tempx[iter]-tempx[iter-1])
            if Tnew < 0 :
                T = 0.9*T
            else :
                T = Tnew
        
        # Inner loop: use the equation of state to calculate K-value to make yi converge
        iter_inner = 0
        dx = np.ones(len(y))

        rho_v = pcsaft.pcsaft_den(T, P, y, pcsaft_args, "vap")
        fugcoef_v = pcsaft.pcsaft_fugcoef(T, rho_v, y, pcsaft_args)
        while (np.linalg.norm(dx)>1e-7) & (iter_inner<max_iter) :
            iter_inner = iter_inner+1
            x = x/sum(x)
            rho_l = pcsaft.pcsaft_den(T, P, x, pcsaft_args, "liq")
            fugcoef_l = pcsaft.pcsaft_fugcoef(T, rho_l, x, pcsaft_args)
            K = fugcoef_l/fugcoef_v   # 聚合物气相逸度系数可能为0, 相平衡常数为Nan
            np.nan_to_num(K,False,0)  # 因此将K中的Nan替换为0
            x =  y/K
            dx = x-x_old
            x_old = x
            
        sum_x = sum(x_old)
        iter = iter+1
        tempx[iter] = sum_x
        tempt[iter] = T
    return T, x, iter_inner, iter


def vle_model(k, t, p, z, args): 
    """
    Return gas-liquid equilibrium model
    
    Parameters
    ----------
    k ：ndarray, shape (2n+1,), n is the number of components
    t : float, temperature(K)
    p : float, pressure(Pa)
    z : ndarray, shape (n, ),  feed composition
    args : A dictionary includes critical temperature "Tc", critical pressure "Pc" and eccentricity factor "w"
        If pc-saft equation is used, segment number "m", segment diameter "s" and 
        dispersion energy "e" are also included. 
        Take the mixed components of propylene, propane and hydrogen as an example:
        args = {"Tc": np.array([364.85, 369.83, 33.19]), 
                "Pc": np.array([4600000, 4248000, 1313000]), 
                "w": np.array([0.137588, 0.152291, -0.215993]),
                "m": np.array([1.9597, 2.002, 0.9846]),
                "s": np.array([3.5356, 3.6184, 2.8263]),
                "e": np.array([207.19, 208.11, 20.893])}

    Returns
    -------
    gas-liquid equilibrium model
    """

    k[k<0] = 0.01    # 迭代过程中,可能出现负数，将其替换为0.01
    k[k>1] = 0.99    # 迭代过程中,可能出现大于1的数，将其替换为0.99

    length = len(k)
    n = int((length-1)/2)  # number of components
    β = k[0]               # Vapor fraction
    x = k[1: n+1]          # Liquid composition
    y = k[n+1: 2*n+1]      # Vapor composition
    x = x/sum(x)           # Normalization
    y = y/sum(y)           # Normalization

    pcsaft_args = {"m":args["m"], "s":args["s"], "e":args["e"]}

    rho_l = pcsaft.pcsaft_den(t, p, x, pcsaft_args, "liq")       # Liquid Phase molar density
    fugcoef_l = pcsaft.pcsaft_fugcoef(t, rho_l, x, pcsaft_args)  # Liquid Phase Fugacities
    rho_v = pcsaft.pcsaft_den(t, p, y, pcsaft_args, "vap")       # Vapor Phase molar density
    fugcoef_v = pcsaft.pcsaft_fugcoef(t, rho_v, y, pcsaft_args)  # Vapor Phase Fugacities

    # VLE model
    model = np.zeros(length)
    model[0] = sum(x) - sum(y)            # 组分归一化
    model[1: n+1] = β*y + (1-β)*x - z     # 物料衡算
    model[n+1: 2*n+1] = y*fugcoef_v - x*fugcoef_l   # 相平衡

    return model

def tp_flash(t, p, z, args, max_iter=100):
    """
    Two-Phase PT Flash
    """

    # Set initial value
    n = len(z)             # Number of components
    x0 = np.zeros(n*2+1)   # Initialize array
    x0[0] = 0.8            # Set initial value of β to 0.5
    x0[1: n+1] = z         # Take feed composition as initial value of liquid composition
    x0[n+1: 2*n+1] = z     # Take feed composition as initial value of vapor composition

    result = optimize.root(vle_model, x0, args=(t, p, z, args), method='hybr', tol=1e-10)
    if result.success :
        β = result.x[0]
        x = result.x[1: n+1]
        x = x/sum(x)
        y = result.x[n+1: 2*n+1]
        y = y/sum(y)
        return β, x, y, result.fun
    else :
        print("Solution failed")
        raise SolutionError(result.message)

def VLE_TPflash(T, P, z, args, tol=1e-7, max_iter=100):
    """
    Two-Phase PT Flash by Successive substitution and the Rachford-Rice equation
    """

    # Get input args
    Tc = args["Tc"]
    Pc = args["Pc"]
    omega = args["w"]
    pcsaft_args = {"m":args["m"], "s":args["s"], "e":args["e"]}

    # Initialize K-value by Wilson equation
    K = Pc/P*np.exp(5.373*(1+omega)*(1-Tc/T))  

    beta, x, y = Rachford_Rice_solver(K,z)        # Call to RR_solver with actual K's

    # Update K-Factors
    rho_l = pcsaft.pcsaft_den(T, P, x, pcsaft_args, "liq")
    fugcoef_l = pcsaft.pcsaft_fugcoef(T, rho_l, x, pcsaft_args)  # Calculate Liquid Phase Fugacities
    rho_v = pcsaft.pcsaft_den(T, P, y, pcsaft_args, "vap")
    fugcoef_v = pcsaft.pcsaft_fugcoef(T, rho_v, y, pcsaft_args)  # Calculate Vapor Phase Fugacities
    # Calculate equilibrium constant (Vapor Phase Fugacities may be 0) 
    K = np.zeros(len(fugcoef_l))
    for i in range(len(fugcoef_l)):
        if fugcoef_v[i] != 0:
            K[i] = fugcoef_l[i]/fugcoef_v[i]
    
    Kb = 1.0    # 1个常数
    u = K/Kb    # 由于K可能为0, 所以没有取对数

    R = beta
    K0 = 1

    error_out = 2*tol   # 外层循环误差
    error_in = 2*tol
    iter = 0
    while error_out > tol and iter < max_iter:
        iter += 1
        inner_iter = 0
        temp_e = np.zeros(max_iter+1)
        temp_r = np.zeros(max_iter+1)
        while abs(error_in) > tol and inner_iter < max_iter:
            p = z/(1-R+K0*R*u)
            # mid(a,b, c) = max(min(a,b),min(max(a,b), c))
            error_in = max(min(R,R-1),min(max(R,R-1), np.sum(p)-Kb*np.sum(u*p)))
            temp_r[inner_iter] = R
            temp_e[inner_iter] = error_in
            if abs(error_in) <= tol:
                break
            if inner_iter <= 1:
                if error_in > 0:
                    R = 0.9*R
                elif error_in < 0:
                    R = 1.1*R
            else:
                R_new = R - error_in/(temp_e[inner_iter]-temp_e[inner_iter-1])*(temp_r[inner_iter]-temp_r[inner_iter-1])
                if R_new >= 1:
                    R = 1.1*R
                elif R_new <= 0:
                    R = 0.9*R
                else:
                    R = R_new
            inner_iter += 1

            # Update K-Factors
        print(": ",p)
        x = p/np.sum(p)
        y = u*p/np.sum(u*p)
        beta = 1-(1-R)*np.sum(p)
        rho_l = pcsaft.pcsaft_den(T, P, x, pcsaft_args, "liq")
        fugcoef_l = pcsaft.pcsaft_fugcoef(T, rho_l, x, pcsaft_args)  # Calculate Liquid Phase Fugacities
        rho_v = pcsaft.pcsaft_den(T, P, y, pcsaft_args, "vap")
        fugcoef_v = pcsaft.pcsaft_fugcoef(T, rho_v, y, pcsaft_args)  # Calculate Vapor Phase Fugacities
        # Calculate equilibrium constant (Vapor Phase Fugacities may be 0) 
        K = np.zeros(len(fugcoef_l))
        for i in range(len(fugcoef_l)):
            if fugcoef_v[i] != 0:
                K[i] = fugcoef_l[i]/fugcoef_v[i]
        u_new = K/Kb    # 由于K可能为0, 所以没有取对数
        error_out = np.max(np.abs(u_new-u))
        u = u_new
    
    return beta, x, y, iter, inner_iter


def Rachford_Rice_solver(K, z, solver_type=0):
    """
    Solver for the Rachford-Rice Equations used in Flash Calculations.
    Given equilibrium constant 'K' and Feed composition 'z', the vapor fraction 'beta', vapor phase composition 'x' 
    and liquid phase composition 'y' are calculated by Newton iterative method.
    Solver types: 0 —— Full RRSolver 
                1 —— Negative RR_Solver
    """
    
    if solver_type == 0:

        # Check if it will flash
        g_0 = Rachford_Rice(0, z, K)[0]  # g for beta=0
        if g_0 < 0:
            print("过冷液体")
            return Rachford_Rice_solver(1.5*K, z, solver_type)  # 增大K值迭代

        if len(np.where(K==0)[0]) == 0:      # 如果K中有0，则有聚合物，必定有液相
            g_1 = Rachford_Rice(1, z, K)[0]  # g for beta=1
            if g_1 > 0:
                print("过热气体:", K)
                return Rachford_Rice_solver(0.75*K, z, solver_type)  # 减小K值迭代
        
        # Set upper and lower bounds for beta
        beta_min = 0
        beta_max = 1
    elif solver_type == 1:
        beta_min = -1/(max(K)-1)
        beta_max = 1/(1-min(K))
    else:
        raise Exception("Invalid solver type")

    # Set a better lower bound for beta
    K_array = K[K>1]
    z_array = z[K>1]
    if len(K_array) > 0:
        beta_array = (K_array*z_array-1)/(K_array-1)
        if max(beta_array) > 0:
            beta_min = max(beta_array)
    
    # Set a better upper bound for beta
    K_array = K[K<1]
    z_array = z[K<1]
    if len(K_array) > 0:
        beta_array = (1-z_array)/(1-K_array)
        if min(beta_array) < 1:
            beta_max = min(beta_array)

    beta = (beta_min + beta_max)/2
    g, dg = Rachford_Rice(beta, z, K)

    delta_beta = 1
    counter = 0
    # There may be dead circulation, which needs to be modified
    while (abs(g) > 1e-7) | (abs(delta_beta) > 1e-7):
        # Update upper or lower bound for beta
        if g > 0:
            beta_min = beta
        else:
            beta_max = beta
    
        delta_beta = -g/dg
        beta_new = beta + delta_beta

        if beta_new >= beta_min and beta_new <= beta_max:
            beta = beta_new
        else:
            delta_beta = (beta_min + beta_max)/2 - beta
            beta = beta + delta_beta

        g, dg = Rachford_Rice(beta, z, K)
        counter = counter+1
    
    x  = z/(1-beta+beta*K)
    y = K*x
    return beta, x, y

def Rachford_Rice(beta, z, K):
    """
    Subroutine that given a vapour fraction, composition and K-factors 
    return the values of g and dg of Rachford-Rice equation. 
    Page 252 of Thermodynamic models:Fundamentals and Computational Aspects
    """

    g = sum(z*(K-1)/(1-beta+beta*K))
    dg = -sum(z*(K-1)**2/(1-beta+beta*K)**2)
    return g, dg


def stability_analysis(T, P, z, K, args, max_iter=100):
    """
    Stability analysis for the Two-Phases PT-Flash
    """

    pcsaft_args = {"m":args["m"], "s":args["s"], "e":args["e"]}
    rho_zl = pcsaft.pcsaft_den(T, P, z, pcsaft_args, "liq")
    fugcoef_zl = pcsaft.pcsaft_fugcoef(T, rho_zl, z, pcsaft_args)  # Assume feed as liquid phase

    rho_zv = pcsaft.pcsaft_den(T, P, z, pcsaft_args, "vap")
    fugcoef_zv = pcsaft.pcsaft_fugcoef(T, rho_zv, z, pcsaft_args)  # Assume feed as liquid phase

    dL = np.log(z) + np.log(fugcoef_zl)   # Assume feed as liquid phase
    dV = np.log(z) + np.log(fugcoef_zv)   # Assume feed as vapor phase

    WL=z/K   # Liquid-Like initial estimate
    WV=K*z   # Vapour-Like initial estimate

    WL=WL/sum(WL)
    WV=WV/sum(WV)

    dW = np.ones(len(z))
    iter_l = 0
    while (max(abs(dW))>1e-7) & (iter_l<max_iter):
        iter_l = iter_l+1
        rho_l = pcsaft.pcsaft_den(T, P, WL, pcsaft_args, "liq")
        fugcoef_l = pcsaft.pcsaft_fugcoef(T, rho_l, WL, pcsaft_args)
        WL_new = np.exp(dV-np.log(fugcoef_l))
        dW = WL_new-WL
        WL=WL_new/sum(WL_new)
    rho_l = pcsaft.pcsaft_den(T, P, WL, pcsaft_args, "liq")
    fugcoef_l = pcsaft.pcsaft_fugcoef(T, rho_l, WL, pcsaft_args)
    tm_l = 1+sum(WL*(np.log(WL)+np.log(fugcoef_l)-dV-1))

    dW = np.ones(len(z))
    iter_v = 0
    while (max(abs(dW))>1e-7) & (iter_v<max_iter):
        iter_v=iter_v+1
        rho_v = pcsaft.pcsaft_den(T, P, WV, pcsaft_args, "vap")
        fugcoef_v = pcsaft.pcsaft_fugcoef(T, rho_v, WV, pcsaft_args)
        WV_new = np.exp(dL-np.log(fugcoef_v))
        dW = WV_new-WV
        WV = WV_new/sum(WV_new)
    rho_v = pcsaft.pcsaft_den(T, P, WV, pcsaft_args, "vap")
    fugcoef_v = pcsaft.pcsaft_fugcoef(T, rho_v, WV, pcsaft_args)
    tm_v = 1+sum(WV*(np.log(WV)+np.log(fugcoef_v)-dL-1)) 
    return tm_l, tm_v, fugcoef_l, fugcoef_v, fugcoef_zl, fugcoef_zv


# Test gibbs minization procedure
if __name__ == "__main__":
    # Flash test 1: "C3H6", "C3H8", "H2", "N2"
    T = 341.15
    P = 3000000
    z = np.array([0.82, 0.16, 0.005, 0.015])
    args = {"Tc": np.array([364.85, 369.83, 33.19, 126.2]), 
            "Pc": np.array([4600000, 4248000, 1313000, 3400000]), 
            "w": np.array([0.137588, 0.152291, -0.215993, 0.0377215]),
            "m": np.array([1.9597, 2.002, 0.9846, 1.2053]),
            "s": np.array([3.5356, 3.6184, 2.8263, 3.313]),
            "e": np.array([207.19, 208.11, 20.893, 90.96])}

    beta, x, y, f, _ = VLE_TPflash(T, P, z, args)
    print("Vapor fraction: ", beta)
    print("Liquid composition: ", x)
    print("Vapor composition: ", y)
    print("f: ", f)

    print("\n",80*"-","\n")

    # Flash test 2: "C3H6", "C3H8", "H2", "N2", "PP", "TICL4", "TEA"
    T = 223.15
    P = 100000
    PP_MWN = 119735
    PP_DPN = 2845.37
    z = np.array([0.7614716, 0.1435221, 0.0784851, 0.00720122, 0.00931595/PP_DPN, 6.25573E-07, 3.44657E-06])
    z = z/sum(z)
    args = {"Tc": np.array([364.85, 369.83, 33.19, 126.2, 2000, 638, 678.15]), 
            "Pc": np.array([4600000, 4248000, 1313000, 3400000, 5000000, 4660950, 8930000]), 
            "w": np.array([0.137588, 0.152291, -0.215993, 0.0377215, 0, 0.283732, 0.841783]),
            "m": np.array([1.9597, 2.002, 0.9846, 1.2053, 0.02528*PP_MWN, 20, 20]),
            "s": np.array([3.5356, 3.6184, 2.8263, 3.313, 4.1473, 3.54, 3.54]),
            "e": np.array([207.19, 208.11, 20.893, 90.96, 298.6392, 210, 210])}

    beta, x, y, f, _ = VLE_TPflash(T, P, z, args)
    print(sum(x), sum(y))
    print("Vapor fraction: ", beta)
    print("Liquid composition: ", x)
    x[4] = x[4]*PP_DPN
    beta = beta/(beta+(1-beta)*sum(x))
    print("新气相分数: ", beta)
    x = x/sum(x)
    print("新液相组成: ", x)
    print("Vapor composition: ", y)
    print("f: ", f)