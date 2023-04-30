"""
The flash algorithm includes: 
    Calculation of bubble point and dew point and TP flash、

Reference: 
Michelsen & Mollerup - Thermodynamic models:Fundamentals and Computational Aspects

"""

import numpy as np
import math
from polymer_model.property_method.property_method_wrapper import PropertyMethod
from polymer_model.utility import InputError

def bubble_pressure(T, x, eos, tol=1e-7, max_iter=100):
    """
    Given the temperature and liquid phase composition, calculate the bubble point pressure
    
    Parameters
    ----------
    T ：float, Given temperature(K)
    x : ndarray, shape (n,), liquid phase composition in molar fractions
    eos : PropertyMethod
    max_iter : Maximum number of iterations

    Returns
    -------
    P : float, bubble point pressure(Pa)
    y :  ndarray,shape (n,), vapour phase composition
    iter_inner : Number of inner iterations until convergence
    iter  : Number of outer iterations until convergence
    """

    # Get input args
    Tc = eos.property_args["Tc"]
    Pc = eos.property_args["Pc"]
    omega = eos.property_args["w"]
    
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
    while (abs(sum_y-1)>tol) & (iter<max_iter) :
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

        rho_l = eos.molar_density(T, P, x, "liq")
        fugcoef_l = eos.fugacity_coefficients(T, P, x, "liq")
        while (np.linalg.norm(dy)>tol) & (iter_inner<max_iter) :
            iter_inner = iter_inner+1
            y = y/sum(y)
            rho_v = eos.molar_density(T, P, y, "vap")
            fugcoef_v = eos.fugacity_coefficients(T, P, y, "vap")

            # Calculate equilibrium constant (Vapor Phase Fugacities may be 0) 
            K = np.zeros(len(fugcoef_l))
            for i in range(len(fugcoef_l)):
                if fugcoef_v[i] != 0:
                    K[i] = fugcoef_l[i]/fugcoef_v[i]

            y = K*x
            dy = y-y_old
            y_old = y
        
        sum_y = sum(y_old)
        iter = iter+1
        tempy[iter] = sum_y
        tempp[iter] = P

    state = True
    if iter_inner >= max_iter or iter >= max_iter:
        state = False
    if abs(rho_l-rho_v) < 1:
        state = False
    return P, y, state

def bubble_temperature(P, x, eos, tol=1e-7, max_iter=100):
    """
    Given the pressure and liquid phase composition, calculate the bubble point temperature
    
    Parameters
    ----------
    P ：float, Given pressure(Pa)
    x : ndarray, shape (n,), liquid phase composition in molar fractions
    eos : PropertyMethod

    Returns
    -------
    T : float, bubble point temperature(K)
    y :  ndarray,shape (n,), vapour phase composition
    iter_inner : Number of inner iterations until convergence
    iter  : Number of outer iterations until convergence
    """

    # Get input args
    Tc = eos.property_args["Tc"]
    Pc = eos.property_args["Pc"]
    omega = eos.property_args["w"]

    # Initialize Temperature and K-values
    T = 350
    dT = 1
    iter = 0
    K = np.zeros(len(x))
    while (abs(dT)>tol) & (iter<max_iter):
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
    while ((abs(sum_y-1)>tol) and (iter<max_iter)):
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
        
        rho_l = eos.molar_density(T, P, x, "liq")
        fugcoef_l = eos.fugacity_coefficients(T, P, x, "liq")
        while (np.linalg.norm(dy)>tol) & (iter_inner<max_iter) :
            iter_inner = iter_inner+1
            y = y/sum(y)
            rho_v = eos.molar_density(T, P, y, "vap")
            fugcoef_v = eos.fugacity_coefficients(T, P, y, "vap")
            # Calculate equilibrium constant (Vapor Phase Fugacities may be 0) 
            K = np.zeros(len(fugcoef_l))
            for i in range(len(fugcoef_l)):
                if fugcoef_v[i] != 0:
                    K[i] = fugcoef_l[i]/fugcoef_v[i]

            y = K*x
            dy = y-y_old
            y_old = y
            
        sum_y = sum(y_old)
        iter = iter+1
        tempy[iter] = sum_y
        tempt[iter] = T

    state = True
    if iter_inner >= max_iter or iter >= max_iter:
        state = False
    if abs(rho_l-rho_v) < 1:
        state = False

    return T, y, state

def dew_pressure(T, y, eos, tol=1e-7, max_iter=100):
    """"
    Given the temperature and vapor phase composition, calculate the dew point pressure
    
    Parameters
    ----------
    T ：float, Given temperature(K)
    y : ndarray, shape (n,), vapor phase composition in molar fractions
    eos : PropertyMethod

    Returns
    -------
    P : float, dew point pressure(Pa)
    x :  ndarray,shape (n,), liquid phase composition
    iter_inner : Number of inner iterations until convergence
    iter  : Number of outer iterations until convergence
    """

    # Get input args
    Tc = eos.property_args["Tc"]
    Pc = eos.property_args["Pc"]
    omega = eos.property_args["w"]
    
    # Initialize pressure and K-values (through Wilson equation)
    P = 1
    dP = 1
    iter = 0
    K = np.zeros(len(y))
    while (abs(dP)>tol) & (iter<max_iter):
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
    while (abs(sum_x-1)>tol) & (iter<max_iter) :
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

        rho_v = eos.molar_density(T, P, y, "vap")
        fugcoef_v = eos.fugacity_coefficients(T, P, y, "vap")
        while (np.linalg.norm(dx)>tol) & (iter_inner<max_iter) :
            iter_inner = iter_inner+1
            x = x/sum(x)
            rho_l = eos.molar_density(T, P, x, "liq")
            fugcoef_l = eos.fugacity_coefficients(T, P, x, "liq")

            # Calculate equilibrium constant (Vapor Phase Fugacities may be 0) 
            K = np.zeros(len(fugcoef_l))
            for i in range(len(fugcoef_l)):
                if fugcoef_v[i] != 0:
                    K[i] = fugcoef_l[i]/fugcoef_v[i]
            x = y/K
            dx = x-x_old
            x_old = x
        
        sum_x = sum(x_old)
        iter = iter+1
        tempx[iter] = sum_x
        tempp[iter] = P
    
    state = True
    if iter_inner >= max_iter or iter >= max_iter:
        state = False
    if abs(rho_l-rho_v) < 1:
        state = False

    return P, x, state

def dew_temperature(P, y, eos, tol=1e-7, max_iter=100):
    """
    Given the pressure and vapor phase composition, calculate the dew point temperature.
    
    Parameters
    ----------
    P ：float, Given pressure(Pa)
    y : ndarray, shape (n,), vapor phase composition in molar fractions
    eos : PropertyMethod

    Returns
    -------
    T : float, dew point temperature(K)
    x :  ndarray,shape (n,), liquid phase composition
    iter_inner : Number of inner iterations until convergence
    iter  : Number of outer iterations until convergence
    """

    # Get input args
    Tc = eos.property_args["Tc"]
    Pc = eos.property_args["Pc"]
    omega = eos.property_args["w"]

    # Initialize Temperature and K-values
    T = 350
    dT = 1
    iter = 0
    K = np.zeros(len(y))
    while (abs(dT)>tol) & (iter<max_iter):
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
    while (abs(sum_x-1)>tol) & (iter<max_iter) :
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

        rho_v = eos.molar_density(T, P, y, "vap")
        fugcoef_v = eos.fugacity_coefficients(T, P, y, "vap")
        while (np.linalg.norm(dx)>tol) & (iter_inner<max_iter) :
            iter_inner = iter_inner+1
            x = x/sum(x)
            rho_l = eos.molar_density(T, P, x, "liq")
            fugcoef_l = eos.fugacity_coefficients(T, P, x, "liq")

            # Calculate equilibrium constant (Vapor Phase Fugacities may be 0) 
            K = np.zeros(len(fugcoef_l))
            for i in range(len(fugcoef_l)):
                if fugcoef_v[i] != 0:
                    K[i] = fugcoef_l[i]/fugcoef_v[i]     

            x =  y/K
            dx = x-x_old
            x_old = x
            
        sum_x = sum(x_old)
        iter = iter+1
        tempx[iter] = sum_x
        tempt[iter] = T

    state = True
    if iter_inner >= max_iter or iter >= max_iter:
        state = False
    if abs(rho_l-rho_v) < 1:
        state = False
         
    return T, x, state


def TP_flash(T, P, z, eos, tol=1e-7, max_iter=100):
    """ 先计算泡露点, 判断是否处于两相区, 再考虑闪蒸 """

    iter = 0

    # 先计算泡露点, 判断是否处于两相区, 再考虑闪蒸
    Pb, _, state1 = bubble_pressure(T, z, eos, tol, max_iter)     # 泡点压力
    Tb, _, state2 = bubble_temperature(P, z, eos, tol, max_iter)  # 泡点温度
    Pd, _, state3 = dew_pressure(T, z, eos, tol, max_iter)        # 露点压力
    Td, _, state4 = dew_temperature(P, z, eos, tol, max_iter)     # 露点温度
    
    # 考虑纯物质的情况(P=Pb=Pd)
    if abs(Pb-Pd)<0.1 and abs(P-Pb)<0.1:
        raise InputError("Please input valid phases of the stream")
    else:
        # 高于泡点压力，处于液相
        if (state1 and P >= Pb) or (state2 and T <= Tb):
            x = z
            y = z
            beta = 0
        # 低于露点压力，处于气相
        elif (state3 and P <= Pd) or (state4 and T >= Td) :
            x = z
            y = z
            beta = 1
        else:
            beta, x, y, iter = VLE_TPflash(T, P, z, eos, tol, max_iter)
    
    return beta, x, y, iter


# 闪蒸算法
def VLE_TPflash(T, P, z, eos, tol=1e-7, max_iter=100):
    """
    闪蒸算法: 先用第一套算法计算, 未收敛, 则启动第二套算法
    """
    beta, x, y, iter = VLE_TPflash_first(T, P, z, eos, tol, max_iter)
    if not(beta>=0 or beta<=1) or iter == max_iter:
        beta, x, y, iter = VLE_TPflash_second(T, P, z, eos, tol, max_iter)
    return beta, x, y, iter


def VLE_TPflash_first(T, P, z, eos, tol=1e-7, max_iter=100):
    """
    Two-Phase PT Flash by Successive substitution and the Rachford-Rice equation
    使用RR算法, beta和x、y分层迭代, 该算法较为稳定, 但如果混合物中的聚合物含量较高,可能不收敛
    优先使用该算法

    """

    # Get input args
    Tc = eos.property_args["Tc"]
    Pc = eos.property_args["Pc"]
    omega = eos.property_args["w"]

    # Initialize K-value by Wilson equation
    K = Pc/P*np.exp(5.373*(1+omega)*(1-Tc/T))  

    # STABILITY ANALYSIS
    # 稳定性分析待完善
    stability = 1

    # Initialization of variables
    iter_total = 0   # Counter of the total number of iterations
    delta_K = 1     # Initialization of the convergence criterium

    K_array = np.ones((max_iter,len(z)))
    while (delta_K > tol) & (iter_total<max_iter):

        if stability == 1:
            beta, x, y = Rachford_Rice_solver(K, z, 0, tol)   # Call to RR_solver with actual K's
        else:
            beta, x, y = Rachford_Rice_solver(K, z, 1, tol)
        
        # Update K-Factors
        rho_l = eos.molar_density(T, P, x, "liq")
        fugcoef_l = eos.fugacity_coefficients(T, P, x, "liq")  # Calculate Liquid Phase Fugacities
        rho_v = eos.molar_density(T, P, y, "vap")
        fugcoef_v = eos.fugacity_coefficients(T, P, y, "vap")  # Calculate Vapor Phase Fugacities
        # Calculate equilibrium constant (Vapor Phase Fugacities may be 0) 
        K = np.zeros(len(fugcoef_l))
        for i in range(len(fugcoef_l)):
            if fugcoef_v[i] != 0:
                K[i] = fugcoef_l[i]/fugcoef_v[i]

        K_array[iter_total,:] = K;       
        if iter_total>=1:
            delta_K=max(abs(K_array[iter_total,:]-K_array[iter_total-1,:])); # Update Convergence criterium
        iter_total=iter_total+1;                        # Iteration Counter
    
    return beta, x, y, iter_total


def VLE_TPflash_second(T, P, z, eos, tol=1e-7, max_iter=100):
    """
    Two-Phase PT Flash by Successive substitution and the Rachford-Rice equation
    使用RR算法, beta和x、y同时迭代, 测试发现在聚合物含量高时,该算法计算结果和Aspen相近,
    但是迭代次数多, 作为第二算法(第一算法不收敛时启用)
    """

    # Get input args
    Tc = eos.property_args["Tc"]
    Pc = eos.property_args["Pc"]
    omega = eos.property_args["w"]

    # Initialize K-value by Wilson equation
    K = Pc/P*np.exp(5.373*(1+omega)*(1-Tc/T))  

    # Initialize x, y, beta
    x = z/K
    x = x/np.sum(x)
    y = K*z
    y = y/np.sum(y)
    beta = 0.9

    # Start iteration
    iter = 0
    while iter < max_iter:
        iter += 1
        rho_l = eos.molar_density(T, P, x, "liq")
        fugcoef_l = eos.fugacity_coefficients(T, P, x, "liq")  # Calculate Liquid Phase Fugacities
        rho_v = eos.molar_density(T, P, y, "vap")
        fugcoef_v = eos.fugacity_coefficients(T, P, y, "vap")  # Calculate Vapor Phase Fugacities
        # Calculate equilibrium constant (Vapor Phase Fugacities may be 0) 
        K = np.zeros(len(fugcoef_l))
        for i in range(len(fugcoef_l)):
            if fugcoef_v[i] != 0:
                K[i] = fugcoef_l[i]/fugcoef_v[i]
        g, dg = Rachford_Rice(beta, z, K)

        # Calculate error
        error = abs(g)
        if error <= tol:
            break

        delta_beta = -g/dg
        # beta_new = beta + delta_beta
        # beta_new = beta + (max_iter-iter)/max_iter*delta_beta
        beta_new = beta + (1-iter/max_iter)**2*delta_beta

        if beta_new >= 1:
            beta_new = (beta+1)/2
        elif beta_new <= 0:
            beta_new = beta*0.9
        beta = beta_new

        x  = z/(1-beta+beta*K)
        y = K*x
        x = x/np.sum(x)
        y = y/np.sum(y)

    return beta, x, y, iter

# RR公式的求解器
def Rachford_Rice_solver(K, z, solver_type=0, tol=1e-7, max_iter=100):
    """
    Solver for the Rachford-Rice Equations used in Flash Calculations.
    Given equilibrium constant 'K' and Feed composition 'z', the vapor fraction 'beta', vapor phase composition 'x' 
    and liquid phase composition 'y' are calculated by Newton iterative method.
    Solver types: 0 —— Full RRSolver 
                1 —— Negative RR_Solver
    """
    
    if solver_type == 0:

        # Check if it will flash
        index = 0
        flag = False
        while True:
            g_0 = Rachford_Rice(1e-8, z, K)[0]  # g for beta=0
            if g_0 >= 0 and flag:
                break
            while g_0 < 0:
                # print("过冷液体")
                K[index % len(K)] = 1.1*K[index % len(K)]
                g_0 = Rachford_Rice(1e-8, z, K)[0]  # g for beta=0
                index += 1
            
            if len(np.where(K==0)[0]) == 0:      # 如果K中有0，则有不挥发物质，必定有液相
                g_1 = Rachford_Rice(1-1e-8, z, K)[0]  # g for beta = 1
                if g_1 <= 0:
                    break  
                while g_1 > 0:
                    # print("过热气体:")
                    K[index % len(K)] = 0.9*K[index % len(K)]
                    g_1 = Rachford_Rice(1-1e-8, z, K)[0]  # g for beta=0
                    index += 1
                flag = True
            else:
                break
        
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
    while (abs(g) > tol) or (abs(delta_beta) > tol):
        # Update upper or lower bound for beta
        if g > 0:
            beta_min = beta
        else:
            beta_max = beta
    
        delta_beta = -g/dg
        beta_new = beta + delta_beta
        # beta_new = beta + (max_iter-iter)/max_iter*delta_beta
        
        # if beta_new >= beta_min and beta_new <= beta_max:
        #     beta = beta_new
        # else:
        #     delta_beta = (beta_min + beta_max)/2 - beta
        #     beta = beta + delta_beta
        if beta_new < beta_min:
            delta_beta = (beta + beta_min)/2 - beta
            beta = beta + delta_beta
        elif beta_new > beta_max:
            delta_beta = (beta + beta_max)/2 - beta
            beta = beta + delta_beta
        else:
            beta = beta_new

        g, dg = Rachford_Rice(beta, z, K)
        counter = counter+1
    x  = z/(1-beta+beta*K)
    y = K*x
    return beta, x, y


# RR公式及其导数
def Rachford_Rice(beta, z, K):
    """
    Subroutine that given a vapour fraction, composition and K-factors 
    return the values of g and dg of Rachford-Rice equation. 
    Page 252 of Thermodynamic models:Fundamentals and Computational Aspects
    """

    g = sum(z*(K-1)/(1-beta+beta*K))
    dg = -sum(z*(K-1)**2/(1-beta+beta*K)**2)
    return g, dg


# 稳定性分析
def stability_analysis(T, P, z, K, eos, max_iter=100):
    """
    Stability analysis for the Two-Phases PT-Flash
    """

    rho_zl = eos.molar_density(T, P, z, "liq")
    fugcoef_zl = eos.fugacity_coefficients(T, P, z, "liq")  # Assume feed as liquid phase

    rho_zv = eos.molar_density(T, P, z, "vap")
    fugcoef_zv = eos.fugacity_coefficients(T, P, z, "vap")  # Assume feed as liquid phase

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
        rho_l = eos.molar_density(T, P, WL, "liq")
        fugcoef_l = eos.fugacity_coefficients(T, P, WL, "liq")
        WL_new = np.exp(dV-np.log(fugcoef_l))
        dW = WL_new-WL
        WL=WL_new/sum(WL_new)
    rho_l = eos.molar_density(T, P, WL, "liq")
    fugcoef_l = eos.fugacity_coefficients(T, P, WL, "liq")
    tm_l = 1+sum(WL*(np.log(WL)+np.log(fugcoef_l)-dV-1))

    dW = np.ones(len(z))
    iter_v = 0
    while (max(abs(dW))>1e-7) & (iter_v<max_iter):
        iter_v=iter_v+1
        rho_v = eos.molar_density(T, P, WV, "vap")
        fugcoef_v = eos.fugacity_coefficients(T, P, WV, "vap")
        WV_new = np.exp(dL-np.log(fugcoef_v))
        dW = WV_new-WV
        WV = WV_new/sum(WV_new)
    rho_v = eos.molar_density(T, P, WV, "vap")
    fugcoef_v = eos.fugacity_coefficients(T, P, WV, "vap")
    tm_v = 1+sum(WV*(np.log(WV)+np.log(fugcoef_v)-dL-1)) 
    return tm_l, tm_v, fugcoef_l, fugcoef_v, fugcoef_zl, fugcoef_zv


if __name__ == "__main__":
    # Flash test 1: "C3H6", "C3H8", "H2", "N2"
    print("测试程序: ")

    T = 341.15
    P = 3000000
    z = np.array([0.82, 0.16, 0.005, 0.015])
    args = {"Tc": np.array([364.85, 369.83, 33.19, 126.2]), 
            "Pc": np.array([4600000, 4248000, 1313000, 3400000]), 
            "w": np.array([0.137588, 0.152291, -0.215993, 0.0377215]),
            "m": np.array([1.9597, 2.002, 0.9846, 1.2053]),
            "s": np.array([3.5356, 3.6184, 2.8263, 3.313]),
            "e": np.array([207.19, 208.11, 20.893, 90.96])}
    method = "PC-SAFT"
    eos = PropertyMethod(method, args)
    
    beta, x, y, iter = VLE_TPflash(T, P, z, eos)
    print("Vapor fraction: ", beta)
    print("Liquid composition: ", x)
    print("Vapor composition: ", y)
    print("iter: ", iter)

    print("\n",80*"-","\n")

    # Flash test 2: "C3H6", "C3H8", "H2", "N2", "PP", "TICL4", "TEA"
    T = 338.15
    P = 500000
    PP_DPN = 535.4335 
    PP_MWN = 42.0806*PP_DPN
    z = np.array([0.2098286, 0.00408862, 0.0005843, 0.000137371, 0.7851512/PP_DPN,
                    3.22413E-05, 0.000177694])
    # z = np.array([0.0251157, 0.000495333, 1.33715E-06, 7.41334E-07, 0.9741264/PP_DPN, 4.00013E-05, 0.000220463])

    z = z/np.sum(z)
    print("z", z)
    args = {"Tc": np.array([364.85, 369.83, 33.19, 126.2, 2000, 638, 678.15]), 
            "Pc": np.array([4600000, 4248000, 1313000, 3400000, 5000000, 4660950, 8930000]), 
            "w": np.array([0.137588, 0.152291, -0.215993, 0.0377215, 0, 0.283732, 0.841783]),
            "m": np.array([1.9597, 2.002, 0.8285, 1.2053, 0.02528*PP_MWN, 20, 20]),
            "s": np.array([3.5356, 3.6184, 2.9729, 3.313, 4.1473, 3.5356, 3.5356]),
            "e": np.array([207.19, 208.11, 12.53, 90.96, 298.6392, 207.19, 207.19])}
    method = "PC-SAFT"
    eos = PropertyMethod(method, args)

    beta, x, y, iter = VLE_TPflash(T, P, z, eos)
    print("Vapor fraction: ", beta)
    print("Liquid composition: ", x)
    x[4] = x[4]*PP_DPN
    beta = beta/(beta+(1-beta)*sum(x))
    print("新气相分数：", beta)
    x = x/sum(x)
    print("新液相组成：", x)

    print("Vapor composition: ", y)
    print("iter: ", iter)

    print(50*"*")

    T = 333.2672
    P = 3000000
    z = np.array([0.8847128, 0.0166278, 0.0903637, 0.0082909, 7.20234E-07, 3.98896E-06])

    pyargs = {
    'Tc': np.array([364.85, 369.83,  33.19, 126.2 , 638.  , 678.15]), 
    'Pc': np.array([4600000., 4248000., 1313000., 3400000., 4660950., 8930000.]), 
    'w': np.array([ 0.137588 ,  0.152291 , -0.215993 ,  0.0377215,  0.283732 ,0.841783 ]), 
    'm': np.array([ 1.9597,  2.002 ,  0.8285,  1.2053, 20.    , 20.    ]), 
    's': np.array([3.5356, 3.6184, 2.9729, 3.313 , 3.5356, 3.5356]), 
    'e': np.array([207.19, 208.11,  12.53,  90.96, 207.19, 207.19])}

    method = "PC-SAFT"
    eos = PropertyMethod(method, pyargs)
    beta, x, y, iter = VLE_TPflash_first(T, P, z, eos)
    print("beta:", beta)
    print(iter)

    # 露点压力计算测试
    print(50*"*")
    T = 100+273.15
    P = 100000
    z = np.array([1])
    pyargs = {
        "w": np.array([0.344861])  ,
        "Tc": np.array([647.096]) ,
        "Pc": np.array([22064000]),
        "mw": np.array([18.01528]),
        "m": np.array([1.0656]),   
        "s": np.array([3.0007]),
        "e": np.array([366.51]),
        # "e_assoc": np.array([2500.7]),
        # "vol_a": np.array([0.034868])
    }
    method = "PC-SAFT"
    eos = PropertyMethod(method, pyargs)
    print(bubble_pressure(T, z, eos))
    print(bubble_temperature(P, z, eos))
    print(dew_pressure(T, z, eos))
    print(dew_temperature(P, z, eos))

    T = 338.15
    P = 500000
    z = np.array([8.14135733e-01, 1.65508367e-01, 5.78011953e-03, 1.17596099e-03,
       1.33856025e-02, 2.17440371e-06, 1.20425133e-05, 0.00000000e+00])

    pyargs = {'Tc': np.array([ 364.85 ,  369.83 ,   33.19 ,  126.2  , 2000.   ,  638.   ,
    678.15 ,  647.096]), 'Pc': np.array([ 4600000.,  4248000.,  1313000.,  3400000.,  5000000.,  4660950.,
    8930000., 22064000.]), 'w': np.array([ 0.137588 ,  0.152291 , -0.215993 ,  0.0377215,  0.       ,
    0.283732 ,  0.841783 ,  0.344861 ]), 'm': np.array([1.95970000e+00, 2.00200000e+00, 8.28500000e-01, 1.20530000e+00,
    2.15394661e+03, 2.00000000e+01, 2.00000000e+01, 1.06560000e+00]), 's': np.array([3.5356, 3.6184, 2.9729, 3.313 , 4.1473, 3.5356, 3.5356, 3.0007]), 'e': np.array([207.19  , 208.11  ,  12.53  ,  90.96  , 298.6392, 207.19  ,
    207.19  , 366.51  ]), 'e_assoc': np.array([   0. ,    0. ,    0. ,    0. ,    0. ,    0. ,    0. , 2500.7]), 'vol_a': np.array([0.      , 0.      , 0.      , 0.      , 0.      , 0.      ,
    0.      , 0.034868])}
    method = "PC-SAFT"
    eos = PropertyMethod(method, pyargs)
    DPN = 2024.781

    z[4] = z[4]/DPN
    z = z/sum(z)

    beta, x, y, iter = VLE_TPflash(T, P, z, eos, max_iter=300)
    x[4] = x[4]*DPN
    beta = beta/(beta+(1-beta)*sum(x))
    print("新气相分数：", beta)
    x = x/sum(x)
    print("新液相组成：", x)

    print("Vapor composition: ", y)
    print("iter: ", iter)







