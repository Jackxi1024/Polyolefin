"""
求解器程序
此程序功能是求解非线性方程组F(X) = X, 可选择直接迭代法或者wegstein法,
该求解程序可用于求解CSTR模型、撕裂流股等。

此程序并没有实现F(X)=0的通用算法(例如Newton法、Broyden法等), 因为这类问题可以调用
Scipy包来计算, 工具包更加成熟。

"""

import numpy as np
from polymer_model.utility import InputError

def solve(fun, x0, args=(), method='wegstein', max_iter=1000, tol=1e-6, w_min=0, w_max=5):
    """
    Find a root of the vector equations: fun(x) = x

    Parameters
    ----------
    fun : callable
        A vector function to find a root of.
    x0 : ndarray
        Initial guess.
    args : tuple, optional
        Extra arguments passed to the objective function
    method : str, optional
        Type of solver. Should be one of

        'direct'    : direct iteration method, x(k+1) = fun(x(k))        
        'wegstein'  : Wegstein method, see reference[1]
    max_iter : int
        Maximum number of iterations
    tol : float
        Tolerance for termination.
    w_min: float
        Lower bound of relaxation factor in wegstein method
    w_max: float
        Upper bound of relaxation factor in wegstein method

    Returns
    -------
    sol : dict
        The solution represented as a dict. The keyword:
        "success" : a Boolean flag indicating if the algorithm exited successfully
        "x" : the solution array
        "message" : which describes the cause of the termination

    References
    ----------
    [1] 


    Examples
    --------
    The following functions define a system of nonlinear equations

    >>> def fun(x):
    ...     return [x[0]  + 0.5 * (x[0] - x[1])**3 - 1.0]

    A solution can be obtained as follows.

    >>> from solver import solve
    >>> sol = solve(fun, [0], method='wegstein')
    >>> sol.x

    """
    
    if not isinstance(args, tuple):
        args = (args,)

    if method == "direct":
        return direct_iterate(fun, x0, args, max_iter, tol)
    elif method == "wegstein":
        return wegstein(fun, x0, args, max_iter, tol, w_min, w_max)
    else:
        raise InputError("Please select direct method or wegstein method")


def direct_iterate(fun, x0, args=(), max_iter=1000, tol=1e-7):
    """ Direct Iteration Method """
    x = x0
    iter = 0
    while iter < max_iter:
        iter += 1
        f = fun(x, *args)
        step = f - x

        error = np.max(np.abs(f-x))
        if error <= tol:
            sol = {"success": True, "x": x, "iter": iter}
            return sol
        else:
            x = x + step
    sol = {"success": False, "x": x, "iter": iter}
    return sol


def wegstein(fun, x0, args=(), max_iter=1000, tol=1e-7, w_min=0, w_max=5):
    """ Wegstein Method """

    x1 = x0
    f1 = fun(x1, *args)
    while np.any(f1<0):
        for i in np.where(f1<0):
            x1[i] = 0.5*x1[i]
        f1 = fun(x1, *args)    
        
    x2 = f1
    f2 = fun(x2, *args)
    while np.any(f2<0):
        for i in np.where(f2<0):
            x2[i] = 0.5*x2[i]
        f2 = fun(x2, *args)  

    iter = 0
    while iter < max_iter:
        # 计算误差
        error = np.max(np.abs(f2-x2))
        if error <= tol:
            sol = {"success": True, "x": x2, "iter": iter}
            return sol

        iter += 1
        s = np.array([0 if x2[i] == x1[i] else (f2[i]-f1[i])/(x2[i]-x1[i]) for i in range(0, x0.size)])
        s[s==1] = 0
        w = 1/1-s
        np.clip(w, w_min, w_max)
        step = w*(f2-x2)
        
        # 更新x1, x2, f1, f2
        x1 = x2
        f1 = f2
        x2 = x2 + step
        f2 = fun(x2, *args)
        while np.any(f2<0):
            for i in np.where(f2<0):
                x2[i] = 0.5*x2[i]
            f2 = fun(x2, *args)  
            
    sol = {"success": False, "x": x2, "iter": iter}
    return sol


if __name__ == "__main__":

    # 测试求解器, 求解方程组
    # x1 = a*sin(x1) + b*cos(x2)   
    # x2 = a*cos(x1) - b*sin(x2)

    from math import sin, cos

    def fun(x, a, b):
        f1 = a*sin(x[0]) + b*cos(x[1])   
        f2 = a*cos(x[0]) - b*sin(x[1])
        return np.array([f1, f2]) 

    args =(0.7, 0.2)
    x0 = np.array([-0.1,0.1])
    sol = solve(fun, x0, args)
    print(sol)
    x = sol["x"]
    print("x: ", x)
    print(type(x))
    y = fun(sol["x"], *args)
    print("y: ", y)

    print()