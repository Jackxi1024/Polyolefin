"""
公共程序
"""

R = 8.314   # 气体常数

class InputError(Exception):
    """Exception raised for errors in the input.
    """
    def __init__(self, message):
        self.message = message

class SolutionError(Exception):
    # Exception raised when a solver does not return a value.
    def __init__(self, message):
        self.message = message

        
# 希腊字母
α, β, γ, δ, ε, ζ, η, κ, μ, ρ, σ, φ, ψ, λ, ξ, π,  = [0]*16
τ = 1 





