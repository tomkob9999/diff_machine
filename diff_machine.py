#
# diff_machine
#
# Description: Calculates difference equation based on the order, target row and initial values specified
# Version: 1.0.2
# Author: Tomio Kobayashi
# Last Update: 2024/9/5

import numpy as np
import math

class diff_machine:
    def __init__(self):
        self.memo = {}
        self.memo_found = 0
        
    def clean_memo(self):
        self.memo = {}

    def calc(x1, x2, exp=False):
        return x1-x2 if not exp else x1/x2
    
    def get_diff(self, ar, i, order, order_exp=[], enable_memo=True):
        if enable_memo and (i, order) in self.memo:
            self.memo_found += 1
            return self.memo[(i, order)]
        
        if order == 1:
            return diff_machine.calc(ar[i-1], ar[i-2], order in order_exp)
        else:
            ret = diff_machine.calc(self.get_diff(ar, i, order-1, order_exp, enable_memo), self.get_diff(ar, i-1, order-1, order_exp, enable_memo), order in order_exp)
            self.memo[(i, order)] = ret
            return ret
    
    def diff_coef(n, order):
        if order==0:
            return n
        return diff_machine.diff_coef(n * order, order-1)

    # Returns array
    # pure - contains only constants and differences in the calculations
    # mix - contains both difference fields and value of x in the calculation
    
    # TIPS:
    #     For y=b^x*a
    #       y(0)->a
    #       y(1)/y(0)->b
    # order=[1], init={a, b/a}
    def solve_pure_array(target, init, order_exp=[]):
        order = len(init)-1
        ar = np.zeros(target+1)
        dd = diff_machine()
        for k, v in init.items():
            ar[k] = v
        for i in range(order+1, target+1, 1):
            order_cum = 0
            for ii in range(order, 0, -1):
                if ii+1 in order_exp:
                    if order_cum == 0:
                        order_cum = 1
                    order_cum *= dd.get_diff(ar, i, ii, order_exp)
                else:
                    order_cum += dd.get_diff(ar, i, ii, order_exp)
            if 1 in order_exp:
                ar[i] = ar[i-1] * order_cum
            else:
                ar[i] = ar[i-1] + order_cum
        return ar
    
    # Returns value
    def solve_pure(target, init, order_exp=[]):
        return diff_machine.solve_pure_array(target, init, order_exp)[-1]
    
    # Returns array
    def solve_mix_array(target, init, coefs, unit=1):
        diff_coefs = []
        for i in range(len(coefs)):
            diff_coefs.append(diff_machine.diff_coef(coefs[i], (i+1)))

        order = len(init)-1
        ar = np.zeros(target+1)
        dd = diff_machine()
        print("order", order)
        fact1 = math.factorial(order+2)
        fact2 = math.factorial(order+1)
        for k, v in init.items():
            ar[k] = v
        for i in range(order+1, target+1, 1):
            order_cum = 0
            for ii in range(order, 0, -1):
                order_cum += dd.get_diff(ar, i, ii)
            ar[i] = ar[i-1] + order_cum + coefs[0]*(i-1)*fact1 + (coefs[1]-(order**2-order))*fact2
        return ar

    
    # Returns array
    def solve_mix(target, init, coefs, unit=1):
        return diff_machine.solve_mix_array(target, init, coefs, unit)[-1]

# Notations throughout
#
# Difference variables are expressed like differential equation notations
# y' = Y(x)-y(x-1)
# y'' = Y'(x)-y'(x-1) 
# y''' = Y''(x)-y''(x-1) 
# ...
# How to set initial values for pure form
# for i=0 to (order-1), y(i) = original_funct()
#
# How to set initial values for mix form
# for i=0 to (order-3), y(i) = original_funct()
# How to set coef values for mix form - internally factorials are used to derive coefficients of derivatives
# for i=order to order-1, y(i) = coef in the original polynomial
#
# Closed form: y=3*x^2+7*x+11 (1 step=1)
# Difference equation: y'=y''+y''', y(0)=11, y(1)=21, y(2)=34
# res = diff_machine.solve_pure(4, {0:11, 1:21, 2:34})
# print("res", res)
# Same as above except step=0.01
# # res = diff_machine.solve_pure(21, {0:0.011, 1:0.021, 2:0.034})
#
# Closed form: y=x^2 (1 step=1)
# Difference equation: y'=y''+y''', y(0)=1, y(1)=1, y(2)=4
# res = diff_machine.solve_pure(9, {0:0, 1:1, 2:4})
# ar = diff_machine.solve_pure_array(4, {0:0, 1:1, 2:4})
#
# Closed form: y=4x^3+3x^2
# Mix Difference Equation: y''=4*6*(x-dx)+6
res = diff_machine.solve_pure_array(4, {0:0, 1:7, 2:44, 3:135})
print("res", res)
# res = diff_machine.solve_mix_array(10, {0:0, 1:7}, [0, 3, 4], unit=1)
res = diff_machine.solve_mix_array(4, {0:0, 1:7}, [4, 3], unit=1)
print("res", res)

# Closed form: y=x^4+3x^3
res = diff_machine.solve_pure_array(5, {0:0, 1:4, 2:40, 3:162, 4:448})
print("res", res)
res = diff_machine.solve_mix_array(5, {0:0, 1:4, 2:40}, [1, 3], unit=1)
print("res", res)

# Closed form: y=2^x (1 step=1)
# Difference equation: y(x)=y(x-1)**2/y(x-2), y(0)=1, y(1)=2
# ar = diff_machine.solve_pure(40, {0:1, 1:2}, order_exp=[1])
# print("ar", ar)
# ar = diff_machine.solve_pure_array(100, {0:1, 1:1.01}, order_exp=[1])
# print("ar", ar)
# ar = diff_machine.solve_pure(100, {0:1, 1:1.01}, order_exp=[1])
# print("ar", ar)
#
# 
# Closed form: y(x)=2^(x-1)+(x-1)
# Difference Equation: Unknown but order calculation is as 1->minus, 2->minus, 3->division
# ar = diff_machine.solve_pure_array(10, {0:1, 1:3, 2:6, 3:11}, order_exp=[3])
# print("ar", ar)
# ar = diff_machine.solve_pure(10, {0:1, 1:3, 2:6, 3:11}, order_exp=[3])
# print("ar", ar)