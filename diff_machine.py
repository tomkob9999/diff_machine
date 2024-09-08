#
# diff_machine
#
# Description: Calculates difference equation based on the order, target row and initial values specified
#
# Version: 1.0.5
# Author: Tomio Kobayashi
# Last Update: 2024/9/8

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
    
    # TIPS:
    #     For y=b^x*a
    #       y(0)->a
    #       y(1)/y(0)->b
    # order=[1], init={a, b/a}
    def solve_array(target, init, order_exp=[]):
        order = len(init)-1
        ar = np.zeros(target+1)
        dd = diff_machine()
        for k, v in init.items():
            ar[k] = v
        for i in range(order+1, target+1, 1):
            order_cum = 0
            for ii in range(order, 0, -1):
#                 print(dd.get_diff(ar, i, ii, order_exp))
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
    def solve(target, init, order_exp=[]):
        return diff_machine.solve_array(target, init, order_exp)[-1]
    
# Notations throughout
#
# Difference variables are expressed like differential equation notations
# y' = y(x)-y(x-1)
# y'' = y'(x)-y'(x-1) 
# y''' = y''(x)-y''(x-1) 
# ...
# How to set initial values for pure form
# for i=0 to (order-1), y(i) = original_funct()
#
# Closed form: y=3*x^2+7*x+11 (1 step=1)
# Difference equation: y'=y''+y''', y(0)=11, y(1)=21, y(2)=34
res = diff_machine.solve(4, {0:11, 1:21, 2:34})
print("res", res)
# Same as above except step=0.01
res = diff_machine.solve(21, {0:0.011, 1:0.021, 2:0.034})
print("res", res)
#
# Closed form: y=x^2 (1 step=1)
# Difference equation: y'=y''+y''', y(0)=1, y(1)=1, y(2)=4
# res = diff_machine.solve_pure(9, {0:0, 1:1, 2:4})
ar = diff_machine.solve_array(4, {0:0, 1:1, 2:4})
print("ar", ar)
#
# Closed form: y=4x^3+3x^2
res = diff_machine.solve_array(4, {0:0, 1:7, 2:44, 3:135})
# res = diff_machine.solve_array(5, {0:0, 1:0.07, 2:0.44, 3:1.35})
# res = diff_machine.solve_array(5, {0:0, 1:0.034, 2:0.152, 3:0.378})
print("res", res)

# Closed form: y=5x^5+4x^4+3x^3+2x^2+1
# res = diff_machine.solve_array(6, {0:0, 1:15, 2:258, 3:1641, 4:6372, 5:18555})
res = diff_machine.solve_array(6, {0:0, 1:0.12345, 2:0.312, 3:0.60555, 4:1.0656, 5:1.78125})
print("res", res)

# Closed form: y=2^x (1 step=1)
# Difference equation: y(x)=y(x-1)**2/y(x-2), y(0)=1, y(1)=2
ar = diff_machine.solve(40, {0:1, 1:2}, order_exp=[1])
print("ar", ar)
ar = diff_machine.solve_array(100, {0:1, 1:1.01}, order_exp=[1])
print("ar", ar)
ar = diff_machine.solve(100, {0:1, 1:1.01}, order_exp=[1])
print("ar", ar)
#
# Closed form: 4*2^(x*3) (1 step=1)
ar = diff_machine.solve(10, {0:4, 1:32}, order_exp=[1])
print("ar", ar)
ar = diff_machine.solve_array(10, {0:4, 1:32}, order_exp=[1])
print("ar", ar)
# 
# Closed form: y(x)=2^(x-1)+(x-1)
ar = diff_machine.solve_array(10, {0:1, 1:3, 2:6, 3:11}, order_exp=[3])
print("ar", ar)
ar = diff_machine.solve(10, {0:1, 1:3, 2:6, 3:11}, order_exp=[3])
print("ar", ar)