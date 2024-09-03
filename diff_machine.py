#
# diff_machine
#
# Description: Calculates difference equation based on the order, target row and initial values specified
# Version: 1.0.0
# Author: Tomio Kobayashi
# Last Update: 2024/9/3

import numpy as np

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
        for k, v in init.items():
            ar[k] = v
        for i in range(order+1, target+1, 1):
            order_cum = 0
            for ii in range(order, 0, -1):
                order_cum += dd.get_diff(ar, i, ii)
            ar[i] = ar[i-1] + order_cum + (np.sum([diff_coefs[m] * (i*unit-unit) ** (m-1) if m != 0 else 0 for m in range(len(diff_coefs))]))*unit**(len(coefs)-1)
        return ar

    
    # Returns array
    def solve_mix(target, init, coefs, unit=1):
        return diff_machine.solve_mix_array(target, init, coefs, unit)[-1]

    
# # res = diff_machine.solve_pure(21, {0:11, 1:21, 2:34})
# # res = diff_machine.solve_pure(21, {0:0.011, 1:0.021, 2:0.034})
# res = diff_machine.solve_pure(9, {0:0, 1:1, 2:4})
# # res = diff_machine.solve_pure(5, {0:0, 1:1, 2:8, 3:27})
# # res = diff_machine.solve_pure(88, {0:0, 1:1, 2:16, 3:81, 4:256})
# print("res", res)
# ar = diff_machine.solve_pure_array(4, {0:0, 1:1, 2:4})
# ar = diff_machine.solve_pure_array(40, {0:1, 1:2}, order_exp=[1])
# print("ar", ar)
# ar = diff_machine.solve_pure(40, {0:1, 1:2}, order_exp=[1])
# print("ar", ar)
# ar = diff_machine.solve_pure_array(100, {0:1, 1:1.01}, order_exp=[1])
# print("ar", ar)
# ar = diff_machine.solve_pure(100, {0:1, 1:1.01}, order_exp=[1])
# print("ar", ar)
ar = diff_machine.solve_pure_array(10, {0:1, 1:3, 2:6, 3:11}, order_exp=[3])
print("ar", ar)
ar = diff_machine.solve_pure(10, {0:1, 1:3, 2:6, 3:11}, order_exp=[3])
print("ar", ar)

# diff_machine.solve_mix_array(10, {0:0, 1:7}, [0, 3, 4], unit=1)
# res = diff_machine.solve_mix_array(11, {0:0, 1:0.034}, [0, 3, 4], unit=0.1)
# print("res", res)
# res = diff_machine.solve_mix(11, {0:0, 1:0.034}, [0, 3, 4], unit=0.1)
# print("res", res)