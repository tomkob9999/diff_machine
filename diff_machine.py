#
# diff_machine
#
# Description: Calculates difference equation based on the order, target row and initial values specified
#
# Version: 1.1.2
# Author: Tomio Kobayashi
# Last Update: 2024/9/9

import numpy as np
import math

# Virtual Array that keeps only a fixed number of array from the highest
class varr:
    def __init__(self, size):
        self.size = size
        self.ar = [0]*size
        self.max = size-1
        self.curr = 0
    def get(self, i):
        return self.ar[self.size - (self.max - i) - 1]
    def get_curr(self):
        return self.ar[min(self.curr, self.size-1)]
    def set(self, i, y):
        if i > self.curr:
            self.curr = i
        if i > self.max:
            if i == self.max+1:
                for ii in range(self.size-1):
                    self.ar[ii] = self.ar[ii+1]
                self.ar[self.size-1] = y
                self.max += 1
            else:
                raise Exception("Tried to increment by more than 1!")
        else:
            self.ar[self.size - (self.max - i) - 1] = y
            
    
class diff_machine:
    def __init__(self, order):
        self.memo = {}
        self.order = order
        
    def clean_memo(self):
        self.memo = {}

    def calc(x1, x2, exp=False):
        return x1-x2 if not exp else x1/x2

    def get_diff(self, ar, i, order, order_exp=[], force=0, print_force=False, enable_memo=True):
        if enable_memo and i in self.memo and order in self.memo[i]:
            pass
        elif order == 1:
            if i not in self.memo:
                self.memo[i] = {}
            self.memo[i][order] =  diff_machine.calc(ar[i-1], ar[i-2], order in order_exp)
        else:
            if i not in self.memo:
                self.memo[i] = {}
            if force > 0 and self.order == order:
                self.memo[i][order] = force
            else:
                self.memo[i][order] = diff_machine.calc(self.get_diff(ar, i, order-1, order_exp, enable_memo), self.get_diff(ar, i-1, order-1, order_exp, enable_memo), order in order_exp)
            if print_force and self.order == order and order+2 == i:
                print("force", self.memo[i][order])
                
        return self.memo[i][order]

    def get_diff2(self, ar, i, order, order_exp=[], force=0, print_force=False, enable_memo=True):
        if enable_memo and i in self.memo and order in self.memo[i]:
            pass
        elif order == 1:
            if i not in self.memo:
                self.memo[i] = {}
            self.memo[i][order] = diff_machine.calc(ar.get(i-1), ar.get(i-2), order in order_exp)
        else:
            if i not in self.memo:
                self.memo[i] = {}
            if force > 0 and self.order == order:
                self.memo[i][order] = force
            else:
                self.memo[i][order] = diff_machine.calc(self.get_diff2(ar, i, order-1, order_exp, enable_memo), self.get_diff2(ar, i-1, order-1, order_exp, enable_memo), order in order_exp)
            if print_force and self.order == order and order+2 == i:
                print("force", self.memo[i][order])
                
        return self.memo[i][order]

    # Returns array
    # contains only constants and differences in the calculations
    
    # TIPS:
    #     For y=b^x*a
    #       y(0)->a
    #       y(1)/y(0)->b
    # order=[1], init={a, b/a}
    def solve_array(target, init, order_exp=[], force=0, print_force=False):
        order = len(init)-1
        ar = np.zeros(target+1)
        dd = diff_machine(order)
#         for k, v in init.items():
        for k, v in enumerate(init):
            ar[k] = v
        for i in range(order+1, target+1, 1):
            order_cum = 0
            for ii in range(order, 0, -1):
                if ii+1 in order_exp:
                    if order_cum == 0:
                        order_cum = 1
                    order_cum *= dd.get_diff(ar, i, ii, order_exp, force, print_force)
                else:
                    order_cum += dd.get_diff(ar, i, ii, order_exp, force, print_force)
            if 1 in order_exp:
                ar[i] = ar[i-1] * order_cum
            else:
                ar[i] = ar[i-1] + order_cum
            if i - ((order+2)*3) in dd.memo:
                del dd.memo[i - ((order+2)*3)]
        return ar
    
    

    def solve(target, init, order_exp=[], force=0, print_force=False):
        order = len(init)-1
        dd = diff_machine(order)
        va = varr((order+2)*2)
#         for k, v in init.items():
        for k, v in enumerate(init):
            va.set(k, v)
        z = 0
        for i in range(order+1, target+1, 1):
            order_cum = 0
            for ii in range(order, 0, -1):
                if ii+1 in order_exp:
                    if order_cum == 0:
                        order_cum = 1
                    order_cum *= dd.get_diff2(va, i, ii, order_exp, force, print_force)
                else:
                    order_cum += dd.get_diff2(va, i, ii, order_exp, force, print_force)
            if 1 in order_exp:
                va.set(i, va.get(i-1) * order_cum)
            else:
                va.set(i, va.get(i-1) + order_cum)
            if i - ((order+2)*3) in dd.memo:
                del dd.memo[i - ((order+2)*3)]
        return va.get_curr()

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
# res = diff_machine.solve(4, {[11, 21, 34])
# print("res", res)
# res = diff_machine.solve_array(4, [11, 21, 34])
# print("res", res)
# res = diff_machine.solve_array2(4, [11, 21, 34])
# print("res", res)
# # Same as above except step=0.01
# res = diff_machine.solve(10000, [0.011, 0.021, 0.034])
# print("res", res)
# #
# # Closed form: y=x^2 (1 step=1)
# Difference equation: y'=y''+y''', y(0)=1, y(1)=1, y(2)=4
# ar = diff_machine.solve_array(50, [0, 1, 4])
# print("ar", ar)
# ar = diff_machine.solve_array(50, [-100, 1, 4], force=2.0)
# print("ar", ar)
ar = diff_machine.solve(10, [0, 1, 4])
print("ar", ar)
# ar = diff_machine.solve(5, [-1, 1, 4], force=2.0)
# ar = diff_machine.solve(10, [-1, 1, 4], force=2.2)
ar = diff_machine.solve(10, [-1, 1, 4], force=2.2)
# ar = diff_machine.solve(50, [1, 4], force=2.0)
print("ar", ar)
# #
# # Closed form: y=4x^3+3x^2
# res = diff_machine.solve_array(4, [0, 7, 44, 135])
# # res = diff_machine.solve_array(5, [0, 0.07, 0.44, 1.35])
# # res = diff_machine.solve_array(5, [0, 0.034, 0.152, 0.378])
# print("res", res)

# # Closed form: y=5x^5+4x^4+3x^3+2x^2+1
import time
start_time = time.time()
res = diff_machine.solve_array(1000, [0, 0.12345, 0.312, 0.60555, 1.0656, 1.78125])
air_time = time.time() - start_time
print("res", res[-1])
print(f"Execution Time: {air_time:.6f} seconds")
start_time = time.time()
res = diff_machine.solve_array(1000, [-100, 0.12345, 0.312, 0.60555, 1.0656, 1.78125], force=0.006)
# res = diff_machine.solve_array(10000, [0, 0.12345, 0.312, 0.60555, 1.0656, 1.78125])
print("res", res[-1])
air_time = time.time() - start_time
print(f"Execution Time: {air_time:.6f} seconds")
start_time = time.time()
res = diff_machine.solve(1000, [0, 0.12345, 0.312, 0.60555, 1.0656, 1.78125])
# res = diff_machine.solve_array(10000, [0, 0.12345, 0.312, 0.60555, 1.0656, 1.78125])
print("res", res)
air_time = time.time() - start_time
print(f"Execution Time: {air_time:.6f} seconds")
start_time = time.time()
res = diff_machine.solve(1000, [-100, 0.12345, 0.312, 0.60555, 1.0656, 1.78125], force=0.006)
# res = diff_machine.solve_array(10000, [0, 0.12345, 0.312, 0.60555, 1.0656, 1.78125])
print("res", res)
air_time = time.time() - start_time
print(f"Execution Time: {air_time:.6f} seconds")

# # Closed form: y=2^x (1 step=1)
# # Difference equation: y(x)=y(x-1)**2/y(x-2), y(0)=1, y(1)=2
# ar = diff_machine.solve(40, [1, 2], order_exp=[1])
# print("ar", ar)
# ar = diff_machine.solve_array(100, [1, 1.01], order_exp=[1])
# print("ar", ar)
# ar = diff_machine.solve(100, [1, 1.01], order_exp=[1])
# print("ar", ar)
# #
# # Closed form: 4*2^(x*3) (1 step=1)
# ar = diff_machine.solve(10, [4, 32}], order_exp=[1])
# print("ar", ar)
# ar = diff_machine.solve_array(10, [4, 32], order_exp=[1])
# print("ar", ar)
# # 
# # Closed form: y(x)=2^(x-1)+(x-1)
# ar = diff_machine.solve_array(10, [1, 3, 6, 11], order_exp=[3])
# print("ar", ar)
# ar = diff_machine.solve(10, [1, 3, 6, 11], order_exp=[3])
# print("ar", ar)