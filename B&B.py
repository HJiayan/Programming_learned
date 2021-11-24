# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 09:27:43 2021

@author: IMGDer
"""

from scipy.optimize import linprog
import numpy as np
import sys
from queue import Queue

class BB(object):
    def __init__(self,c,A,b,Aeq,beq,bounds):
        
        self.c,self.A_eq,self.b_eq,self.bounds = c,Aeq,beq,bounds
        self.Q = Queue()
        self.opt = None
        self.lower_bound = -sys.maxsize
        self.upper_bound = sys.maxsize
        self.best_x = None 
        #解初始问题
        res = linprog(c,A,b,Aeq,beq,bounds)
        
        if not res.success:
            raise ValueError("No feasible solution")
        self.Q.put((res,A,b))
    
    def branch_and_bound(self):
        
        while not self.Q.empty():
            res,A,b = self.Q.get(block=False)   
            #Bound.分支的解是原问题的上界，是分支的下界（最小化问题中）
            if res.fun >= self.upper_bound:
                continue
            
            X = res.x
            variable_list = []
            for i in X:
                if not i.is_integer():
                    variable_list.append(i)
            #if all integer
            if len(variable_list) == 0:
                if res.fun <= self.upper_bound:
                    self.opt = res.fun
                    self.upper_bound = res.fun  #原问题上界更新，即原问题不会大于这个解
                    self.best_x = res.x
            else:
                #find 1st fraction
                index = 0
                for i in X:
                    if not i.is_integer():
                        variable = i
                        break
                    index += 1

                #Branch
                new_constrain_1 = [0] * len(A[0])
                new_constrain_2 = [0] * len(A[0])
                new_constrain_1[index] = 1
                new_constrain_2[index] = -1
                new_A_1 = np.insert(A, len(A), values = new_constrain_1, axis=0)
                new_A_2 = np.insert(A, len(A), values = new_constrain_2, axis=0)
                new_b_1 = np.insert(b, len(b), values = np.floor(variable), axis=0)
                new_b_2 = np.insert(b, len(b), values = -np.ceil(variable), axis=0)
                
                res_1 = linprog(self.c,new_A_1,new_b_1,self.A_eq,self.b_eq,self.bounds)
                res_2 = linprog(self.c,new_A_2,new_b_2,self.A_eq,self.b_eq,self.bounds)
                
                if res_1.success and res_2.success:
                    self.Q.put((res_1,new_A_1,new_b_1))
                    self.Q.put((res_2,new_A_2,new_b_2))
                elif res_1.success and not res_2.success:
                    self.Q.put((res_1,new_A_1,new_b_1))
                elif not res_1.success and res_2.success:
                    self.Q.put((res_2,new_A_2,new_b_2))

def test1():
    """ 此测试的真实最优解为 [4, 2] """
    c = np.array([-40, -90])
    A = np.array([[9, 7], [7, 20]])
    b = np.array([56, 70])
    Aeq = None
    beq = None
    bounds = [(0, None), (0, None)]
 
    solver = BB(c, A, b, Aeq, beq, bounds)
    solver.branch_and_bound()
 
    print("Test 1's result:", solver.opt, solver.best_x)
    #print("Test 1's true optimal x: [4, 2]\n")

def test2():
    c = np.array([-4,-5])
    A = np.array([[1,4],[3,-4]])
    b = np.array([10,6])
    Aeq = None
    beq = None
    bounds = [(0, None), (0, None)]
    
    solver = BB(c, A, b, Aeq, beq, bounds)
    solver.branch_and_bound()
 
    print("Test 1's result:", solver.opt, solver.best_x)
    print("Test 1's true optimal x: [4, 2]\n")
            
if __name__ == '__main__':
    test1()
    #test2()                                                            
                
                
                
                
                    
            
            
            
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        