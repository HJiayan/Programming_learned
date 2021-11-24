# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 14:01:12 2021

@author: GDer
"""

import numpy as np
import sys

class BC():
    def __init__(self,c,A,b,Agq,bgq,Aeq,beq,Bounds):
        self.c = c
        self.A = A
        self.b = b
        self.Agq = Agq
        self.bgq = bgq
        self.Aeq = Aeq
        self.beq = beq
        self.Bounds = Bounds
        self.frac = 3
    
    def standardize(self):
        #初始不等式
        standard_c = self.c
        standard_A = self.A
        
        #标准化b，先b，再bgq,在beq
        if np.any(self.bgq):
            standard_b = np.append(self.b,self.bgq)
        else:
            standard_b = self.b
            
        if np.any(self.beq):
            standard_b = np.append(standard_b,self.beq)
        else:
            standard_b = standard_b
        
        #标准化A,先A，再Agq,再Aeq
        if np.any(self.Agq):
            standard_A = np.insert(self.A,self.A.shape[0],self.Agq,axis=0)
        else:
            standard_A = self.A
        
        identity_matrix = np.identity(standard_A.shape[0])
        standard_A = np.insert(standard_A,standard_A.shape[1],identity_matrix,axis=1)
        for i in range(self.Agq.shape[0]):
            for j in range(self.Agq.shape[0]):
                standard_A[-(i+1)][-(j+1)] = -1
        
        if np.any(self.Aeq):
            line_Aeq = self.Aeq.shape[0]
            column_Aeq = self.Aeq.shape[1]
            new_array = np.zeros((line_Aeq,standard_A.shape[1]-column_Aeq))
            standard_Aeq = np.append(self.Aeq,new_array,axis=1)
            standard_A   = np.append(standard_A,standard_Aeq,axis=0)
        else:
            standard_A = standard_A 
        #标准化c
        standard_c = np.append(self.c,np.zeros(standard_A.shape[1]-self.c.shape[0]))
        
        return standard_A,standard_b,standard_c
        
    def judge(self):
        standardization = self.standardize()
        
        standard_A,standard_b,standard_c\
                = standardization[0].astype(float),\
                  standardization[1].astype(float),\
                  standardization[2].astype(float)

        #确定基变量
        Basic_variable = []
        for i in range(standard_A.shape[1]):
            for j in range(standard_A.shape[0]):
                unit_vector = np.zeros(standard_A.shape[0])
                unit_vector[j] = 1
                if (unit_vector == standard_A[:,i]).all():
                    Basic_variable.append(i)
                    break
        
        num_artifical_variable = 0
        status = 0
        artifical_variable = []
        #如果基变量个数等于方程式直接算
        if len(Basic_variable) == standard_A.shape[0]:    
            C_B = []
            sigma = []
            #确定基变量系数
            for i in range(len(Basic_variable)):
                C_B.append(standard_c[Basic_variable[i]])
            C_B = np.array(C_B)
            
            #计算sigma
            for i in range(len(standard_c)):
                num = standard_c[i] - np.dot(C_B,standard_A[:,i])
                sigma.append(num)
            
            return standard_A,standard_b,standard_c,Basic_variable,C_B,sigma,num_artifical_variable,artifical_variable,status
        #如果没有单位矩阵，添加人工变量
        else:
            status = 1
            #在A中添加人工变量列
            standard_A_copy = standard_A.copy()
            
            for i in range(standard_A_copy.shape[0]):
                unit_vector = np.zeros(standard_A_copy.shape[0])
                unit_vector[i] = 1
                kk = 1
                for j in range(standard_A_copy.shape[1]):
                    if (unit_vector == standard_A_copy[:,j]).all():
                        kk += 1
                if kk <= 1:
                    standard_A = np.insert(standard_A,standard_A.shape[1],unit_vector,axis=1)            
                    num_artifical_variable +=1
                    artifical_variable.append(standard_A.shape[1] - 1)
                    
            #在c中添加人工变量系数
            standard_c = np.zeros(standard_A.shape[1])
            for i in range(num_artifical_variable):
                standard_c[-(i+1)] = -1

            #确定基变量
            Basic_variable = []
            for i in range(standard_A.shape[1]):
                for j in range(standard_A.shape[0]):
                    unit_vector = np.zeros(standard_A.shape[0])
                    unit_vector[j] = 1
                    if (unit_vector == standard_A[:,i]).all():
                        Basic_variable.append(i)
                        break      
            #print(standard_A)
            #print('Basic_variable is ' + str(Basic_variable))
            #确定基变量系数
            C_B = []
            for i in range(len(Basic_variable)):
                C_B.append(standard_c[Basic_variable[i]])
            C_B = np.array(C_B)

            #计算sigma
            sigma = []
            for i in range(len(standard_c)):
                num = standard_c[i] - np.dot(C_B,standard_A[:,i])
                sigma.append(num)

            return standard_A,standard_b,standard_c,Basic_variable,C_B,sigma,num_artifical_variable,artifical_variable,status

    def first_iterate(self):
        print('**********************************************')
        print('-------------------第一阶段--------------------')
        print('**********************************************')
        parameters = self.judge()
        standard_A,standard_b,standard_c,Basic_variable,C_B,sigma,num_artifical_variable,artifical_variable \
            = parameters[0].astype(float),\
              parameters[1].astype(float),\
              parameters[2].astype(float),\
              parameters[3],\
              parameters[4],\
              parameters[5],\
              parameters[6],\
              parameters[7]
        #松弛变量矩阵
        print('inital A is \n' + str(standard_A))
        print('Basic_variable is ' + str(Basic_variable))
        print('C_B is' + str(C_B))
        print('sigma is ' + str(sigma))
        print('standard_c is' + str(standard_c))
        print('')
        #print('人工变量是' + str(artifical_variable))
        #迭代
        while max(sigma) > 0 :
            #找进基变量
            print('sigma is ' + str(sigma))
            entering_variable = sigma.index(max(sigma))
            theta = []
            print("entering_variable is " + str(entering_variable))
            #找出基变量
            print("standard_b is " + str(standard_b))
            for i in range(standard_A.shape[0]):
                if standard_A[i][entering_variable] != 0 and standard_A[i][entering_variable] > 0:
                    theta.append(np.round(standard_b[i]/standard_A[i][entering_variable],self.frac))
                else:
                    theta.append(sys.maxsize)
            print("theta is " + str(theta))
            theta_copy = theta.copy()
            for i in theta_copy:
                if i <= 0 :
                    theta_copy.remove(i)

            leaving_variable_index = theta.index(min(theta_copy))
            if len(artifical_variable) != 0:
                for i in range(len(theta)):
                    if (min(theta_copy) == theta[i]) and (Basic_variable[i] in artifical_variable):
                        leaving_variable_index = i

            leaving_variable = Basic_variable[leaving_variable_index]
            print("leaving_variable is " + str(leaving_variable))
            #归一
            #找到出基变量的行
            leaving_variable_line = Basic_variable.index(leaving_variable)
            entering_variable_column = entering_variable
            num = standard_A[leaving_variable_line][entering_variable_column]
            for i in range(standard_A.shape[1]):
                standard_A[leaving_variable_line][i] = np.round(standard_A[leaving_variable_line][i]/num,self.frac)
            
            
            #更新b
            standard_b[leaving_variable_line] = np.round(standard_b[leaving_variable_line]/num,self.frac)
            
            #变为单位向量
            for i in range(standard_A.shape[0]):
                if i != leaving_variable_line:
                    num = np.round(standard_A[i][entering_variable_column]/ \
                        standard_A[leaving_variable_line][entering_variable_column],self.frac)
                        
                    for j in range(standard_A.shape[1]):
                        standard_A[i][j] = np.round(standard_A[i][j]-num*standard_A[leaving_variable_line][j],self.frac)
            #更新b
                    standard_b[i] = np.round(standard_b[i]- num * standard_b[leaving_variable_line],self.frac)
            #更新基变量
            Basic_variable[leaving_variable_index] = entering_variable
            #更新C_B
            C_B = []
            sigma = []
            for i in range(len(Basic_variable)):
                C_B.append(standard_c[Basic_variable[i]])
            C_B = np.array(C_B)
            #更新sigma
            for i in range(len(standard_c)):
                num = standard_c[i] - np.dot(C_B,standard_A[:,i])
                sigma.append(num)                
            #break
            print('A is \n' +str(standard_A))
            print('Basic_variable is'+ str(Basic_variable))
            print('b is '+ str(standard_b))
            print('C_B is' + str(C_B))
            print('')
            
        for i in range(self.A.shape[1]+self.A.shape[0]+self.Aeq.shape[0],standard_A.shape[1]+1):
            if i in Basic_variable:
                raise ValueError('无可行解')
        
        return standard_A,standard_b,standard_c,Basic_variable,C_B,sigma,num_artifical_variable
    
    def second_iterate(self):

        parameters = self.first_iterate()
        standard_A,standard_b,standard_c,Basic_variable,num_artifical_variable \
            = parameters[0].astype(float),\
              parameters[1].astype(float),\
              parameters[2].astype(float),\
              parameters[3],\
              parameters[-1]
        
        standard_A = standard_A[:,:-num_artifical_variable]
        standard_c = np.zeros(standard_A.shape[1])
        for i in range(self.c.shape[0]):
            standard_c[i] = self.c[i]
        
        C_B = []
        for i in range(len(Basic_variable)):
            C_B.append(standard_c[Basic_variable[i]])
        C_B = np.array(C_B)
        #计算sigma
        sigma = []
        for i in range(len(standard_c)):
            num = np.round(standard_c[i] - np.dot(C_B,standard_A[:,i]),self.frac)
            sigma.append(num)    

        print('**********************************************')
        print('-------------------第二阶段--------------------')
        print('**********************************************')        
        #松弛变量矩阵
        print('inital A is \n' + str(standard_A))
        print('Basic_variable is ' + str(Basic_variable))
        print('C_B is' + str(C_B))
        print('sigma is ' + str(sigma))
        print('standard_b is ' + str(standard_b))
        print('')
        #迭代

        while max(sigma) > 0 :
            #找进基变量
            entering_variable = sigma.index(max(sigma))
            theta = []
            print("entering_variable is " + str(entering_variable))
            #找出基变量
            for i in range(standard_A.shape[0]):
                if standard_A[i][entering_variable] != 0 and standard_A[i][entering_variable] > 0:
                    theta.append(np.round(standard_b[i]/standard_A[i][entering_variable],self.frac))
                else:
                    theta.append(sys.maxsize)
            print("theta is " + str(theta))
            theta_copy = theta.copy()
            for i in theta_copy:
                if i <= 0 :
                    theta_copy.remove(i)
            
            
            leaving_variable_index = theta.index(min(theta_copy))
            leaving_variable = Basic_variable[leaving_variable_index]
            print("leaving_variable is " + str(leaving_variable))
            #归一
            #找到出基变量的行
            leaving_variable_line = Basic_variable.index(leaving_variable)
            entering_variable_column = entering_variable
            num = standard_A[leaving_variable_line][entering_variable_column]
            for i in range(standard_A.shape[1]):
                standard_A[leaving_variable_line][i] = np.round(standard_A[leaving_variable_line][i]/num,self.frac)
            
            
            
            #更新b
            standard_b[leaving_variable_line] = np.round(standard_b[leaving_variable_line]/num,self.frac)
            
            #变为单位向量
            for i in range(standard_A.shape[0]):
                if i != leaving_variable_line:
                    num = np.round(standard_A[i][entering_variable_column]/ \
                        standard_A[leaving_variable_line][entering_variable_column],self.frac)
                    for j in range(standard_A.shape[1]):
                        standard_A[i][j] = np.round(standard_A[i][j]-num*standard_A[leaving_variable_line][j],self.frac)
            #更新b
                    standard_b[i] = np.round(standard_b[i]- num * standard_b[leaving_variable_line],self.frac)
            #更新基变量
            Basic_variable[leaving_variable_index] = entering_variable
            #更新C_B
            C_B = []
            sigma = []
            for i in range(len(Basic_variable)):
                C_B.append(standard_c[Basic_variable[i]])
            C_B = np.array(C_B)
            #更新sigma
            for i in range(len(standard_c)):
                num = standard_c[i] - np.dot(C_B,standard_A[:,i])
                sigma.append(num)
            #break
            print('A is \n' +str(standard_A))
            print('Basic_variable is'+ str(Basic_variable))
            print('b is '+ str(standard_b))
            print('C_B is' + str(C_B))
            print('')
        return standard_A,standard_b,standard_c,Basic_variable,C_B,sigma
    
    def solve(self):
        judgement = self.judge()
        if judgement[-1] == 0:
            solve_1 = self.first_iterate()
            return solve_1[0],solve_1[1],solve_1[2],solve_1[3],solve_1[4],solve_1[5]
        
        else:
            solve_2 = self.second_iterate()
            return solve_2[0],solve_2[1],solve_2[2],solve_2[3],solve_2[4],solve_2[5]
            
            
    def Generate_Cut(self):
        Dual_Simplex_parameters = self.solve()
        print('**********************************************')
        print('-------------------第三阶段--------------------')
        print('**********************************************')          
        
        
        standard_A,standard_b,standard_c,Basic_variable,C_B,sigma = \
            Dual_Simplex_parameters[0].astype(float),\
            Dual_Simplex_parameters[1].astype(float),\
            Dual_Simplex_parameters[2].astype(float),\
            Dual_Simplex_parameters[3],\
            Dual_Simplex_parameters[4],\
            Dual_Simplex_parameters[5]
        
        

        frac_list = []
        frac_num = 0
        for j in range(standard_b.shape[0]):
            if standard_b[j] - np.floor(standard_b[j]) > 1e-3 or np.ceil(standard_b[j]) - standard_b[j] > 1e-3:
                frac_num += 1
        aa = 1
        while frac_num != 0 and aa <=3:
            #print('**********************************************')
            #print('-------------------第四阶段--------------------')
            #print('**********************************************')   
            #找到最大小数部分
            for i in range(standard_b.shape[0]):
                frac = standard_b[i] - np.floor(standard_b[i])
                frac_list.append(frac)

            most_frac_variable_index = frac_list.index(max(frac_list))
            #most_frac_variable = Basic_variable[most_frac_variable_index]
            the_line = most_frac_variable_index

            print(the_line)
            #加入到约束中
            Cut = np.zeros(standard_A.shape[1])
            for i in range(Cut.shape[0]):
                if i not in Basic_variable:
                    Cut[i] = -(standard_A[the_line][i] - np.floor(standard_A[the_line][i]))
            standard_A = np.insert(standard_A,standard_A.shape[0],Cut,axis =0 )
            
            #加松弛变量
            relaxation_variable = np.zeros(standard_A.shape[0])
            relaxation_variable[-1] = 1
            standard_A = np.insert(standard_A,standard_A.shape[1],relaxation_variable,axis = 1)
            standard_b = np.insert(standard_b,standard_b.shape[0],\
                                   -(standard_b[the_line] - np.floor(standard_b[the_line])), axis = 0)
            sigma.append(0)
            
            standard_c = np.insert(standard_c,standard_c.shape[0],0,axis=0)
            
            Basic_variable.append(standard_A.shape[1]-1)
            
            #找出基变量
            leaving_variable = Basic_variable[-1]
            leaving_variable_line = Basic_variable.index(leaving_variable)
            #print('出基变量是' + str(leaving_variable))
            
            print(standard_A)
            while min(standard_b) < 0:
                print('基变量是' + str(Basic_variable))
                print('simga 是' + str(sigma))
                print('b是 ' + str(standard_b))
                
                #找进基变量
                theta = []
                for i in range(len(sigma)):
                    if standard_A[leaving_variable_line][i] != 0 and standard_A[leaving_variable_line][i] < 0:
                        theta.append(np.round(sigma[i] / standard_A[leaving_variable_line][i],self.frac))
                    else:
                        theta.append(123456)
                print('theta 是' + str(theta))
                entering_variable_column = theta.index(min(theta))
                print('进基变量是' + str(entering_variable_column))
                
                #找出基变量
                for i in range(standard_b.shape[0]):
                    if standard_b[i] < 0 :
                        leaving_variable_line = i
                        leaving_variable = Basic_variable[i]
                print('出基变量是' + str(leaving_variable))
                #归一
                num = standard_A[leaving_variable_line][entering_variable_column]
                for i in range(standard_A.shape[1]):
                    standard_A[leaving_variable_line][i] = np.round(standard_A[leaving_variable_line][i]/num,self.frac)
                
                #更新b
                standard_b[leaving_variable_line] = np.round(standard_b[leaving_variable_line]/num,self.frac)
                
                #变为单位向量
                for i in range(standard_A.shape[0]):
                    if i != leaving_variable_line:
                        num = standard_A[i][entering_variable_column] / standard_A[leaving_variable_line][entering_variable_column]
                        for j in range(standard_A.shape[1]):
                            standard_A[i][j] = np.round(standard_A[i][j] - num *standard_A[leaving_variable_line][j],self.frac)
                    
                    #更新b
                        standard_b[i] = np.round(standard_b[i] - num*standard_b[leaving_variable_line])
                
                #更新基变量
                Basic_variable[leaving_variable_line] = entering_variable_column
                #更新C_B
                C_B = []
                sigma = []
                for i in range(len(Basic_variable)):
                    C_B.append(standard_c[Basic_variable[i]])
                C_B = np.array(C_B)
                
                #更新sigma
                for i in range(len(standard_c)):
                    num = standard_c[i] - np.dot(C_B,standard_A[:,i])
                    sigma.append(num)
                print('')
                print('A是'+str(standard_A))        
                print('基变量是' + str(Basic_variable))
                print('b是' + str(standard_b))
                print('')
            break
            frac_list = []
            frac_num = 0
            for j in range(standard_b.shape[0]):
                if standard_b[j] - np.floor(standard_b[j]) > 1e-3 or np.ceil(standard_b[j]) - standard_b[j] > 1e-3:
                    frac_num += 1
            aa += 1
            print(standard_b)
            
        
        
def test1():
    """ 此测试的真实最优解为x^T = [36/11,,40/11,0,0,75/11]^T """
    c = np.array([3,-1])
    A = np.array([[3, -2],[2,1]])
    b = np.array([3,5])
    Agq = np.array([[5,4]])
    bgq = np.array([10])
    Aeq = np.array([])
    beq = np.array([])
    bounds = [(0, None), (0, None),(0,None)]
 
    solver = BC(c, A, b,Agq,bgq, Aeq, beq, bounds)
    #solver.standardize()
    solver.solve()
    solver.Generate_Cut()

if __name__ == '__main__':
    test1() 

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    