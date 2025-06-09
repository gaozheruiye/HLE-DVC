from py_ecc.bls12_381 import G1, G2
from py_ecc.bls12_381 import  pairing, add, multiply

import random
import numpy as np
from polynomials import evaluate_univariate_polynomial_at_x, univariate_polynomial_division, bivariate_polynomial_division_x, mod_inverse
import copy
from utils import gen_variable_power_list, gen_binary_list, gen_rho_plus_1_dim_array, get_power_cycle
from Lagrange_interpolation import univariate_lagrange_interpolation_coefficients
from poly_utils import PrimeField
# from utils import get_power_cycle





# M: 证明者的人数
# l: 子向量长度
# rho：log M
# 下标k：M；下标i：l；下标c:rho
class Hybrid_mul_polynomial_commitment_scheme:
    def __init__(self, M, n, N, rho, omega_n_s, modulus, vector):
        ''' 构造函数，初始化对象，并且包含了setup函数的功能'''
        self.M = M
        self.n = n
        self.N = N
        self.rho = rho
        self.omega_n_s = omega_n_s
        self.modulus = modulus
        self.vector = vector
        self.subvector = [vector[i * n : (i + 1) * n] for i in range(M)]
        # self.P = G1
        print("self.vector:", self.vector)
        print("self.subvector:", self.subvector)

        # 秘密陷门，应当销毁 前rho个是tau_s,最后一个是tau_x
        # self.alpha_s_x = [random.randint(2, modulus) for _ in range(rho +1)]
        self.alpha_s_x = [3,4,5,6,7]
        print("self.alpha_s_x:", self.alpha_s_x)


        # 产生一个二进制二维列表，大小为（M,rho）。里面是[0,M)的二进制分解
        self.binary_list = gen_binary_list(rho)
        print("self.binary_list:", self.binary_list)

        # 产生一个rho+1 维度的list，其中前rho个度数为2，最后一个为n
        self.pp_PST = gen_rho_plus_1_dim_array(rho,n)
        for index, value in np.ndenumerate(self.pp_PST):
            # print("index: ",index)
            # print("value: ",value)
            temp = 1
            for i,v in enumerate(index):
                temp *= pow(self.alpha_s_x[i], v, self.modulus)
            self.pp_PST[index] = multiply(G1, temp% self.modulus)
        # pp_PST_pie指的是[eq(s,k)x^i]_{k\in [0,M),i\in [0,n)};size (2,2,n)
        self.pp_PST_pie = gen_rho_plus_1_dim_array(rho,n)
        # print(self.pp_PST_pie)
        for index, value in np.ndenumerate(self.pp_PST_pie):
            # print("index: ",index) # 指的是pp_PST_pie内部的index
            # print("value: ",value)
            temp = 1
            for i in range(len(index)-1):
                if index[i]==0:
                    temp *= (1-self.alpha_s_x[i])
                else:
                    temp *= self.alpha_s_x[i]
            temp *= pow(self.alpha_s_x[-1], index[-1], self.modulus)
            temp = temp % self.modulus
            self.pp_PST_pie[index] = multiply(G1, temp% self.modulus)
        # print("self.pp_PST_pie",self.pp_PST_pie)

        # pp_PST_pie_pie ;size (2,2,rho,n)
        shape = tuple([2] * rho + [rho]+ [n])
        self.pp_PST_pie_pie =  np.full(shape, 0, dtype=object)
        for index, value in np.ndenumerate(self.pp_PST_pie_pie):
            # self.pp_PST_pie_pie[index] = self.pp_PST_pie[*index[:-2],index[-1]] * mod_inverse(self.alpha_s_x[index[-2]]-(1-index[index[-2]]), self. modulus) % self.modulus
            self.pp_PST_pie_pie[index] = multiply(self.pp_PST_pie[*index[:-2],index[-1]],
                                                  mod_inverse(self.alpha_s_x[index[-2]]-(1-index[index[-2]]), self. modulus) % self.modulus)

    def dist_commit(self):
        '''
        分布式产生承诺
        '''
        # 每个证明者先产生子多项式f(x)
        self.coefficients = [[0]*self.n]*self.M
        self.C_v = G1

        for k in range(self.M):
            index_k = self.binary_list[k]
            # self.coefficients size: (M,n)
            self.coefficients[k] = univariate_lagrange_interpolation_coefficients(self.omega_n_s, self.subvector[k], self.modulus)[::-1]
            # 注意！！！！这里coefficients[k][0]是常数项，从小到大排序
            # print("self.coefficients[k]:",self.coefficients[k])
            # print("求值",evaluate_univariate_polynomial_at_x(self.coefficients[k][::-1],self.omega_n_s[0],self.modulus))
            # 这个库好像没有给零点，所以需要先设它为G1，再把它减掉
            C_v_k = G1
            # print("!!!!!!!!!!",add(multiply(G1,2),multiply(G1,mod_inverse(-1,self.modulus))))
            # print("G1",G1)
            for i in range(self.n):
                # C_v_k +=  self.coefficients[k][i] * self.pp_PST_pie[*index_k,i]
                C_v_k = add(C_v_k,multiply(self.pp_PST_pie[*index_k,i],self.coefficients[k][i]))
            C_v_k = add(C_v_k,multiply(G1,mod_inverse(-1,self.modulus)))
            # print("C_v_k:",C_v_k)
            self.C_v = add(C_v_k,self.C_v)
        self.C_v = add(self.C_v,multiply(G1,mod_inverse(-1,self.modulus)))
        print("self.coefficients:",self.coefficients)
        print("承诺C_v：", self.C_v)

    def genAux(self):
        '''
        产生AuxData; size: (M,rho)
        '''
        self.Aux = np.full((M,rho), 0, dtype=object)
        # print("self.Aux.shape",self.Aux.shape)
        for k in range(M):
            index_k_binary = self.binary_list[k]
            for c in range(rho):
                self.Aux[k, c] = G1
                for i in range(n):
                    self.Aux[k,c] = add(self.Aux[k,c], multiply(self.pp_PST_pie_pie[*index_k_binary,c,i],self.coefficients[k][i]))
                self.Aux[k,c] =add(self.Aux[k,c],multiply(G1,mod_inverse(-1,self.modulus)))


    def genPartialProof(self): #显然每个机器的视角是等价的，工作量
        # 仅仅以P_0为视角 Q_{0,0} = {2,3}, Q_{0,1}={1}
        # partialProof_P_0 = [add(self.Aux[2,0], self.Aux[3,0]), self.Aux[1,1]]
        tmp0 = add(self.Aux[8,0], self.Aux[9,0])
        tmp0 = add(tmp0, self.Aux[10,0])
        tmp0 = add(tmp0, self.Aux[11,0])
        tmp0 = add(tmp0, self.Aux[12,0])
        tmp0 = add(tmp0, self.Aux[13,0])
        tmp0 = add(tmp0, self.Aux[14,0])
        tmp0 = add(tmp0, self.Aux[15,0])
        tmp1 = add(self.Aux[4,1], self.Aux[5,1])
        tmp1 = add(tmp1, self.Aux[6,1])
        tmp1 = add(tmp1, self.Aux[7,1])

        partialProof_P_0 = [tmp0,tmp1,add(self.Aux[2,2], self.Aux[3,2]), self.Aux[1,3]]
        # print("partialProof_P_0:",partialProof_P_0)
        return partialProof_P_0

    def genAllPartialProof(self):
        list =[[]*2*n]*self.rho
        for i in range(self.rho):
            for j in range(self.M):
                list[i][j] = self.Aux[j,i]
        for i in range(4):
            for j in range(8):
                list[i][j+16] = list[i][2*j]+list[i][2*j+1]
            for j in range(4):
                list[i][j+24] = list[i][2*j+16]+list[i][2*j+1+16]
            list[i][29] = list[i][26]+list[i][25]
            list[i][30] = list[i][27]+list[i][28]
            list[i][31] = list[i][30]+list[i][29]



    def prove(self, k, i):
        index_k_binary_P_0 = self.binary_list[0]
        value = self.subvector[k][i]
        # 仅仅以P_0为视角
        # 需要计算多项式q(x) =  (witness_univariate_polynomial - value)/(x-\omega_n^)
        numerator_polynomial = self.coefficients[k]
        # 因为原来的代码是高位在前，没办法只能先倒置一下，等会再倒置回来
        numerator_polynomial = numerator_polynomial[::-1]
        numerator_polynomial[-1] = numerator_polynomial[-1] - value
        denominator_polynomial = [1, -omega_n_s[i]]
        print("numerator_polynomial",numerator_polynomial)
        proof_polynomial, remainder_polynomial = univariate_polynomial_division(numerator_polynomial, denominator_polynomial, self.modulus)
        print("证明多项式:", proof_polynomial)
        assert remainder_polynomial == [], "证明多项式生成失败"
        proof_polynomial = proof_polynomial[::-1] + [0]
        pi_d_rho_i = G1
        for i in range(n):
            pi_d_rho_i =add(pi_d_rho_i, multiply(self.pp_PST_pie[*index_k_binary_P_0][i],proof_polynomial[i]))
        # commit_proof_polynomial = self.univariate_commit(proof_polynomial)
        pi_d_rho_i = add(pi_d_rho_i, multiply(G1,mod_inverse(-1,self.modulus)))
        print("pi_d_rho_i:",pi_d_rho_i)
        return pi_d_rho_i


    def verify(self, partialProof,pi_d_rho_i,value,k,i):
        # left = pairing(self.C_v,)
        left = self.C_v
        right = G1
        index_binary_P_0 = self.binary_list[0] #不是P_0的话就把这个换了

        for c in range(self.rho):
            # right += (self.alpha_s_x[c]-index_binary_P_0[c])*partialProof[c]
            # right += (self.alpha_s_x[c]-index_binary_P_0[c])*partialProof[c]
            right =add(right, multiply(partialProof[c], self.alpha_s_x[c]-index_binary_P_0[c]))
        addition1 = multiply(pi_d_rho_i,(self.alpha_s_x[self.rho]-self.omega_n_s[i])% self.modulus)
        addition2 = multiply(self.pp_PST_pie[*index_binary_P_0,0],value)
        temp = add(addition1,addition2)
        # right += (self.alpha_s_x[self.rho]-self.omega_n_s[i])*pi_d_rho_i + value * self.pp_PST_pie[*index_binary_P_0,0]
        right = add(right, temp)

        right = add(right, multiply(G1,mod_inverse(-1,self.modulus)))

        print("left:",left)
        print("right:",right)
        assert left == right, "Verification failed."
        print("Verification passed")


if __name__ == "__main__":
    # M is fixed to 16
    rho = 4
    M = 16
    n = 4
    N = M * n
    modulus = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
    # modulus = 257
    f = PrimeField(modulus)
    #单位根
    omega_n = f.exp(7, (modulus-1)//n)
    omega_n_s = get_power_cycle(omega_n, modulus)
    # print ("omega_n_s:",omega_n_s)
    vector = [random.randint(1, modulus) for _ in range(N)]
    # vector = [i+1 for i in range(16)]
    # print(vector_long)
    example = Hybrid_mul_polynomial_commitment_scheme(M, n, N, rho,omega_n_s, modulus, vector)
    example.dist_commit()
    example.genAux()
    partialProof_P_0 = example.genPartialProof()
    k = 0
    i = 1
    pi_d_rho_i = example.prove(k,i)
    value = vector[k*n+i]
    print("value:", value)
    example.verify(partialProof_P_0, pi_d_rho_i, value,k,i)










# # M: 证明者的人数
# # l: 子向量长度
# # rho：log M
# # 下标k：M；下标i：l；下标c:rho
# class Hybrid_mul_polynomial_commitment_scheme:
#     def __init__(self, M, n, N, rho, omega_n_s, modulus, vector):
#         ''' 构造函数，初始化对象，并且包含了setup函数的功能'''
#         self.M = M
#         self.n = n
#         self.N = N
#         self.rho = rho
#         self.omega_n_s = omega_n_s
#         self.modulus = modulus
#         self.vector = vector
#         self.subvector = [vector[i * n : (i + 1) * n] for i in range(M)]
#         print("self.vector:", self.vector)
#         print("self.subvector:", self.subvector)

#         # 秘密陷门，应当销毁 前rho个是tau_s,最后一个是tau_x
#         # self.alpha_s_x = [random.randint(2, modulus) for _ in range(rho +1)]
#         self.alpha_s_x = [3,4,5]
#         print("self.alpha_s_x:", self.alpha_s_x)


#         # 产生一个二进制二维列表，大小为（M,rho）。里面是[0,M)的二进制分解
#         self.binary_list = gen_binary_list(rho)
#         print("self.binary_list:", self.binary_list)

#         # 产生一个rho+1 维度的list，其中前rho个度数为2，最后一个为n
#         self.pp_PST = gen_rho_plus_1_dim_array(rho,n)
#         for index, value in np.ndenumerate(self.pp_PST):
#             # print("index: ",index)
#             # print("value: ",value)
#             temp = 1
#             for i,v in enumerate(index):
#                 temp *= pow(self.alpha_s_x[i], v, self.modulus)
#             self.pp_PST[index] = temp % self.modulus
#         # pp_PST_pie指的是[eq(s,k)x^i]_{k\in [0,M),i\in [0,n)};size (2,2,n)
#         self.pp_PST_pie = gen_rho_plus_1_dim_array(rho,n)
#         # print(self.pp_PST_pie)
#         for index, value in np.ndenumerate(self.pp_PST_pie):
#             print("index: ",index) # 指的是pp_PST_pie内部的index
#             # print("value: ",value)
#             temp = 1
#             for i in range(len(index)-1):
#                 if index[i]==0:
#                     temp *= (1-self.alpha_s_x[i])
#                 else:
#                     temp *= self.alpha_s_x[i]
#             temp *= pow(self.alpha_s_x[-1], index[-1], self.modulus)
#             temp = temp % self.modulus
#             print("temp: ",temp) # 指的是pp_PST_pie内部的index
#             print("index: ",index) # 指的是pp_PST_pie内部的index
#             self.pp_PST_pie[index] = temp
#         print("self.pp_PST_pie",self.pp_PST_pie)

#         # pp_PST_pie_pie ;size (2,2,rho,n)
#         shape = tuple([2] * rho + [rho]+ [n])
#         self.pp_PST_pie_pie =  np.full(shape, 0)
#         for index, value in np.ndenumerate(self.pp_PST_pie_pie):
#             # print("index: ",index) # 指的是pp_PST_pie_pie内部的index
#             # if index[index[-2]]==0:
#             # print("index[:-3],index[-1]",index[:-2],index[-1])
#             # print("self.pp_PST_pie",self.pp_PST_pie)
#             # print("self.pp_PST_pie[index[:-3],index[-1]]",self.pp_PST_pie[*index[:-2],index[-1]])
#             self.pp_PST_pie_pie[index] = self.pp_PST_pie[*index[:-2],index[-1]] * mod_inverse(self.alpha_s_x[index[-2]]-(1-index[index[-2]]), self. modulus) % self.modulus
#             # print("self.pp_PST_pie_pie[index]",self.pp_PST_pie_pie[index])

#         # for index, value in np.ndenumerate(self.pp_PST):
#         #     print("index: ",index)
#         #     print("value: ",value)

#         # 产生PST多项式承诺的public parameter：
#         ## 产生阶数列表，长度为(M*l,rho+1),这个实际上不怎么用得到，因为有更快捷的方法（基于拉格朗日插值法的pp）


#     def dist_commit(self):
#         '''
#         分布式产生承诺
#         '''
#         # 每个证明者先产生子多项式f(x)
#         self.coefficients = [[0]*n]*M
#         self.C_v = 0
#         # print("self.pp_PST",self.pp_PST)

#         for k in range(self.M):
#             index_k = self.binary_list[k]
#             # print("k:",k)
#             # print("self.subvector[k]:",self.subvector[k])
#             # self.coefficients size: (M,n)
#             self.coefficients[k] = univariate_lagrange_interpolation_coefficients(omega_n_s, self.subvector[k], self.modulus)[::-1]
#             # 注意！！！！这里coefficients[k][0]是常数项，从小到大排序
#             # print("self.coefficients[k]:",self.coefficients[k])
#             # print("求值",evaluate_univariate_polynomial_at_x(self.coefficients[k][::-1],self.omega_n_s[0],self.modulus))
#             C_v_k = 0
#             for i in range(self.n):
#                 C_v_k +=  self.coefficients[k][i] * self.pp_PST_pie[*index_k,i]
#                 print("self.coefficients[k][i] * self.pp_PST_pie[*index_k,i]",self.coefficients[k][i] * self.pp_PST_pie[*index_k,i])
#             # print("C_v_k:",C_v_k)
#             self.C_v += C_v_k
#             print("承诺C_v_",k,":", C_v_k %self.modulus)
#         self.C_v = self.C_v  % self.modulus
#         print("self.coefficients:",self.coefficients)
#         print("承诺C_v：", self.C_v)

#     def genAux(self):
#         '''
#         产生AuxData; size: (M,rho)
#         '''
#         self.Aux = np.full((M,rho), 0)
#         # print("self.Aux.shape",self.Aux.shape)
#         for k in range(M):
#             index_k_binary = self.binary_list[k]
#             for c in range(rho):
#                 self.Aux[k, c] = 0
#                 for i in range(n):
#                     self.Aux[k,c] += self.coefficients[k][i] * self.pp_PST_pie_pie[*index_k_binary,c,i]
#                 self.Aux[k,c] =self.Aux[k,c] % self.modulus
#                 # print("k:",k)
#                 # print("c:",c)
#                 # print("self.Aux[k,c]:",self.Aux[k,c])

#     def genPartialProof(self): #显然每个机器的视角是等价的，工作量
#         # 仅仅以P_0为视角 Q_{0,0} = {2,3}, Q_{0,1}={1}
#         partialProof_P_0 = [self.Aux[2,0] + self.Aux[3,0], self.Aux[1,1]]
#         print("partialProof_P_0:",partialProof_P_0)
#         return partialProof_P_0

#     def prove(self, k, i):
#         index_k_binary_P_0 = self.binary_list[0]
#         value = self.subvector[k][i]
#         # 仅仅以P_0为视角
#         # 需要计算多项式q(x) =  (witness_univariate_polynomial - value)/(x-\omega_n^)
#         numerator_polynomial = self.coefficients[k]
#         # 因为原来的代码是高位在前，没办法只能先倒置一下，等会再倒置回来
#         numerator_polynomial = numerator_polynomial[::-1]
#         numerator_polynomial[-1] = numerator_polynomial[-1] - value
#         denominator_polynomial = [1, -omega_n_s[i]]
#         print("numerator_polynomial",numerator_polynomial)
#         proof_polynomial, remainder_polynomial = univariate_polynomial_division(numerator_polynomial, denominator_polynomial, self.modulus)
#         print("证明多项式:", proof_polynomial)
#         assert remainder_polynomial == [], "证明多项式生成失败"
#         proof_polynomial = proof_polynomial[::-1] + [0]
#         pi_d_rho_i = 0
#         for i in range(n):
#             pi_d_rho_i += self.pp_PST_pie[*index_k_binary_P_0][i] * proof_polynomial[i]
#         # commit_proof_polynomial = self.univariate_commit(proof_polynomial)
#         pi_d_rho_i = pi_d_rho_i % self.modulus
#         print("pi_d_rho_i:",pi_d_rho_i)
#         return pi_d_rho_i % self.modulus


#     def verify(self, partialProof,pi_d_rho_i,value,k,i):
#         left = self.C_v
#         right = 0
#         index_binary_P_0 = self.binary_list[0] #不是P_0的话就把这个换了
#         print("index_binary_P_0:",index_binary_P_0)

#         for c in range(self.rho):
#             right += (self.alpha_s_x[c]-index_binary_P_0[c])*partialProof[c]
#         right += (self.alpha_s_x[self.rho]-self.omega_n_s[i])*pi_d_rho_i + value * self.pp_PST_pie[*index_binary_P_0,0]
#         print(pi_d_rho_i)

#         right = right % self.modulus
#         print("right:",right)
#         assert left == right, "Verification failed."
#         print("Verification passed")


# if __name__ == "__main__":

#     rho = 2
#     M = 4
#     n = 4
#     N = M * n
#     modulus = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
#     modulus = 257
#     f = PrimeField(modulus)
#     #单位根
#     omega_n = f.exp(7, (modulus-1)//n)
#     omega_n_s = get_power_cycle(omega_n, modulus)
#     # print ("omega_n_s:",omega_n_s)
#     # vector = [random.randint(1, modulus) for _ in range(N)]
#     vector = [i+1 for i in range(16)]
#     # print(vector_long)
#     example = Hybrid_mul_polynomial_commitment_scheme(M, n, N, rho,omega_n_s, modulus, vector)
#     example.dist_commit()
#     example.genAux()
#     partialProof_P_0 = example.genPartialProof()
#     k = 0
#     i = 1
#     commit_proof_polynomial = example.prove(k,i)
#     pi_d_rho_i = example.prove(k,i)
#     value = vector[k*n+i]
#     print("value:", value)
#     example.verify(partialProof_P_0, pi_d_rho_i, value,k,i)
