from py_ecc.bls12_381 import  G1, G2
from py_ecc.bls12_381 import  pairing, add, multiply
import PCS_group
import random
import numpy as np
from polynomials import evaluate_univariate_polynomial_at_x, univariate_polynomial_division, bivariate_polynomial_division_x, mod_inverse
import copy
from utils import gen_variable_power_list, gen_binary_list, gen_rho_plus_1_dim_array, get_power_cycle
from Lagrange_interpolation import univariate_lagrange_interpolation_coefficients
from poly_utils import PrimeField

mod = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

'''
M: machine number， 修改M需要修改 PCS_group.py
n: subvector length
N: vector length
在本实验中，批量打开（在拥有 partial Proo的情况下）单机聚合的代码与aSVC一致，因此可以通过运行下面的代码测算时间开销
t https://github.com/sunblaze-ucb/eVSS
在本实验中，多机聚合的代码与Hyperproofs一致，因此可以通过运行下面的代码测算时间开销
https://github.com/hyperproofs/hyperproofs-go.git
'''

#
rho = 4
M = 16
n = 4
N = M * n
modulus = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
f = PrimeField(modulus)
#单位根
omega_n = f.exp(7, (modulus-1)//n)
omega_n_s = get_power_cycle(omega_n, modulus)
# print ("omega_n_s:",omega_n_s)
vector = [random.randint(1, modulus) for _ in range(N)]
# vector = [i+1 for i in range(16)]
# print(vector_long)
example = PCS_group.Hybrid_mul_polynomial_commitment_scheme(M, n, N, rho,omega_n_s, modulus, vector)
example.dist_commit()
example.genAux()
partialProof_P_0 = example.genPartialProof()
k = 0
i = 1
pi_d_rho_i = example.prove(k,i)
value = vector[k*n+i]
print("value:", value)
example.verify(partialProof_P_0, pi_d_rho_i, value,k,i)



# def G1Add():
#     # 定义 G1 群的生成元（基点）
#     P = G1  # 这是 G1 群的生成元
#     Q = add(P, P)  # 计算 2P
#     print(f"2P = {Q}")

# def G1Mul():
#     P = G1  # 生成元
#     # 选择一个标量 k
#     k = 123456789

#     # 计算 k * P（即 G1 群上的点乘）
#     Q = multiply(P, k)

#     print(f"k * P = {Q}")

# def Pair():
#     P = G1  # 生成元
#     Q = G2  # 生成元
#     # # 选择一个标量 k
#     # k = 123456789

#     # # 计算 k * P（即 G1 群上的点乘）
#     # Q = multiply(P, k)

#     # print(f"k * P = {Q}")
#     # # 计算 e(G1, G2)
#     result = pairing(Q, P)

#     print(f"e(G2, G1) = {result}")



# if __name__ == "__main__":
#     G1Add()
#     G1Mul()
#     Pair()
