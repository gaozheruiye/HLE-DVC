from py_ecc.bls12_381 import  G1, G2
from py_ecc.bls12_381 import  pairing, add, multiply
import time

# from py_ecc.bls12_381 import (
#     pairing,
#     G1,
#     G2,
#     multiply,
# )

# # 取两个点
# P1 = multiply(G1, 3)
# Q1 = multiply(G2, 5)

# P2 = multiply(G1, 6)
# Q2 = multiply(G2, 5)

# # 分别做pairing
# e1 = pairing(Q1,P1)
# e2 = pairing(Q1,P1)
# e3 = pairing(Q2,P2)

# # 在GT中相乘
# product = e1 * e2 *pairing(G2,multiply(G1,0)) # FQ12 元素可以直接用 * 运算符

# print("product",product)
# print("e3",e3)





def G1Add():
    # 定义 G1 群的生成元（基点）
    P = G1  # 这是 G1 群的生成元
    Q = add(P, P)  # 计算 2P
    # print(f"2P = {Q}")

def G1Mul():
    P = G1  # 生成元
    # 选择一个标量 k
    k = 123456789

    # 计算 k * P（即 G1 群上的点乘）
    Q = multiply(P, k)

    # print(f"k * P = {Q}")

def Pair():
    P = G1  # 生成元
    Q = G2  # 生成元
    # # 选择一个标量 k
    # k = 123456789

    # # 计算 k * P（即 G1 群上的点乘）
    # Q = multiply(P, k)

    # print(f"k * P = {Q}")
    # # 计算 e(G1, G2)
    result = pairing(Q, P)

    # print(f"e(G2, G1) = {result}")



if __name__ == "__main__":
    add_time = 0
    pair_time = 0
    mul_time = 0
    for i in range(10):
        add_start_time = time.time()
        G1Add()
        add_end_time = time.time()
        add_time += add_end_time - add_start_time

    print("add_time",add_time/10)
    for i in range(10):

        mul_start_time = time.time()
        G1Mul()
        mul_end_time = time.time()
        mul_time += mul_end_time - mul_start_time
    print("mul_time",mul_time/10)

    for i in range(10):
        pair_start_time = time.time()
        Pair()
        pair_end_time = time.time()
        pair_time += pair_end_time - pair_start_time
    print("pair_time",pair_time/10)





















# import galois

# # 定义有限域 GF(17)
# GF = galois.GF(17)

# # 定义多项式，例如 a(x) = x^2 - 3x + 2，b(x) = x - 1
# # 注意系数也要转为 GF 元素
# a = galois.Poly([1, -3, 2], field=GF)  # 对应 x^2 - 3x + 2
# b = galois.Poly([1, -1], field=GF)     # 对应 x - 1

# # 执行多项式除法
# q, r = divmod(a, b)

# print("Quotient:", q)
# print("Remainder:", r)
