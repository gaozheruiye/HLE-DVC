import time
import random
from utils import get_power_cycle
from fft import fft
from util_G1G2 import G1, G2, pairing, add, multiply, field_element
from poly_utils import PrimeField
modulus = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001


def G1Add():
    # 定义 G1 群的生成元（基点）
    P = G1  # 这是 G1 群的生成元
    Q = add(P, P)  # 计算 2P
    # print(f"2P = {Q}")

def G1Mul():
    P = G1  # 生成元
    # 选择一个标量 k
    k = random.randint(1, modulus)

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
    g1_add_time = 0
    g1_mul_time = 0
    g2_add_time = 0
    g2_mul_time = 0
    pair_time = 0
    fft_time = 0
    num = 100
    P = G1

    add_start_time = time.time()
    for i in range(num):
        P = add(P, G1)
    add_end_time = time.time()
    g1_add_time = add_end_time - add_start_time
    print("G1 g1_add_time",g1_add_time/num)


    k = random.randint(1, modulus)
    mul_start_time = time.time()
    for i in range(num):
        P= multiply(P,k)
    mul_end_time = time.time()
    g1_mul_time = mul_end_time - mul_start_time
    print("G1 g1_mul_time",g1_mul_time/num)

    pair_start_time = time.time()
    for i in range(num):
        a=pairing(G2, P)
    pair_end_time = time.time()
    pair_time += pair_end_time - pair_start_time
    print("pair_time",pair_time/num)

    P = G2
    add_start_time = time.time()
    for i in range(num):
        P = add(P, G2)
    add_end_time = time.time()
    g2_add_time = add_end_time - add_start_time
    print("G2 g2_add_time",g2_add_time/num)


    k = random.randint(1, modulus)
    mul_start_time = time.time()
    for i in range(num):
        P= multiply(P,k)
    mul_end_time = time.time()
    g2_mul_time = mul_end_time - mul_start_time
    print("G2 g2_mul_time",g2_mul_time/num)

    # 测试fft,p是vector length
    n = pow(2,20)
    f = PrimeField(modulus)
    omega_n = f.exp(7, (modulus-1)//n)
    omega_n_s = get_power_cycle(omega_n, modulus)

    vector = [random.randint(1, modulus) for i in range(pow(2,20))]
    fft_start_time = time.time()
    fft(vector, modulus, omega_n_s[1], inv=True)
    fft_end_time = time.time()
    fft_time = fft_end_time - fft_start_time
    print("fft_mul_time",fft_time)

