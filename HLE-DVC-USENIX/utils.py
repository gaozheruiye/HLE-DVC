from typing import List
import numpy as np






# Get the set of powers of R, until but not including when the powers
# loop back to 1
# 获得子群，就是omega的子群
def get_power_cycle(r, modulus):
    o = [1, r]
    while o[-1] != 1:
        o.append((o[-1] * r) % modulus)
    return o[:-1]



def is_a_power_of_2(x):
    return True if x==1 else False if x%2 else is_a_power_of_2(x//2)

# #输入一个变量power列表返回下标，这个下标代表着对应公开参数在长向量中的位置
# def get_index(list, rho, l):
#     index = 0
#     for i in range(rho):
#         index += 2**(rho-i-1)*list[i]
#     index *= l
#     index += list[rho]-1
#     print(index)

def gen_binary_list(list_len:int)->List[List[int]]:
    out=[[0]*list_len for _ in range(1<<list_len)]
    for l in range(list_len):
        for j in range(1<<list_len):
            out[j][list_len-l-1]=(j>>l)&0x1
    return out

def gen_variable_power_list(n:int, rho:int,l:int)->List[List[int]]:
    list_len = n * l
    assert n == 1<<rho
    output = [[0]*(rho+1) for i in range(l*n) ]
    # print(output)

    out=[[0]*rho for _ in range(n)]
    for i in range(rho):
        for j in range(n):
            out[j][rho-i-1]=(j>>i)&0x1

    for i in range(n):
        for j in range(l):
            output[i*l+j][:-1] = out[i]
            output[i*l+j][-1] = j

    return output

def gen_rho_plus_1_dim_array(rho, n, size=2, fill=0):
    """
    创建一个rho+1维数组，前rho维长度为2，最后一个维度为n，元素为fill
    """
    shape = tuple([size] * rho + [n])
    return np.full(shape, fill, dtype=object)

if __name__ == "__main__":
    # get_index((1,0,1,4),3,4)

    # 示例：创建一个3维，形状为(2, 2, 2)的数组
    rho = 4
    n =4
    arr = gen_rho_plus_1_dim_array(rho, n )
    print(arr)