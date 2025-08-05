from typing import List

def gen_list(list_len:int)->List[List[int]]:
    out=[[0]*list_len for _ in range(1<<list_len)]
    for l in range(list_len):
        for j in range(1<<list_len):
            out[j][list_len-l-1]=(j>>l)&0x1
    return out


def gen_list_new(n:int, rho:int,l:int)->List[List[int]]:
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
            # print(output[i*l+j])
            output[i*l+j][-1] = j
            # print(output[i*l+j])
            # print(out[i])
    # print
    return output

# def num_2_list(val:int,len:int):
#     out=[0]*len
#     for i in range(len):
#         out[len-i-1]=val>>i & 0x1
#     return out

rho = 4
n = 16
l = 3
print(gen_list_new(n, rho,l))
# print([[0]*(rho+1) for i in range(l*n) ])

# print(num_2_list(5,3))