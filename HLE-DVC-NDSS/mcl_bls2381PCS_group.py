# from py_ecc.bls12_381 import G1, G2
# from py_ecc.bls12_381 import  pairing, add, multiply

from util_G1G2 import G1, G2, pairing, add, multiply, field_element

import random
import numpy as np
from polynomials import evaluate_univariate_polynomial_at_x, univariate_polynomial_division, mod_inverse, univariate_minus_polynomials,univariate_multiply_polynomials
import copy
from utils import gen_variable_power_list, gen_binary_list, gen_rho_plus_1_dim_array, get_power_cycle
from Lagrange_interpolation import univariate_lagrange_interpolation_coefficients
from poly_utils import PrimeField
# from utils import get_power_cycle
import math
from fft import fft
from aux_tree import build_sum_tree,TreeNode
import mcl
# from mcl import *
# import mcl.hook

BN254 = 0
BLS12_381 = 5
MCLBN_FR_UNIT_SIZE = 4
MCLBN_FP_UNIT_SIZE = 6

FP_SIZE = MCLBN_FP_UNIT_SIZE
FR_SIZE = MCLBN_FR_UNIT_SIZE + 2
G1_SIZE = MCLBN_FP_UNIT_SIZE * 3
G2_SIZE = MCLBN_FP_UNIT_SIZE * 6
GT_SIZE = MCLBN_FP_UNIT_SIZE * 12

SEC_SIZE = FR_SIZE * 2
PUB_SIZE = G1_SIZE + G2_SIZE
G1_CIPHER_SIZE = G1_SIZE * 2
G2_CIPHER_SIZE = G2_SIZE * 2
GT_CIPHER_SIZE = GT_SIZE * 4

MCLBN_COMPILED_TIME_VAR = (MCLBN_FR_UNIT_SIZE * 10) + MCLBN_FP_UNIT_SIZE


assert mcl.hook.mcl.mclBn_init(BLS12_381, MCLBN_COMPILED_TIME_VAR)==0

G1_STR = b"1 3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569"
G2_STR = b"1 352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160 3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758 1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905 927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582"




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


        # 秘密陷门，应当销毁 前rho个是tau_s,最后一个是tau_x
        self.alpha_s_x = [random.randint(2, modulus) for _ in range(self.rho +1)]
        # 产生一个二进制二维列表，大小为（M,rho）。里面是[0,M)的二进制分解
        self.binary_list = gen_binary_list(self.rho)
        # print("self.binary_list:", self.binary_list)

        # 产生一个rho+1 维度的list，其中前rho个度数为2，最后一个为n
        # print("rho:",self.rho)
        self.pp_PST = gen_rho_plus_1_dim_array(self.rho,self.n)
        # print("self.pp_PST:",self.pp_PST)
        for index, value in np.ndenumerate(self.pp_PST):
            temp = 1
            for i,v in enumerate(index):
                temp *= pow(self.alpha_s_x[i], v, self.modulus)
            self.pp_PST[index] = multiply(G1, temp% self.modulus)
        # pp_PST_pie指的是[eq(s,k)x^i]_{k\in [0,M),i\in [0,n)};size (2,...,2,n)
        self.pp_PST_pie = gen_rho_plus_1_dim_array(self.rho,self.n)
        # print(self.pp_PST_pie)
        for index, value in np.ndenumerate(self.pp_PST_pie):
            # print("index: ",index) # 指的是pp_PST_pie内部的index
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

        # pp_PST_pie_pie指的是[eq(s,k)x^i/(s_c-(1-<k>_c))]_{k\in [0,M),c\in[0,\rho),i\in [0,n)} ;size (2,2,rho,n) 这边是不是写错了？
        shape = tuple([2] * self.rho + [self.rho]+ [self.n])
        self.pp_PST_pie_pie =  np.full(shape, 0, dtype=object)
        for index, value in np.ndenumerate(self.pp_PST_pie_pie):
            # 这里的index[-2]指的就是c,index[index[-2]]指的就是k的二进制分解的第c位
            self.pp_PST_pie_pie[index] = multiply(self.pp_PST_pie[*index[:-2],index[-1]],
                                                  mod_inverse(self.alpha_s_x[index[-2]]-(1-index[index[-2]]), self. modulus) % self.modulus)
        print("Setup done")

    def dist_commit(self):
        '''
        分布式产生承诺
        '''
        # 每个证明者先产生子多项式f(x)
        # self.coefficients = [[0]*self.n]*self.M
        self.coefficients = [[0]*self.n for k in range(self.M)]
        self.C_v = multiply(G1,0)

        for k in range(self.M):
            index_k = self.binary_list[k]
            # self.coefficients size: (M,n)
            # self.coefficients[k] = univariate_lagrange_interpolation_coefficients(self.omega_n_s, self.subvector[k], self.modulus)[::-1] # 老版本使用拉格朗日插值，新版本使用fft
            self.coefficients[k] = fft(self.subvector[k], self.modulus, self.omega_n_s[1], inv=True)
            # 注意！！！！这里coefficients[k][0]是常数项，从小到大排序
            C_v_k = multiply(G1,0)
            for i in range(self.n):
                C_v_k = add(C_v_k,multiply(self.pp_PST_pie[*index_k,i],self.coefficients[k][i]))
            self.C_v = add(C_v_k,self.C_v)
        print("承诺C_v：", self.C_v)

    def genAux(self):
        '''
        产生AuxData; size: (M,rho)
        '''
        self.Aux = np.full((self.M,self.rho), 0, dtype=object)
        for k in range(self.M):
            index_k_binary = self.binary_list[k]
            for c in range(self.rho):
                self.Aux[k, c] = G1
                for i in range(self.n):
                    self.Aux[k,c] = add(self.Aux[k,c], multiply(self.pp_PST_pie_pie[*index_k_binary,c,i],self.coefficients[k][i]))
                self.Aux[k,c] =add(self.Aux[k,c],multiply(G1,mod_inverse(-1,self.modulus)))

    def genPartialProof(self): #显然每个机器的视角是等价的，工作量
        # 仅仅以P_0为视角 Q_{0,0} = {2,3}, Q_{0,1}={1}, n = 4
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
        print("partialProof_P_0 Done:",partialProof_P_0)
        return partialProof_P_0

    def genAllPartialProof(self):
        self.tree_list = [[0]*(2*self.M) for i in range(self.rho)]
        leaves = [[0]*self.M for i in range(self.rho)]
        for c in range(self.rho):
            for k in range(self.M):
                leaves[c][k] = self.Aux[k,c]
        # 构造rho 棵Aux-tree
        for c in range(self.rho):
            self.tree_list[c] = build_sum_tree(leaves[c])

        self.partialProofList = [[multiply(G1,0)]*self.rho for i in range(self.M)]
        for k in range(self.M):
            index_k_binary = copy.copy(self.binary_list[k])
            for c in range(self.rho):
                index = 0
                # if c == 0:
                #     index = 1-index_k_binary[0]
                # else:
                for q in range(c):
                    index += index_k_binary[q]*pow(2,c-q)
                index += 1-index_k_binary[c]
                self.partialProofList[k][c] = self.tree_list[c][self.rho - c -1][index].value + self.partialProofList[k][c]
        print("partialProof_P_0:",self.partialProofList[0])
        print("genAllPartialProof Done")

    def prove(self, d, i):
        index_d_binary_P_0 = self.binary_list[0]
        value = self.subvector[d][i]
        # 仅仅以P_0为视角
        # 需要计算多项式q(x) =  (witness_univariate_polynomial - value)/(x-\omega_n^)
        numerator_polynomial = self.coefficients[d]
        # 因为原来的代码是高位在前，没办法只能先倒置一下，等会再倒置回来
        numerator_polynomial = numerator_polynomial[::-1]
        numerator_polynomial[-1] = numerator_polynomial[-1] - value
        denominator_polynomial = [1, -self.omega_n_s[i]]
        # print("numerator_polynomial",numerator_polynomial)
        proof_polynomial, remainder_polynomial = univariate_polynomial_division(numerator_polynomial, denominator_polynomial, self.modulus)
        # print("证明多项式:", proof_polynomial)
        assert remainder_polynomial == [], "证明多项式生成失败"
        print("证明多项式生成成功！")
        # assert remainder_polynomial == [], "证明多项式生成失败"
        proof_polynomial = proof_polynomial[::-1] + [0]
        pi_d_rho_i = multiply(G1,0)
        for i in range(self.n):
            pi_d_rho_i =add(pi_d_rho_i, multiply(self.pp_PST_pie[*index_d_binary_P_0][i],proof_polynomial[i]))
        # commit_proof_polynomial = self.univariate_commit(proof_polynomial)
        # pi_d_rho_i = add(pi_d_rho_i, multiply(G1,mod_inverse(-1,self.modulus)))
        print("Prove done")
        print("pi_d_rho_i:",pi_d_rho_i)
        return value, pi_d_rho_i

    def verify(self, partialProof,pi_d_rho_i,value,d,i):
        index_binary_P_d = self.binary_list[d]
        verification_key_d = [0] * (self.rho+1)
        for c in range(self.rho):
            verification_key_d[c] = multiply(G2,self.alpha_s_x[c]-index_binary_P_d[c])
        verification_key_d[self.rho] = multiply(G2,(self.alpha_s_x[self.rho]-self.omega_n_s[i])% self.modulus)
        right = pairing(G2,multiply(G1,0))
        # 从这里开始计时
        print("value:",value)
        left =pairing( G2, add(self.C_v,multiply(self.pp_PST_pie[*index_binary_P_d,0], (-value)%self.modulus)))
        for c in range(self.rho):
            right = right * pairing( verification_key_d[c],partialProof[c])
        right = right * pairing( verification_key_d[self.rho],pi_d_rho_i)
        assert left == right, "Verification failed."
        print("Verification passed")

    def BatchProve(self, d, I):
        """
        I是下标列表
        """
        # 把子向量取出来
        I_list = [0]*len(I)
        for index, value  in enumerate(I):
            I_list[index] = self.subvector[d][value]
        # print("I_list", I_list)
        # 先计算A_I(x)=\prod_{i\in I} (x-\omega_n^i)
        A_I = [1]
        for i in I:
            A_I = univariate_multiply_polynomials(A_I,[1, -self.omega_n_s[i]],self.modulus)

        x_points = [self.omega_n_s[i] for i in I]
        # print("x_points",x_points)
        y_points = I_list
        r_x = univariate_lagrange_interpolation_coefficients(x_points, y_points, self.modulus)[::-1]
        # numerator_polynomial = f_k(x)
        numerator_polynomial = self.coefficients[d]
        # 因为原来的代码是高位在前，没办法只能先倒置一下，等会再倒置回来
        numerator_polynomial = numerator_polynomial[::-1]
        numerator_polynomial = univariate_minus_polynomials(numerator_polynomial, r_x[::-1], self.modulus)

        denominator_polynomial = A_I
        proof_polynomial, remainder_polynomial = univariate_polynomial_division(numerator_polynomial, denominator_polynomial, self.modulus)

        # print("证明多项式:", proof_polynomial)
        # print("remainder_polynomial:", remainder_polynomial)
        assert remainder_polynomial == [], "证明多项式生成失败"

        # 再换回来，变成低的在前面
        proof_polynomial = proof_polynomial[::-1]

        pi_d_rho_I_batch = multiply(G1,0)

        index_k_binary = self.binary_list[d]
        for i in range(self.n-len(I)):
            pi_d_rho_I_batch =add(pi_d_rho_I_batch, multiply(self.pp_PST_pie[*index_k_binary][i],proof_polynomial[i]))

        # print("pi_d_rho_I_batch:",pi_d_rho_I_batch)
        print("BatchProve Done")
        return I_list, pi_d_rho_I_batch
        # return I_list, pi_d_rho_I_batch

    def ProveAll(self):
        #  https://github.com/sunblaze-ucb/eVSS
        for i in range(3*int(math.log2(self.n))):
            index_k_binary_P_0 = self.binary_list[0]
            value = self.subvector[0][i]
            # 需要计算多项式q(x) =  (witness_univariate_polynomial - value)/(x-\omega_n^)
            numerator_polynomial = self.coefficients[0]
            # 因为原来的代码是高位在前，没办法只能先倒置一下，等会再倒置回来
            numerator_polynomial = numerator_polynomial[::-1]
            numerator_polynomial[-1] = numerator_polynomial[-1] - value
            denominator_polynomial = [1, -self.omega_n_s[i]]
            # print("numerator_polynomial",numerator_polynomial)
            proof_polynomial, remainder_polynomial = univariate_polynomial_division(numerator_polynomial, denominator_polynomial, self.modulus)
            # print("证明多项式:", proof_polynomial)
            assert remainder_polynomial == [], "证明多项式生成失败"
            print("证明多项式生成成功")
            proof_polynomial = proof_polynomial[::-1] + [0]
            pi_d_rho_i = G1
            for i in range(self.n):
                pi_d_rho_i =add(pi_d_rho_i, multiply(self.pp_PST_pie[*index_k_binary_P_0][i],proof_polynomial[i]))
            # commit_proof_polynomial = self.univariate_commit(proof_polynomial)
            pi_d_rho_i = add(pi_d_rho_i, multiply(G1,mod_inverse(-1,self.modulus)))
        print("ProveAll Done")

    def BatchVerify(self, partialProof, pi_d_rho_I_batch, value_list, d, I):
        """
        I是下标列表
        """
        # 先计算A_I(x)=\prod_{i\in I} (x-\omega_n^i)
        A_I = [1]
        for i in I:
            A_I = univariate_multiply_polynomials(A_I,[1, -self.omega_n_s[i]],self.modulus)
        # 计算verification key, 通常来说这个是在setup里面被计算的，为了方便，我们写在这里。
        # 当然它有更合理的生成方式，不过不重要。
        # 计算A_I的PST承诺
        A_I = A_I[::-1] # A_I:[-1,1]
        PST_A_I = multiply(G2,0)

        # index_0_binary = self.binary_list[0] # [0, 0, 0, 0]!!!!!!!!注意这里不要动，就是这样的，我们要取到与s无关的量
        for i in range(len(A_I)):
            # PST_A_I =add(PST_A_I, multiply(self.pp_PST[*index_0_binary][i],A_I[i]))
            PST_A_I =add(PST_A_I, multiply(multiply(G2, pow(self.alpha_s_x[self.rho],i, self.modulus)),A_I[i]))

        # 用户k的二进制展开\langle k \rangle
        index_k_binary = self.binary_list[d]
        # 新建验证密钥，长度为rho+1
        verification_key_d = [0] * (self.rho + 1)
        for c in range(self.rho):
            verification_key_d[c] = multiply(G2,self.alpha_s_x[c]-index_k_binary[c])
        # verification_key_d[self.rho] = multiply(G2,(self.alpha_s_x[self.rho]-self.omega_n_s[i])% self.modulus)

        # 正式开始
        # 再计算r(x),使得r(\omega_n^i)  = v_{d,i},r_x是高位在前
        x_points = [self.omega_n_s[i] for i in I]
        # print("x_points",x_points)
        y_points = value_list
        r_x = univariate_lagrange_interpolation_coefficients(x_points, y_points, self.modulus)[::-1]

        # 计算r_x的PST承诺
        PST_r_x = multiply(G1,0)
        index_k_binary = self.binary_list[d]
        for i in range(len(r_x)):
            PST_r_x =add(PST_r_x, multiply(self.pp_PST_pie[*index_k_binary][i],r_x[i]))
            # PST_r_x += self.pp_PST_pie[*index_k_binary][i] * r_x[i]

        # print("A_I:",A_I)
        # print("PST_A_I:",PST_A_I)
        # print("r_x:",r_x)
        # print("PST_r_x:",PST_r_x)

        right = pairing(G2,multiply(G1,0)) # 初始化，置零
        left =pairing( G2, add(self.C_v,multiply(PST_r_x, -1%self.modulus)))
        for c in range(self.rho):
            right = right * pairing( verification_key_d[c],partialProof[c])
        right = right * pairing( PST_A_I, pi_d_rho_I_batch)

        print("left:",left)
        print("right:",right)

        print("Verification right:",right)
        assert left == right, "Batch Verification failed."
        print("Batch Verification passed")

    def aggregate(self, Single_proof_list, I):
        # setup 阶段，可复用， 因此不算计算时间
        # $A_I^{\prime}(\omega^i)=\sum_{j \in I} \prod_{k \in I,k\neq j}(\omega^i-\omega^k)=\prod_{k \in I,k\neq i}(\omega^i-\omega^k)$
        # c_i = A_I^{\prime}(\omega^i)
        c_list =[1] * len(I)

        for index,i in enumerate(I):
            # 循环1，计算每个c_i
            for j in I:
                if j==i:
                    continue
                else:
                    c_list[index] *= (self.omega_n_s[i]-self.omega_n_s[j])
                    c_list[index] = c_list[index] % self.modulus
                c_list[index] = mod_inverse(c_list[index],self.modulus) % self.modulus


        # 正式开始聚合证明
        pi_d_rho_I_Agg = multiply(G1,0)
        for i in range(len(I)):
            pi_d_rho_I_Agg = add(pi_d_rho_I_Agg,multiply(Single_proof_list[i],c_list[i]))
            # pi_d_rho_I_Agg += (c_list[i] * Single_proof_list[i])%self.modulus

        # print("理论值",((self.omega_n_s[2]-self.omega_n_s[3])*Single_proof_list[0]+(self.omega_n_s[3]-self.omega_n_s[2]) *Single_proof_list[1])%self.modulus)
        return pi_d_rho_I_Agg

    def UpdateAuxTree(self,u,j,delta_v,upk_aux_u_j):
        # UpdateSetup(self)
        for c in range(self.rho):
            self.Aux[u,c] = add(self.Aux[u,c],multiply(upk_aux_u_j[0],delta_v))
        # 更新AUX和Aux树本质上开销是一样的。

    def UpdateProof(self,u,delta_v,pi_u_rho_i, a_i, a_j,i,j):
        w_i_j = add(multiply(a_i,mod_inverse(self.omega_n_s[i]-self.omega_n_s[j],self.modulus)),
                    multiply(a_j,mod_inverse(self.omega_n_s[j]-self.omega_n_s[i],self.modulus)))

        u_i_j = multiply(w_i_j,mod_inverse(self.n* mod_inverse(self.omega_n_s[i],self.modulus),self.modulus))

        new_pi_u_rho_i = add(pi_u_rho_i, multiply(u_i_j,delta_v))

    def UpdateSetup(self,u,j,i):
        '''
        产生update key
        '''
        # 产生upk_commitment_u_j_G1
        vector_temp = [0]*self.n
        vector_temp[j] = 1
        Lagrange_basis_u = fft(vector_temp, self.modulus, self.omega_n_s[1], inv=True)
        # 注意！！！！这里Lagrange_basis_u[0]是常数项，从小到大排序
        upk_commitment_u_j = evaluate_univariate_polynomial_at_x(Lagrange_basis_u,self.alpha_s_x[-1],self.modulus)
        index_u_binary = self.binary_list[u]

        temp = 1
        for i in range(len(index_u_binary)-1):
            if index_u_binary[i]==0:
                temp *= (1-self.alpha_s_x[i])
            else:
                temp *= self.alpha_s_x[i]
        upk_commitment_u_j *= temp
        upk_commitment_u_j_G1 = multiply(G1,upk_commitment_u_j)

        # 产生upk_Aux_u_j_G1; size: rho
        # upk_aux指的是update aux 的 key eq(s,k)L_i(x)/(s_c-(1-<k>_c)) ;size (2,2,rho,n)
        upk_aux_u_j_G1 = [0]*self.rho
        for c in range(self.rho):
            upk_aux_u_j_G1[c] = multiply(upk_commitment_u_j_G1,mod_inverse(self.alpha_s_x[c]-(1-index_u_binary[c]), self. modulus) % self.modulus)

        # 产生P_u update 元素i的 key:
        # a_i = g^{F(tau_s,\tau_x)/(tau_x-\omega^i)}
        a_i = multiply(self.C_v, mod_inverse(self.alpha_s_x[-1]-self.omega_n_s[i], self. modulus) % self.modulus)
        # a_j 一样的
        a_j = multiply(self.C_v, mod_inverse(self.alpha_s_x[-1]-self.omega_n_s[j], self. modulus) % self.modulus)

        shape = tuple([2] * self.rho + [self.rho]+ [self.n])
        self.upk_aux =  np.full(shape, 0, dtype=object)
        for index, _ in np.ndenumerate(self.pp_PST_pie_pie):
            self.upk_aux[index] = multiply(self.pp_PST_pie[*index[:-2],index[-1]],
                                                  mod_inverse(self.alpha_s_x[index[-2]]-(1-index[index[-2]]), self. modulus) % self.modulus)
        print("Update key generation done")


        return  upk_commitment_u_j_G1, upk_aux_u_j_G1, a_i, a_j

def AggregateTest(example,d, partialProof_P):
    # 这里给出的是单机聚合版本，至于跨subvectors 聚合，我们使用的是这个代码库测试数据:Hyperproofs
    # 。
    d = 0 # By default
    I = [2,3]
    Single_proof_list = [0]*len(I)
    value_list = [0]*len(I)
    for i in range(len(I)):
        value_list[i],Single_proof_list[i] = example.prove(d,I[i])
        print("value_list[i],Single_proof_list[i]",value_list[i],Single_proof_list[i])
        example.verify(partialProof_P, Single_proof_list[i], value_list[i], d, I[i])
    pi_d_rho_I_Agg = example.aggregate(Single_proof_list,I) # 这步没问题啊
    # print("!!!!!!!!!!!!!!", pi_d_rho_I_Agg, value_list, d, I)
    # print("pi_d_rho_I_Agg",pi_d_rho_I_Agg)
    example.BatchVerify(partialProof_P, pi_d_rho_I_Agg, value_list, d, I)


def UpdateTest(example, u, j ,delta_v):
    i = 1
    value, pi_u_rho_i = example.prove(u,i) # 这个是可以提前计算的，不算在计算时间中
    upk_commitment_u_j, upk_aux_u_j, a_i, a_j = example.UpdateSetup(u,j,i) # 产生upk，实际上它可以被提前产生
    # upk_aux_u_j 用来更新AUX树，  a_i, a_j 用来更新pi_rho
    example.UpdateAuxTree(u,j,delta_v,upk_aux_u_j)

    example.UpdateProof(u,delta_v,pi_u_rho_i, a_i, a_j,i,j)
    print("Update done!")
    return





if __name__ == "__main__":
    # M is fixed to 16
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
    # vector = [random.randint(1, modulus) for _ in range(N)]
    vector = [i+1 for i in range(N)]
    # vector = [i+1 for i in range(16)]
    # print(vector_long)
    example = Hybrid_mul_polynomial_commitment_scheme(M, n, N, rho,omega_n_s, modulus, vector)
    example.dist_commit()
    print("-----------------------------------------------generate partial proofs----------------------------------------")
    example.genAux()
    partialProof_P_0 = example.genPartialProof()
    example.genAllPartialProof()
    # print("-----------------------------------------------generate single proof----------------------------------------")
    # k = 0
    # d = k # d和k差不多
    # i = 1
    # value, pi_d_rho_i = example.prove(k,i)
    # value = vector[k*n+i]
    # print("value:", value)
    # example.verify(partialProof_P_0, pi_d_rho_i, value,k,i)
    # print("-----------------------------------------------generate Batch proof----------------------------------------")
    # I=[0,1]
    # value_list, pi_d_rho_I_batch = example.BatchProve(k,I)
    # print("pi_d_rho_I_batch", pi_d_rho_I_batch)
    # example.BatchVerify(partialProof_P_0, pi_d_rho_I_batch, value_list, d, I)
    # print("-----------------------------------------------Aggregate proof----------------------------------------")
    # AggregateTest(example,d, partialProof_P_0)
    # print("-----------------------------------------------update proof----------------------------------------")
    UpdateTest(example, u=0,j=0,delta_v=1)

