from py_ecc.bls12_381 import G1, G2
from py_ecc.bls12_381 import  pairing, add, multiply

import random
import numpy as np
from polynomials import evaluate_univariate_polynomial_at_x, univariate_polynomial_division, mod_inverse, univariate_multiply_polynomials, univariate_minus_polynomials
import copy
from utils import gen_variable_power_list, gen_binary_list, gen_rho_plus_1_dim_array, get_power_cycle
from Lagrange_interpolation import univariate_lagrange_interpolation_coefficients
from poly_utils import PrimeField
# from utils import get_power_cycle
import galois
import time
from fft import fft
from aux_tree import build_sum_tree,TreeNode





# M: 证明者的人数
# l: 子向量长度
# rho：log M
# 下标k：M；下标i：l；下标c:rho
# 整个代码中所有的多项式都是[a, b, c], f= a + bx + cx^2这种，就是地位在前的，
# 但是多项式运算库中所有的代码都是高位在前，就是f= ax^2 + bx + c这种，写一半没办法改了
# 所以每次进行多项式运算的时候要先进行倒置，然后运算完了再倒置回来。
# 默认地，始终是低位在前面，即使出现高位在前面，也会立刻翻转成低位在前面
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
        print("self.vector:", self.vector)
        print("self.subvector:", self.subvector[0])
        # start_time = time.time()

        # # self.GF = galois.GF(self.modulus) # 伽罗瓦域，为了多项式运算
        # end_time = time.time()
        # print("galois field is generated, cost time:", end_time - start_time)

        # 秘密陷门，应当销毁 前rho个是tau_s,最后一个是tau_x
        # self.alpha_s_x = [random.randint(2, modulus) for _ in range(rho +1)]
        self.alpha_s_x = [3,4,5,6,7]
        # print("self.alpha_s_x:", self.alpha_s_x)
        # self.alpha_x


        # 产生一个二进制二维列表，大小为（M,rho）。里面是[0,M)的二进制分解
        self.binary_list = gen_binary_list(rho)
        # print("self.binary_list:", self.binary_list)

        # 产生一个rho+1 维度的list，其中前rho个度数为2，最后一个为n
        self.pp_PST = gen_rho_plus_1_dim_array(rho,n)
        for index, value in np.ndenumerate(self.pp_PST):
            # print("index: ",index)
            # print("value: ",value)
            temp = 1
            for i,v in enumerate(index):
                temp *= pow(self.alpha_s_x[i], v, self.modulus)
            self.pp_PST[index] = temp % self.modulus
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
            # print("temp: ",temp) # 指的是pp_PST_pie内部的index
            # print("index: ",index) # 指的是pp_PST_pie内部的index
            self.pp_PST_pie[index] = temp
        # print("self.pp_PST_pie",self.pp_PST_pie)

        # pp_PST_pie_pie ;size (2,2,rho,n)
        shape = tuple([2] * rho + [rho]+ [n])
        self.pp_PST_pie_pie =  np.full(shape, 0, dtype=object)
        for index, value in np.ndenumerate(self.pp_PST_pie_pie):
            self.pp_PST_pie_pie[index] = self.pp_PST_pie[*index[:-2],index[-1]] * mod_inverse(self.alpha_s_x[index[-2]]-(1-index[index[-2]]), self. modulus) % self.modulus

    def dist_commit(self):
        '''
        分布式产生承诺
        '''
        # 每个证明者先产生子多项式f(x)
        # self.coefficients 的size 是(M,n), self.coefficients是每个
        # self.coefficients = [[0]*self.n]*self.M
        self.coefficients = [[0]*self.n for k in range(self.M)]


        self.C_v = 0

        for k in range(self.M):
            index_k = self.binary_list[k]
            # self.coefficients[k] = univariate_lagrange_interpolation_coefficients(self.omega_n_s, self.subvector[k], self.modulus)[::-1] # 老版本使用拉格朗日插值，新版本使用fft
            self.coefficients[k] = fft(self.subvector[k], modulus, self.omega_n_s[1], inv=True)
            # 注意！！！！这里coefficients[k][0]是常数项，从小到大排序
            C_v_k = 0
            for i in range(self.n):
                C_v_k +=  self.coefficients[k][i] * self.pp_PST_pie[*index_k,i]
                # print("self.coefficients[k][i] * self.pp_PST_pie[*index_k,i]",self.coefficients[k][i] * self.pp_PST_pie[*index_k,i])
            # print("C_v_k:",C_v_k)
            self.C_v += C_v_k
            # print("承诺C_v_",k,":", C_v_k %self.modulus)
        self.C_v = self.C_v  % self.modulus
        # print("self.coefficients:",self.coefficients)
        print("承诺C_v:", self.C_v)

    def genAux(self):
        '''
        产生AuxData; size: (M,rho)
        '''
        self.Aux = np.full((self.M,self.rho), 0, dtype=object)
        # print("self.Aux.shape",self.Aux.shape)
        for k in range(self.M):
            index_k_binary = self.binary_list[k]
            for c in range(self.rho):
                self.Aux[k, c] = 0
                for i in range(self.n):
                    self.Aux[k,c] += self.coefficients[k][i] * self.pp_PST_pie_pie[*index_k_binary,c,i]
                self.Aux[k,c] =self.Aux[k,c] % self.modulus
                # print("k:",k)
                # print("c:",c)
                # print("self.Aux[k,c]:",self.Aux[k,c])

    def genPartialProof(self): #显然每个机器的视角是等价的，工作量
        # 仅仅以P_0为视角 Q_{0,0} = {2,3}, Q_{0,1}={1}
        # partialProof_P_0 = [self.Aux[2,0] + self.Aux[3,0], self.Aux[1,1]]
        # print("partialProof_P_0:",partialProof_P_0)
        # return partialProof_P_0

        # partialProof_P_0 = [self.Aux[2,0] + self.Aux[3,0], self.Aux[1,1]]
        # print("partialProof_P_0:",partialProof_P_0)
        # return partialProof_P_0
        tmp0 = self.Aux[8,0]+ self.Aux[9,0]
        tmp0 = tmp0+ self.Aux[10,0]
        tmp0 = tmp0+ self.Aux[11,0]
        tmp0 = tmp0+ self.Aux[12,0]
        tmp0 = tmp0+ self.Aux[13,0]
        tmp0 = tmp0+ self.Aux[14,0]
        tmp0 = tmp0+ self.Aux[15,0]

        tmp1 = self.Aux[4,1]+ self.Aux[5,1]
        tmp1 = tmp1+ self.Aux[6,1]
        tmp1 = tmp1+ self.Aux[7,1]

        partialProof_P_0 = [tmp0,tmp1,self.Aux[2,2]+ self.Aux[3,2], self.Aux[1,3]]
        print("partialProof_P_0 Done:",partialProof_P_0)
        return partialProof_P_0

    def genAllPartialProof(self):
        tree_list = [[0]*(2*self.M) for i in range(self.rho)]
        leaves = [[0]*self.M for i in range(self.rho)]
        for c in range(self.rho):
            for k in range(self.M):
                leaves[c][k] = self.Aux[k,c]
        # 构造rho 棵Aux-tree
        for c in range(self.rho):
            tree_list[c] = build_sum_tree(leaves[c])

        self.partialProofList = [[0]*self.rho for i in range(self.M)]
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
                self.partialProofList[k][c] += tree_list[c][self.rho - c -1][index].value
        print("partialProof_P_0:",self.partialProofList[0])
        print("genAllPartialProof Done")

    def prove(self, d, i):
        index_k_binary = self.binary_list[0]
        value = self.subvector[d][i]
        # 仅仅以P_0为视角
        # 需要计算多项式q(x) =  (witness_univariate_polynomial - value)/(x-\omega_n^)
        numerator_polynomial = self.coefficients[d]
        # 因为原来的代码是高位在前，没办法只能先倒置一下，等会再倒置回来
        numerator_polynomial = numerator_polynomial[::-1]
        numerator_polynomial[-1] = numerator_polynomial[-1] - value
        denominator_polynomial = [1, -self.omega_n_s[i]]
        # print("numerator_polynomial",numerator_polynomial)
        # proof_polynomial, remainder_polynomial = univariate_polynomial_division(numerator_polynomial, denominator_polynomial, self.modulus) 这是用老库做的
        proof_polynomial, remainder_polynomial = univariate_polynomial_division(numerator_polynomial, denominator_polynomial, self.modulus)
        # print("证明多项式:", proof_polynomial)
        assert remainder_polynomial == [], "证明多项式生成失败"
        proof_polynomial = proof_polynomial[::-1] + [0]
        pi_d_rho_i = 0
        for i in range(self.n):
            pi_d_rho_i += self.pp_PST_pie[*index_k_binary][i] * proof_polynomial[i]
        # commit_proof_polynomial = self.univariate_commit(proof_polynomial)
        pi_d_rho_i = pi_d_rho_i % self.modulus
        # print("pi_d_rho_i:",pi_d_rho_i)
        print("Prove Done")
        return value, pi_d_rho_i

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
        # 测试一下运算是否正确
        # print("A_I", A_I)

        # 再计算r(x),使得r(\omega_n^i)  = v_{d,i},r_x是高位在前

        x_points = [self.omega_n_s[i] for i in I]
        # print("x_points",x_points)
        y_points = I_list
        r_x = univariate_lagrange_interpolation_coefficients(x_points, y_points, self.modulus)[::-1]
        # 测试一下运算是否正确
        # print("r_x",r_x)
        # print("r_x(\omega_n^i)",evaluate_univariate_polynomial_at_x(r_x[::-1], self.omega_n_s[0], self.modulus))
        # print("r_x(\omega_n^i)",evaluate_univariate_polynomial_at_x(r_x[::-1], self.omega_n_s[1], self.modulus))
        # print("r_x(\omega_n^i)",evaluate_univariate_polynomial_at_x(r_x[::-1], self.omega_n_s[2], self.modulus))
        # print("r_x(\omega_n^i)",evaluate_univariate_polynomial_at_x(r_x[::-1], self.omega_n_s[3], self.modulus))

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
        pi_d_rho_I_batch = 0
        index_k_binary = self.binary_list[d]
        for i in range(self.n-len(I)):
            pi_d_rho_I_batch += self.pp_PST_pie[*index_k_binary][i] * proof_polynomial[i]
        # commit_proof_polynomial = self.univariate_commit(proof_polynomial)
        pi_d_rho_I_batch = pi_d_rho_I_batch % self.modulus
        # print("pi_d_rho_I_batch:",pi_d_rho_I_batch)
        print("BatchProve Done")
        return I_list, pi_d_rho_I_batch
        # return I_list, pi_d_rho_I_batch

    def verify(self, partialProof,pi_d_rho_i,value,d,i):
        left = self.C_v
        right = 0
        # 机器的下标
        index_binary_P = self.binary_list[d]
        # print("index_binary_P:",index_binary_P)

        for c in range(self.rho):
            right += (self.alpha_s_x[c]-index_binary_P[c])*partialProof[c]
        right += (self.alpha_s_x[self.rho]-self.omega_n_s[i])*pi_d_rho_i + value * self.pp_PST_pie[*index_binary_P,0]
        # print(pi_d_rho_i)
        # print("self.alpha_s_x[self.rho]",self.alpha_s_x[self.rho])
        # print("self.omega_n_s[i])", self.omega_n_s[i])
        # print("self.alpha_s_x[self.rho]-self.omega_n_s[i])",self.alpha_s_x[self.rho]-self.omega_n_s[i])
        # print("value * self.pp_PST_pie[*index_binary_P,0]",value * self.pp_PST_pie[*index_binary_P,0])

        right = right % self.modulus
        # print("Verification right:",right)
        assert left == right, "Verification failed."
        print("Verification passed")

    def BatchVerify(self, partialProof, pi_d_rho_I_batch, value_list, d, I):
        """
        I是下标列表
        """
        # 先计算A_I(x)=\prod_{i\in I} (x-\omega_n^i)
        A_I = [1]
        for i in I:
            A_I = univariate_multiply_polynomials(A_I,[1, -self.omega_n_s[i]],self.modulus)
        # 再计算r(x),使得r(\omega_n^i)  = v_{d,i},r_x是高位在前
        x_points = [self.omega_n_s[i] for i in I]
        # print("x_points",x_points)
        y_points = value_list
        r_x = univariate_lagrange_interpolation_coefficients(x_points, y_points, self.modulus)[::-1]
        # 计算A_I的PST承诺
        A_I = A_I[::-1] # A_I:[-1,1]
        PST_A_I = 0
        index_k_binary = self.binary_list[0] # [0, 0, 0, 0]!!!!!!!!注意这里不要动，就是这样的，我们要取到与s无关的量
        for i in range(len(A_I)):
            PST_A_I += self.pp_PST[*index_k_binary][i] * A_I[i]
        PST_A_I = PST_A_I % self.modulus

        # 计算r_x的PST承诺
        PST_r_x = 0
        index_k_binary = self.binary_list[d]
        for i in range(len(r_x)):
            PST_r_x += self.pp_PST_pie[*index_k_binary][i] * r_x[i]
        PST_r_x = PST_r_x  % self.modulus

        # print("A_I:",A_I)
        # print("PST_A_I:",PST_A_I)
        # print("r_x:",r_x)
        # print("PST_r_x:",PST_r_x)

        left = self.C_v
        right = 0
        index_k_binary = self.binary_list[d]
        # print("index_k_binary:",index_k_binary)
        for c in range(self.rho):
            right += (self.alpha_s_x[c]-index_k_binary[c])*partialProof[c]
        right += PST_A_I*pi_d_rho_I_batch + PST_r_x

        right = right % self.modulus
        # print("Verification right:",right)
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
        #             print("index",index)
        #             print("i",i)
        #             print("j",j)
        #             print("c_list[index]",c_list[index])
        # print("c_list",c_list)

        # 正式开始聚合证明
        pi_d_rho_I_Agg = 0
        for i in range(len(I)):
            pi_d_rho_I_Agg += (c_list[i] * Single_proof_list[i])%self.modulus

        # print("理论值",((self.omega_n_s[2]-self.omega_n_s[3])*Single_proof_list[0]+(self.omega_n_s[3]-self.omega_n_s[2]) *Single_proof_list[1])%self.modulus)

        pi_d_rho_I_Agg = pi_d_rho_I_Agg %self.modulus
        return pi_d_rho_I_Agg

    def update_same_machine(self,pi_u_rho_i, u,i, j):

        new_pi_u_rho_j = 0
        return

    def update_different_machine(self,pi_d_rho_i,u, d):
        # 其实就是改变partial proof

        self.partialProofList[k][c]


def AggregateTest(example, d, partialProof_P):
    d = 0 # By default
    I = [2,3]
    Single_proof_list = [0]*len(I)
    value_list = [0]*len(I)
    for i in range(len(I)):
        value_list[i],Single_proof_list[i] = example.prove(d,I[i])
        print("value_list[i],Single_proof_list[i]",value_list[i],Single_proof_list[i])
        example.verify(partialProof_P, Single_proof_list[i], value_list[i], d, I[i])
    pi_d_rho_I_Agg = example.aggregate(Single_proof_list,I) # 这步没问题啊
    print("!!!!!!!!!!!!!!", pi_d_rho_I_Agg, value_list, d, I)
    # print("pi_d_rho_I_Agg",pi_d_rho_I_Agg)
    example.BatchVerify(partialProof_P, pi_d_rho_I_Agg, value_list, d, I)


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
    print ("omega_n_s:",omega_n_s)
    # vector = [random.randint(1, modulus) for _ in range(N)]
    vector = [i+1 for i in range(N)]
    # print(vector_long)
    example = Hybrid_mul_polynomial_commitment_scheme(M, n, N, rho,omega_n_s, modulus, vector)
    example.dist_commit()
    print("-----------------------------------------------generate partial proofs----------------------------------------")
    example.genAux()
    partialProof_P_0 = example.genPartialProof()
    example.genAllPartialProof()
    print("-----------------------------------------------generate single proof----------------------------------------")
    k = 0
    d = k # d和k差不多
    i = 1
    value, pi_d_rho_i = example.prove(k,i)
    print("pi_d_rho_i",pi_d_rho_i)
    print("value",value)
    example.verify(partialProof_P_0, pi_d_rho_i, value,k,i)
    print("-----------------------------------------------generate Batch proof----------------------------------------")
    I=[2,3]
    value_list, pi_d_rho_I_batch = example.BatchProve(k,I)
    print("pi_d_rho_I_batch", pi_d_rho_I_batch)
    example.BatchVerify(partialProof_P_0, pi_d_rho_I_batch, value_list, d, I)
    print("!!!!!!!!!!!!!!", pi_d_rho_I_batch, value_list, d, I)
    print("-----------------------------------------------Aggregate proof----------------------------------------")
    AggregateTest(example,d, partialProof_P_0)
    print("-----------------------------------------------Update proof----------------------------------------")



