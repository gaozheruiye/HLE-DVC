import mcl
from mcl import *
import mcl.hook
import mcl_bls2381PCS_group
from mcl_bls2381PCS_group import AggregateTest, UpdateTest
import random
import numpy as np


from utils import gen_variable_power_list, gen_binary_list, gen_rho_plus_1_dim_array, get_power_cycle
from Lagrange_interpolation import univariate_lagrange_interpolation_coefficients
from poly_utils import PrimeField
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


rho = 4
M = pow(2,rho)
# 调整M后需要调整genPartialProof()函数或者使用example.genAllPartialProof()然后把对应machine的取出来

# After adjusting M, you need to either modify the genPartialProof() function accordingly, or use example.genAllPartialProof() and extract the one corresponding to the target machine.
n = 32
N = M * n
modulus = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
f = PrimeField(modulus)
#单位根
omega_n = f.exp(7, (modulus-1)//n)
omega_n_s = get_power_cycle(omega_n, modulus)
# print ("omega_n_s:",omega_n_s)
# vector = [random.randint(1, modulus) for _ in range(N)]
vector = [random.randint(1, modulus) for i in range(N)]
# vector = [i+1 for i in range(16)]
# print(vector_long)
example = mcl_bls2381PCS_group.Hybrid_mul_polynomial_commitment_scheme(M, n, N, rho,omega_n_s, modulus, vector)
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
value = vector[k*n+i]
print("value:", value)
example.verify(partialProof_P_0, pi_d_rho_i, value,k,i)
example.ProveAll()
print("-----------------------------------------------generate Batch proof----------------------------------------")
I=[0,1]
value_list, pi_d_rho_I_batch = example.BatchProve(k,I)
print("pi_d_rho_I_batch", pi_d_rho_I_batch)
example.BatchVerify(partialProof_P_0, pi_d_rho_I_batch, value_list, d, I)
print("-----------------------------------------------Aggregate proof----------------------------------------")
AggregateTest(example,d, partialProof_P_0)
print("-----------------------------------------------update proof----------------------------------------")
UpdateTest(example,u=0,j=0,delta_v=5)# the updated subvector and the element index and v' - v