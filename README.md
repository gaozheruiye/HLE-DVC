# The implementation of the HLE-DVC
USENIX 2026

This artifact accompanies the USENIX Security paper titled “Distributed Vector Commitments and Their Applications”. 
We introduce a new notion—distributed VC (DVC), which allows multiple machines, 
each holding only a subvector of the input vector, to collectively commit to the entire vector and generate position proofs in a distributed manner. 
To the best of our knowledge, there is no prior work on DVCs and no existing work can trivially derive an efficient DVC scheme. We propose the first DVC scheme, HLE-DVC, conduct the experiments and open-source the code.

The artifact provides the experimental data and the source code of the implementation of HLE-DVC for reproducing the experimental results reported in the paper. 
In particular, the artifact supports reproducing the performance, scalability, and correctness evaluations presented in Sections 5.


The project can run on both Windows and Linux environments (Ubuntu version 20.04 or higher).
We evaluate the performance of the HLE-DVC scheme
using mcl libaray (Go language) on a Tencent cloud server
(16-core CPU, 64 GB memory, 5 Mbps bandwidth).
This implementation has no specific hardware requirements. However, larger-scale experiments require more memory.





This artifact provides the implementation of HLE-DVC based on the Python binding of the mcl library, 
which is originally implemented in Go. While the underlying library is written in Go, our project itself is implemented in Python.

The main.py in this folder can be executed directly, while the core code and key functions are implemented in the PCS_group.py.
Most of the remaining files are utility modules used to support polynomial operations.
  
The implementation is organized as follows:

* [`HLE-DVC-USENIX/`](HLE-DVC-USENIX): Implements most of the algorithms in HLE-DVC.
* [`AggCross/`](AggCross):  Benchmarks the performance of the cross-subvector aggregation in HLE-DVC.
* [`HLE-DVC-USENIX/main.py`](HLE-DVC-USENIX/main.py): An example, which present the complete proof generation process of HLE-DVC. Moreover, it includes update and single-subvector aggregation.
* [`HLE-DVC-USENIX/mcl_bls2381PCS_group.py`](HLE-DVC-USENIX/mcl_bls2381PCS_group.py): Contains implementations of most of the algorithms in HLE-DVC.
* [`HLE-DVC-USENIX/poly_utils.py`](HLE-DVC-USENIX/poly_utils.py) and [`HLE-DVC-USENIX/polynomials.py`](HLE-DVC-USENIX/polynomials.py): Supports finite field and polynomial computations.


**Configurable Parameters**

n: Length of each subvector

M: Number of machines

N: Length of the full vector

modulus: Modulus of the finite field. This parameter can be modified if the BLS12-381 curve is not used.




**Mapping Between Algorithms in the Paper and Code Functions**
(Note: Some algorithm inputs and outputs are internal attributes of classes and are not explicitly listed.)

|  algorithms of HLE-DVC   | functions in the code  |
|  ----  | ----  |
| Setup  | example = mcl_bls2381PCS_group.Hybrid_mul_polynomial_commitment_scheme(M, n, N, rho, omega_n_s, modulus, vector) <br>  omega_n_s is the multiplicative subgroup of the finite field|
| DistCommit  | example.dist_commit() |
| genAux  | example.genAux() |
| GenPartialProof  | partialProof_P_0 = example.genPartialProof() |
| GenAllPartialProof  | AllPartialProof = example.genAllPartialProof() |
| Prove  | value, pi_d_rho_i = example.prove(k,i) <br> k is the machine index and i is the elment index in the subvector|
| Verify  | example.verify(partialProof_P_0, pi_d_rho_i, value,k,i)    |
| ProveAll  | We measured the time overhead of ProveAll using https://github.com/sunblaze-ucb/eVSS |
| BatchProve  | value_list, pi_d_rho_I_batch = example.BatchProve(k,I)  <br> I is the index set of the subvector |
| Aggregate(Single)  | The whole process of Aggregation (including verification) is presented in AggregateTest(example,d, partialProof_P_0)  |
| Aggregate(Cross)  | The whole aggregation process (including verification) is presented in the folder "AggCross"|
| UpdComProof  | UpdateTest(example,u=0,j=0,delta_v=5) |




## Reproducing Experimental Results (Tables 4, 5, and 6)

This section provides instructions on how to reproduce the results reported in Tables 4, 5, and 6 of the Experiments section in the paper. We clarify which results are obtained via theoretical analysis and which are measured empirically using the provided code.

### Table 4: Storage Overhead of the HLE-DVC Scheme

The results in Table 4 report the **storage overhead** of the HLE-DVC scheme.  
These values are obtained via **theoretical calculation**, based on the construction of HLE-DVC and the sizes of group elements and proofs.  
No additional code execution is required to reproduce Table 4.

---

### Table 5: Time Overhead of the HLE-DVC Scheme

Table 5 reports the **time overhead** of HLE-DVC for different algorithms.

To reproduce the results in Table 5, set the parameters in the code as follows:
- `n = 2^30`
- `rho = 4`

The time overhead of **ProveAll** is measured using the implementation from the following repository:
- https://github.com/sunblaze-ucb/eVSS

The time overhead of **Agg(Cross)** is measured using the implementation in the folder "AggCross".

The reported values in Table 5 are obtained by running the corresponding algorithms with the above settings and averaging over multiple runs.

---

### Table 6: Comparison with Other Schemes

Table 6 compares HLE-DVC with KZG-DPCS-DVC, Hyperproofs, and Hyperproofs-DVC.

- **Storage size** results in Table 6 are obtained via **theoretical analysis** for all schemes.

- **KZG-DPCS-DVC**:
  - The time overheads are obtained through **theoretical analysis**, based on the costs of the underlying cryptographic operations.

- **Hyperproofs**:
  - The time overheads are obtained by setting:
    - `n = 1`
    - `rho = 34`
  - As pointed out in Appendix E.3 of the paper, Hyperproofs is essentially equivalent to the **MLE version of HLE-DVC** under these parameters.

- **Hyperproofs-DVC**:
  - Hyperproofs-DVC is a variant proposed in the *extension* section of the appendix and will be discussed further in future work.
  - The algorithms and time overheads of `DistCommit`, `GenAux`, and `GenPartialProof` are **identical to those of HLE-DVC**.
  - The `Verify` algorithm has the same cost as Hyperproofs.
  - The time overhead of `ProveAll` is computed as:
    ```
    (Hyperproofs ProveAll time) / M
    + time of GenAux
    + time of GenAllPartialProof
    ```

This combination of theoretical analysis and experimental measurement yields the results reported in Table 6.

### A small demo: How to generate a single position proof

```python
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

rho = 4
M = pow(2,rho)
n = 32
N = M * n

modulus = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
f = PrimeField(modulus)

omega_n = f.exp(7, (modulus - 1) // n)
omega_n_s = get_power_cycle(omega_n, modulus)

vector = [random.randint(1, modulus) for i in range(N)]

example = mcl_bls2381PCS_group.Hybrid_mul_polynomial_commitment_scheme(
    M, n, N, rho, omega_n_s, modulus, vector
)

example.dist_commit()

print("-----------------------------------------------generate partial proofs----------------------------------------")
example.genAux()

AllPartialProof = example.genAllPartialProof()
partialProof_P_0 = AllPartialProof[0]

print("-----------------------------------------------generate single proof----------------------------------------")
k = 0
i = 1

value, pi_d_rho_i = example.prove(k, i)
value = vector[k * n + i]

print("value:", value)
example.verify(partialProof_P_0, pi_d_rho_i, value, k, i)
print("-----------------------------------------------generate Batch proof----------------------------------------")
I=[0,1]
value_list, pi_d_rho_I_batch = example.BatchProve(k,I)
example.BatchVerify(partialProof_P_0, pi_d_rho_I_batch, value_list, d, I)
print("-----------------------------------------------Aggregate proof----------------------------------------")
AggregateTest(example,d, partialProof_P_0)
print("-----------------------------------------------update proof----------------------------------------")
UpdateTest(example,u=0,j=0,delta_v=5)
```

## The expected outputs:

```text
Setup time: 0.6898074999917299
Setup done
DistCommit time： 0.00542835624946747
Commitment C_v： <class 'mcl.structures.G1.G1'> 1 3690804181930197328350565422360803301969100863518835010754805841405630555150403607237889634527387431206600518712706 1055797299926905174841154417832477937643209798334699546578187888479526096414429831375883901145828008956946992067582
-----------------------------------------------generate partial proofs----------------------------------------
GenAux time： 0.020340199997008312
GenAux Done
GenAllPartialProof time： 0.0004296000115573406
GenAllPartialProof Done
-----------------------------------------------generate single proof----------------------------------------
witness polynomial q_{d,rho,i} generated successfully!
GenAllPartialProof time： 0.0051096000242978334
Prove done
value: 34474471475560831495890525465817897714957660052858741296885640342666713953404
Verification time: 0.0063587999902665615
Verification passed
-----------------------------------------------generate Batch proof----------------------------------------
BatchProve time: 0.005104599986225367
BatchProve Done
Batch Verification time: 0.009578299941495061
Batch Verification passed
-----------------------------------------------Aggregate proof----------------------------------------
witness polynomial q_{d,rho,i} generated successfully!
GenAllPartialProof time： 0.005799600039608777
Prove done
Verification time: 0.006426200037822127
Verification passed
witness polynomial q_{d,rho,i} generated successfully!
GenAllPartialProof time： 0.0052423999877646565
Prove done
Verification time: 0.0065484000369906425
Verification passed
Aggregate time: 0.0003554999129846692
Batch Verification time: 0.00981870002578944
Batch Verification passed
-----------------------------------------------update proof----------------------------------------
witness polynomial q_{d,rho,i} generated successfully!
GenAllPartialProof time： 0.00555010000243783
Prove done
Update key generation done
UpdateTime: 0.0010937000624835491
Update done!
```
