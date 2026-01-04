# The implementation of the HLE-DVC
USENIX 2026

The project can run on both Windows and Linux environments (Ubuntu version 20.04 or higher).

This project is built on the Python binding of the mcl library, 
which is originally implemented in Go. While the underlying library is written in Go, our project itself is implemented in Python.


The main.py in this folder can be executed directly, while the core code and key functions are implemented in the PCS_group.py.
Most of the remaining files are utility modules used to support polynomial operations.

Configurable Parameters

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
| Aggregate  | The whole process of Aggregation (including verification) is presented in AggregateTest(example,d, partialProof_P_0)  |
| UpdComProof  | UpdateTest(example,u=0,j=0,delta_v=5) |


### A small demo: How to generate a single position proof

```python
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
```

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
  - The algorithms and time overheads of `DistCommit`, `GenAux`, and `GenPartialProof` are **identical to those of HLE-DVC**, and only the `Setup` phase differs from HLE-DVC.
  - The `Verify` algorithm has the same cost as Hyperproofs.
  - The time overhead of `ProveAll` is computed as:
    ```
    (Hyperproofs ProveAll time) / M
    + time of GenAux
    + time of GenAllPartialProof
    ```

This combination of theoretical analysis and experimental measurement yields the results reported in Table 6.

