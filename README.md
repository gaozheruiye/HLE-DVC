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
| Setup  | example = mcl_bls2381PCS_group.Hybrid_mul_polynomial_commitment_scheme(M, n, N, rho, omega_n_s, modulus, vector) |
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
