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
| 单元格  | 单元格 |
