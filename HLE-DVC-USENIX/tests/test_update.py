import io
import random
import unittest
from contextlib import redirect_stdout

from mcl_bls2381PCS_group import Hybrid_mul_polynomial_commitment_scheme, UpdateTest
from poly_utils import PrimeField
from utils import get_power_cycle


class UpdateTests(unittest.TestCase):
    def _new_example(self):
        rho = 2
        m = 2 ** rho
        n = 4
        modulus = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
        field = PrimeField(modulus)
        omega_n = field.exp(7, (modulus - 1) // n)
        omega_n_s = get_power_cycle(omega_n, modulus)
        random.seed(123)
        vector = [random.randint(1, modulus) for _ in range(m * n)]
        return Hybrid_mul_polynomial_commitment_scheme(m, n, m * n, rho, omega_n_s, modulus, vector)

    def test_update_keeps_position_proofs_valid(self):
        with redirect_stdout(io.StringIO()):
            example = self._new_example()
            example.dist_commit()
            example.genAux()
            partial_before = example.genAllPartialProof()

            other_value, other_pi = example.prove(1, 0)
            example.verify(partial_before[1], other_pi, other_value, 1, 0)

            partial_after, updated_proofs = UpdateTest(example, u=2, j=1, delta_v=7)

            example.verify(partial_after[1], other_pi, other_value, 1, 0)
            for i, updated_pi in updated_proofs.items():
                value, recomputed_pi = example.prove(2, i)
                example.verify(partial_after[2], updated_pi, value, 2, i)
                self.assertEqual(updated_pi, recomputed_pi)


if __name__ == "__main__":
    unittest.main()
