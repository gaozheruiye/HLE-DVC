import argparse
import csv
import io
import random
import time
from contextlib import redirect_stdout

from mcl_bls2381PCS_group import Hybrid_mul_polynomial_commitment_scheme
from poly_utils import PrimeField
from utils import get_power_cycle


MODULUS = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001


def quiet_call(func, *args, **kwargs):
    with redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


def time_call(rows, name, func, *args, **kwargs):
    start = time.perf_counter()
    result = quiet_call(func, *args, **kwargs)
    elapsed = time.perf_counter() - start
    rows.append((name, elapsed))
    return result


def make_example(rho, n, seed):
    m = 2 ** rho
    field = PrimeField(MODULUS)
    omega_n = field.exp(7, (MODULUS - 1) // n)
    omega_n_s = get_power_cycle(omega_n, MODULUS)
    random.seed(seed)
    vector = [random.randint(1, MODULUS) for _ in range(m * n)]
    return Hybrid_mul_polynomial_commitment_scheme(m, n, m * n, rho, omega_n_s, MODULUS, vector)


def run_benchmark(args):
    rows = []
    m = 2 ** args.rho
    if args.n <= 0 or args.n & (args.n - 1) != 0:
        raise ValueError("--n must be a power of two")
    if args.batch_size > args.n:
        raise ValueError("--batch-size must be <= --n")
    if args.d >= m or args.update_u >= m:
        raise ValueError("--d and --update-u must be smaller than M=2^rho")
    if args.i >= args.n or args.update_j >= args.n:
        raise ValueError("--i and --update-j must be smaller than n")

    example = time_call(rows, "Setup", make_example, args.rho, args.n, args.seed)
    time_call(rows, "DistCommit", example.dist_commit)
    time_call(rows, "GenAux", example.genAux)
    partial_proofs = time_call(rows, "GenAllPartialProof", example.genAllPartialProof)
    partial_proof_d = partial_proofs[args.d]

    value, proof = time_call(rows, "Prove", example.prove, args.d, args.i)
    time_call(rows, "Verify", example.verify, partial_proof_d, proof, value, args.d, args.i)

    batch_indices = list(range(args.batch_size))
    value_list, batch_proof = time_call(rows, "BatchProve", example.BatchProve, args.d, batch_indices)
    time_call(
        rows,
        "BatchVerify",
        example.BatchVerify,
        partial_proof_d,
        batch_proof,
        value_list,
        args.d,
        batch_indices,
    )

    single_values = []
    single_proofs = []
    for index in batch_indices:
        single_value, single_proof = quiet_call(example.prove, args.d, index)
        single_values.append(single_value)
        single_proofs.append(single_proof)
    aggregate_proof = time_call(rows, "Agg.Prove(single)", example.aggregate, single_proofs, batch_indices)
    time_call(
        rows,
        "Agg.Verify(single)",
        example.BatchVerify,
        partial_proof_d,
        aggregate_proof,
        single_values,
        args.d,
        batch_indices,
    )

    update_i = args.update_i
    if update_i is None:
        update_i = 1 if args.update_j != 1 else 0
    old_value_i, old_proof_i = quiet_call(example.prove, args.update_u, update_i)
    old_value_j, old_proof_j = quiet_call(example.prove, args.update_u, args.update_j)

    update_start = time.perf_counter()
    upk_commitment, upk_aux, upk_position, lagrange_basis_j = time_call(
        rows, "UpdateSetup", example.UpdateSetup, args.update_u, args.update_j
    )
    time_call(rows, "UpdateCommitment", example.UpdateCommitment, upk_commitment, args.delta)
    partial_after = time_call(rows, "UpdateAuxTree", example.UpdateAuxTree, args.update_u, args.update_j, args.delta, upk_aux)
    time_call(
        rows,
        "UpdateVectorAndPolynomial",
        example.UpdateVectorAndPolynomial,
        args.update_u,
        args.update_j,
        args.delta,
        lagrange_basis_j,
    )
    updated_proof_i = time_call(
        rows,
        "UpdateProof(i!=j)" if update_i != args.update_j else "UpdateProof(i=j)",
        example.UpdateProof,
        args.update_u,
        args.delta,
        old_proof_i,
        upk_position,
        update_i,
        args.update_j,
    )
    if update_i != args.update_j:
        time_call(
            rows,
            "UpdateProof(i=j)",
            example.UpdateProof,
            args.update_u,
            args.delta,
            old_proof_j,
            upk_position,
            args.update_j,
            args.update_j,
        )
    rows.append(("Update(total phases)", time.perf_counter() - update_start))

    time_call(
        rows,
        "Verify(updated)",
        example.verify,
        partial_after[args.update_u],
        updated_proof_i,
        example.subvector[args.update_u][update_i],
        args.update_u,
        update_i,
    )

    return rows, {
        "rho": args.rho,
        "M": m,
        "n": args.n,
        "N": m * args.n,
        "batch_size": args.batch_size,
        "d": args.d,
        "i": args.i,
        "update_u": args.update_u,
        "update_j": args.update_j,
        "update_i": update_i,
        "delta": args.delta,
        "seed": args.seed,
    }


def print_table(rows, params):
    print("HLE-DVC benchmark")
    print(
        "Parameters: "
        f"rho={params['rho']}, M={params['M']}, n={params['n']}, N={params['N']}, "
        f"|I|={params['batch_size']}, d={params['d']}, i={params['i']}, "
        f"update=(u={params['update_u']}, j={params['update_j']}, i={params['update_i']}), "
        f"delta={params['delta']}, seed={params['seed']}"
    )
    print()
    print(f"{'Function':<30} {'seconds':>12} {'ms':>12}")
    print("-" * 56)
    for name, elapsed in rows:
        print(f"{name:<30} {elapsed:12.6f} {elapsed * 1000:12.3f}")


def write_csv(path, rows, params):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["parameter", "value"])
        for key, value in params.items():
            writer.writerow([key, value])
        writer.writerow([])
        writer.writerow(["function", "seconds", "milliseconds"])
        for name, elapsed in rows:
            writer.writerow([name, f"{elapsed:.9f}", f"{elapsed * 1000:.6f}"])


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark HLE-DVC algorithm runtimes.")
    parser.add_argument("--rho", type=int, default=4, help="log2(M), default: 4")
    parser.add_argument("--n", type=int, default=32, help="subvector length, default: 32")
    parser.add_argument("--batch-size", type=int, default=3, help="batch/aggregation size, default: 3")
    parser.add_argument("--d", type=int, default=0, help="subvector index used for prove/batch/aggregate")
    parser.add_argument("--i", type=int, default=1, help="position index used for single prove/verify")
    parser.add_argument("--update-u", type=int, default=0, help="updated subvector index")
    parser.add_argument("--update-j", type=int, default=0, help="updated position index")
    parser.add_argument("--update-i", type=int, default=None, help="position proof updated for timing; default picks i != j")
    parser.add_argument("--delta", type=int, default=5, help="update delta")
    parser.add_argument("--seed", type=int, default=1234, help="random seed")
    parser.add_argument("--csv", default=None, help="optional CSV output path")
    return parser.parse_args()


if __name__ == "__main__":
    cli_args = parse_args()
    benchmark_rows, benchmark_params = run_benchmark(cli_args)
    print_table(benchmark_rows, benchmark_params)
    if cli_args.csv:
        write_csv(cli_args.csv, benchmark_rows, benchmark_params)
