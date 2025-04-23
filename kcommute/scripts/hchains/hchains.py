from typing import Dict
import argparse
import json

import cirq
import openfermion

from kcommute import commute, get_si_sets, r_hat_measurement_count
from kcommute.hamlib_interface import read_openfermion_hdf5



def rhat_vs_k(ham, k_min, step) -> Dict[int, float]:
    """Compute Rhat for k values in [k_min, k_min + step, ..., num_qubits]"""
    nq = openfermion.utils.count_qubits(ham)
    qubits = cirq.LineQubit.range(nq)
    rhats = {}
    ks = range(k_min, openfermion.utils.count_qubits(ham) + 1, step)
    for k in ks:
        print("On k =", k)
        blocks = commute.compute_blocks(qubits, k)
        # print("Blocks are:")
        # print(blocks)
        groupings = get_si_sets(ham, blocks=blocks)
        print(f"Finished grouping, there are {len(groupings)} groups.")
        # print("Groups are:")
        # print(groupings)
        rhats[k] = r_hat_measurement_count(groupings)
        print("rhat =", rhats[k])
    return rhats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str)
    parser.add_argument("output_file", type=str)
    parser.add_argument("k_min", type=int)
    parser.add_argument("step", type=int)
    parser.add_argument("threshold", type=float)
    args = parser.parse_args()

    n_hydrogens = int(args.input_file.split("_")[1][1:])   
    hamiltonian = read_openfermion_hdf5(args.input_file, f"ham_JW-{2 * n_hydrogens}" if n_hydrogens >= 40 else "ham_JW")
    nq = openfermion.utils.count_qubits(hamiltonian)
    print(f"Hamiltonian acts on {nq} qubit(s) has {len(hamiltonian.terms)} term(s).")

    hamiltonian.compress(float(args.threshold))
    print(f"Hamiltonian compressed to threshold {args.threshold} acts on {nq} qubit(s) has {len(hamiltonian.terms)} term(s).")

    k_rhat_dict = rhat_vs_k(hamiltonian, args.k_min, args.step)

    output_dict = {}
    output_dict["input_file"] = args.input_file
    output_dict["nqubits"] = nq
    output_dict["nterms"] = len(hamiltonian.terms)
    output_dict["rhat_vs_k"] = k_rhat_dict
    with open(args.output_file, "w", encoding="utf8") as f:
        json.dump(output_dict, f)


if __name__ == "__main__":
    main()
