from typing import Tuple, List

import argparse

from math import sqrt

import pandas as pd
import matplotlib.pyplot as plt

import cirq
import openfermion as of

from kcommute.sorted_insertion import get_si_sets
from kcommute.shot_metrics import r_hat_measurement_count
from kcommute.commute import compute_blocks
from kcommute.hamlib_interface import read_openfermion_hdf5, print_hdf5_structure

def get_rhats_ngroups(
    ham: of.QubitOperator, step: int, verbose: bool=False
) -> Tuple[List[int], List[int], List[float]]:
    """Get the rhat and ngroups values from the Hamiltonian
    for different k values."""

    nq = of.utils.count_qubits(ham)
    qubits = cirq.LineQubit.range(nq)
    print(f"H acts on {nq} qubits")
    print(f"H has {len(ham.terms)} terms")
    kvals = list(range(1, nq + 1, step))
    rhats = []
    ngroups = []

    for k in kvals:
        blocks = compute_blocks(qubits, k)
        if verbose:
            print("On k =", k)
            # print("Blocks are:")
            # print(blocks)
        
        groups = get_si_sets(ham, blocks=blocks)
        ngroups.append(len(groups))
        rhats.append(r_hat_measurement_count(groups))
        if verbose:
            print(f"Finished grouping, there are {len(groups)} groups.")
            # print("Groups are:")
            # print(groups)
            print("rhat =", rhats[-1])

    return kvals, ngroups, rhats

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, help="HDF5 file.")
    parser.add_argument('key', type=str, help="HDF5 key.")
    parser.add_argument('step', type=int, help="Step in k value.", default=1)
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()

    print("Reading in H")
    ham = read_openfermion_hdf5(args.filename, args.key)
    print("Read in H")
    kvals, ngroups, rhats = get_rhats_ngroups(ham, args.step, verbose=args.verbose)
    df = pd.DataFrame({"ngroups": ngroups, "rhat": rhats})
    df.index.name = "i"
    df.to_csv(args.key + ".csv")

    plt.plot(kvals, rhats, label="$\\hat{R}$")
    plt.plot(kvals, ngroups, label="# groups")

    plt.xlabel("$k$")
    plt.legend();
    plt.savefig(args.key + ".png")

if __name__ == "__main__":
    main()
