from typing import List
import numpy as np
import cirq
import openfermion as of
from quimb.tensor.tensor_1d import MatrixProductOperator, MatrixProductState
from quimb.tensor.tensor_1d_compress import tensor_network_1d_compress_direct

def pauli_string_to_mpo(pstring: cirq.PauliString, qs: List[cirq.Qid]) -> MatrixProductOperator:
    """Convert a Pauli string to a matrix product operator."""

    # Make a list of matrices for each operator in the string.
    ps_dense = pstring.dense(qs)
    matrices: List[np.ndarray] = []
    for pauli_int in ps_dense.pauli_mask:
        if pauli_int == 0:
            matrices.append(np.eye(2))
        elif pauli_int == 1:
            matrices.append(cirq.unitary(cirq.X))
        elif pauli_int == 2:
            matrices.append(cirq.unitary(cirq.Y))
        else: # pauli_int == 3
            matrices.append(cirq.unitary(cirq.Z))
    # Convert the matrices into tensors. We have a bond dim chi=1 for a Pauli string MPO.
    tensors: List[np.ndarray] = []
    for i, m in enumerate(matrices):
        if i == 0:
            tensors.append(m.reshape((2, 2, 1)))
        elif i == len(matrices) - 1:
            tensors.append(m.reshape((1, 2, 2)))
        else:
            tensors.append(m.reshape((1, 2, 2, 1)))
    return pstring.coefficient * MatrixProductOperator(tensors, shape="ludr")


def pauli_sum_to_mpo(psum: cirq.PauliSum, qs: List[cirq.Qid], max_bond: int) -> MatrixProductOperator:
    """Convert a Pauli sum to an MPO."""

    for i, p in enumerate(psum):
        if i == 0:
            mpo = pauli_string_to_mpo(p, qs)
        else:
            mpo += pauli_string_to_mpo(p, qs)
            tensor_network_1d_compress_direct(mpo, max_bond=max_bond, inplace=True)
    return mpo


def groups_of_to_mpos(
    groups_of: List[List[of.QubitOperator]], qs: List[cirq.Qid], max_bond: int
) -> List[MatrixProductOperator]:
    """Convert a list of lists of QubitOperators to a list of MPOs,
    where each MPO is the sum of a list of QubitOperators."""

    # Sum each list of QubitOperators into a single PauliSum,
    # then convert that to an MPO.
    group_mpos = []
    for group in groups_of:
        group_psum = cirq.PauliSum()
        for qubop in group:
            group_psum += of.transforms.qubit_operator_to_pauli_sum(qubop)
        group_mpos.append(pauli_sum_to_mpo(group_psum, qs, max_bond))
    return group_mpos


def mpo_mps_exepctation(mpo: MatrixProductOperator, mps: MatrixProductState, power: int=1) -> complex:
    """Get the expectation of an operator given the state."""

    assert power > 0

    mpo_times_mps = mps.copy()
    for i in range(power):
        mpo_times_mps = mpo.apply(mpo_times_mps)
    return mps.H @ mpo_times_mps


def mpo_mps_variance(mpo: MatrixProductOperator, mps: MatrixProductState) -> complex:
    """Get the variance of the operator mpo w.r.t. the state mps."""

    mpo_expectation = mpo_mps_exepctation(mpo, mps, power=1)
    mpo_squared_expectation = mpo_mps_exepctation(mpo, mps, power=2)
    return mpo_squared_expectation - mpo_expectation ** 2
