from typing import List
import numpy as np
from math import sqrt
from openfermion import qubit_operator_to_pauli_sum
from openfermion import QubitOperator
from quimb.tensor.tensor_1d import MatrixProductOperator, MatrixProductState
from .tensor_nets import mpo_mps_variance

def grouped_measurement_count(groupings, epsilon, circuit, simulator):
    # Based on the M_g metric in section 2 of 
    # https://quantum-journal.org/papers/q-2021-01-20-385/
    
    # Sum up the operators in each group and convert to PauliSums.
    group_operators = [sum(group) for group in groupings]
    group_pauli_sums = [qubit_operator_to_pauli_sum(op) for op in group_operators]
    group_pauli_squares = [op ** 2 for op in group_pauli_sums]
    # For each summed operator, get a variance from the sim result.
    exp_values = np.array(simulator.simulate_expectation_values(circuit, group_pauli_sums))
    squared_exp_values = np.array(simulator.simulate_expectation_values(circuit, group_pauli_squares))
    variances = squared_exp_values - exp_values ** 2
    # Now compute n_i and M_g
    n_i = (1.0 / epsilon ** 2) * np.sqrt(variances) * np.sum(np.sqrt(variances))
    m_g = np.sum(n_i)
    return m_g


def r_hat_measurement_count(groupings):
    r_numerator = 0
    r_denominator = 0
    for group in groupings:
        if isinstance(group, QubitOperator):
            a_ij = sum(list(group.terms.values()))
            r_numerator += abs(a_ij)
            r_denominator += sqrt(abs(a_ij) ** 2)
        else:
            # group_sum = 0
            # group_sum_squares = 0
            """
            for op in group:
                a_ij = sum(list(op.terms.values()))
                group_sum += abs(a_ij)
                group_sum_squares += abs(a_ij) ** 2
            """
            a_ij = np.array([list(op.terms.values())[0] for op in group])
            group_sum = np.sum(np.abs(a_ij))
            group_sum_squares = np.sum(np.abs(a_ij) ** 2)
            r_numerator += group_sum
            r_denominator += sqrt(group_sum_squares)
    return (r_numerator / r_denominator) ** 2


def get_variance(qubop,psi,nq=None):
    '''Returns variance <O^2> - <O>^2 w.r.t. given state'''

    # Ensure is instance
    assert isinstance(qubop, (QubitOperator,np.ndarray) ) or issparse(qubop)

    if isinstance(qubop,QubitOperator):
        if nq is None:
            nq = count_qubits(qubop)

        op = get_sparse_operator(qubop, nq)
        opsq = get_sparse_operator(qubop * qubop, nq)

    elif isinstance(qubop,(np.ndarray)) or issparse(qubop):
        opsq = qubop @ qubop
        op   = qubop

    mean        = psi.conj().T @ (op @ psi)
    mean_sq     = mean**2
    opsq_expect = psi.conj().T @ (opsq @ psi)

    return opsq_expect - mean_sq


def get_shotcounts_from_opsum(ops_to_sum,psi,epsilon,nq=None):
    """

    This allows for "overlapping" partitions, where the same Pauli string
    belongs to multiple partitions.

    Args:
        ops_to_sum (list of {QubitOperators or numpy matrix or sparse matrix}): the
            Hamiltonian = the sum of the list's members
        psi (numpy vector): state
        epsilon (float): desired error

    Returns:
        shot counts (float)
    """

    if nq is None:
        nq = int( np.log2(len(psi)) )

    temp = 0
    for op in ops_to_sum:
        # op can be QubitOperator or numpy matrix or sparse matrix
        var = get_variance(op,psi,nq)

        temp += np.sqrt(var)

    temp = temp**2
    shotcounts = temp.real / epsilon**2

    return shotcounts.real


def get_shotcounts_mpo_mps(mpos: List[MatrixProductOperator], mps: MatrixProductState, eps: float) -> float:
    """Get the shot count for the list of operators mpos w.r.t. the state
    mps, given the desired accuracy eps.
    
    Args:
    mpos - List of MPO's. Each represents a k-commuting group.
    The sum of the MPO's should equal the observable.
    mps - The state we are measuring the expectation for.
    eps - Desired accuracy. Should be a positive, floating point value."""

    temp = 0
    for mpo in mpos:
        var = mpo_mps_variance(mpo, mps)
        temp += np.sqrt(var)
    
    temp = temp ** 2
    shotcount = temp.real / eps**2
    return shotcount.real
