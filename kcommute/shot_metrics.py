import numpy as np
from math import sqrt
from openfermion import qubit_operator_to_pauli_sum
from openfermion import QubitOperator

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
