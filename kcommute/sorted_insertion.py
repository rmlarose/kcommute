# Sorted insertion (SI) algorithm for measurement reduction

"""
Algorithm description:
https://quantum-journal.org/papers/q-2021-01-20-385/
Crawford et al. Quantum 5, 385 (2021)

Note that this algo is deterministic.
"""

from openfermion import QubitOperator
from openfermion.transforms import qubit_operator_to_pauli_sum

import kcommute.commute


def get_terms_ordered_by_abscoeff(op):
    '''Returns terms of QubitOperator, ordered by abs(coeff)

    Args:
        op (QubitOperator)

    Returns:
        list of tuples
    '''

    # Ensure is instance
    assert isinstance(op, QubitOperator)
    
    # Order the terms by absolute val of coefficient
    terms = sorted(op.terms.items(), key=lambda x: abs(x[1]), reverse=True)
    # terms = [t[0] for t in terms]
    terms = [QubitOperator(t[0], t[1]) for t in terms]

    # Return terms
    return terms


def get_si_sets(op, blocks, verbosity: int = 0):
    '''Returns grouping from sorted insertion algo.

    TODO: Add docstring and type function.
    '''
    nterms = len(op.terms)

    # Basic assertions
    assert isinstance(op, QubitOperator)

    qo_to_ps = lambda qo: next(iter(qubit_operator_to_pauli_sum(qo)))
    comm_func = lambda ps1, ps2: kcommute.commute.commutes(qo_to_ps(ps1), qo_to_ps(ps2), blocks=blocks)

    # Commuting sets (as list datatype)
    commuting_sets = []

    # Order the terms by absolute val of coefficient, or shuffle randomly.
    terms_ord = get_terms_ordered_by_abscoeff(op)
    
    # Remove any identity operators from the list of terms.
    for op in terms_ord:
        if () in op.terms.keys():
            # This is just the identity operator. Remove it. 
            # Remember when have one operator per term now.
            terms_ord.remove(op)

    # Loop over terms
    for i, pstring in enumerate(terms_ord):
        if verbosity > 0:
            print(f"Status: On Pauli string {i} / {nterms}", end="\r")
        if verbosity > 1:
            print(f"There are currently {len(commuting_sets)} group(s) of terms.", end="\r")
        if verbosity > 2:
            print(f"The groups are:")
            for j, group in enumerate(commuting_sets):
                print(f"Group {j + 1}: {group}")
        found_commuting_set = False

        # Loop over existing commuting sets
        for commset in commuting_sets:
            """
            comm_checks = [comm_func(pstring,pstring2) for pstring2 in commset]

            # All must be true
            if all(comm_checks):
                # Add to set
                commset.append(pstring)
                found_commuting_set = True
                break
            """

            all_strings_in_commset_commute = True
            for pstring2 in commset:
                if not comm_func(pstring, pstring2):
                    all_strings_in_commset_commute = False
                    break
            if all_strings_in_commset_commute:
                found_commuting_set = True
                commset.append(pstring)
                break
            else:
                found_commuting_set = False
            
        if not found_commuting_set:
            # Create new commuting set
            commuting_sets.append([pstring])

    return commuting_sets
