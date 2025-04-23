#!/usr/bin/env python3

"""Defines k-qubit wise commuting."""

import math
from typing import Iterable, List

import cirq


def restrict_to(
    pauli: cirq.PauliString, qubits: Iterable[cirq.Qid]
) -> cirq.PauliString:
    return cirq.PauliString(p.on(q) for q, p in pauli.items() if q in qubits)


def compute_blocks(qubits: Iterable[cirq.Qid], k: int) -> List[List[cirq.Qid]]:
    return [qubits[k * i : k * (i + 1)] for i in range(math.ceil(len(qubits) / k))]


def commutes(pauli1: cirq.PauliString, pauli2: cirq.PauliString, blocks: List[List[cirq.Qid]]) -> bool:
    """Returns True if pauli1 k-commutes with pauli2, else False.

    Args:
        pauli1: A Pauli string.
        pauli2: A Pauli string.
        blocks: The block partitioning.
    """
    for block in blocks:
        if not cirq.commutes(restrict_to(pauli1, block), restrict_to(pauli2, block)):
            return False
    return True


if __name__ == "__main__":
    a, b, c, d = cirq.LineQubit.range(4)
    pauli1 = cirq.Z.on(a) * cirq.Z.on(b) * cirq.X.on(c) * cirq.X.on(d)
    pauli2 = cirq.X.on(a) * cirq.X.on(b) * cirq.X.on(c) * cirq.X.on(d)

    assert not commutes(pauli1, pauli2, k=1)
    assert commutes(pauli1, pauli2, k=2)
    assert commutes(pauli2, pauli2, k=3)
    assert commutes(pauli1, pauli2, k=4)

    a, b, c, d, e = cirq.LineQubit.range(5)
    pauli1 = cirq.Z.on(a) * cirq.I.on(b) * cirq.X.on(c) * cirq.I.on(d) * cirq.X.on(e)
    pauli2 = cirq.Z.on(a) * cirq.Z.on(b) * cirq.X.on(c) * cirq.X.on(d) * cirq.X.on(e)

    assert commutes(pauli1, pauli2, k=1)
    assert commutes(pauli1, pauli2, k=2)
    assert commutes(pauli1, pauli2, k=3)
    assert commutes(pauli1, pauli2, k=4)
    assert commutes(pauli1, pauli2, k=5)
