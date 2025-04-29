#!/bin/bash

fname="BH_D-1_d-4.hdf5"

keys=(
    bh_graph-1D-grid-pbc-qubitnodes_Lx-4_U-2_enc-unary_d-4
    bh_graph-1D-grid-pbc-qubitnodes_Lx-4_U-2_enc-stdbinary_d-4
    bh_graph-1D-grid-nonpbc-qubitnodes_Lx-10_U-2_enc-stdbinary_d-4
    bh_graph-1D-grid-nonpbc-qubitnodes_Lx-10_U-2_enc-unary_d-4
    bh_graph-1D-grid-nonpbc-qubitnodes_Lx-12_U-2_enc-stdbinary_d-4
    bh_graph-1D-grid-nonpbc-qubitnodes_Lx-12_U-2_enc-unary_d-4
    bh_graph-1D-grid-nonpbc-qubitnodes_Lx-16_U-10_enc-stdbinary_d-4
    bh_graph-1D-grid-nonpbc-qubitnodes_Lx-16_U-10_enc-unary_d-4
    bh_graph-1D-grid-nonpbc-qubitnodes_Lx-22_U-40_enc-stdbinary_d-4
    bh_graph-1D-grid-nonpbc-qubitnodes_Lx-22_U-40_enc-unary_d-4
)

for i in ${!keys[@]}; do
    echo "On $i / ${#keys[@]} ${keys[$i]}"
    python bose_hubbard.py $fname ${keys[$i]} 1
done