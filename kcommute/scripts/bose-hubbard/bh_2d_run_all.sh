#!/bin/sh

fname="BH_D-2_d-4.hdf5"

keys=(
    bh_graph-2D-grid-nonpbc-qubitnodes_Lx-5_Ly-5_U-100_enc-unary_d-4
    bh_graph-2D-grid-nonpbc-qubitnodes_Lx-7_Ly-7_U-100_enc-unary_d-4
    bh_graph-2D-grid-nonpbc-qubitnodes_Lx-10_Ly-10_U-100_enc-unary_d-4
)

for i in ${!keys[@]}; do
    echo "On $i / ${#keys[@]} ${keys[$i]}"
    python bose_hubbard.py $fname ${keys[$i]} 1
done
