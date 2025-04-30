# Copyright (C) 2020-2022 Intel Corporation
# SPDX-License-Identifier: Apache-2.0

import sys
sys.path.append('../../utility')
import kcommute.hamlib_interface as hamlib_snippets

import os
from copy import deepcopy
import h5py

# import numpy as np
import math
import networkx as nx
import mat2qubit as m2q
from itertools import product
from openfermion import count_qubits

import datetime





def print_timestamp(start,flush=True):
    now = datetime.datetime.now()
    print("Elapsed: ",str(now - start) ,flush=flush)

# ***********************
# Bose-Hubbard
# ***********************
def bose_hubbard_bosop(thegraph,ssid_shift=0):

    # H = -t \sum b_i^ b_j + (U/2) \sum n_i(n_i-1)  {- \mu \sum n_i}

    bhOp = m2q.qSymbOp()

    for (I,J) in thegraph.edges():
        bhOp += m2q.qSymbOp("[ad_{0} a_{1}] ++ [ad_{1} a_{0}]".format(I,J))
       
    for I in thegraph.nodes():
        #print(f"I: {I}" )
        bhOp += (m2q.qSymbOp("(U/2) [n_{}]".format(I)))*(m2q.qSymbOp("[n_{}] ++ (-1) []".format(I)))

    return bhOp


# def parse_filename_vars(fname):
#     '''Function to return vars from file'''

#     fname = fname[:-4] # remove ".dat"
#     varpairs = fname.split('_')
#     ftype = varpairs[0]
#     varpairs = dict([ a.split('-') for a in varpairs[1:]  ])

#     return ftype,varpairs


def generate_all_bh(graphs_fname, outdir, start_keystr=None): #, max_nqub=math.inf, min_nqub=0):

    # Start time
    start = datetime.datetime.now()

    encodings = ('stdbinary','gray','unary')
    dvals = (4,8) # Doing only 4 & 8
    # (Uvals specified below)

    ctr = 0 
    still_skipping = False # Determines whether to skip until start_keystr
    if start_keystr is not None:
        still_skipping = True

    # for fname in os.listdir(graphdir):
    for graph_keystr in hamlib_snippets.get_hdf5_keys(graphs_fname):
        
        print("--",graph_keystr,flush=True)
        if still_skipping:
            if graph_keystr == start_keystr:
                still_skipping = False
            else:
                print('skipping')
                continue

        # Get properties for this graph instance
        graph_prop = hamlib_snippets.process_keystring(graph_keystr)


        # Get graph
        thegraph = hamlib_snippets.read_graph_hdf5(graphs_fname,graph_keystr[1:-1])
        nnodes = thegraph.number_of_nodes()

        # # Calc number of nodes
        # nnodes = 1
        # nnodes *= graph_prop['Lx']
        # if 'Ly' in graph_prop:
        #     nnodes *= graph_prop['Ly']
        # if 'Lz' in graph_prop:
        #     nnodes *= graph_prop['Lz']
        # nnodes = int(nnodes)
        # # Ignore if too many nodes
        maxnodes = math.inf
        # if nnodes>maxnodes:
        #     print(f"Skipping {graph_keystr} because nnodes={nnodes} > maxnodes={maxnodes}")
        #     continue



        # Ignore if too many nodes
        # maxnodes = 40 # ****************************************** <------
        # if nnodes>maxnodes:
        #     print(f"Skipping {graph_keystr} because nnodes={nnodes} > maxnodes={maxnodes}")
        #     continue

        
        # Create symbolic Ham (uses symbolic vars)
        symb_bh = bose_hubbard_bosop(thegraph)
        
        # Uvals are dependent on dimensionality
        if '1D' in graph_prop['probclass']:
            dimstr = '1D'
            Uvals = (2,10,20,30,40)
        elif '2D' in graph_prop['probclass']:
            dimstr = '2D'
            Uvals = (10,30,50,70,100)
        elif '3D' in graph_prop['probclass']:
            dimstr = '3D'
            Uvals = (20,40,60,90,120)
        else:
            raise Exception(f"graph_type should contain one of ('1D','2D','3D'). graph_type={graph_type}")

        # Save all BH models into one HDF5 file (Uvals are still in symbolic form)
        fname_symbops = f"{outdir}/all-bh-symbolic-unencoded.hdf5"
        hamlib_snippets.save_mat2qubit_hdf5(symb_bh,fname_symbops,graph_keystr)


        for enc,d,U in product(encodings,dvals,Uvals):
            
            qub_per_site = math.log2(d) if enc in ('stdbinary','gray') else d

            bh = deepcopy(symb_bh)
            vals = {'U':U}
            bh.scalar_subs(vals) 
            # print(bh)

            # Convert to mat2qubit compositeOperator
            ssid_order = [str(i) for i in range(nnodes)] # ['0','1','2']
            dlev_obj = m2q.symbop_to_dlevcompositeop(bh, ssid_order,d,enc)

            # Convert to QubitOperator
            qub_op = dlev_obj.opToPauli()
            nqub = count_qubits(qub_op)

            # Categorize based on dimensionality & d (levels-per-boson)
            bh_encoded_fname = f"{outdir}/BH_D-{dimstr[0]}_d-{d}.hdf5"
            haminst_name = "bh_" + graph_keystr[1:-1] + f"_U-{U}_enc-{enc}_d-{d}"
            # print(haminst_name)
            hamlib_snippets.save_openfermion_hdf5(qub_op,bh_encoded_fname,haminst_name)



    
    

if __name__=="__main__":
    
    print("===== Starting bh.py =====")
    # graphdir = "../condmat_graphs/"
    outdir   = "."

    # Grids file
    fname_grids = "Grids.hdf5"

    keystr_start = None
    if len(sys.argv)>1:
        keystr_start = sys.argv[1]
        print(f"*** Starting from keystr={keystr_start}")

    generate_all_bh(fname_grids,outdir, keystr_start)




    




