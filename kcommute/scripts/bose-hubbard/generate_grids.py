#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append('../../utility')


#!/usr/bin/env python
# -*- coding: utf-8 -*-
import networkx as nwx
import numpy as np
import math
import itertools
import sys
import h5py
from kcommute.hamlib_interface import save_graph_hdf5

nribbons = 30  # number of ribbons of every width

def _save_to_files(Gp, Gnp, grid_pos, prefix, postfix):
    ''' Save graph edges to files as specified in the document.
    Args:
        Gp:      graph with periodic boundary conditions
        Gnp:     graph with non periodic boundary conditions
        grid_pos:   dict of grid positions associated with each node in graph
        prefix:  prefix to filename
        postfix  postfix to filename
    '''
    print(f'size: {postfix}')
    
    entry_name = f"{prefix}-nonpbc-qubitnodes_{postfix}"
    entry_name_pbc = f"{prefix}-pbc-qubitnodes_{postfix}"

    save_graph_hdf5(Gp, "Grids.hdf5", entry_name_pbc, overwrite=True, grid_pos=grid_pos)
    save_graph_hdf5(Gnp, "Grids.hdf5", entry_name, overwrite=True, grid_pos=grid_pos)
    
    
    

    
        


def snaking_2d_square_grids(Gp_grid):
    node_array = np.array(list(Gp_grid.nodes))
    num_col = max(node_array[:,0])

    nodes_per_col = []
    for i in range(0,num_col+1):    
        nodes_per_col = nodes_per_col + [len(np.where(node_array[:,0] == i)[0])]
    
    num_nodes = Gp_grid.number_of_nodes()
    orig = list(range(num_nodes))
    new = []

    i = 0
    j = 0

    # define new labeling
    while i < num_nodes:
        n = nodes_per_col[j] 
        c = orig[i:i+n]
        if j % 2:
            c.reverse()
        new = new + c

        i += n
        j += 1

    # do relabeling
    mapping = {orig[i]:new[i] for i in range(0, num_nodes)}
    Gp_int = nwx.convert_node_labels_to_integers(Gp_grid) 
    Gp = nwx.relabel_nodes(Gp_int, mapping)
    
    # get grid positions
    nodes = list(Gp.nodes)
    gridnodes = list(Gp_grid.nodes)
    grid_pos = {str(nodes[i]):gridnodes[i] for i in range(0,num_nodes)}
    
    return Gp, grid_pos


def snaking_2d_triag_grids(Gp_grid):
    node_array = np.array(list(Gp_grid.nodes))
    num_col = max(node_array[:,1])

    nodes_per_col = []
    for i in range(0,num_col+1):    
        nodes_per_col = nodes_per_col + [len(np.where(node_array[:,1] == i)[0])]
    
    num_nodes = Gp_grid.number_of_nodes()
    orig = list(range(num_nodes))
    new = []

    i = 0
    j = 0
    while i < num_nodes:
        n = nodes_per_col[j] 
        c = orig[i:i+n]
        if j % 2:
            c.reverse()
        new = new + c

        i += n
        j += 1

    mapping = {orig[i]:new[i] for i in range(0, num_nodes)}
    Gp_int = nwx.convert_node_labels_to_integers(Gp_grid)
    Gp = nwx.relabel_nodes(Gp_int, mapping)
    
    # get grid positions
    nodes = list(Gp.nodes)
    gridnodes = list(Gp_grid.nodes)
    grid_pos = {str(nodes[i]):gridnodes[i] for i in range(0,num_nodes)}
    
    return Gp, grid_pos

        
def snaking_2d_hex_grids(Gp_grid):
    node_array = np.array(list(Gp_grid.nodes))
    num_col = max(node_array[:,0])

    nodes_per_col = []
    for i in range(0,num_col+1):    
        nodes_per_col = nodes_per_col + [len(np.where(node_array[:,0] == i)[0])]

    num_nodes = Gp_grid.number_of_nodes()
    orig = list(range(num_nodes))
    new = []

    i = 0
    j = 0
    while i < num_nodes:
        n = nodes_per_col[j] 
        c = orig[i:i+n]
        if j % 2 == 0:
            c.reverse()
        new = new + c

        i += n
        j += 1

    mapping = {orig[i]:new[i] for i in range(0, num_nodes)}
    Gp_int = nwx.convert_node_labels_to_integers(Gp_grid)
    Gp = nwx.relabel_nodes(Gp_int, mapping)
    
    # get grid positions
    nodes = list(Gp.nodes)
    gridnodes = list(Gp_grid.nodes)
    grid_pos = {str(nodes[i]):gridnodes[i] for i in range(0,num_nodes)}
    
    return Gp, grid_pos

        
def snaking_3d_square_grids(Gp_grid,n):
    orig = list(range(n**3))
    new = []

    for j in range(0,n):    
        d = orig[j*(n**2):(j+1)*(n**2)]
    
        if j % 2:
            d.reverse()
    
        for i in range(0,n):
            c = d[i*(n):(i+1)*n]
            if i % 2:
                c.reverse()
            new = new + c    

    mapping = {orig[i]:new[i] for i in range(0, n**3)}
    Gp_int = nwx.convert_node_labels_to_integers(Gp_grid)
    Gp = nwx.relabel_nodes(Gp_int, mapping)
    
    # get grid positions
    nodes = list(Gp.nodes)
    gridnodes = list(Gp_grid.nodes)
    grid_pos = {str(nodes[i]):gridnodes[i] for i in range(0,n**3)}
    
    return Gp, grid_pos




def generate_1d_grids(nmax=np.inf,nmin=0):
    print('Generating 1D grids ...')
    r1 = range(2, 10)
    r2 = range(10, 50, 2)
    r3 = range(50, 200, 10)
    r4 = range(200, 1001, 100)
    prefix = 'graph-1D-grid'
    postfix = lambda x : f'Lx-{x}'
    for i in itertools.chain(r1, r2, r3, r4):
        if i>nmax or i<nmin:
            continue        
            
        Gp_grid = nwx.grid_graph(dim=(i,), periodic=True)  
        Gnp_grid = nwx.grid_graph(dim=(i,), periodic=False)
        
        Gp = nwx.convert_node_labels_to_integers(Gp_grid)
        Gnp = nwx.convert_node_labels_to_integers(Gnp_grid)
        
        # get grid positions        
        nodes = list(Gnp.nodes)
        gridnodes = list(Gnp_grid.nodes)
        grid_pos = {str(nodes[i]):gridnodes[i] for i in range(0,i)}
        
        _save_to_files(Gp, Gnp, grid_pos, prefix, postfix(i))

        
        
def generate_2d_square_grids(): #(nmax=np.inf,nmin=0):
    print('Generating 2D square grids ...')
    prefix = 'graph-2D-grid'
    postfix = lambda x, y : f'Lx-{x}_Ly-{y}'
    # square boundaries
    for i in [5, 10, 15, 20, 25, 30]:
        # if i**2>nmax or i**2<nmin:
        #     continue
            
        Gp_grid = nwx.grid_graph(dim=(i, i), periodic=True)        
        Gnp_grid = nwx.grid_graph(dim=(i, i), periodic=False)
        
        Gp, grid_pos_p = snaking_2d_square_grids(Gp_grid)
        Gnp, grid_pos = snaking_2d_square_grids(Gnp_grid)
        
        _save_to_files(Gp, Gnp, grid_pos, prefix, postfix(i, i))
        
    # ribbon grids
    #for x in range(2, 6):
    #    maxval = int(np.ceil(1000 / x))
    #    for y in np.linspace(x, maxval, nribbons, dtype=int):
#   #          print(x,y)
#
#            # if x*y>nmax or x*y<nmin or x==y:
#            #     continue
#            Gp_grid = nwx.grid_graph(dim=(x, y), periodic=True)
#            Gnp_grid = nwx.grid_graph(dim=(x, y), periodic=False)
#            
#            Gp, grid_pos_p = snaking_2d_square_grids(Gp_grid)
#            Gnp, grid_pos = snaking_2d_square_grids(Gnp_grid)            
#            
#            _save_to_files(Gp, Gnp, grid_pos, prefix, postfix(x, y))

 
            
def generate_2d_triag_grids(): #(nmax=np.inf,nmin=0):
    print('Generating 2D triangular grids ...')
    prefix = 'graph-2D-triag'
    postfix = lambda x, y : f'Lx-{x}_Ly-{y}'
    # square boundaries
    for i in range(5, 47): # start at 5 to allow periodic; 46x46 in fact is ~1000 nodes
        # if i**2>1.5*nmax: # Eliminate most, before creating any graph
        #     continue
        Gp_grid = nwx.triangular_lattice_graph(m=i, n=i, periodic=True)
        Gp, grid_pos_p = snaking_2d_triag_grids(Gp_grid)
        
        
        #
        nn = Gp.number_of_nodes()
        print(f"nnodes: {nn}")
        if nn>1000:
            break
        #


        Gnp_grid = nwx.triangular_lattice_graph(m=i, n=i, periodic=False)
        Gnp, grid_pos = snaking_2d_triag_grids(Gnp_grid)
        
        _save_to_files(Gp, Gnp, grid_pos,prefix, postfix(i, i))
        
    # # ribbon grids
    # for x in range(3, 6):
    #     maxyval = int(np.ceil(1000 / x))
    #     for y in np.linspace(x+2, maxyval, nribbons, dtype=int):
    #         # if x*y>1.5*nmax: # Eliminate most, before creating any graph
    #         #     continue
    #         Gp_grid = nwx.triangular_lattice_graph(m=x, n=y, periodic=True)
    #         Gp, grid_pos_p = snaking_2d_triag_grids(Gp_grid)
    #         #
    #         nn = Gp.number_of_nodes()
    #         if nn>nmax or nn<nmin:
    #             continue
    #         #
            
    #         Gnp_grid = nwx.triangular_lattice_graph(m=x, n=y, periodic=False)
    #         Gnp, grid_pos = snaking_2d_triag_grids(Gnp_grid)
            
    #         _save_to_files(Gp, Gnp, grid_pos, prefix, postfix(x, y))


def generate_2d_hex_grids(): #(nmax=np.inf,nmin=0):
    print('Generating 2D hexagonal grids ...')
    prefix = 'graph-2D-hex'
    postfix = lambda x, y : f'Lx-{x}_Ly-{y}'
    # square boundaries
    for i in range(2, 33, 2):  # even n required for periodic
        # if i**2>1.5*nmax: # Eliminate most, before creating any graph
        #     continue
        Gp_grid = nwx.hexagonal_lattice_graph(m=i, n=i, periodic=True)
        Gp, grid_pos_p = snaking_2d_hex_grids(Gp_grid)
        
        #
        nn = Gp.number_of_nodes()
        print(f"nnodes: {nn}")
        if nn>1000:
            break
        #
        Gnp_grid = nwx.hexagonal_lattice_graph(m=i, n=i, periodic=False)
        Gnp, grid_pos = snaking_2d_hex_grids(Gnp_grid)
        
        _save_to_files(Gp, Gnp, grid_pos, prefix, postfix(i, i))
        
#     # ribbon grids
#     def round_up_to_even(f): # only even n is allowed
#         return math.ceil(f / 2.) * 2
#     for x in range(2, 6):
#         maxval = int(np.ceil(1000 / x))
#         for y in np.linspace(x, maxval, nribbons):
#             y = round_up_to_even(y)
# #             print(x,y)
#             if x*y>1.5*nmax or x == y: # Eliminate most, before creating any graph
#                 continue
#             Gp_grid = nwx.hexagonal_lattice_graph(m=x, n=y, periodic=True)
#             Gp, grid_pos_p = snaking_2d_hex_grids(Gp_grid)
#             #
#             nn = Gp.number_of_nodes()
#             if nn>nmax or nn<nmin:
#                 continue
#             #
#             Gnp_grid = nwx.hexagonal_lattice_graph(m=x, n=y, periodic=False)
#             Gnp, grid_pos = snaking_2d_hex_grids(Gnp_grid)
            
#             _save_to_files(Gp, Gnp, grid_pos, prefix, postfix(x, y))


def generate_3d_square_grids(nmax=np.inf,nmin=0):
    print('Generating 3D square grids ...')
    prefix = 'graph-3D-grid'
    postfix = lambda x, y, z : f'Lx-{x}_Ly-{y}_Lz-{z}'
    # cube grids
    for i in range(2, 11):
        if i**3>nmax or i**3<nmin:
            continue
        Gp_grid = nwx.grid_graph(dim=(i, i, i), periodic=True)
        Gnp_grid = nwx.grid_graph(dim=(i, i, i), periodic=False)
        
        Gp, grid_pos_p = snaking_3d_square_grids(Gp_grid,i)
        Gnp, grid_pos = snaking_3d_square_grids(Gnp_grid,i)
                
        _save_to_files(Gp, Gnp, grid_pos, prefix, postfix(i, i, i))  



def main(graph_dataset_name=None):
    
    if graph_dataset_name is None:
        generate_1d_grids()
        generate_2d_square_grids()
        #generate_2d_triag_grids()
        #generate_2d_hex_grids()
        #generate_3d_square_grids()
        return

    # if graph_dataset_name=="30u":
    #     maxn = 30
    #     generate_1d_grids(maxn)
    #     generate_2d_square_grids(maxn)
    #     generate_2d_triag_grids(maxn)
    #     generate_2d_hex_grids(maxn)
    #     generate_3d_square_grids(maxn)





if __name__ == "__main__":
    if len(sys.argv)>1:
        main(sys.argv[1])
    else:
        main()





# nmax = np.inf
# generate_1d_grids(nmax)
# generate_2d_square_grids(nmax)
# generate_2d_triag_grids(nmax)
# generate_2d_hex_grids(nmax)
# generate_3d_square_grids(nmax)