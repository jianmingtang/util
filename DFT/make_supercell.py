#!/usr/bin/env python3

"""
    Generate supercell coordinates from a base file.
    bc_pv: base cell primitive vectors
    sc_pv: supercell primitive vectors
    sc_basis: supercell atomic basis vectors
    base_coords: position data in the base primitive cell
"""

import math
import numpy
import sys
import argparse


def txt_to_num(x):
    """
    Input: a list contains 7 or 4 items. (See QE output)
    Output: a list of three lists (1 string, 3 float, 3 int)
    0: atom type
    1-3: atom coordinates
    4-6: atom move flag
    """
    xs = x.split()
    if len(xs) == 7:
        flags = list(map(int, xs[4:7]))
    elif len(xs) == 4:
        flags = [1,1,1]
    else:
        sys.stderr.write('Error: The coordiante file is bad.\n\n')
        exit(1)
    return xs[0], list(map(float, xs[1:4])), flags

def read_coords(fs):
    """
    Assuming the coordinate date files always hava exactly
        one header and one footer line
    """
    sys.stderr.write('Processing ' + fs.name + '\n')
    rawdata = fs.readlines()
    fs.close()
    sys.stderr.write('Number of layers/atoms: %i\n\n' % (len(rawdata)-2))
    return list(map(txt_to_num, rawdata[1:-1]))

def print_atomic_position(pos, rs=[0,0]):
    """
    rs: in-plane shift of the atoms
    """
    r = pos[1]
    sys.stdout.write('   {0}{1:17.9f}{2:15.9f}{3:15.9f}'.\
        format(pos[0], r[0]+rs[0], r[1]+rs[1], r[2]))
    if pos[2] == [0,0,0]:
        print(' '*6 + '0  0  0')
    else:
        print()

def make_simple_sc(bc_pv, sc_basis, base_coords):
    """
    Use the base coordinates to fill the supercell.
    The number of atomic layers remains the same.
    """
    for vb in sc_basis:
        # in-plane shift of the atoms
        rs = vb[0] * bc_pv[0] + vb[1] * bc_pv[1]
        for pos in base_coords:
            print_atomic_position(pos, rs)

def make_merged_sc(bc_pv, sc_basis, base_coords, sub_coords, Ns):
    """
    Assume the substructure has the same supercell structure as the target
    Use the base coordinates to fill the supercell.
    The number of atomic layers remains the same.
    """
    # Number of atoms per sc layer:
    na_sc = len(sc_basis)
    # Number of base layers:
    nl_base = len(base_coords) 
    # Number of substructure layers:
    nl_sub = (len(sub_coords) - Ns*2) // na_sc
    if nl_sub >= nl_base:
        print('Nothing to merge!')
        exit(0)

    # start the index base at the first "cell"
    sub_atom = Ns * 2

    for i in range(Ns):
        pos = sub_coords[i]
        print_atomic_position(pos)

    # Find out the size of half of the substructure
    moved_layer = 0
    for i in range(nl_sub):
        pos = sub_coords[sub_atom+i]
        if pos[2] == [0,0,0]:
            break
        else:
            moved_layer += 1

    # Work out the shifts
    zs = pos[1][2] - base_coords[moved_layer][1][2]
    zs += base_coords[nl_base-moved_layer-1][1][2]
    zs -= sub_coords[sub_atom+nl_sub-moved_layer-1][1][2]

    for i in range(Ns):
        pos = sub_coords[Ns+i]
        pos[1][2] += zs
        print_atomic_position(pos)

    # For each "cell" in the supercell
    for vb in sc_basis:
        rb = vb[0] * bc_pv[0] + vb[1] * bc_pv[1]

        # print out the first half of the substructure
        for i in range(moved_layer):
            pos = sub_coords[sub_atom+i]
            print_atomic_position(pos)

        # fill in the base coordinates
        zs = sub_coords[sub_atom+moved_layer][1][2]
        zs -= base_coords[moved_layer][1][2]
        for i in range(moved_layer, nl_base-moved_layer):
            pos = base_coords[i]
            pos[1][2] += zs
            print_atomic_position(pos, rb)

        # print out the second half of the substructure
        zs = pos[1][2] - sub_coords[sub_atom+nl_sub-moved_layer-1][1][2]
        for i in range(nl_sub-moved_layer, nl_sub):
            pos = sub_coords[sub_atom+i]
            pos[1][2] += zs
            print_atomic_position(pos)

        # move the index base to the next "cell"
        sub_atom += nl_sub


# Main Program

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                description='Supercell Coordinates Generator')
    parser.add_argument('basefile', type=argparse.FileType('r'),
                help='Base coordinates file')
    parser.add_argument('-t', '--type', required=True,
                help='Supercell type: Hr3xr3, H2x2, S2x2, S4x4;\
                The first letter stands for a hex or square supercell.')
    parser.add_argument('-m', '--merge', dest='subfile',
                type=argparse.FileType('r'),
                help='Merge a substructure into the supercell.')
    parser.add_argument('-n', metavar='overlayers', type=int, default=0,
                help='Number of atoms in the overlayer (each side)')
    args = parser.parse_args()

    # base_coords: atomic coordinates in the base primitive cell
    base_coords = read_coords(args.basefile)

    sc_type = args.type[0]
    sc_size = args.type[1:]

    # Define base-cell primitive vectors (bc_pv)
    if sc_type == 'H': 
        bc_pv = numpy.array([[1.,0.],[-0.5,math.sqrt(3.)/2]])
    elif sc_type == 'S':
        bc_pv = numpy.array([[1.,0.],[0.,1.]])
    else:
        sys.stderr.write('Supercell type ' + args.type + ' is undefined.\n\n')
        exit(1)

    # Define supercell primitive vectors (sc_pv) and
    #       atomic basis vectors (sc_basis)
    #       in the basis of base primitive vectors (bc_pv)
    if sc_size == 'r3xr3':
        sc_pv = numpy.array([[2.,1.],[1.,2.]])
        sc_basis = numpy.array([[0.,0.],[1.,1.],[2.,2.]])
    elif sc_size == '2x2':
        sc_pv = numpy.array([[2.,0.],[0.,2.]])
        sc_basis = numpy.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.]])
    else:
        sys.stderr.write('Supercell size ' + args.type + ' is undefined.\n\n')
        exit(1)

    print('CELL_PARAMETERS alat')
    for va in sc_pv:
        # Expand the supercell vectors in Cartesian coordinates
        sc = va[0] * bc_pv[0] + va[1] * bc_pv[1]
        print('{0:16.9f}{1:15.9f}{2:15.9f}'.format(sc[0], sc[1], 0.))
    # Placeholder for the z direction
    print('{0:16.9f}{1:15.9f}{2:15.9f}\n'.format(0., 0., 0.))

    print('ATOMIC_POSITIONS alat')
    if args.subfile != None:
        sub_coords = read_coords(args.subfile)
        make_merged_sc(bc_pv, sc_basis, base_coords, sub_coords, args.n)
    else:
        make_simple_sc(bc_pv, sc_basis, base_coords)
