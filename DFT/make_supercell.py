#!/usr/bin/env python

"""
    Generate supercell coordinates from a base file.
"""

import math
import numpy
import sys
import argparse


def txt_to_num(x):
    return x.split()[0], map(float, x.split()[1:4]), map(int, x.split()[4:7])

def read_coords(fs):
    """
    Assuming the coordinate date files always hava exactly
        one header and one footer line
    """
    sys.stderr.write('Processing ' + fs.name + '\n\n')
    rawdata = fs.readlines()
    fs.close()
    return map(txt_to_num, rawdata[1:-1])

def make_simple_sc(ca, nb, base_coords):
    for i in range(len(nb)):
        nr = nb[i][0] * ca[0] + nb[i][1] * ca[1]
        for j in range(len(base_coords)):
            r = base_coords[j][1]
            f = base_coords[j][2]
            sys.stdout.write('   {0}{1:17.9f}{2:15.9f}{3:15.9f}'.\
                format(base_coords[j][0], r[0]+nr[0], r[1]+nr[1], r[2]))
            if f[2] == 0:
                print '      0  0  0'
            else:
                print

# Main Program

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                description='Supercell Coordinates Generator')
    parser.add_argument('coordfile', type=argparse.FileType('r'),
                help='Base coordinates file')
    parser.add_argument('-t', '--type', required=True,
                help='e.g.: Hr3xr3, H2x2, S2x2, S4x4;\
                The first letter stands for a hex or square supercell.')
    args = parser.parse_args()

    base_coords = read_coords(args.coordfile)

    sc_base = args.type[0]
    sc_size = args.type[1:]

    # Define base primitive vectors
    if sc_base == 'H': 
        ca = numpy.array([[1.,0.],[-0.5,math.sqrt(3.)/2]])
    elif sc_base == 'S':
        ca = numpy.array([[1.,0.],[0.,1.]])
    else:
        sys.stderr.write('Supercell type ' + args.type + ' is undefined.\n\n')
        exit(1)

    # Define cell vectors a's and shift vectors b's
    #        in units of the base primitive vectors
    if sc_size == 'r3xr3':
        na = numpy.array([[2.,1.],[1.,2.]])
        nb = numpy.array([[0.,0.],[1.,1.],[2.,2.]])
    elif sc_size == '2x2':
        na = numpy.array([[2.,0.],[0.,2.]])
        nb = numpy.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.]])
    else:
        sys.stderr.write('Supercell size ' + args.type + ' is undefined.\n\n')
        exit(1)

    print 'CELL_PARAMETERS alat'
    c1 = na[0][0] * ca[0] + na[0][1] * ca[1]
    c2 = na[1][0] * ca[0] + na[1][1] * ca[1]
    print '{0:16.9f}{1:15.9f}{2:15.9f}'.format(c1[0], c1[1], 0.)
    print '{0:16.9f}{1:15.9f}{2:15.9f}'.format(c2[0], c2[1], 0.)
    print '{0:16.9f}{1:15.9f}{2:15.9f}\n'.format(0., 0., 0.)

    print 'ATOMIC_POSITIONS alat'
    make_simple_sc(ca, nb, base_coords)
