##!/usr/bin/env python3
"""filter_map.py: Applies beam and low- and high-pass filtering to an 
image FITS file (e.g. from generateMockMaps).
If used as main, expects as CLI args fits file name (including path),
output file name prefix, and path to save directory.
"""
import os, sys
sys.path.append('../ThumbStack')
# from flat_map import FlatMap
import numpy as np
import fitsio
from scipy.special import erf
from functools import partial

import flat_map
from flat_map import *

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Apply beam and LPF/HPF to a flatmap file.')
    
    parser.add_argument('--lpf_map',
                        help='Full Path to the LPF flatmap FITS file')
    parser.add_argument('--hpf_map',
                        help='Full Path to the HPF flatmap FITS file')
   parser.add_argument('--outpath',
                        help='path to directory to save output files')
   parser.add_argument('--filter_type', default='tauring',
                        help='Choose thumbstack filter (can be string or list of strings)'
   parser.add_argument('--filter_width', default=100.,
                        help='LPF/HPF filter width for erf function filtering (i.e. `theo` filter_type)')
    
    args = parser.parse_args()

    return args

def initialize():
    u = ts.universe.UnivMariana()
    massConversion = ts.mass_conversion.MassConversionKravtsov14()

    boxmap, boxmask = make_box()

    return u, massConversion, boxmap, boxmask

def make_box():

    return boxmap, boxmask

def get_catalog():

    return catalog

def setup_maps():

    return lpf_map, hpf_map

def run_thumbstack():

def main(argv):
    args = parse_args()
    
    initialize()


if __name__ == '__main__':
    main(sys.argv)
