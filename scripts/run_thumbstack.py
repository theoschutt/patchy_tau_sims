#!/usr/bin/env python3
"""run_thumbstack.py: Run thumbstack on a given set of filtered maps.
"""
import os, sys
import numpy as np
import fitsio
from scipy.special import erf
from functools import partial

sys.path.append('../../ThumbStack')
from headers import *
from flat_map import FlatMap
from universe import UnivMariana
from catalog import Catalog
from mass_conversion import MassConversionKravtsov14
from thumbstack import ThumbStack

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description='Run thumbstack on a given set of filtered maps.')
    parser.add_argument('--workdir',
        default='/home/theo/Documents/research/CMB/patchy_tau_sims',
        help='parent directory of input/output/figures directories')
    parser.add_argument('--lpfpath',
        help='Full path to the LPF flatmap FITS file')
    parser.add_argument('--test',
        default=False, action='store_const', const=True,
        help='Test script using first 10 objects in catalog')
    parser.add_argument('--hpfpath',
        help='Full path to the HPF flatmap FITS file')
    parser.add_argument('--catpath',
        help='Full path to input catalog file')
    parser.add_argument('--pix_scale',
        default=0.5,
        type=float,
        help='Pixel scale in arcmin/px')
    parser.add_argument('--side_length',
        default=10.,
        type=float,
        help='Name of catalog (used in path)')
    parser.add_argument('--ra_min',
        default=200.,
        help='Min RA of map')
    parser.add_argument('--dec_min',
        default=10.,
        help='Min DEC of map')
    parser.add_argument('--catname',
        help='Name of catalog (used in path)')
    parser.add_argument('--tsname',
        help='Name for thumbstack output files')
    parser.add_argument('--filtertype',
        default='tauring',
        help='Choose thumbstack filter (can be string or list of strings)')
    parser.add_argument('--esttype',
        default=['tau_ti_uniformweight', 'tau_sgn_uniformweight'],
        help=('Choose thumbstack estimator and weight'
              ' (can be string or list of strings)'))
    parser.add_argument('--dobootstrap',
        default=True,
        help='Boolean whether to compute bootstrap covariance')
    parser.add_argument('--dostackmap',
        default=False,
        help='Boolean whether to compute stacked map')
    parser.add_argument('--equalsignedweights',
        default=False, action='store_const', const=True,
        help='Whether to equalize number of +/- weights')
    parser.add_argument('--nproc',
        default=1,
        help='Number of processes for parallel computing')
    parser.add_argument('--outpath',
        help='Non-default path to directory to save output files')
 
    args = parser.parse_args()

    return args

def initialize():
    u = UnivMariana()
    massConversion = MassConversionKravtsov14()

    return u, massConversion

def make_basemap(sizeX, sizeY, pixel_scale):
    # sizeX, sizeY: map dimensions in degrees

    # number of pixels for the flat map
    nX = int(sizeX * 60. / pixel_scale)
    nY = int(sizeY * 60. / pixel_scale)

    # basic map object
    baseMap = FlatMap(nX=nX, nY=nY,
        sizeX=sizeX*np.pi/180., sizeY=sizeY*np.pi/180.)
    
    return baseMap

def make_box(side_length=10., ra_min=200., dec_min=10., pix_scale=0.5):
    """Generate empty square map.
    side_length, ra_min, and dec_min in degrees
    pix_scale is in arcmin/px
    """
    ra_max = ra_min + side_length
    dec_max = dec_min + side_length
    # convention for defining box corners: [[dec_min, ra_max],[dec_max, ra_min]]
    box = np.array([[dec_min, ra_max],[dec_max, ra_min]]) * utils.degree 
    shape, wcs = enmap.geometry(pos=box, res=pix_scale * utils.arcmin,
        proj='car')

    # create a mask that keeps the whole area
    boxmask = enmap.ones(shape, wcs=wcs)

    # while we're here, let's make the base flat_map
    nX = int(side_length * 60. / pix_scale)
    nY = int(side_length * 60. / pix_scale)
    basemap = FlatMap(nX=nX, nY=nY, sizeX=side_length*np.pi/180.,
        sizeY=side_length*np.pi/180.)

    return shape, wcs, boxmask, basemap

def setup_maps(lpf_path, hpf_path, side_length=10., ra_min=200.,
    dec_min=10., pix_scale=0.5):
    """Thumbstack needs the input maps as enmaps
    """
    shape, wcs, boxmask, basemap = make_box(
        side_length, ra_min, dec_min, pix_scale)
    lpf_map = basemap.copy() # flatMap loads in place
    hpf_map = basemap.copy()
    lpf_map.read(lpf_path)
    hpf_map.read(hpf_path)

    assert lpf_map.data.shape == shape
    assert hpf_map.data.shape == shape
    lpf_enmap = enmap.enmap(lpf_map.data, wcs=wcs)
    hpf_enmap = enmap.enmap(hpf_map.data, wcs=wcs)

    return lpf_enmap, hpf_enmap, boxmask

def main():
    args = parse_args()

    # Write text file logging what command line args were used
    outpath = os.path.join(args.workdir, "output/thumbstack/"+args.tsname)
    if not os.path.exists(outpath):
       os.makedirs(outpath)
    log_fn = os.path.join(outpath, '%s_args.log'%args.tsname)
    print('Writing argument log file:', log_fn)
    arg_dict = vars(args)
    with open(log_fn, 'w') as f:
        for arg in arg_dict:
            f.write('%s: %s\n'%(str(arg), str(arg_dict[arg])))

    u, massConv = initialize()

    if args.test:
        nObj = 11 #ensures unequal pos/neg weights
        do_bootstrap = False
        do_stacked_map = True
    else:
        nObj = None
        do_stacked_map = args.dostackmap
        do_bootstrap = args.dobootstrap

    galcat = Catalog(
        u,
        massConv,
        name=args.catname,
        nObj=nObj,
        pathInCatalog=args.catpath,
        workDir=args.workdir)

    lpf_enmap, hpf_enmap, boxmask = setup_maps(
        args.lpfpath,
        args.hpfpath,
        args.side_length,
        args.ra_min,
        args.dec_min,
        args.pix_scale
    )

    # run thumbstack
    ts = ThumbStack(
        u,
        galcat,
        hpf_enmap,
        boxmask,
        cmbHit=None,
        cmbMap2=lpf_enmap,
        name=args.tsname,
        save=True,
        nProc=args.nproc,
        filterTypes=args.filtertype,
        estimatorTypes=args.esttype,
        doBootstrap=do_bootstrap,
        equalSignedWeights=args.equalsignedweights,
        workDir=args.workdir,
        runEndToEnd=True,
        test=args.test,
        doStackedMap=do_stacked_map
    )

if __name__ == '__main__':
    main()
