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
    
    parser.add_argument('--flatmap',
                        help='Full Path to the unfiltered flatmap FITS file')
    parser.add_argument('--renorm', default=0.000245, # from Will's fit
                        help='target peak value for the tau Gaussian profile (before applying beam)')
    parser.add_argument('--tau_fwhm', default=5., # from Will's fit
                        help='tau Gaussian profile FWHM in arcmin used to make unfiltered tau map (before applying beam)')    
    parser.add_argument('--outpath',
                        help='path to directory to save output files')
    parser.add_argument('--beam_fwhm', default=1.6,
                        help='FWHM in arcmin of beam to apply to map')
    parser.add_argument('--filter_type', default='theo',
                        help='Choose lpf/hpf functions (`theo` or `will`')
    parser.add_argument('--lpf_loc', default=1000,
                        help='LPF center location for erf function filtering (i.e. `theo` filter_type)')
    parser.add_argument('--hpf_loc', default=1500,
                        help='HPF center location for erf function filtering (i.e. `theo` filter_type)')
    parser.add_argument('--filter_width', default=100.,
                        help='LPF/HPF filter width for erf function filtering (i.e. `theo` filter_type)')
    
    args = parser.parse_args()

    return args

def fbeam(ell, real_fwhm):
    real_sigma = real_fwhm / np.sqrt(8.*np.log(2.))
    return np.exp(-0.5*(ell * real_sigma)**2)

def high_pass(ell, loc=1500, width=100.):
    half_width = width/2
    return erf(np.sqrt(2)/half_width*(ell-loc))/2 + 0.5

def low_pass(ell, loc=1000, width=100.):
    half_width = width/2
    return -1 * erf(np.sqrt(2)/half_width*(ell-loc))/2 + 0.5

def hpf_will(ell):
    return np.sin((ell - 2350) * np.pi / 300.)

def lpf_will(ell):
    return np.cos((ell - 2000) * np.pi / 300.)

def make_basemap(sizeX=11.4, sizeY=11.4, pixel_scale=0.5):
    # sizeX, sizeY: map dimensions in degrees

    # number of pixels for the flat map
    nX = int(sizeX * 60. / pixel_scale)
    nY = int(sizeY * 60. / pixel_scale)

    # basic map object
    baseMap = FlatMap(nX=nX, nY=nY, sizeX=sizeX*np.pi/180., sizeY=sizeY*np.pi/180.)
    
    return baseMap

def renorm_image(image_data, target_peak=0.000245, tau_fwhm=5.0):
    # default profile_fwhm and target_peak from Will's best fit
    sigma = tau_fwhm / np.sqrt(8.*np.log(2))
    gauss_norm = 1. / 2. / np.pi / sigma**2
    norm = target_peak / gauss_norm
    return norm * image_data

def make_flatmap(image_fits, renorm=None, tau_fwhm=5.0):
    print('Creating flatmap from image file:', image_fits)
    flatmap = make_basemap()
    flatmap.read(image_fits)
    print(flatmap.name)
    flatmap.name = '10kgal_gauss_peak1'
    print(flatmap.name)
    if renorm is not None:
        renormed_image = renorm_image(flatmap.data, target_peak=renorm, tau_fwhm=tau_fwhm)
        flatmap.data = renormed_image
        flatmap.name += '_renorm%.2e'%renorm
    return flatmap

def apply_beam(flatmap, fwhm=1.6):
    print('Applying beam with FWHM=%s arcmin.'%str(fwhm))
    # convert fwhm from arcmin to radians
    fwhm_rad = fwhm * np.pi / 180. / 60.
    # v1
    beamfn = partial(fbeam, real_fwhm=fwhm_rad)
    flatmap.dataFourier = flatmap.filterFourierIsotropic(fW=beamfn)
    # v2
    # flatmap.dataFourier *= fbeam(flatmap.l, fwhm_rad)
    flatmap.data = flatmap.inverseFourier()
    flatmap.name += '_beam%s'%str(fwhm)
    
    return flatmap

def apply_filtering(flatmap, filter_type='theo', lpf_loc=1000, hpf_loc=1500, width=100.):
    if filter_type == 'theo':
        lpf = partial(low_pass, loc=lpf_loc, width=width)
        hpf = partial(high_pass, loc=hpf_loc, width=width)
    elif filter_type == 'will':
        lpf = lambda l: lpf_will(l)
        hpf = lambda l: hpf_will(l)
    else:
        raise ValueError('Invalid filter type:', filter_type)

    # make LPF map
    print('Applying low-pass filter.')
    filtFourier_lo = flatmap.filterFourierIsotropic(fW=lpf)
    filtData_lo = flatmap.inverseFourier(filtFourier_lo)
    lpf_map = flatmap.copy()
    lpf_map.data = filtData_lo
    lpf_map.dataFourier = filtFourier_lo
    
    # make HPF map
    print('Applying high-pass filter.')
    filtFourier_high = flatmap.filterFourierIsotropic(fW=hpf)
    filtData_high = flatmap.inverseFourier(filtFourier_high)
    hpf_map = flatmap.copy()
    hpf_map.data = filtData_high
    hpf_map.dataFourier = filtFourier_high
    
    if filter_type == 'theo':
        lpf_map.name += '_lpf%it'%lpf_loc
        hpf_map.name += '_hpf%it'%hpf_loc
    elif filter_type == 'will':
        lpf_map.name += '_lpf2000w'
        hpf_map.name += '_hpf2350w'
    
    return lpf_map, hpf_map

def save_flatmap(flatmap, path=None, save_diagnostics=True):
    if save_diagnostics:
        ell_max = np.max(flatmap.l.flatten())
        ps_path = os.path.join(path, '%s_powspec.png'%flatmap.name)
        print('Saving power spectrum plot:', ps_path)
        flatmap.powerSpectrum(nBins=int(ell_max/200), lRange=[1., ell_max],
            plot=True, save=True, path=ps_path)
        map_path = os.path.join(path, '%s_map.png'%flatmap.name)
        print('Saving map image:', map_path)
        flatmap.plot(save=True, title=flatmap.name, path=map_path)
    fm_path = os.path.join(path, '%s_flatmap.fits'%flatmap.name)
    print('Saving flatmap:', fm_path)
    flatmap.write(fm_path)

def main(argv):
    args = parse_args()
    
    flatmap = make_flatmap(args.flatmap, args.renorm, args.tau_fwhm)
    save_flatmap(flatmap, path=args.outpath)
    
    beamed_map = apply_beam(flatmap, fwhm=args.beam_fwhm)
    save_flatmap(beamed_map, path=args.outpath)
    
    lpf_map, hpf_map = apply_filtering(beamed_map, args.filter_type, args.lpf_loc,
                                       args.hpf_loc, args.filter_width)
    save_flatmap(lpf_map, path=args.outpath)
    save_flatmap(hpf_map, path=args.outpath)


if __name__ == '__main__':
    main(sys.argv)
