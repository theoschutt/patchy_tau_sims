"""filter_map.py: Applies beam and low- and high-pass filtering to an
image FITS file (e.g. from generateMockMaps).
If used as main, expects as CLI args fits file name (including path),
output file name prefix, and path to save directory.
"""
import os, sys
sys.path.append('../../ThumbStack')
import numpy as np
import fitsio
from scipy.special import erf
from functools import partial

from flat_map import FlatMap

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Apply beam and LPF/HPF to a flatmap or image file.')
    parser.add_argument('--image_fits',
                        default=None,
                        help='Full Path to the unfiltered image FITS file')
    parser.add_argument('--flatmap_fits',
                        default=None,
                        help='Full Path to the unfiltered flatmap FITS file')
    parser.add_argument('--flatmap_name',
                        help='Name stem for the flatmap to be made from the image FITS')
    parser.add_argument('--renorm',
                        default=None, #0.000245, # from Will's fit
                        type=float,
                        help='target peak value for the tau Gaussian profile (before applying beam)')
    parser.add_argument('--sub_const',
                        default=None, type=float,
                        help='Constant to subtract from the entire map (before applying beam). For pure tau maps.')
    parser.add_argument('--halfconst_map',
                        default=None, type=float,
                        help='Constant to use when making a half-n-half constant map (i.e. 50%% +const, 50%% -const map) for systematics studies.')
    parser.add_argument('--tau_fwhm',
                        default=5., # from Will's fit
                        type=float,
                        help='tau Gaussian profile FWHM in arcmin that used to make unfiltered tau map (before applying beam). Needed when renorming. [default: 5.]')
    parser.add_argument('--outpath',
                        default=None,
                        help='path to directory to save output files')
    parser.add_argument('--save_image', default=True,
                        action='store_const', const=True,
                        help='in addition to flatmaps, save image fits files')
    parser.add_argument('--beam_fwhm',
                        default=None,
                        type=float,
                        help='FWHM in arcmin of beam to apply to map')
    parser.add_argument('--sensitivity',
                        default=None,
                        type=float,
                        help='Noise in uK*arcmin to apply to map')
    parser.add_argument('--filter_type', default='will',
            help='Choose lpf/hpf functions (`theo` or `will`) [default: will]')
    parser.add_argument('--lpf_loc',
                        default=1000,
                        type=int,
                        help='LPF center location for erf function filtering (i.e. `theo` filter_type)')
    parser.add_argument('--hpf_loc',
                        default=1500,
                        type=int,
                        help='HPF center location for erf function filtering (i.e. `theo` filter_type)')
    parser.add_argument('--filter_width',
                        default=100.,
                        type=float,
                        help='LPF/HPF filter width for erf function filtering (i.e. `theo` filter_type)')

    args = parser.parse_args()

    return args

def fbeam(ell, real_fwhm):
    """real_fwhm [rad]"""
    real_sigma = real_fwhm / np.sqrt(8.*np.log(2.))
    return np.exp(-0.5*(ell * real_sigma)**2)

def detnoise(ell, sens_amin, real_fwhm):
    sens_rad = sens_amin * (np.pi/180.)/60.
    return sens_rad**2 / fbeam(ell, real_fwhm)**2

def high_pass(ell, loc=1500, width=100.):
    half_width = width/2
    return erf(np.sqrt(2)/half_width*(ell-loc))/2 + 0.5

def low_pass(ell, loc=1000, width=100.):
    half_width = width/2
    return -1 * erf(np.sqrt(2)/half_width*(ell-loc))/2 + 0.5

def hpf_will(ell):
    # do piecewise filtering
    low_ell = ell[ell<2350]
    mid_ell = ell[np.logical_and(2350<ell, ell<2500)]
    hi_ell = ell[ell>2500]
    low_filt = np.zeros_like(low_ell)
    mid_filt = np.sin((mid_ell-2350) * np.pi/300.)
    hi_filt = np.ones_like(hi_ell)
    tot_filt = np.concatenate([low_filt, mid_filt, hi_filt])
    return tot_filt

def lpf_will(ell):
    # do piecewise filtering
    low_ell = ell[ell<2000]
    mid_ell = ell[np.logical_and(2000<ell, ell<2150)]
    hi_ell = ell[ell>2150]
    low_filt = np.ones_like(low_ell)
    mid_filt = np.cos((mid_ell-2000) * np.pi/300.)
    hi_filt = np.zeros_like(hi_ell)
    tot_filt = np.concatenate([low_filt, mid_filt, hi_filt])
    return tot_filt

def make_basemap(sizeX=10., sizeY=10., pixel_scale=0.5):
    # sizeX, sizeY: map dimensions in degrees

    # number of pixels for the flat map
    nX = int(sizeX * 60. / pixel_scale)
    nY = int(sizeY * 60. / pixel_scale)

    # basic map object
    baseMap = FlatMap(nX=nX, nY=nY, sizeX=sizeX*np.pi/180., sizeY=sizeY*np.pi/180.)

    return baseMap

def renorm_map(flatmap, target_peak=0.000245, tau_fwhm=5.0):
    # default profile_fwhm and target_peak from Will's best fit
    renormmap = flatmap.copy()
    sigma = tau_fwhm / np.sqrt(8.*np.log(2))
    gauss_norm = 1. / 2. / np.pi / sigma**2
    norm = target_peak / gauss_norm
    renormmap.data *= norm
    renormmap.dataFourier = renormmap.fourier()
    renormmap.name += '_renorm%.2e'%target_peak
    return renormmap

def subtract_const(flatmap, const):
    """Create constant map with tau profiles added.
    For tau signal map only."""
    constmap = flatmap.copy()
    constmap.data = - const * (np.ones_like(flatmap.data) - flatmap.data)
    constmap.dataFourier = constmap.fourier(constmap.data)
    constmap.name += '_c-%s'%str(const)
    return constmap

def make_const_map(flatmap, const):
    # need constant map for T_large for pure tau measurement
    # normal LPF removes the mean over the map, which is usually
    # zero anyway, but for the signal-only maps, is nonzero.
    constmap = flatmap.copy()
    constmap.data = - const * np.ones_like(flatmap.data)
    constmap.dataFourier = constmap.fourier(constmap.data)
    constmap.name += '_CONST%s'%str(const)
    return constmap

def make_halfconst_map(flatmap, const=1.):
    #constmap_v = flatmap.copy()
    #constmap_h = flatmap.copy()
    constmap_4 = flatmap.copy()
    x, y = flatmap.data.shape
    print(x,y)

    # constmap_v.data[:x//2,:] = const * np.ones_like(flatmap.data[:x//2,:])
    # constmap_v.data[x//2:,:] = - const * np.ones_like(flatmap.data[x//2:,:])
    # constmap_v.dataFourier = constmap_v.fourier(constmap_v.data)
    # constmap_v.name += '_HALF-N-HALF_vert_CONST%s'%str(const)
    # constmap_h.data[:,:x//2] = const * np.ones_like(flatmap.data[:,:x//2])
    # constmap_h.data[:,x//2:] = - const * np.ones_like(flatmap.data[:,x//2:])
    # constmap_h.dataFourier = constmap_h.fourier(constmap_h.data)
    # constmap_h.name += '_HALF-N-HALF_horz_CONST%s'%str(const)

    constmap_4.data[:x//2,:x//2] = const * np.ones_like(flatmap.data[:x//2,:x//2])
    constmap_4.data[x//2:,x//2:] = const * np.ones_like(flatmap.data[x//2:,x//2:])
    constmap_4.data[:x//2,x//2:] = -const * np.ones_like(flatmap.data[:x//2,x//2:])
    constmap_4.data[x//2:,:x//2] = -const * np.ones_like(flatmap.data[x//2:,:x//2])
    constmap_4.dataFourier = constmap_4.fourier(constmap_4.data)
    constmap_4.name += '_HALF-N-HALF_quad_CONST%s'%str(const)

    # make screened quadrant map
    corr_map = flatmap.copy()
    corr_map.data = - constmap_4.data * corr_map.data
    corr_map.dataFourier = corr_map.fourier()
    corr_map.name += '_HALF-N-HALF_quad_CONST%s_weighted'%str(const)
    #return constmap_v, constmap_h, constmap_4
    return constmap_4, corr_map

def load_flatmap(flatmap_fits):
    print('Loading flatmap file:', flatmap_fits)
    flatmap = make_basemap()
    flatmap.read(flatmap_fits)
    return flatmap

def make_flatmap(image_fits, name='test'):
    print('Creating flatmap with name %s from image file:'%name, image_fits)
    flatmap = make_basemap()
    image_data = fitsio.read(image_fits)
    flatmap.data = image_data
    flatmap.dataFourier = flatmap.fourier(image_data)
    flatmap.name = name
    return flatmap

def apply_beam(flatmap, fwhm):
    print('Applying beam with FWHM=%s arcmin.'%str(fwhm))
    # convert fwhm from arcmin to radians
    beamedmap = flatmap.copy()
    fwhm_rad = fwhm * np.pi / 180. / 60.
    beamfn = partial(fbeam, real_fwhm=fwhm_rad)
    beamedmap.dataFourier = beamedmap.filterFourierIsotropic(fW=beamfn)
    beamedmap.data = beamedmap.inverseFourier()
    beamedmap.name += '_beam%s'%str(fwhm)

    return beamedmap

def add_noise(flatmap, sens_amin, real_fwhm):
    print('Adding detector noise with sensitivity: %.2f [muK*arcmin]'%sens_amin)
    noisemap = flatmap.copy()
    noisefn = partial(detnoise, sens_amin=sens_amin, real_fwhm=real_fwhm)
    noisemap.dataFourier = noisemap.filterFourierIsotropic(fW=noisefn)
    noisemap.data = noisemap.inverseFourier()
    noisemap.name += '_noise%.2f'%sens_amin

    return noisemap

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
        lpf_map.name += '_lpf2075w'
        hpf_map.name += '_hpf2425w'

    return lpf_map, hpf_map

def save_flatmap(flatmap, path=None, save_image=True, save_diagnostics=True):
    print('--------------\nSaving to directory:', path)
    if save_diagnostics:
        ell_max = np.max(flatmap.l.flatten())
        ps_path = os.path.join(path, '%s_powspec.png'%flatmap.name)
        flatmap.powerSpectrum(nBins=int(ell_max/200), lRange=[1., ell_max],
            plot=True, save=True, path=ps_path)

        map_path = os.path.join(path, '%s_map.png'%flatmap.name)
        flatmap.plot(save=True, title=flatmap.name, path=map_path)

    fm_path = os.path.join(path, '%s_flatmap.fits'%flatmap.name)
    flatmap.write(fm_path)

    if save_image:
        im_path = os.path.join(path, '%s_image.fits'%flatmap.name)
        print('Saving image fits:', os.path.basename(im_path))
        fitsio.write(im_path, flatmap.data)

    return fm_path

def main():
    args = parse_args()

    if args.outpath is None:
        if args.flatmap_fits is not None:
            args.outpath = os.path.split(args.flatmap_fits)[0]
        elif args.image_fits is not None:
            args.outpath = os.path.split(args.image_fits)[0]

    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)

    if args.flatmap_fits is not None:
        flatmap = load_flatmap(args.flatmap_fits)
        save_flatmap(flatmap, path=args.outpath, save_image=args.save_image)
    elif args.image_fits is not None:
        flatmap = make_flatmap(args.image_fits, name=args.flatmap_name)
        save_flatmap(flatmap, path=args.outpath, save_image=args.save_image)
    else:
        raise ValueError('Must specify an input image or flatmap FITS file.')

    # Write text file logging what command line args were used
    if args.flatmap_name is None:
        args.flatmap_name = flatmap.name
    log_fn = os.path.join(args.outpath, 'args_%s.log'%args.flatmap_name)
    print('Writing argument log file:', log_fn)
    arg_dict = vars(args)
    with open(log_fn, 'w') as f:
        for arg in arg_dict:
            f.write('%s: %s\n'%(str(arg), str(arg_dict[arg])))

    if args.renorm is not None:
        flatmap = renorm_map(flatmap, target_peak=args.renorm,
            tau_fwhm=args.tau_fwhm)
        save_flatmap(flatmap, path=args.outpath, save_image=args.save_image)

    if args.halfconst_map is not None:
        quad_map, flatmap = make_halfconst_map(flatmap, const=args.halfconst_map)
        # hc_map_v, hc_map_h, hc_map_4 = make_halfconst_map(flatmap, const=args.halfconst_map)
        #save_flatmap(hc_map_v, path=args.outpath, save_image=args.save_image)
        #save_flatmap(hc_map_h, path=args.outpath, save_image=args.save_image)
        save_flatmap(quad_map, path=args.outpath, save_image=args.save_image)
        save_flatmap(flatmap, path=args.outpath, save_image=args.save_image)

    if args.sub_const is not None:
        flatmap = subtract_const(flatmap, args.sub_const)
        save_flatmap(flatmap, path=args.outpath, save_image=args.save_image)
        constmap = make_const_map(flatmap, args.sub_const)
        save_flatmap(constmap, path=args.outpath, save_image=args.save_image)

    if args.beam_fwhm is not None:
        flatmap = apply_beam(flatmap, fwhm=args.beam_fwhm)
        save_flatmap(flatmap, path=args.outpath, save_image=args.save_image)

    if args.sensitivity is not None:
        flatmap = add_noise(flatmap, sens_amin=args.sensitivity, real_fwhm=args.beam_fwhm)
        save_flatmap(flatmap, path=args.outpath, save_image=args.save_image)

    lpf_map, hpf_map = apply_filtering(flatmap, args.filter_type, args.lpf_loc,
                                       args.hpf_loc, args.filter_width)
    save_flatmap(lpf_map, path=args.outpath, save_image=args.save_image)
    save_flatmap(hpf_map, path=args.outpath, save_image=args.save_image)

if __name__ == '__main__':
    main()
