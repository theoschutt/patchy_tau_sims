import os, sys
sys.path.append('../../ThumbStack')
sys.path.append('../../LensQuEst')
import numpy as np
import fitsio
#import flat_map
from flat_map import FlatMap
from pn_2d import interp1d
from cmb import StageIVCMB
#from filter_map import apply_filtering_and_save

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Create map from noise curve or expt parameters.')
    parser.add_argument(
        '--noise_curve_file',
        default=None,
        help='Full Path to the input noise curve data file'
    )
    parser.add_argument(
        '--experiment',
        default='s4',
        help='choice of experiment. Sets parameters such as noise, beam, filename tags'
    )
    parser.add_argument(
        '--name',
        default='test',
        help='name stem for the flatmap to be made'
    )
    parser.add_argument(
        '--outpath',
        default='../output/noise_maps/',
        help='output path where head directory (named --name) will be created'
    )
    parser.add_argument(
        '--map_side_length',
        default=10., type=float,
        help='side length of map in degrees. Will make a square map.'
    )
    parser.add_argument(
        '--pixel_scale',
        default=0.5,
        help='pixel scale of map in arcmin'
    )
    parser.add_argument(
        '--ell_min',
        default=40.,
        help='minimum ell to be used in creating the power spectrum to generate the map'
    )
    parser.add_argument(
        '--ell_max',
        default=10000.,
        help='minimum ell to be used in creating the power spectrum to generate the map'
    )
    parser.add_argument(
        '--save_image',
        default=True,
        action='store_const',
        const=True,
        help='in addition to flatmaps, save image fits files'
    )
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

def make_map(sizeX=10., sizeY=10., pixel_scale=0.5):
    # map dimensions in degrees

    # number of pixels for the flat map, let's do 0.5' pixels
    nX = int(sizeX * 60. / pixel_scale)
    nY = int(sizeY * 60. / pixel_scale)

    # basic map object
    baseMap = FlatMap(nX=nX, nY=nY, sizeX=sizeX*np.pi/180., sizeY=sizeY*np.pi/180.)

    return baseMap

def make_cmb(lMin=30., lMax=10000., nBins=150, beam=1., noise=1.):
    # Adjust the lMin and lMax to the assumptions of the analysis.
    # multipoles to include in the lensing reconstruction
    # lMin = 30.; lMax = 2.5e4

    # ell bins for power spectra
    # nBins = 150  # number of bins

    # CMB S4 specs:
    #     beam: 1'
    #     noise: 1 uK*arcmin
    # SO specs:
    #     beam: 1.3'
    #     noise: 10 uK*arcmin
    # AdvACT specs:
    #     beam: 1.4' update: 1.6'
    #     noise: 15 uK*arcmin
    cmb = StageIVCMB(beam=beam, noise=noise, lMin=lMin, lMaxT=lMax, lMaxP=lMax, atm=False)

    return cmb

def gen_screened_map(tau_fits, name):
    basemap = make_map(10., 10., 0.5)
    cmb = make_cmb()
    dat = fitsio.read(tau_fits)
    forCtot = lambda l: cmb.funlensedTT(l)# * cmb.fbeam(l)**2
    cmbmap = gen_map_from_fn(forCtot, cmb, basemap, name)
    cmbmap.data *= (1-dat)
    cmbmap.dataFourier = cmbmap.fourier(cmbmap.data)

    return cmbmap

def gen_cmb_fg_map(fg_fits, cmb_fits=None, name='test'):
    cmbmap = make_map(10., 10., 0.5)
    dat = fitsio.read(fg_fits)

    if cmb_fits is None:
        cmb = make_cmb()
        forCtot = lambda l: cmb.funlensedTT(l)
        cmbmap = gen_map_from_fn(forCtot, cmb, basemap, name)
    else:
        cmbdat = fitsio.read(cmb_fits)
    cmbmap.data = cmbdat
    cmbmap.data += dat
    cmbmap.dataFourier = cmbmap.fourier(cmbmap.data)
    cmbmap.name = name

    return cmbmap

def gen_map_from_fn(powspec_fn, cmb, flat_map, name, lMin=30., lMax=10000., seed=42):
    print('Generating map')
    # reinterpolate: gain factor 10 in speed
    L = np.logspace(np.log10(lMin/2.), np.log10(2.*lMax), 1001, 10.)
    F = np.array(list(map(powspec_fn, L)))
    cmb.fCtotal = interp1d(L, F, kind='linear', bounds_error=False, fill_value=0.)

    # make map data
    cmbtotFourier = flat_map.genGRF(cmb.fCtotal, seed=seed, test=False)
    cmbtot = flat_map.inverseFourier(cmbtotFourier)

    # make new flat_map with map data
    this_cmb_map = flat_map.copy()
    this_cmb_map.name = name
    this_cmb_map.data = cmbtot
    this_cmb_map.dataFourier = cmbtotFourier

    return this_cmb_map

def gen_map_from_curve(ell, f, cmb, flat_map, name, seed=42):
    print('Generating map')
    # reinterpolate: gain factor 10 in speed
    # first interpolate lensed CMB
    forCtotal = lambda l: cmb.flensedTT(l)
    F = np.array(list(map(forCtotal, ell)))
    # then add on interpolated noise curve
    F += f
    # add the beam
    F = F * cmb.fbeam(ell)**2
    cmb.fCtotal = interp1d(ell, F, kind='linear', bounds_error=False, fill_value=0.)

    # make map data
    cmbtotFourier = flat_map.genGRF(cmb.fCtotal, seed=seed, test=False)
    cmbtot = flat_map.inverseFourier(cmbtotFourier)

    # make new flat_map with map data
    this_cmb_map = flat_map.copy()
    this_cmb_map.name = name
    this_cmb_map.data = cmbtot
    this_cmb_map.dataFourier = cmbtotFourier

    return this_cmb_map

def save_map(flatmap, cmb, path, name):
    this_cmb_path = os.path.join(path, name)
    if not os.path.exists(this_cmb_path):
        os.makedirs(this_cmb_path)

    # plot and save diagnostic plots
    print("plot and save CMB map")
    map_path = os.path.join(this_cmb_path, '%s_map.png'%name)
    flatmap.plot(save=True, path=map_path)
    print("plot and save the power spectrum")
    ps_path = os.path.join(this_cmb_path, '%s_powspec.png'%name)
    lCen, Cl, sCl = flatmap.powerSpectrum(theory=[cmb.funlensedTT,
        cmb.flensedTT, cmb.fdetectorNoise], nBins=150, plot=True, 
        save=True, path=ps_path)

    # save CMB map
    fm_path = os.path.join(this_cmb_path, '%s_flatmap.fits'%name)
    print('Saving flatmap:', fm_path)
    flatmap.write(fm_path)

    # and then make and save filtered versions of the map
    # apply_filtering_and_save(fm_path, lpf_loc=1000, hpf_loc=1500,
    #     half_width=50., path=this_cmb_path, save_diagnostics=True)

def gen_and_save_maps_from_file(file, name, path='output/cmb_maps/10x10_noise_maps'):
    cmb = make_cmb()
    fm = make_map()
    ncurve = np.genfromtxt(file)
    ell = ncurve[:,0]
    f = ncurve[:,1]
    nmap = gen_map_from_curve(ell, f, cmb, fm, name)
    save_map(nmap, cmb, path, name)

def gen_and_save_maps_from_npz(npz, name='advACT_nilc', path='output/cmb_maps/10x10_noise_maps'):
    cmb = make_cmb(beam=1.6)
    fm = make_map()
    dat = np.load(npz)
    ell = dat['ells']
    f = dat['cl_tt']
    nmap = gen_map_from_curve(ell, f, cmb, fm, name)
    save_map(nmap, cmb, path, name)

def gen_and_save_spt_maps(file, name='spt_2023ilc', path='output/cmb_maps/10x10_noise_maps/'):
    cmb = make_cmb(beam=1.2)
    fm = make_map()
    sptdict = np.load(file, allow_pickle = 1, encoding = 'latin1').item()
    ell = sptdict['el']
    f = sptdict['cl_residual']['TT']
    nmap = gen_map_from_curve(ell, f, cmb, fm, name)
    save_map(nmap, cmb, path, name)

def make_s4_noise_maps(path='output/cmb_maps/10x10_noise_maps/s4_l30-10k_v1'):

    # specs and file paths for 10kgal-related maps
    # path = 'output/cmb_maps/10kgal_noise_maps/s4_l30-25k_v2'
    if not os.path.exists(path):
        os.makedirs(path)
    cmb = make_cmb()
    baseMap = make_map()

    # Make series of CMB noise maps
    # 1. CMB, no lensing, atmosphere, detnoise, or foregrounds
    # 2. lensed CMB
    # 3. lensed CMB + foregrounds
    # 4. lensed CMB + foregrounds + S4 noise
    # 5. lensed CMB + foregrounds + S4 noise + atmosphere

    # make power spectra for maps
    forCtotal_1 = lambda l: cmb.funlensedTT(l) * cmb.fbeam(l)**2
    forCtotal_2 = lambda l: cmb.flensedTT(l) * cmb.fbeam(l)**2
    forCtotal_3 = lambda l: (cmb.ftotal(l) - cmb.fdetectorNoise(l)) * cmb.fbeam(l)**2
    forCtotal_4 = lambda l: cmb.ftotal(l) * cmb.fbeam(l)**2
    forCtotal_5 = lambda l: (cmb.ftotal(l) + cmb.fatmosphericNoiseTT(l)) * cmb.fbeam(l)**2

    cmb_names = ['cmb_beam_nolens', 'cmb_beam_lens', 'cmb_beam_lens+fg',
                 'cmb_beam_lens+fg+det', 'cmb_beam_lens+fg+det+atm']
    ps_list = [forCtotal_1, forCtotal_2, forCtotal_3, forCtotal_4, forCtotal_5]

    lMin = 30.
    lMax = 10000.

    # interpolate power spectrum and make map
    for forCtotal, name in zip(ps_list, cmb_names):
        print('Generating map:', name)
        # reinterpolate: gain factor 10 in speed
        L = np.logspace(np.log10(lMin/2.), np.log10(2.*lMax), 1001, 10.)
        F = np.array(list(map(forCtotal, L)))
        cmb.fCtotal = interp1d(L, F, kind='linear', bounds_error=False, fill_value=0.)

        # make map data
        cmbtotFourier = baseMap.genGRF(cmb.fCtotal, seed=42, test=False)
        cmbtot = baseMap.inverseFourier(cmbtotFourier)

        # make new flat_map with map data
        this_cmb = baseMap.copy()
        this_cmb.name = name
        this_cmb.data = cmbtot
        this_cmb.dataFourier = cmbtotFourier

        this_cmb_path = os.path.join(path, name)
        if not os.path.exists(this_cmb_path):
            os.makedirs(this_cmb_path)

        # plot and save diagnostic plots
        print("plot and save CMB map")
        map_path = os.path.join(this_cmb_path, '%s_map.png'%name)
        this_cmb.plot(save=True, path=map_path)
        print("plot and save the power spectrum")
        ps_path = os.path.join(this_cmb_path, '%s_powspec.png'%name)
        lCen, Cl, sCl = this_cmb.powerSpectrum(theory=[cmb.funlensedTT,
            cmb.flensedTT, cmb.fdetectorNoise], nBins=150, plot=True,
            save=True, path=ps_path)

        # save CMB map
        fm_path = os.path.join(this_cmb_path, '%s_flatmap.fits'%name)
        print('Saving flatmap:', fm_path)
        this_cmb.write(fm_path)

        # and then make and save filtered versions of the map
        # apply_filtering_and_save(fm_path, lpf_loc=1000, hpf_loc=1500,
        #     half_width=50., path=this_cmb_path, save_diagnostics=True)

def main(argv):
    args = parse_args()

