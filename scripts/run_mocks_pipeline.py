#!/usr/bin/env python3
import os,sys
import fitsio

sys.path.append('../../ThumbStack')
from catalog import Catalog
from thumbstack import ThumbStack

from make_noise_maps import make_map, make_cmb, gen_map_from_fn, gen_cmb_fg_map, save_flatmap
from filter_tau_map import apply_filtering
from run_thumbstack import initialize, setup_maps

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description='Run thumbstack on a given set of filtered maps.')
    parser.add_argument('--nmocks',
        type=int,
        help='Number of mocks to generate and process.')
    parser.add_argument('--test',
        default=False, action='store_const', const=True,
        help='Test script using first 10 objects in catalog')
    parser.add_argument('--t_large_min',
        default=None,
        type=float,
        help='Minimum abs value temperature for T_large in stack')
    parser.add_argument('--equalsignedweights',
        default=False, action='store_const', const=True,
        help='Whether to equalize number of +/- weights')
    parser.add_argument('--nproc',
        default=1,
        type=int,
        help=('Number of processes to split seeds over for parallel computing.'
              'Each seed runs on one core/thread.'))
    args = parser.parse_args()

    return args

def write_log(args, outfile):
    # Write text file logging what command line args were used
    log_fn =  outfile + '_args.log'
    print('Writing argument log file:', log_fn)
    arg_dict = vars(args)
    with open(log_fn, 'w') as f:
        for arg in arg_dict:
            f.write('%s: %s\n'%(str(arg), str(arg_dict[arg])))

def setup_fixed_dirs():
    FG_DIR = '/home/theo/Documents/research/CMB/patchy_tau_sims/output/fg_maps'

    kappa_fits = os.path.join(
        FG_DIR,
        'kappamap_0.53arcminGauss_correlatedwithTaumap_cmass_10x10_image.fits'
    )

    fg_fits = os.path.join(
        FG_DIR,
        'tSZmap_0.88arcminGauss_multiplicativefact1.3_correlatedwithTaumap_cmass_10x10_image.fits'
    )

    fg_hpf_path = os.path.join(
        FG_DIR,
        'tSZ_sig0.88_x1.3_corr-w-tau',
        'tSZ_sig0.88_x1.3_corr-w-tau_beam1.6_flatmap.fits'
    )

    return kappa_fits, fg_fits, fg_hpf_path

def setup_run_dir(nmocks):
    runname = 'tsz+lens_x_%imocks'%nmocks
    rundir = os.path.join(
        '/home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs',
        runName
    )

    return runname, rundir

def setup_seed_dir(rundir, runname, seed):
    # make map dir with seed in name
    seedname = '%s_seed%s_cmb'%(runname, seed)
    seeddir = os.path.join(rundir, seedname)
    if not os.path.exists(seeddir):
        os.makedirs(seeddir)

    return seedname, seeddir

def gen_unfiltered_maps(seed, seedname, seeddir, kappa_fits, fg_fits):
    """generates and saves CMB GRF map, same GRF lensed by the `kappa_fits` kappa map, and the unlensed GRF plus the `fg_fits` FG map"""

    # make base flatmap for all map generation
    basemap = make_map()

    # make unlensed CMB GRF map
    cmb = make_cmb()
    Ctot = lambda l: cmb.funlensedTT(l)
    cmbname =  seedname+'_cmb':
    cmbmap = gen_map_from_fn(Ctot, cmb, basemap, cmbname, seed=seed)

    # make lensed CMB GRF
    lensname =  cmbname+'+lens':
    kappaData = fitsio.read(kappa_fits)
    kappaFourier = cmbmap.fourier(kappaData)
    cmb_lens_map = cmbmap.doLensing(kappaFourier)
    cmb_lens_map.name = lensname

    # make unlensed CMB GRF + tsz map
    fgname =  cmbname+'+tsz'
    cmb_fg_map = gen_cmb_fg_map(tsz_fits, fgname, cmb_fm=cmbmap)

    # save in seed dir
    if test:
        # includes saving diagnostics: png image and power spectrum plot
        save_flatmap(cmbmap, cmb, seeddir, cmbname)
        save_flatmap(cmb_lens_map, cmb, seeddir, lensname)
        save_flatmap(cmb_fg_map, cmb, seeddir, fgname)
    else:
        cmbmap.write(seeddir)
        cmb_lens_map.write(seeddir)
        cmb_fg_map.write(seeddir)

    return cmbmap, cmb_lens_map, cmb_fg_map

def gen_filtered_maps(cmbmap, cmb_lens_map, cmb_fg_map):
    """applies LPF, HPF filtering, saves these new maps in the seeddir.
    @return the list of directory paths for these maps"""

    map_paths = []
    for cmap in [cmbmap, cmb_lens_map, cmb_fg_map]:
        lpfmap, hpfmap = apply_filtering(cmap, filter_type='will')
        lpfpath = save_flatmap(lpfmap, path=seeddir, save_diagnostics=test)
        map_paths.append(lpfpath)
        hpfpath = save_flatmap(hpfmap, path=seeddir, save_diagnostics=test)
        map_paths.append(hpfpath)

    return map_paths

# def filter_cmb_maps(cmbmap):
#     lpfmap, hpfmap = apply_filtering(cmap, filter_type='will')
#         lpfpath = save_flatmap(lpfmap, path=seeddir, save_diagnostics=test)
#         map_paths.append(lpfpath)
#         hpfpath = save_flatmap(hpfmap, path=seeddir, save_diagnostics=test)
#         map_paths.append(hpfpath)
#
#     return lpfpath, hpfpath

def setup_for_ts():
    u, massConv = initialize()

    galcat = Catalog(
        u,
        massConv,
        name='cmass_m_10x10_v2',
        nObj=nObj,
        workDir='/home/theo/Documents/research/CMB/patchy_tau_sims'
    )

    return u, galcat

# def run_ts_lens(lpfpath, hpfpath, tsname):
#
#     lpf_enmap, hpf_enmap, boxmask = setup_maps(lpfpath, hpfpath)
#
#     # run thumbstack, no bootstrap,
#     # outdirs: seeddir/output/thumbstack, seeddir/figures/thumbstack
#     _ = ThumbStack(
#         u,
#         galcat,
#         hpf_enmap,
#         boxmask,
#         cmbHit=None,
#         cmbMap2=lpf_enmap,
#         name='%imocks_seed%s'%(nMocks, seed) + '_cmb-lpf_cmblens-hpf',
#         save=True,
#         nProc=args.nproc,
#         filterTypes='tauring',
#         estimatorTypes=['tau_ti_uniformweight', 'tau_sgn_uniformweight'],
#         doBootstrap=False,
#         tLargeMin=args.t_large_min,
#         equalSignedWeights=args.equalsignedweights,
#         workDir=seeddir,
#         runEndToEnd=True,
#         test=args.test,
#         doStackedMap=False
#     )
#
# def run_ts_tsz():
#     map_pairs = [(map_paths[0], map_paths[1]), # CMB LPF    x CMB HPF
#                  (map_paths[0], map_paths[1]), # CMB LPF    x CMB+lens HPF
#                  (map_paths[0], fg_hpf_path),  # CMB LPF    x fg HPF
#                  (map_paths[2], fg_hpf_path)]  # CMB+fg LPF x fg HPF
#
#     # suffixes for thumbstack naming
#     map_suffix = ['_cmb-lpf_cmb-hpf',
#                   '_cmb-lpf_%s-hpf'%fg_type,
#                   '_cmb+%s-lpf_%s-hpf'%(fg_type, fg_type)]
#
#     for i, pair in enumerate(map_pairs):
#         lpf_enmap, hpf_enmap, boxmask = setup_maps(pair[0], pair[1])
#
#         # run thumbstack, no bootstrap,
#         # outdirs: seeddir/output/thumbstack, seeddir/figures/thumbstack
#         _ = ThumbStack(
#             u,
#             galcat,
#             hpf_enmap,
#             boxmask,
#             cmbHit=None,
#             cmbMap2=lpf_enmap,
#             name='%imocks_seed%s'%(nMocks, seed) + map_suffix[i],
#             save=True,
#             nProc=args.nproc,
#             filterTypes='tauring',
#             estimatorTypes=['tau_ti_uniformweight', 'tau_sgn_uniformweight'],
#             doBootstrap=False,
#             tLargeMin=args.t_large_min,
#             equalSignedWeights=args.equalsignedweights,
#             workDir=seeddir,
#             runEndToEnd=True,
#             test=args.test,
#             doStackedMap=False
#         )

def end2end(args, kappa_fits, tsz_fits, tsz_hpf_path, runname, rundir, seed):
    """runs pipeline for a single seed"""

    seedname, seeddir = setup_seed_dir(rundir, runname, seed)

    cmbmap, cmb_lens_map, cmb_tsz_map = gen_unfiltered_maps(seed, seedname, seeddir, kappa_fits, tsz_fits)

    map_paths = gen_filtered_maps(cmbmap, cmb_lens_map, cmb_tsz_map)

    # setup for 4 thumbstack runs
    map_pairs = [(map_paths[0], map_paths[1]), # CMB LPF    x CMB HPF
                 (map_paths[0], map_paths[1]), # CMB LPF    x CMB+lens HPF
                 (map_paths[0], tsz_hpf_path),  # CMB LPF    x fg HPF
                 (map_paths[2], tsz_hpf_path)]  # CMB+fg LPF x fg HPF

    # suffixes for thumbstack naming
    map_suffix = ['_cmb-lpf_cmb-hpf',
                  '_cmb-lpf_cmb+lens-hpf',
                  '_cmb-lpf_tsz-hpf',
                  '_cmb+tsz-lpf_tsz-hpf']

    u, galcat = setup_for_ts()

    for i, pair, suffix in enumerate(zip(map_pairs, map_suffix)):
        tsname = seedname + suffix
        lpf_enmap, hpf_enmap, boxmask = setup_maps(pair[0], pair[1])

        # run thumbstack, no bootstrap,
        # outdirs: seeddir/output/thumbstack, seeddir/figures/thumbstack
        _ = ThumbStack(
            u,
            galcat,
            hpf_enmap,
            boxmask,
            cmbHit=None,
            cmbMap2=lpf_enmap,
            name=tsname,
            save=True,
            nProc=1, # we're parallelizing over seeds
                     # so each seed needs to use only one thread
            filterTypes='tauring',
            estimatorTypes=['tau_ti_uniformweight', 'tau_sgn_uniformweight'],
            doBootstrap=False,
            tLargeMin=args.t_large_min,
            equalSignedWeights=args.equalsignedweights,
            workDir=seeddir,
            runEndToEnd=True,
            test=args.test,
            doStackedMap=False
        )

def main():
    args = parse_args()
    # same for all seeds
    kappa_fits, tsz_fits, tsz_hpf_path = setup_fixed_dirs()
    runname, rundir = setup_run_dir(nmocks=args.nmocks)
    write_log(args, os.path.join(rundir, runname))

    end2end(args, kappa_fits, tsz_fits, tsz_hpf_path, runname, rundir, seed) #joblib this
