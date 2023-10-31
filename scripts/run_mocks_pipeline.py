#!/usr/bin/env python3
import os,sys
import pickle
import sharedmem
import fitsio

sys.path.append('../../ThumbStack')
from catalog import Catalog
from thumbstack import ThumbStack

from make_noise_maps import make_map, make_cmb, gen_map_from_fn, gen_cmb_fg_map
from filter_tau_map import apply_beam, apply_filtering, save_flatmap
from run_thumbstack import initialize, setup_maps

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description='Make, filter, and run thumbstack on mock CMB maps.')
    parser.add_argument('--nmocks',
        type=int,
        help='Number of mocks to generate and process.')
    parser.add_argument('--njobs',
        type=int,
        help='Number of jobs to spawn in multiprocessing.')
    parser.add_argument('--start_seed',
        type=int,
        help='Value to use as the first seed in sequence.')
    parser.add_argument('--version',
        type=str,
        default='01',
        help='Version tag for directory naming.')

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
    fg_dir = '/home/theo/Documents/research/CMB/patchy_tau_sims/output/fg_maps'

    kappa_fits = os.path.join(
        fg_dir,
        'kappamap_0.53arcminGauss_correlatedwithTaumap_cmass_10x10_image.fits'
    )

    # fg_fits = os.path.join(
    #     fg_dir,
    #     'tSZmap_0.88arcminGauss_multiplicativefact1.3_correlatedwithTaumap_cmass_10x10_image.fits'
    # )
    fg_fits = os.path.join(
        fg_dir,
        'tSZmap_2.46arcminGauss_correlatedwithTaumap_cmass_10x10_image.fits'
    )


    # fg_hpf_path = os.path.join(
    #     fg_dir,
    #     'tSZ_sig0.88_x1.3_corr-w-tau',
    #     'tSZ_sig0.88_x1.3_corr-w-tau_beam1.6_hpf2425w_flatmap.fits'
    # )


    fg_hpf_path = os.path.join(
        fg_dir,
        'tSZ_sig2.46_corr-w-tau',
        'tSZ_sig2.46_corr-w-tau_beam1.6_hpf2425w_lmax5950_flatmap.fits'
    )

    # advact_spec_dir = '/home/theo/Documents/research/CMB/patchy_tau_sims/data/AdvACT_NILC_cls_fullRes_TT'

    return kappa_fits, fg_fits, fg_hpf_path

def setup_run_dir(args):
    runname = 'lo-tsz_final_%imocks'%args.nmocks

    if args.equalsignedweights:
        runname += '_eqsgn'
    if args.t_large_min is not None:
        runname += '_tmin%f'%args.t_large_min
    if args.test:
        runname += '_test'
    runname += '_'+args.version

    rundir = os.path.join(
        '/home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs',
        runname
    )

    if not os.path.exists(rundir):
        os.makedirs(rundir)

    return runname, rundir

def setup_seed_dir(rundir, runname, seed):
    # make map dir with seed in name
    seedname = '%s_seed%s'%(runname, seed)
    seeddir = os.path.join(rundir, seedname)
    if not os.path.exists(seeddir):
        os.makedirs(seeddir)

    return seedname, seeddir

def gen_unfiltered_maps(seed, seedname, seeddir, kappa_fits, fg_fits):
    """generates and saves CMB GRF map, same GRF lensed by the `kappa_fits` kappa map,
    and the unlensed GRF plus the `fg_fits` FG map"""

    # make base flatmap for all map generation
    basemap = make_map()
    maplist = []
    # make unlensed CMB GRF map
    cmb = make_cmb()
    Ctot = lambda l: cmb.funlensedTT(l)
    cmbname = seedname + '_cmb'
    cmbmap = gen_map_from_fn(Ctot, cmb, basemap, cmbname, seed=seed)
    # maplist.append(cmbmap)

    # make lensed CMB GRF
    # lensname = cmbname + '+lens'
    # kappaData = fitsio.read(kappa_fits)
    # kappaFourier = cmbmap.fourier(kappaData)
    # cmb_lens_map = basemap.copy()
    # cmb_lens_map.data = cmbmap.doLensing(kappaFourier=kappaFourier)
    # cmb_lens_map.dataFourier = cmb_lens_map.fourier()
    # cmb_lens_map.name = lensname
    # maplist.append(cmb_lens_map)

    # make unlensed CMB GRF + tsz map
    fgname =  cmbname+'+tsz'
    cmb_fg_map = gen_cmb_fg_map(fg_fits, fgname, cmb_fm=cmbmap)
    maplist.append(cmb_fg_map)

    # make AdvACT-like CMB (with noise) + tSZ map
    # advact_name = seedname + '_advact-cmb'
    # advact_fgname = advact_name + '+tsz'
    # ell = np.load(os.path.join(advact_spec_dir, 'ells.npy'))
    # cl_tt = np.load(os.path.join(advact_spec_dir, 'cl_tt.npy'))

    # advact_map = gen_map_from_curve(ell, cl_tt, cmb, basemap,
    #     advact_name, seed=seed)
    # advact_fg_map = gen_cmb_fg_map(fg_fits, advact_fgname, cmb_fm=advact_map)
    # maplist.append(advact_map)
    # maplist.append(advact_fg_map)

    # save in seed dir
    # save_flatmap(cmbmap, seeddir, save_image=False)
    # save_flatmap(cmb_lens_map, seeddir, save_image=False)
    save_flatmap(cmb_fg_map, seeddir, save_image=False)
    # save_flatmap(advact_map, seeddir, save_image=False)
    # save_flatmap(advact_fg_map, seeddir, save_image=False)

    return maplist

def gen_filtered_maps(maplist, seeddir):
    """applies LPF, HPF filtering, saves these new maps in the seeddir.
    @return the list of directory paths for these maps"""

    map_paths = []
    for cmap in maplist:
        beamedmap = apply_beam(cmap, 1.6)
        lpfmap, hpfmap = apply_filtering(beamedmap, filter_type='will')
        lpfpath = save_flatmap(lpfmap, path=seeddir, save_image=False)
        map_paths.append(lpfpath)
        # hpfpath = save_flatmap(hpfmap, path=seeddir, save_image=False)
        # map_paths.append(hpfpath)

    return map_paths

def setup_for_ts(test):
    u, massConv = initialize()

    if test:
        nObj = 10
    else:
        nObj = None
    galcat = Catalog(
        u,
        massConv,
        name='cmass_m_10x10_v2',
        nObj=nObj,
        workDir='/home/theo/Documents/research/CMB/patchy_tau_sims'
    )

    return u, galcat

def end2end(args, u, galcat, kappa_fits, tsz_fits, tsz_hpf_path, runname, rundir, seed):
    """runs pipeline for a single seed"""

    seedname, seeddir = setup_seed_dir(rundir, runname, seed)

    # TODO: add branch to skip map generation, just return
    # paths to relevant maps
    maplist = gen_unfiltered_maps(seed, seedname, seeddir, kappa_fits, tsz_fits)

    map_paths = gen_filtered_maps(maplist, seeddir)
    print('map_paths:\n', map_paths)
    #map_paths scheme:
    # [CMB LPF, CMB HPF,
    # lensCMB LPF, lensCMB HPF,
    # CMB+tSZ LPF, CMB+tSZ HPF]

    # setup for 4 thumbstack runs
    # map_pairs = [(map_paths[0], map_paths[1]), # CMB LPF    x CMB HPF
    #              (map_paths[0], map_paths[1]), # CMB LPF    x CMB+lens HPF
    #              (map_paths[0], tsz_hpf_path),  # CMB LPF    x fg HPF
    #              (map_paths[2], tsz_hpf_path)]  # CMB+fg LPF x fg HPF
    # not sure about these indices above... let's redo lensing and CMBxfg
    # map_pairs = [(map_paths[0], map_paths[1]), # CMB LPF x CMB HPF
    #              (map_paths[0], map_paths[3]), # CMB LPF x lensCMB HPF
    #              (map_paths[2], map_paths[3]), # lensCMB LPF x lensCMB HPF
    #              (map_paths[4], tsz_hpf_path)] # CMB+fg LPF x fg HPF
    # [cmbmap, cmb_lens_map, advact_map, advact_fg_map]
    # map_pairs = [(map_paths[0], map_paths[1]), # CMB LPF x CMB HPF
    #              (map_paths[2], map_paths[3]), # lensCMB LPF x lensCMB HPF
    #              (map_paths[4], map_paths[5]), # ACT LPF x ACT HPF
    #              (map_paths[4], tsz_hpf_path), # ACT LPF x fg HPF
    #              (map_paths[6], tsz_hpf_path)] # ACT+fg LPF x fg HPF
    # for hi-tsz+lens_final_100nmocks
    # map_pairs = [(map_paths[0], map_paths[1]), # CMB LPF x CMB HPF
    #              (map_paths[2], map_paths[3]), # lensCMB LPF x lensCMB HPF
    #              (map_paths[4], tsz_hpf_path)] # CMB+fg LPF x fg HPF
    # for low-tsz run
    map_pairs = [(map_paths[0], tsz_hpf_path)] # CMB+fg LPF x fg HPF

    # suffixes for thumbstack naming
    # map_suffix = ['_cmb-lpf_cmb-hpf',
    #               '_cmb-lpf_cmb+lens-hpf',
    #               '_cmb-lpf_tsz-hpf',
    #               '_cmb+tsz-lpf_tsz-hpf']
    # map_suffix = ['_cmb-lpf_cmb-hpf',
    #               '_cmb-lpf_cmb+lens-hpf',
    #               '_cmb+lens-lpf_cmb+lens-hpf',
    #               '_cmb+tsz-lpf_tsz-hpf']
#     map_suffix = ['_cmb-lpf_cmb-hpf',
#                   '_lens-cmb-lpf_lens-cmb-hpf',
#                   '_advact-lpf_advact-hpf',
#                   '_advact-lpf_tsz-hpf',
#                   '_advact+tsz-lpf_tsz-hpf']
    # for hi-tsz+lens_final_100nmocks
#     map_suffix = ['_cmb-lpf_cmb-hpf',
#                   '_cmb+lens-lpf_cmb+lens-hpf',
#                   '_cmb+tsz-lpf_tsz-hpf']
    map_suffix = ['_cmb+tsz-lpf_tsz-hpf']


    # build dictionaries for the nested all_stack_dict
    seed_dict = {}
    pair_dict = {}
    est_dict = {}

    # for each seed 4 correlations are run
    for pair, suffix in zip(map_pairs, map_suffix):
        print('--------------------------------------------------------------------------------')
        tsname = seedname + suffix
        lpf_enmap, hpf_enmap, boxmask = setup_maps(pair[0], pair[1])

        # run thumbstack, no bootstrap,
        # outdirs: seeddir/output/thumbstack, seeddir/figures/thumbstack
        print('--------------------------------------------------------------------------------')
        print('Beginning thumbstack:', tsname)
        print('--------------------------------------------------------------------------------')
        ts = ThumbStack(
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
            test=False,
            doStackedMap=False,
        )

        for est in ['ti', 'sgn']:
            est_dict.update(
                {est:
                    {'stack':ts.stackedProfile['tauring_tau_%s_uniformweight'%est],
                     'sStack':ts.sStackedProfile['tauring_tau_%s_uniformweight'%est]
                    }
                }
            )
            pair_dict.update({suffix[1:]:est_dict})
            seed_dict.update({seed:pair_dict})

    return seed_dict

def main():
    args = parse_args()
    print(args)

    # same for all seeds
    kappa_fits, tsz_fits, tsz_hpf_path = setup_fixed_dirs()
    runname, rundir = setup_run_dir(args)
    write_log(args, os.path.join(rundir, runname))

    u, galcat = setup_for_ts(args.test)

    # we'll store all the stacked profiles in a dict
    # 1st level: seed
    # 2nd level: n thumbstack types (cmb-cmb, cmb-cmb+lens etc)
    # 3rd level: TI, sgn estimators
    # 4th level: stackedProfile, sStackedProfile
    # 3-4th levels already get made in thumbstack processing
    all_stack_dict = {}

    with sharedmem.MapReduce(np=args.njobs) as pool:
        f = lambda seed: end2end(
            args,
            u,
            galcat,
            kappa_fits,
            tsz_fits,
            tsz_hpf_path,
            runname,
            rundir,
            seed
        )

        stacks_list = pool.map(f, list(range(args.start_seed, args.start_seed+args.nmocks)))

    dictfn = os.path.join(rundir, '%s_all-stacked-profiles.pkl'%runname)
    for seed_dict in stacks_list:
        all_stack_dict.update(seed_dict)

    with open(dictfn, 'wb') as f:
        pickle.dump(all_stack_dict, f)

if __name__ == '__main__':
    main()
