#!/usr/bin/env python3
import os,sys
import pickle
import sharedmem
import numpy as np
import fitsio

sys.path.append('../../ThumbStack')
from catalog import Catalog
from thumbstack import ThumbStack

from make_noise_maps import make_map, make_cmb, gen_map_from_fn, gen_cmb_fg_map, gen_map_from_data_powspec
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
    parser.add_argument('--paralleltype',
        type=str,
        help=("What to parallelize over. Options are ['seed', 'catalog']."
              "'seed' good for smallish (<1M entries) catalog."
              "'catalog' good for large catalogs."))
    parser.add_argument('--start_seed',
        type=int,
        help='Value to use as the first seed in sequence.')
    parser.add_argument('--version',
        type=str,
        default='01',
        help='Version tag for directory naming.')
    parser.add_argument('--runtype',
        type=str,
        help="Type of run. Options are ['act', 'spt', 'so', 's4', 'bias']")
    parser.add_argument('--runname',
        type=str,
        help="Run name template for directory naming. Don't include no. of mocks.")
    parser.add_argument('--catname',
        type=str,
        help="Catalog name (i.e. directory at ../output/catalog).")
    parser.add_argument('--cattype',
        type=str,
        help="Catalog type. Options are ['radec', 'cmass']")
    parser.add_argument('--gen_maps',
        default=False, action='store_const', const=True,
        help='Generate maps instead of loading presaved ones.')
    parser.add_argument('--test',
        default=False, action='store_const', const=True,
        help='Test script using first 10 objects in catalog')
    parser.add_argument('--t_large_min',
        default=None,
        type=float,
        help='Minimum abs value temperature for T_large in stack. [default: None]')
    parser.add_argument('--equalsignedweights',
        default=False, action='store_const', const=True,
        help='Whether to equalize number of +/- weights')
    parser.add_argument('--outdir',
        default='../output/multi_mock_runs',
        type=str,
        help="Where to save the parent run directory. [default: '../output/multi_mock_runs']")
    args = parser.parse_args()

    return args

def write_log(args, outfile):
    # Write text file logging what command line args were used
    log_fn =  outfile + '_args.log'
    print('Writing argument log file:', log_fn)
    arg_dict = vars(args)
    with open(log_fn, 'w') as f:
        for arg in arg_dict:
            f.write(f'{str(arg)}: {arg_dict[arg]}\n')

def setup_fixed_dirs():
    fg_dir = '../output/fg_maps'

    kappa_fits = os.path.join(
        fg_dir,
        'kappamap_0.53arcminGauss_correlatedwithTaumap_cmass_10x10_image.fits'
    )

    fg_fits = os.path.join(
        fg_dir,
        'tSZmap_0.88arcminGauss_multiplicativefact1.3_correlatedwithTaumap_cmass_10x10_image.fits'
    )
    # fg_fits = os.path.join(
    #     fg_dir,
    #     'tSZmap_2.46arcminGauss_correlatedwithTaumap_cmass_10x10_image.fits'
    # )

    # fg_hpf_path = os.path.join(
    #     fg_dir,
    #     'tSZ_sig0.88_x1.3_corr-w-tau',
    #     'tSZ_sig0.88_x1.3_corr-w-tau_beam1.6_hpf2425w_flatmap.fits'
    # )
    # fg_hpf_path = os.path.join(
    #     fg_dir,
    #     'tSZ_sig2.46_corr-w-tau',
    #     'tSZ_sig2.46_corr-w-tau_beam1.6_hpf2425w_lmax5950_flatmap.fits'
    # )

    # advact_spec_dir = '/home/theo/Documents/research/CMB/patchy_tau_sims/data/AdvACT_NILC_cls_fullRes_TT'

    return kappa_fits, fg_fits #, fg_hpf_path

def setup_run_dir(args):
    runname = f'{args.runname}_{args.nmocks}mocks'

    if args.equalsignedweights:
        runname += '_eqsgn'
    if args.t_large_min is not None:
        runname += f'_tmin{args.t_large_min:.1f}'
    if args.test:
        runname += '_test'
    runname += '_'+args.version

    rundir = os.path.join(
        args.outdir,
        runname
    )

    if not os.path.exists(rundir):
        print('Creating directory:', rundir)
        os.makedirs(rundir)

    return runname, rundir

def setup_seed_dir(rundir, runname, seed):
    # make map dir with seed in name
    seedname = f'{runname}_seed{seed}'
    seeddir = os.path.join(rundir, seedname)
    if not os.path.exists(seeddir):
        print('Creating directory:', seeddir)
        os.makedirs(seeddir)

    return seedname, seeddir

def gen_advact_map(seed, seedname, seeddir):
    """generates and saves a CMB GRF map with the AdvACT NILC power spectrum.
    NILC PS includes CMB and beam already."""

    # make base flatmap for all map generation
    basemap = make_map()
    cmb = make_cmb()
 
    print('Loading AdvACT NILC Spectrum.')
    advact_spec_dir = '/home/groups/roodman/schutt20/cmb/patchy_tau_sims/data/AdvACT_NILC_cls_fullRes_TT'
    ell = np.load(os.path.join(advact_spec_dir, 'ells.npy'))
    cl_tt = np.load(os.path.join(advact_spec_dir, 'cl_tt.npy'))

    advact_name = seedname + '_advact-cmb'
    advact_map = gen_map_from_data_powspec(ell, cl_tt, cmb, basemap,
        advact_name, seed=seed)

    # save in seed dir
    save_flatmap(advact_map, seeddir, save_image=False)

    return [advact_map]

def gen_spt_map(seed, seedname, seeddir):
    """generates and saves a map with lensed CMB + SPT noise power spectrum.
    Beam is not applied here."""

    basemap = make_map()
    cmb = make_cmb()
 
    print('Loading SPT-3G spectrum.')
    fn = ('/home/groups/roodman/schutt20/cmb/patchy_tau_sims/data/noisecurve_subset/'
          'spt_proposal_2023_ilc_cmb_90-150-220_TT-EE_fsky1500_5years.npy')
    sptdict = np.load(fn, allow_pickle=True).tolist()
    ell = sptdict['el']
    cl_tt = sptdict['cl_residual']['TT']
 
    # add lensed CMB
    cl_cmb = cmb.flensedTT(ell)
    cl_tot = cl_tt + cl_cmb

    spt_name = seedname + '_spt-cmb'
    spt_map = gen_map_from_data_powspec(ell, cl_tot, cmb, basemap,
        spt_name, seed=seed)

    # save in seed dir
    save_flatmap(spt_map, seeddir, save_image=False)

    return [spt_map]

def gen_so_map(seed, seedname, seeddir):
    """generates and saves a map with lensed CMB + SO noise power spectrum.
    Beam is not applied here."""
  
    basemap = make_map()
    cmb = make_cmb()
    print('Loading SO spectrum.')
    so_fn = ('/home/groups/roodman/schutt20/cmb/patchy_tau_sims/data/noisecurve_subset/'
             'SOV3_T_default1-4-2_noisecurves_deproj0_SENS2_mask_16000_ell_TT_yy.txt')
    so_ps = np.loadtxt(so_fn)
    ell = so_ps[:,0]
    cl_tt = so_ps[:,1]

    # add lensed CMB
    cl_cmb = cmb.flensedTT(ell)
    cl_tot = cl_tt + cl_cmb

    so_name = seedname + '_so-cmb'
    so_map = gen_map_from_data_powspec(ell, cl_tot, cmb, basemap,
        so_name, seed=seed)

    save_flatmap(so_map, seeddir, save_image=False)

    return [so_map]

def gen_s4_map(seed, seedname, seeddir):
    """generates and saves a map with lensed CMB + SO noise power spectrum.
    Beam is not applied here."""

    basemap = make_map()
    cmb = make_cmb()
    # make power spectra for map
    # fdetectorNoise is debeamed
    print('Generating CMB-S4 spectrum.')
    lMin = 30.
    lMax = 10000.
    ell = np.logspace(np.log10(lMin/2.), np.log10(2.*lMax), 1001, 10.)
    cl_tot = cmb.flensedTT(ell) + cmb.fdetectorNoise(ell)

    s4_name = seedname + '_s4-cmb'
    s4_map = gen_map_from_data_powspec(ell, cl_tot, cmb, basemap, s4_name, seed=seed)

    save_flatmap(s4_map, seeddir, save_image=False)

    return [s4_map]

def gen_unfiltered_maps(seed, seedname, seeddir, kappa_fits, fg_fits):
    """generates and saves CMB GRF map, same GRF lensed by the `kappa_fits` kappa map,
    and the unlensed GRF plus the `fg_fits` FG map"""

    # make base flatmap for all map generation
    basemap = make_map()
    maplist = []

    # make unlensed CMB GRF map
    cmb = make_cmb()
    cmbname = seedname + '_cmb'
    cmbmap = gen_map_from_fn(cmb.funlensedTT, cmb, basemap, cmbname, seed=seed)
    maplist.append(cmbmap)

    # make lensed CMB GRF
    lensname = cmbname + '+lens'
    kappaData = fitsio.read(kappa_fits)
    kappaFourier = cmbmap.fourier(kappaData)
    cmb_lens_map = basemap.copy()
    cmb_lens_map.data = cmbmap.doLensing(kappaFourier=kappaFourier)
    cmb_lens_map.dataFourier = cmb_lens_map.fourier()
    cmb_lens_map.name = lensname
    maplist.append(cmb_lens_map)

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
    save_flatmap(cmbmap, seeddir, save_image=False)
    save_flatmap(cmb_lens_map, seeddir, save_image=False)
    save_flatmap(cmb_fg_map, seeddir, save_image=False)
    # save_flatmap(advact_map, seeddir, save_image=False)
    # save_flatmap(advact_fg_map, seeddir, save_image=False)

    return maplist

def gen_filtered_maps(maplist, seeddir, beam=None):
    """applies LPF, HPF filtering, saves these new maps in the seeddir.
    @return the list of directory paths for these maps"""

    map_paths = []
    for cmap in maplist:
        if beam is not None:
            beamedmap = apply_beam(cmap, beam)
        else:
            beamedmap = cmap
        lpfmap, hpfmap = apply_filtering(beamedmap, filter_type='will')
        lpfpath = save_flatmap(lpfmap, path=seeddir, save_image=False)
        map_paths.append(lpfpath)
        hpfpath = save_flatmap(hpfmap, path=seeddir, save_image=False)
        map_paths.append(hpfpath)

    return map_paths

def get_map_paths(expt, seed):
    """returns the paths to maps already saved in previous runs."""
    # ../output/multi_mock_runs/act-nilc_x_unwise-b+g_10x10_v2_128mocks_eqsgn/act-nilc_x_unwise-b+g_10x10_v2_32mocks_eqsgn_0062/act-nilc_x_unwise-b+g_10x10_v2_32mocks_eqsgn_0062_seed2000/act-nilc_x_unwise-b+g_10x10_v2_32mocks_eqsgn_0062_seed2000_advact-cmb_hpf2425w_lmax5950_flatmap.fits
    # tag = 'act-nilc_x_unwise-b+g_10x10_v2'
    # suffix = '32mocks_eqsgn'
    # in_rundir = f'/home/groups/roodman/schutt20/cmb/patchy_tau_sims/output/multi_mock_runs/{tag}_128mocks_eqsgn'
    # map_seeddir = f'{tag}_{suffix}_{int(start_seed/32):04}/{tag}_{suffix}_{int(start_seed/32):04}_seed{seed}'
    # in_cmbname = f'{tag}_{suffix}_{int(start_seed/32):04}_seed{seed}_advact-cmb'

    long_names = {
        'act': 'advact-nilc',
        'spt': 'spt3g-ilc',
        'so': 'so-sens2-dp0',
        's4': 's4-nofg'
    }
    beams = {
        'act': '1.6',
        'spt': '1.6',
        'so': '1.4',
        's4': '1.0'
    }
 
    scratch_dir = "/scratch/users/schutt20/cmb/patchy_tau_sims/output"
    expt_prefix = f"{expt}_10x10_seed2000-2127_maps/{long_names[expt]}_seed{seed}_beam{beams[expt]}"
    # lpfname = in_cmbname + '_lpf2075w'
    # hpfname = in_cmbname + '_hpf2425w_lmax5950'
    # lpf_path = os.path.join(in_rundir, map_seeddir, f'{lpfname}_flatmap.fits')
    # hpf_path = os.path.join(in_rundir, map_seeddir, f'{hpfname}_flatmap.fits')
    lpfname = expt_prefix + '_lpf2075w'
    hpfname = expt_prefix + '_hpf2425w_lmax5950'
    lpf_path = os.path.join(scratch_dir, f'{lpfname}_flatmap.fits')
    hpf_path = os.path.join(scratch_dir, f'{hpfname}_flatmap.fits')

    return [lpf_path, hpf_path]

def setup_for_ts(catname, cattype, test):
    u, massConv = initialize()

    if test:
        nObj = 10
    else:
        nObj = None
    galcat = Catalog(
        u,
        massConv,
        catType=cattype,
        name=catname,
        nObj=nObj,
        workDir='..'
    )

    return u, galcat

def run_parallel_on_seeds(args, u, galcat, kappa_fits, tsz_fits, runname, rundir):
    """Run thumbstack while parallelizing over seeds (1 seed per core). Useful for when the
    catalog is small (<1M entries) and running on one core will finish in a reasonable amount
    of time (~hours).
    """
    # we'll store all the stacked profiles in a dict
    # 1st level: seed
    # 2nd level: n thumbstack types (cmb-cmb, cmb-cmb+lens etc)
    # 3rd level: TI, sgn estimators
    # 4th level: stackedProfile, sStackedProfile
    # 3-4th levels already get made in thumbstack processing
    all_stacks_dict = {}

    with sharedmem.MapReduce(np=args.njobs) as pool:
        f = lambda seed: end2end(
            args,
            u,
            galcat,
            kappa_fits,
            tsz_fits,
            runname,
            rundir,
            nproc=1,  # we're parallelizing over seeds
                      # so each seed needs to use exactly one thread
            seed=seed
        )

        stacks_list = pool.map(f, list(range(args.start_seed, args.start_seed+args.nmocks)))

    for seed_dict in stacks_list:
        all_stacks_dict.update(seed_dict)

    dictfn = os.path.join(rundir, f'{runname}_all-stacked-profiles.pkl')
    with open(dictfn, 'wb') as f:
        pickle.dump(all_stacks_dict, f)

def run_parallel_on_catalog(args, u, galcat, kappa_fits, tsz_fits, runname, rundir):
    """Run thumbstack while parallelizing over the catalog. Useful for when the catalog is too
    large for each core to hold the full catalog and/or we want to save at shorter time intervals
    (i.e. finish one seed in ~1/n_cores time).
    """
    all_stacks_dict = {}
    for seed in range(args.start_seed, args.start_seed+args.nmocks): 
        seed_dict = end2end(
            args,
            u,
            galcat,
            kappa_fits,
            tsz_fits,
            runname,
            rundir,
            nproc=args.njobs,
            seed=seed
        )

        all_stacks_dict.update(seed_dict)

    dictfn = os.path.join(rundir, f'{runname}_all-stacked-profiles.pkl')
    with open(dictfn, 'wb') as f:
        pickle.dump(all_stacks_dict, f)

def end2end(args, u, galcat, kappa_fits, tsz_fits, runname, rundir, nproc, seed):
    """runs pipeline for a single seed"""

    seedname, seeddir = setup_seed_dir(rundir, runname, seed)

    # set up input maps
    if args.gen_maps:
        # for bias runs
        if args.runtype == 'bias':
            maplist = gen_unfiltered_maps(seed, seedname, seeddir, kappa_fits, tsz_fits)
            map_paths = gen_filtered_maps(maplist, seeddir, beam=1.6)

	# for sim covmat runs
        elif args.runtype == 'act':
            maplist = gen_advact_map(seed, seedname, seeddir)
            map_paths = gen_filtered_maps(maplist, seeddir, beam=None)
        elif args.runtype == 'spt':
            maplist = gen_spt_map(seed, seedname, seeddir)
            map_paths = gen_filtered_maps(maplist, seeddir, beam=1.6)
        elif args.runtype == 'so':
            maplist = gen_so_map(seed, seedname, seeddir)
            map_paths = gen_filtered_maps(maplist, seeddir, beam=1.4)
        elif args.runtype == 's4':
            maplist = gen_s4_map(seed, seedname, seeddir)
            map_paths = gen_filtered_maps(maplist, seeddir, beam=1.0)
        else:
            raise ValueError('Invalid runtype:', args.runtype) 
    else:
        map_paths = get_map_paths(args.runtype, seed)
    print('map_paths:\n', map_paths)

    # set up paths to maps and directory naming suffixes
    # for 1000 bias mocks run
    if args.runtype == 'bias':
        map_pairs = [(map_paths[0], map_paths[1]), # CMB LPF x CMB HPF
                     (map_paths[2], map_paths[3]), # lensCMB LPF x lensCMB HPF
                     (map_paths[4], map_paths[5])] # CMB+fg LPF x CMB+fg HPF
        map_suffix = ['_cmb-lpf_cmb-hpf',
                      '_cmb+lens-lpf_cmb+lens-hpf',
                      '_cmb+tsz-lpf_tsz-hpf']
    # for unwise and LSST covariance sim runs
    else:
        map_pairs = [(map_paths[0], map_paths[1])] # CMB LPF x CMB HPF
        map_suffix = ['_cmb-lpf_cmb-hpf']

    # build dictionaries for the nested all_stack_dict
    seed_dict = {}
    pair_dict = {}
    est_dict = {}

    # get one bootstrap covariance matrix from the first seed for comparison
    # do_bootstrap = (seed == 2000)
    do_bootstrap = True

    # for each seed N=len(map_pairs) correlations are run
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
            nProc=nproc,
            filterTypes='tauring',
            estimatorTypes=['tau_ti_uniformweight', 'tau_sgn_uniformweight'],
            doBootstrap=do_bootstrap,
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
                    {'stack':ts.stackedProfile[f'tauring_tau_{est}_uniformweight'],
                     'sStack':ts.sStackedProfile[f'tauring_tau_{est}_uniformweight']
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
    kappa_fits, tsz_fits = setup_fixed_dirs()
    runname, rundir = setup_run_dir(args)
    write_log(args, os.path.join(rundir, runname))

    u, galcat = setup_for_ts(args.catname, args.cattype, args.test)

    if args.paralleltype == 'seed':
        run_parallel_on_seeds(args, u, galcat, kappa_fits, tsz_fits, runname, rundir)
    elif args.paralleltype == 'catalog':
        run_parallel_on_catalog(args, u, galcat, kappa_fits, tsz_fits, runname, rundir)

if __name__ == '__main__':
    main()
