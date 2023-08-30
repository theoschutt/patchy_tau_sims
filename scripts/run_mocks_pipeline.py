#!/usr/bin/env python3
import os,sys
import argparse

sys.path.append('../../ThumbStack')
from catalog import Catalog
from thumbstack import ThumbStack

from make_noise_maps import make_map, make_cmb, gen_map_from_fn, gen_cmb_fg_map, save_flatmap
from filter_tau_map import apply_filtering
from run_thumbstack import initialize, setup_maps

# info same for all seeds
FG_DIR = '/home/theo/Documents/research/CMB/patchy_tau_sims/output/fg_maps'
fg_type = args.fg_type
if fg_type == 'tsz':
    fg_fits = os.path.join(
        FG_DIR,
        'tSZmap_0.88arcminGauss_multiplicativefact1.3_correlatedwithTaumap_cmass_10x10_image.fits'
    )
    fg_hpf_path = os.path.join(
        FG_DIR,
        'tSZ_sig0.88_x1.3_corr-w-tau',
        'tSZ_sig0.88_x1.3_corr-w-tau_beam1.6_flatmap.fits'
    )

test = args.test
nmocks = args.nmocks
runname = '%s_x_%imocks'%(fg_type, nmocks)
rundir = os.path.join(
    '/home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs',
    runName
)
# TODO: write log at rundir level

# make map dir with seed in name
seed = 
seedname = '%imocks_seed%s_cmb'%(nMocks, seed)
newdir = os.path.join(rundir, seedname)
if not os.path.exists(newdir):
    os.makedirs(newdir)

# make base flatmap for all map generation
basemap = make_map()

# make unlensed CMB GRF map, save in map dir
cmb = make_cmb()
Ctot = lambda l: cmb.funlensedTT(l)
cmbmap = gen_map_from_fn(Ctot, cmb, basemap, seedname, seed=seed)

if test:
    # includes saving diagnostics: png image and power spectrum plot
    save_flatmap(cmbmap, cmb, newdir, seedname)
else:
    cmbmap.write(newdir)

# make unlensed CMB GRF + fg map, save in map dir
fgname = seedname+'+%s'%fg_type
cmb_fg_map = gen_cmb_fg_map(fg_fits, fgname, cmb_fm=cmbmap)

if test:
    # includes saving diagnostics: png image and power spectrum plot
    save_flatmap(cmb_fg_map, cmb, newdir, fgname)
else:
    cmb_fg_map.write(newdir)

# filter GRF, GRF+fg map, save
map_paths = [] # order: cmb_lpf, cmb_hpf, cmb+fg_lpf, cmb+fg_hpf
for cmap in [cmbmap, cmb_fg_map]:
    lpfmap, hpfmap = apply_filtering(cmap, filter_type='will')
    lpfpath = save_flatmap(lpfmap, path=newdir, save_diagnostics=test)
    map_paths.append(lpfpath)
    hpfpath = save_flatmap(hpfmap, path=newdir, save_diagnostics=test)
    map_paths.append(hpfpath)

# setup for thumbstack
u, massConv = initialize()

galcat = Catalog(
    u,
    massConv,
    name='cmass_m_10x10_v2',
    nObj=nObj,
    workDir='/home/theo/Documents/research/CMB/patchy_tau_sims'
)

map_pairs = [(map_paths[0], map_paths[1]), # CMB LPF    x CMB HPF
             (map_paths[0], fg_hpf_path),  # CMB LPF    x fg HPF
             (map_paths[2], fg_hpf_path)]  # CMB+fg LPF x fg HPF

# suffixes for thumbstack naming
map_suffix = ['_cmb-lpf_cmb-hpf',
              '_cmb-lpf_%s-hpf'%fg_type,
              '_cmb+%s-lpf_%s-hpf'%(fg_type, fg_type)]

for i, pair in enumerate(map_pairs):
    lpf_enmap, hpf_enmap, boxmask = setup_maps(pair[0], pair[1])

    # run thumbstack, no bootstrap,
    # outdirs: newdir/output/thumbstack, newdir/figures/thumbstack
    ts = ThumbStack(
        u,
        galcat,
        hpf_enmap,
        boxmask,
        cmbHit=None,
        cmbMap2=lpf_enmap,
        name='%imocks_seed%s'%(nMocks, seed) + map_suffix[i],
        save=True,
        nProc=args.nproc,
        filterTypes='tauring',
        estimatorTypes=['tau_ti_uniformweight', 'tau_sgn_uniformweight'],
        doBootstrap=False,
        tLargeMin=args.t_large_min,
        equalSignedWeights=args.equalsignedweights,
        workDir=newdir,
        runEndToEnd=True,
        test=args.test,
        doStackedMap=False
    )

