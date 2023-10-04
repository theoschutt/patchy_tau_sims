#!/bin/bash/python
"""compute_mock_stacks.py: Computes stacks for all 100 mocks using arbitrary
mix of T_small, T_large_num (numerator), and T_large_den (denominator).
"""
import os, sys
import glob
import numpy as np
from compute_stack import compute_stack, write_log, load_photometry, save_stack

def main():

    # T_small, etc determined by parent directory name
    args = {
        't_small' : 'cmb-lpf_tsz-hpf',
        't_large_num' : 'cmb-lpf_tsz-hpf',
        'tag' : 'cmb_x_tsz_eqsgn_w_detnoise'
    }

    # get the files
    # dir_stem = ('/home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/tsz+lens_100mocks_eqsgn_v2.1/tsz+lens*seed*/output/thumbstack/')
    dir_stem = ('/media/theo/1TB_WD_Passport/2023_ubuntu_backup/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/hi-tsz+lens_100mocks_eqsgn/tsz+lens*seed*/output/thumbstack/')
    outdirs = glob.glob(dir_stem)
    t_small_list = glob.glob(dir_stem + '*_%s/tauring_filtmap.txt'%(args['t_small']))
    t_large_num_list = glob.glob(dir_stem + '*_%s/tauring_filtnoisestddev.txt'%(args['t_large_num']))
    # just going to use this one for all of runs
    t_large_den_fn = ('/home/theo/Documents/research/CMB/patchy_tau_sims/output/thumbstack/cmass_m_10x10_advACT_nilc/tauring_filtnoisestddev.txt')
    # sort
    outdirs.sort()
    t_small_list.sort()
    t_large_num_list.sort()
    print(outdirs[40])
    print(t_small_list[40])
    print(t_large_num_list[40])

    # basically do compute_stack main() for each seed
    for out, ts_fn, tln_fn in zip(outdirs, t_small_list, t_large_num_list):
        # make new dir in the seed's thumbstack dir
        path_out = os.path.join(out, args['tag'])
        if not os.path.exists(path_out):
            os.makedirs(path_out)

        write_log(args, path_out, args['tag'])

        ts, tln, tld = load_photometry(ts_fn, tln_fn, t_large_den_fn)

        for est in ['ti', 'sgn']:
            stack, sStack = compute_stack(ts, tln, tld, est, None, True)
            save_stack(stack, sStack, est, path_out)

if __name__ == '__main__':
    main()
