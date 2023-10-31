#!/bin/bash/python
"""compute_mock_stacks.py: Computes stacks for all 100 mocks using arbitrary
mix of T_small, T_large_num (numerator), and T_large_den (denominator).
"""
import os
import glob
from compute_stack import compute_stack, write_log, load_photometry, save_stack

def main():

    # T_small, etc determined by parent directory name
    # args = {
    #     't_small' : 'cmb-lpf_cmb-hpf',
    #     't_large_num' : 'cmb-lpf_cmb-hpf',
    #     't_large_den' : 'cmb+lens-lpf_cmb+lens-hpf',
    #     'tag' : 'cmb_x_cmb_lensnorm'
    # }
    # args = {
    #     't_small' : 'cmb-lpf_cmb-hpf',
    #     't_large_num' : 'cmb+lens-lpf_cmb+lens-hpf',
    #     't_large_den' : 'cmb+lens-lpf_cmb+lens-hpf',
    #     'tag' : 'cmb+lens_x_cmb_lensnorm'
    # }
    args = {
        't_small' : 'cmb+lens-lpf_cmb+lens-hpf',
        't_large_num' : 'cmb-lpf_cmb-hpf',
        't_large_den' : 'cmb+lens-lpf_cmb+lens-hpf',
        'tag' : 'cmb_x_cmb+lens_lensnorm'
    }
    # args = {
    #     't_small' : 'cmb+tsz-lpf_tsz-hpf',
    #     't_large_num' : 'cmb-lpf_cmb-hpf',
    #     't_large_den' : 'cmb+tsz-lpf_tsz-hpf',
    #     'tag' : 'cmb_x_tsz_tsznorm'
    # }

    # get the files
    # dir_stem = ('/home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/tsz+lens_100mocks_eqsgn_v2.1/tsz+lens*seed*/output/thumbstack/')
    # dir_stem_orig = ('/home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/hi-tsz+lens_final_100mocks_eqsgn*/hi-tsz+lens_final_*/hi-tsz+lens_final_*mocks_eqsgn_*_seed*/output/thumbstack/')
    dir_stem = ('/home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/hi-tsz+lens_final_200mocks_eqsgn*/hi-tsz+lens_final_*/hi-tsz+lens_final_*mocks_eqsgn_*_seed*/output/thumbstack/')
    # dir_stem = ('/home/theo/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/lo-tsz_final_200mocks_eqsgn*/lo-tsz_final_*/lo-tsz_final_*mocks_eqsgn_*_seed*/output/thumbstack/')
    outdirs = glob.glob(dir_stem)
    # dir_stem = ('/media/theo/1TB_WD_Passport/2023_ubuntu_backup/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/hi-tsz+lens_100mocks_eqsgn/tsz+lens*seed*/output/thumbstack/')
    t_small_list = glob.glob(dir_stem + '*_%s/tauring_filtmap.txt'%(args['t_small']))
    t_large_num_list = glob.glob(dir_stem + '*_%s/tauring_filtnoisestddev.txt'%(args['t_large_num']))
    t_large_den_list = glob.glob(dir_stem + '*_%s/tauring_filtnoisestddev.txt'%(args['t_large_den']))
    # sort
    outdirs.sort()
    t_small_list.sort()
    t_large_num_list.sort()
    t_large_den_list.sort()
    print(outdirs[40])
    print(t_small_list[40])
    print(t_large_num_list[40])
    print(t_large_den_list[40])

    # basically do compute_stack main() for each seed
    for out, ts_fn, tln_fn, tld_fn in zip(outdirs, t_small_list, t_large_num_list, t_large_den_list):
        # make new dir in the seed's thumbstack dir
        path_out = os.path.join(out, args['tag'])
        if not os.path.exists(path_out):
            os.makedirs(path_out)

        write_log(args, path_out, args['tag'])

        ts, tln, tld = load_photometry(ts_fn, tln_fn, tld_fn)

        for est in ['ti', 'sgn']:
            stack, sStack = compute_stack(ts, tln, tld, est, None, True)
            save_stack(stack, sStack, est, path_out)

if __name__ == '__main__':
    main()
