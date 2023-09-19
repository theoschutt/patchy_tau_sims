#!/bin/bash/python
import os
import glob
import numpy as np

mocks_dir = ('/home/theo/Documents/research/CMB/patchy_tau_sims/'
             'output/multi_mock_runs/hi-tsz+lens_100mocks_eqsgn/')

# def collate_dicts():
#     # glob the pkl files
#     pkl_list = glob.glob(mocks_dir + '*.pkl')
#     _ = [print(pkl) for pkl in pkl_list]
#     # iterate over pkl files
#     # update top empty dict with the saved dict
#     tot_dict = {}
#     for pkl in pkl_list:
#         with open(pkl, 'rb') as f:
#             run_dict = pickle.load(f)
#             tot_dict.update(run_dict)
#
#     # voila: combined 100 mock dict
#     _ = [print(key) for key in tot_dict.keys()]
#
#     # and save
#     with open(mocks_dir + 'hi-tsz+lens_100mocks_eqsgn_collated-stacked-profiles.pkl') as f:
#         pickle.dump(tot_dict, f)

# collated_dict has 8 sets of thumbstack outputs
# *-lpf_*-hpf x (ti or sgn)
corr_type = ['cmb-lpf_cmb-hpf',
             'cmb-lpf_cmb+lens-hpf',
             'cmb-lpf_tsz-hpf',
             'cmb+tsz-lpf_tsz-hpf']
est_type = ['ti', 'sgn']

for corr in corr_type:
    for est in est_type:
        print(corr, est)
        stack_arr = np.zeros((0,9))
        std_arr = np.zeros((0,9))
        fn_list = glob.glob(
            os.path.join(
                mocks_dir,
                'tsz+lens_*_seed*/output/thumbstack/*_%s'%corr,
                'tauring_tau_%s_uniformweight_measured.txt'%est)
        )
        print(len(fn_list))

        # collate profiles
        for fn in fn_list:
            seed_info = np.genfromtxt(fn)
            seed_stack = seed_info[:,1]
            seed_std = seed_info[:,2]
            stack_arr = np.vstack((stack_arr, seed_stack))
            std_arr = np.vstack((std_arr, seed_std))
        print('stack_arr, std_arr:', stack_arr.shape, std_arr.shape)

        # save collated profiles
        stack_fn = os.path.join(mocks_dir, '%s_%s_collated-stack.txt'%(corr, est))
        std_fn = os.path.join(mocks_dir, '%s_%s_collated-stddev.txt'%(corr, est))
        np.savetxt(stack_fn, stack_arr)
        np.savetxt(std_fn, std_arr)

        # calculate average profile
        tot_stack = np.average(stack_arr, axis=0)
        # error on the mean => include 1/sqrt(N)
        tot_std = np.sqrt(np.sum(std_arr**2, axis=0) / len(fn_list))

        # save average profile
        # to make output file same as thumbstack's
        # we include radial bin value as first col
        r = np.linspace(1., 6., 9)
        tot_stack_info = np.vstack((r, tot_stack, tot_std)).T
        print('tot_stack_info:', tot_stack_info.shape)
        np.savetxt(mocks_dir+'%s_%s_avg-stack-measured.txt'%(corr, est), tot_stack_info)

# def calc_stack(stack_fn, std_fn):
#
#     stack_info = np.genfromtxt(fn)
#     stack_arr = stack_info[:,1]
#     std_arr = stack_info[:,2]
#
#     stack = np.average(stack_arr, axis=0)
#     std = np.sqrt(len(collated_dict)) * np.average(std_arr, axis=0)
#
#     stack_info = np.vstack([r, stack, std]).T
#
#     np.savetxt('hi-tsz+lens_100mocks_%s_%s_measured.txt'%(corr, est), stack_info, fmt='%.5e')
