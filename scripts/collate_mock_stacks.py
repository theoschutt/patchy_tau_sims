#!/bin/bash/python
import os
import glob
import numpy as np

MOCKS_DIR = ('/home/theo/Documents/research/CMB/patchy_tau_sims/'
             'output/multi_mock_runs/tsz+lens_100mocks_eqsgn_v2.1/')

MOCKS_DIR_2 = ('/media/theo/1TB_WD_Passport/2023_ubuntu_backup/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/hi-tsz+lens_100mocks_eqsgn')

def calc_stacked_profile(corr, est, mock_dir):
    print(corr, est)
    stack_arr = np.zeros((0,9))
    std_arr = np.zeros((0,9))
    fn_list = glob.glob(
        os.path.join(
            mock_dir,
            'tsz+lens_*_seed*/output/thumbstack/*%s'%corr,
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
    stack_fn = os.path.join(MOCKS_DIR, '%s_%s_collated-stack.txt'%(corr, est))
    std_fn = os.path.join(MOCKS_DIR, '%s_%s_collated-stddev.txt'%(corr, est))
    np.savetxt(stack_fn, stack_arr)
    np.savetxt(std_fn, std_arr)

    # calculate average profile
    tot_stack = np.average(stack_arr, axis=0)
    # error on the mean => include 1/sqrt(N)
    # tot_std = np.sqrt(np.sum(std_arr**2, axis=0) / len(fn_list))
    tot_std = np.std(stack_arr, axis=0) / np.sqrt(len(fn_list))

    # save average profile
    # to make output file same as thumbstack's
    # we include radial bin value as first col
    r = np.linspace(1., 6., 9)
    tot_stack_info = np.vstack((r, tot_stack, tot_std)).T
    print('tot_stack_info:', tot_stack_info.shape)
    np.savetxt(MOCKS_DIR+'%s_%s_avg-stack-measured_v2.txt'%(corr, est), tot_stack_info)

def calc_diff_profile(corr1, corr2, est):
    print(corr1, corr2, est)
    stack1_arr = np.zeros((0,9))
    fn1_list = glob.glob(
        os.path.join(
            MOCKS_DIR,
            'tsz+lens_*_seed*/output/thumbstack/*%s'%corr1,
            'tauring_tau_%s_uniformweight_measured.txt'%est)
    )
    fn2_list = glob.glob(
        os.path.join(
            MOCKS_DIR_2,
            'tsz+lens_*_seed*/output/thumbstack/*%s'%corr2,
            'tauring_tau_%s_uniformweight_measured.txt'%est)
    )
    fn1_list.sort()
    fn2_list.sort()
    print(len(fn1_list), len(fn2_list))

    # collate profiles for corr1 and corr2
   # stack_pair = (stack1_arr, stack2_arr)
    #std1_arr = np.zeros((0,9))
    #std2_arr = np.zeros((0,9))
    stack1_list, stack2_list = [], []
    std1_list, std2_list = [], []
    stack_pair = (stack1_list, stack2_list)
    std_pair = (std1_list, std2_list)
    fn_pair = (fn1_list, fn2_list)
    for fn_list, stack_arr, std_arr in zip(fn_pair, stack_pair, std_pair):
        #print(stack_arr)
        for fn in fn_list:
            seed_info = np.genfromtxt(fn)
            seed_stack = seed_info[:,1]
            seed_std = seed_info[:,2]
            stack_arr.append(seed_stack)
            std_arr.append(seed_std)
            #stack_arr = np.vstack((stack_arr, seed_stack))
            #std_arr = np.vstack((std_arr, seed_std))
            print(len(stack1_list), len(stack2_list))
        print('stack_arr, std_arr:', stack_arr[-1], std_arr[-1])
    stack1_arr = np.array(stack1_list)
    std1_arr = np.array(std1_list)
    stack2_arr = np.array(stack2_list)
    std2_arr = np.array(std2_list)
    print(stack1_arr.shape)
    # stack1_list = []
    # stack2_list = []
    # std1_list = []
    # std2_list = []

    # take difference
    stack_arr = stack1_arr - stack2_arr
    #print(stack1_arr, stack2_arr, stack_arr)
    # get error of difference for each seed
    # this error doesn't really make sense though since we know the 2 stacks are
    # highly correlated
    # std_arr = np.sqrt(std1_arr**2 + std2_arr**2)

    # save collated profiles
    stack_fn = os.path.join(MOCKS_DIR, '%s-MINUS-%s_%s_collated-stack.txt'%(corr1, corr2, est))
    # std_fn = os.path.join(MOCKS_DIR, '%s-MINUS-%s_%s_collated-stddev.txt'%(corr1, corr2, est))
    np.savetxt(stack_fn, stack_arr)
    # np.savetxt(std_fn, std_arr)

    # calculate average profile
    tot_stack = np.average(stack_arr, axis=0)

    # more relevant error: scatter on mean over seeds
    tot_std = np.std(stack_arr, axis=0) / np.sqrt(len(stack_arr))

    # save average profile
    # to make output file same as thumbstack's
    # we include radial bin value as first col
    r = np.linspace(1., 6., 9)
    tot_stack_info = np.vstack((r, tot_stack, tot_std)).T
    print('tot_stack_info:', tot_stack_info.shape)
    np.savetxt(MOCKS_DIR+'%s-MINUS-%s_%s_avg-stack-measured.txt'%(corr1, corr2, est), tot_stack_info)

# corr_type = ['cmb-lpf_cmb-hpf',
#              'cmb-lpf_cmb+lens-hpf',
#              'cmb-lpf_tsz-hpf',
#              'cmb+tsz-lpf_tsz-hpf']

# corr_type = ['cmb-lpf_cmb-hpf',
#              'cmb+lens-lpf_cmb+lens-hpf',
#              'cmb-lpf_cmb+lens-hpf',
#              'cmb+tsz-lpf_tsz-hpf']
# corr_type = ['cmb+lens_x_cmb+lens_eqsgn_w_detnoise',
#              'cmb_x_cmb_eqsgn_w_detnoise']
corr_type = ['cmb+tsz_x_tsz_eqsgn_w_detnoise',
             'cmb_x_tsz_eqsgn_w_detnoise']


# corr_type = ['cmb_x_cmb_eqsgn_w_detnoise',
#              'cmb_x_cmb+lens_eqsgn_w_detnoise',
#              'cmb_x_tsz_eqsgn_w_detnoise',
#              'cmb+tsz_x_tsz_eqsgn_w_detnoise']

# diff_type = [(corr_type[0], corr_type[1]), # lensing bias
#              (corr_type[3], corr_type[2])] # tsz mean bias
diff_type = [(corr_type[0], corr_type[1])] # lensing bias

est_type = ['ti', 'sgn']

for corr, mock_dir in zip(corr_type, [MOCKS_DIR, MOCKS_DIR_2]):
    for est in est_type:
        calc_stacked_profile(corr, est, mock_dir)

for diff in diff_type:
    for est in est_type:
        calc_diff_profile(diff[0], diff[1], est)

