#!/bin/bash/python
import os
import glob
import numpy as np

MOCKS_DIR = ('/home/theo/Documents/research/CMB/patchy_tau_sims/'
             'output/multi_mock_runs/lo-tsz_final_200mocks_eqsgn/')

# MOCKS_DIR = ('/home/theo/Documents/research/CMB/patchy_tau_sims/'
#              'output/multi_mock_runs/hi-tsz+lens_final_200mocks_eqsgn/')
# MOCKS_DIR_2 = ('/media/theo/1TB_WD_Passport/2023_ubuntu_backup/Documents/research/CMB/patchy_tau_sims/output/multi_mock_runs/hi-tsz+lens_100mocks_eqsgn')
# TAG = 'hi-tsz+lens_final'
TAG = 'lo-tsz_final'

def calc_stacked_profile(corr, est):
    print(corr, est)
    stack_arr = np.zeros((0,9))
    std_arr = np.zeros((0,9))
    fn_list = glob.glob(
        os.path.join(
            MOCKS_DIR,
            '%s_*/%s_*_seed*/output/thumbstack/*%s'%(TAG, TAG, corr),
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
    np.savetxt(MOCKS_DIR+'%s_%s_avg-stack-measured.txt'%(corr, est), tot_stack_info)

def calc_diff_profile(corr1, corr2, est):
    print(corr1, corr2, est)
    stack1_arr = np.zeros((0,9))
    fn1_list = glob.glob(
        os.path.join(
            MOCKS_DIR,
            '%s_*/%s_*_seed*/output/thumbstack/*%s'%(TAG, TAG, corr1),
            'tauring_tau_%s_uniformweight_measured.txt'%est)
    )
    fn2_list = glob.glob(
        os.path.join(
            MOCKS_DIR,
            '%s_*/%s_*_seed*/output/thumbstack/*%s'%(TAG, TAG, corr2),
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

def calc_lens_diff_profile(est):
    """calculate the exact mean lensing bias for the TI estimator.
    This is CMBxCMB + lensedCMBxlensedCMB - lensedCMBxCMB - CMBxlensedCMB.
    """
    cmbcmb_fn = glob.glob(
        os.path.join(
            MOCKS_DIR,
            'hi-tsz+lens_final_*/hi-tsz+lens_*_seed*/output/thumbstack/cmb_x_cmb_lensnorm',
            'tauring_tau_%s_uniformweight_measured.txt'%est)
    )
    cmblens_cmb_fn = glob.glob(
        os.path.join(
            MOCKS_DIR,
            'hi-tsz+lens_final_*/hi-tsz+lens_*_seed*/output/thumbstack/cmb+lens_x_cmb_lensnorm',
            'tauring_tau_%s_uniformweight_measured.txt'%est)
    )
    cmb_cmblens_fn = glob.glob(
        os.path.join(
            MOCKS_DIR,
            'hi-tsz+lens_final_*/hi-tsz+lens_*_seed*/output/thumbstack/cmb_x_cmb+lens_lensnorm',
            'tauring_tau_%s_uniformweight_measured.txt'%est)
    )
    cmblens_cmblens_fn = glob.glob(
        os.path.join(
            MOCKS_DIR,
            'hi-tsz+lens_final_*/hi-tsz+lens_*_seed*/output/thumbstack/*cmb+lens-lpf_cmb+lens-hpf',
            'tauring_tau_%s_uniformweight_measured.txt'%est)
    )

    cmbcmb_fn.sort()
    cmblens_cmb_fn.sort()
    cmb_cmblens_fn.sort()
    cmblens_cmblens_fn.sort()
    print(len(cmbcmb_fn), len(cmblens_cmb_fn),
          len(cmb_cmblens_fn), len(cmblens_cmblens_fn))

    # collate profiles for the 4 correlations
    cc_list, cl_c_list, c_cl_list, cl_cl_list = [], [], [], []
    stack_quad = (cc_list, cl_c_list, c_cl_list, cl_cl_list)
    fn_quad = (cmbcmb_fn, cmblens_cmb_fn, cmb_cmblens_fn, cmblens_cmblens_fn)
    for fn_list, stack_arr in zip(fn_quad, stack_quad):
        #print(stack_arr)
        for fn in fn_list:
            seed_info = np.genfromtxt(fn)
            seed_stack = seed_info[:,1]
            stack_arr.append(seed_stack)
        print(len(cc_list), len(cl_c_list), len(c_cl_list), len(cl_cl_list))
        print('stack_arr:', stack_arr[-1])
    cc_arr = np.array(cc_list)
    cl_c_arr = np.array(cl_c_list)
    c_cl_arr = np.array(c_cl_list)
    cl_cl_arr = np.array(cl_cl_list)
    print('cc_arr shape:', cc_arr.shape)

    # take difference
    stack_arr = cc_arr + cl_cl_arr - cl_c_arr - c_cl_arr

    # save collated profiles
    stack_fn = os.path.join(MOCKS_DIR, 'cc+cl_cl-cl_c-c_cl_%s_collated-stack.txt'%est)
    np.savetxt(stack_fn, stack_arr)

    # calculate average profile
    tot_stack = np.average(stack_arr, axis=0)

    # relevant error: scatter on mean over seeds
    tot_std = np.std(stack_arr, axis=0) / np.sqrt(len(stack_arr))

    # save average profile
    # to make output file same as thumbstack's
    # we include radial bin value as first col
    r = np.linspace(1., 6., 9)
    tot_stack_info = np.vstack((r, tot_stack, tot_std)).T
    print('tot_stack_info:', tot_stack_info.shape)
    tot_stack_fn = os.path.join(MOCKS_DIR, 'cc+cl_cl-cl_c-c_cl_%s_avg-stack-measured.txt'%est)
    np.savetxt(tot_stack_fn, tot_stack_info)

def main():
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
    # corr_type = ['cmb+tsz_x_tsz_eqsgn_w_detnoise',
    #              'cmb_x_tsz_eqsgn_w_detnoise']
    # corr_type = ['cmb_x_cmb_eqsgn_w_detnoise',
    #              'cmb_x_cmb+lens_eqsgn_w_detnoise',
    #              'cmb_x_tsz_eqsgn_w_detnoise',
    #              'cmb+tsz_x_tsz_eqsgn_w_detnoise']
    # for hi-tsz+lens run
    # corr_type = ['cmb+tsz-lpf_tsz-hpf',
    #              'cmb_x_tsz_tsznorm',
    #              'cmb_x_cmb_lensnorm',
    #              'cmb+lens-lpf_cmb+lens-hpf',
    #              'cmb_x_cmb+lens_lensnorm',
    #              'cmb+lens_x_cmb_lensnorm']

    # for lo-tsz run
    corr_type = ['cmb+tsz-lpf_tsz-hpf',
                 'cmb_x_tsz_tsznorm']

    # diff_type = [(corr_type[0], corr_type[1]), # lensing bias
    #              (corr_type[3], corr_type[2])] # tsz mean bias
    # diff_type = [(corr_type[0], corr_type[1])] # lensing bias

    # for hi-tsz+lens run
    # calculating lens diff profiles was done on command line using calc_lens_diff_profile
    # diff_type = [(corr_type[0], corr_type[1]), # exact tsz mean bias for TI, approx for sgn
    #              (corr_type[3], corr_type[5])] # exact but noisy lensing bias for sgn

    # for lo-tsz run
    diff_type = [(corr_type[0], corr_type[1])] # exact tsz mean bias for TI, approx for sgn

    est_type = ['ti', 'sgn']

    for corr in corr_type:
        for est in est_type:
            calc_stacked_profile(corr, est)

    for diff in diff_type:
        for est in est_type:
            calc_diff_profile(diff[0], diff[1], est)

if __name__ == '__main__':
    main()
