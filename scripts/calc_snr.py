#!/bin/bash/python

import os, sys
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description='Calc SNR for different experiments and N_galaxies')
    parser.add_argument(
        '--noise_cov_files',
        nargs='+',
        type=str,
        help='list of paths to noise covariance matrix txt files'
    )
    parser.add_argument(
        '--signal_file',
        type=str,
        help='path to signal radial profile txt data file'
    )
    parser.add_argument(
        '--unwise_only',
        default=False,
        action='store_const',
        const=True,
        help='Whether to only calculate expt SNRs for unWISE galaxy sample [default: False]'
    )
    parser.add_argument(
        '--outfile',
        type=str,
        help='full path and filename where to save output. Do not include file extension.'
    )

    args = parser.parse_args()

    return args

def write_log(args, outfile):
    # Write text file logging what command line args were used
    if not os.path.exists(os.path.split(outfile)[0]):
        os.makedirs(os.path.split(outfile)[0])
    log_fn =  outfile + '_args.log'
    print('Writing argument log file:', log_fn)
    arg_dict = vars(args)
    with open(log_fn, 'w') as f:
        for arg in arg_dict:
            f.write('%s: %s\n'%(str(arg), str(arg_dict[arg])))

def calc_snr_per_object(datavector, cov, N=7659):
    snr_sq = datavector @ np.linalg.solve(cov, datavector)
    snr_per_object = np.sqrt(snr_sq / N)
    print('snr_sq:', snr_sq)
    print('snr per obj (%i objs):'%N, snr_per_object)
    return snr_per_object

def calc_snr_expt(snr_per_object, n_gal_dict, print_snr=True):
    snr_expt_list = []
    for gal_survey in n_gal_dict:
        snr_expt_list.append(snr_per_object * np.sqrt(n_gal_dict[gal_survey]))

    if print_snr:
        for gal_survey, snr in zip(n_gal_dict, snr_expt_list):
            print(gal_survey, n_gal_dict[gal_survey], '--->', '%.3f'%snr)
    return snr_expt_list

def calc_snr_unwise_blue(signal_dv, noise_cov_files):
    """Calcs SNR for AdvACT, SPT-3G, SO, and S4 with the expected galaxy
    number from the unWISE blue sample. To calc n_gal just scales by
    survey area vs AdvACT."""
    # assumes order of cov files is act, spt, so, s4
    # advACT area: 12200 deg2
    # SPT-3G area: 1500 deg2
    # SO area: 15000 deg2 (40% f_sky)
    # S4 area: 20000 deg2 (50% f_sky)
    ngal_list = list(3.8e7, 4.7e6, 4.7e7, 6.2e7)
    snr_expt_list = []
    for cov_file, ngal_expt in zip(noise_cov_files, ngal_list):
        cov = np.genfromtxt(cov_file)
        snr_per_obj = calc_snr_per_object(signal_dv, cov)
        snr_expt_list.append(snr_per_obj * np.sqrt(ngal_expt))
    print("advact, spt, so, s4:", snr_expt_list)
    return snr_expt_list

def print_snr_table():
    """generate table string to print to terminal"""
    raise NotImplementedError()

def write_snr_tables(snr_expt_list, outfile):
    """write latex formatted and more human-readable tables
    for SNR values
    """
    snr_fn =  outfile + '.txt'
    print('Writing SNR file:', snr_fn)
    with open(log_fn, 'w') as f:
        for snr in snr_expt_list:
            f.write('%s: %s\n'%(str(arg), str(arg_dict[arg])))

def main():
    args = parse_args()

    write_log(args, args.outfile)

    signal_dv = np.genfromtxt(args.signal_file)[:,1]
    print(args.signal_file)
    advact_gal_dict = {
        'BOSS': 568776,
        'DESI LRGs Y1': 903000,
        'DESI LRGs+Ext. LRGs Y1': 3418500,
        'DESI LRGs Y5': 3010000,
        'DESI LRGs+Ext. LRGs Y5': 11395000,
        'unWISE': 50750600,
        'unWISE blue': 3.8e7,
        'unWISE green': 2.1e7,
        'DES (approx)': 100000000,
        'LSST': 1624320000
    }

    if args.unwise_only:
        snr_expt_list = calc_snr_unwise_blue(signal_dv, args.noise_cov_files)
        # write_snr_table(snr_expt_list)

    for cov_file in args.noise_cov_files:
        print(cov_file)
        cov = np.genfromtxt(cov_file)
        snr_per_obj = calc_snr_per_object(signal_dv, cov)
        snr_expt = calc_snr_expt(snr_per_obj, advact_gal_dict)

if __name__ == '__main__':
    main()
