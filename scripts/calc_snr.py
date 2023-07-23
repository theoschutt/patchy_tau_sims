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
        '--outpath',
        type=str,
        default='../output/snr/',
        help='where to save output'
    )

    args = parser.parse_args()

    return args

def calc_snr_per_object(datavector, cov, N=7676):
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
            print(gal_survey, n_gal_dict[gal_survey], '%.3f'%snr)
    return snr_expt_list

def print_snr_table():
    """generate table string to print to terminal"""
    raise NotImplementedError()

def write_snr_tables():
    """write latex formatted and more human-readable tables
    for SNR values
    """
    raise NotImplementedError()

def main():
    args = parse_args()
    signal_dv = np.genfromtxt(args.signal_file)[:,1]
    print(args.signal_file)
    advact_gal_dict = {
        'BOSS': 568776,
        'DESI LRGs Y1': 903000,
        'DESI LRGs+Ext. LRGs Y1': 3418500,
        'DESI LRGs Y5': 3010000,
        'DESI LRGs+Ext. LRGs Y5': 11395000,
        'unWISE': 50750600,
        'DES (approx)': 100000000,
        'LSST': 1624320000
    }
    for cov_file in args.noise_cov_files:
        print(cov_file)
        cov = np.genfromtxt(cov_file)
        snr_per_obj = calc_snr_per_object(signal_dv, cov)
        snr_expt = calc_snr_expt(snr_per_obj, advact_gal_dict)

if __name__ == '__main__':
    main()
