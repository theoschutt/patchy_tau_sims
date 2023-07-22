#!/bin/bash/python

import os, sys
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Plot radial dependent curves (signal, noise, etc')
    parser.add_argument(
        '--plot_type',
        type=str,
        default='',
        help='which type of plot to plot'
    )
    parser.add_argument(
        '--data_files',
        nargs='+',
        default=None,
        type=str,
        help='list of curve data filenames'
    )

    args = parser.parse_args()

    return args

def plot_signal_curve(data_txt_list):
    pass

def plot_noise_curves(data_txt_list):
    pass

def plot_noise_ratios(data_txt_list):
    pass

def main():
    args = parse_args()

    if args.plot_type == 'signal':
        plot_signal_curve(args.data_files)

    if args.plot_type == 'noise':
        plot_noise_curves(args.data_files)

    if args.plot_type == 'noise_ratio':
        plot_noise_ratios(args.data_files)

    else:
        raise ValueError('Invalid plot type specified:', args.plot_type)
