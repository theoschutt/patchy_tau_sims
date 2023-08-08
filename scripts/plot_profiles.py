#!/bin/bash/python

import os, sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('default.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Plot radial dependent curves (signal, noise, etc')
    parser.add_argument(
        '--plot_type',
        type=str,
        default='signal',
        help='which type of plot to plot'
    )
    parser.add_argument(
        '--data_files',
        nargs='+',
        default=None,
        type=str,
        help='list of curve data filenames'
    )
    parser.add_argument(
        '--outfile',
        type=str,
        help='path+filename for output. Do not include file extension.'
    )

    args = parser.parse_args()

    return args

def plot_signal_curve(data_txt_list, fn):
    fig, ax = plt.subplots(1,1, figsize=(10,5))

    labels = [r'$\tau$', 'tSZ TI', 'tSZ Sgn', 'kSZ', 'CIB', 'AdvACT spectrum']
    colors = ['k', 'r', 'b', 'orange', 'c']
    lss = ['-', '-', '-', '-', '--']

    rad2arcmin = 180. * 60. / np.pi

    for i, txt in enumerate(data_txt_list):
        print(txt)
        dat = np.genfromtxt(txt)
        r = dat[:,0]
        t = dat[:,1] * rad2arcmin**2
        t_err = dat[:,2] * rad2arcmin**2

        # divide by ring area to get mean tau
        areas = r.copy()
        for j in range(len(r)):
            if j == 0:
                areas[j] = np.pi * r[j]**2
            else:
                areas[j] = np.pi * (r[j]**2 - r[j-1]**2)
        print(areas)
        t /= areas
        t_err /= areas

        ax.errorbar(r + i*0.03, t, t_err, label=labels[i], c=colors[i], ls=lss[i])
        ax.set_ylabel(r'Mean $\tau$')
        ax.set_xlabel('Radius [arcmin]')
        ax.legend()
    plt.savefig(fn+'.pdf', metadata=dict(Subject=data_txt_list), dpi=300, bbox_inches='tight')
    plt.savefig(fn+'.png', dpi=100, bbox_inches='tight')


def plot_noise_curves(data_txt_list):
    pass

def plot_noise_ratios(data_txt_list):
    pass

def main():
    args = parse_args()

    if args.plot_type == 'signal':
        plot_signal_curve(args.data_files, args.outfile)

    elif args.plot_type == 'noise':
        plot_noise_curves(args.data_files)

    elif args.plot_type == 'noise_ratio':
        plot_noise_ratios(args.data_files)

    else:
        raise ValueError('Invalid plot type specified:', args.plot_type)

if __name__ == '__main__':
    main()
