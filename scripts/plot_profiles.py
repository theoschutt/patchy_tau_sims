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
        '--plot_title',
        type=str,
        default='',
        help='which type of plot to plot'
    )
    parser.add_argument(
        '--plot_labels',
        nargs='+',
        type=str,
        help='list of plot labels'
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

def plot_signal_curve(data_txt_list, fn, title, labels):
    fig, ax = plt.subplots(1,1, figsize=(10,5))
    # labels = ['tSZ TI', 'tSZ Sgn']
    # labels = ['kSZ TI', 'kSZ Sgn']
    # labels = ['CIB TI', 'CIB Sgn']
    # labels = ['TI', 'Sgn', r'TI, +=-', r'Sgn, +=-']
    # labels = ['TI', 'Sgn', r'TI, |T_L| > 2\mu$K', r'Sgn, |T_L| > 2\mu$K']
    ax.set_title(title)

#     labels = [r'$\tau$ TI', 'tSZ TI', 'kSZ TI', 'CIB TI', 'AdvACT spectrum']
    # labels = [r'$\tau$ Sgn', 'tSZ Sgn', 'kSZ Sgn', 'CIB Sgn', 'AdvACT spectrum']
    # colors = ['k', 'r', 'b', 'orange', 'c']
    colors = ['k', 'r', 'k', 'r', 'c']
    lss = ['-', '-', '--', '--', '--']

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

        ax.errorbar(r+i*0.03, t, t_err, label=labels[i], c=colors[i], ls=lss[i])
        ax.set_ylabel(r'$\left<\tau\right>$')
        ax.set_xlabel('Radius [arcmin]')
        ax.legend(framealpha=1.)
    plt.savefig(fn+'.pdf', metadata=dict(Subject=str(data_txt_list)), dpi=300, bbox_inches='tight')
    plt.savefig(fn+'.png', dpi=100, bbox_inches='tight')


def plot_noise_curves(data_txt_list, fn):
    pass

def plot_signal_ratios(data_txt_list, fn, title, labels):
    fig, ax = plt.subplots(1,1, figsize=(10,5))

    # labels = [r'$\tau$ TI', 'tSZ TI', 'kSZ TI', 'CIB TI', 'AdvACT spectrum']
    # labels = [r'$\tau$ Sgn', 'tSZ Sgn', 'kSZ Sgn', 'CIB Sgn', 'AdvACT spectrum']
    colors = ['k', 'r', 'b', 'orange', 'c']
    lss = ['-', '-', '-', '-', '--']

    rad2arcmin = 180. * 60. / np.pi

    t_dat = np.genfromtxt(data_txt_list[0])
    r = t_dat[:,0]
    tau = t_dat[:,1] * rad2arcmin**2
    tau_err = t_dat[:,2] * rad2arcmin**2

    # divide by ring area to get mean tau
    areas = r.copy()
    for j in range(len(r)):
        if j == 0:
            areas[j] = np.pi * r[j]**2
        else:
            areas[j] = np.pi * (r[j]**2 - r[j-1]**2)
    tau /= areas
    tau_err /= areas

    for i, txt in enumerate(data_txt_list):
        dat = np.genfromtxt(txt)
        t = dat[:,1] * rad2arcmin**2
        t_err = dat[:,2] * rad2arcmin**2

        t /= areas
        t_err /= areas

        ax.plot(r, t/tau, label=labels[i], c=colors[i], ls=lss[i])
        ax.set_ylabel(r'$\left<\tau\right>$/$\left<\tau_{\rm true}\right>$')
        ax.set_xlabel('Radius [arcmin]')
        ax.legend(framealpha=1.)
    plt.savefig(fn+'.pdf', metadata=dict(Subject=str(data_txt_list)), dpi=300, bbox_inches='tight')
    plt.savefig(fn+'.png', dpi=100, bbox_inches='tight')

def write_log(args, outfile):
    # Write text file logging what command line args were used
    log_fn =  outfile + '_args.log'
    print('Writing argument log file:', log_fn)
    arg_dict = vars(args)
    with open(log_fn, 'w') as f:
        for arg in arg_dict:
            f.write('%s: %s\n'%(str(arg), str(arg_dict[arg])))


def main():
    args = parse_args()

    if not os.path.exists(os.path.split(args.outfile)[0]):
        os.makedirs(os.path.split(args.outfile)[0])

    write_log(args, args.outfile)

    if args.plot_type == 'signal':
        plot_signal_curve(args.data_files, args.outfile, args.plot_title, args.plot_labels)

    elif args.plot_type == 'noise':
        plot_noise_curves(args.data_files, args.outfile)

    elif args.plot_type == 'tau_ratio':
        plot_signal_ratios(args.data_files, args.outfile, args.plot_title, args.plot_labels)

    else:
        raise ValueError('Invalid plot type specified:', args.plot_type)

if __name__ == '__main__':
    main()
