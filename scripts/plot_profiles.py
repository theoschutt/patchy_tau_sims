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
        '--title',
        type=str,
        default='',
        help='Plot title'
    )
    parser.add_argument(
        '--labels',
        nargs='+',
        type=str,
        default=['TI', 'Sgn'],
        help='list of plot labels'
    )
    parser.add_argument(
        '--ls',
        nargs='+',
        type=str,
        default=['-', '-'],
        help="list of linestyles to use, e.g. ['-', '--', '-.', ':']"
    )
    parser.add_argument(
        '--data_files',
        nargs='+',
        default=None,
        type=str,
        help='list of curve data filenames'
    )
    parser.add_argument(
        '--noise_files',
        nargs='+',
        default=None,
        type=str,
        help='list of noise data filenames'
    )
    parser.add_argument(
        '--outfile',
        type=str,
        help='path+filename for output. Do not include file extension.'
    )

    args = parser.parse_args()

    return args

def plot_signal_curve(data_txt_list, fn, title, labels, lss):
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
    colors = ['#17685f', '#b33210', 'k', 'r', 'c']
    # lss = ['-', '-', '--', '--', '--']

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

        ax.errorbar(r+i*0.01, t, t_err, label=labels[i], c=colors[i], ls=lss[i])
        ax.set_ylabel(r'$\left<\tau\right>$')
        ax.set_xlabel('Radius [arcmin]')
        ax.legend(framealpha=1., loc='upper right')
    plt.savefig(fn+'.pdf', metadata=dict(Subject=str(data_txt_list)), dpi=300, bbox_inches='tight')
    plt.savefig(fn+'.png', dpi=100, bbox_inches='tight')

def plot_signal_with_noise_band(data_txt_list, noise_txt_list, fn, title, labels, lss):
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
    darker_colors = ['#003931', '#730000']
    colors = ['#17685f', '#b33210']
    #colors = ['#ffb348', '#946e6d']
    # colors = ['#17685f', '#cb4d4b', 'b', 'r', 'c']
    # lss = ['-', '-', '--', '--', '--']
    formats = ['x', 'x']
    formats = ['o', 'D']

    rad2arcmin = 180. * 60. / np.pi

    # tau_dat = np.genfromtxt('/home/theo/Documents/research/CMB/patchy_tau_sims/output/thumbstack/cmass_m_10s10_tau-fwhm_renorm2.45e-4_v3_c-1.0_willfilt/tauring_tau_ti_uniformweight_measured.txt')
    tau_dat = np.genfromtxt('/home/theo/Documents/research/CMB/patchy_tau_sims/output/thumbstack/cmass_m_10x10_fwhm5_v5/tauring_tau_ti_uniformweight_measured.txt')
    r = tau_dat[:,0]
    plot_r = r.copy()
    plot_r[0] -= 0.1
    plot_r[-1] += 0.1
    tau = tau_dat[:,1] * rad2arcmin**2
    print(tau)
    # divide by ring area to get mean tau
    areas = r.copy()
    for j in range(len(r)):
        if j == 0:
            areas[j] = np.pi * r[j]**2
        else:
            areas[j] = np.pi * (r[j]**2 - r[j-1]**2)
    # these are the actual ring areas used to integrate the
    # aperture filters
    # pre bugfix
    # areas = [
    #     2.25,
    #     6.75,
    #     8.,
    #     8.,
    #     11.5,
    #     18.5,
    #     18.,
    #     16.,
    #     20.5
    # ]
    # post bugfix
#     areas = [
#         2.5,
#         6.75,
#         8.,
#         8.,
#         11.5,
#         18.5,
#         18.,
#         16.,
#         20.5
#     ]
    tau /= areas
    print(tau)

    for i, (txt, noise) in enumerate(zip(data_txt_list, noise_txt_list)):
        print(txt)
        dat = np.genfromtxt(txt)
        r = dat[:,0]
        t = dat[:,1] * rad2arcmin**2 / areas
        t_err = dat[:,2] * rad2arcmin**2 / areas

        noise_dat = np.genfromtxt(noise)
        noise_err = noise_dat[:,2] * rad2arcmin**2 / areas

        ax.fill_between(plot_r, -noise_err/10, noise_err/10, color=colors[i], alpha=0.2, label="%s stat err"%labels[i], zorder=12-i)
        #ax.fill_between(plot_r, t-t_err, t+t_err, color=darker_colors[i], alpha=0.4, label="%s bias"%labels[i], zorder=8-i)

        ax.errorbar(r+i*0.05 - 0.025, t, t_err, label="%s bias"%labels[i], fmt=formats[i], c=colors[i], markeredgecolor=colors[i], lw=1.5, markersize=6, markeredgewidth=1.5, markerfacecolor='white', zorder=1)
        # ax.errorbar(r, t, t_err, label="%s bias"%labels[i], c=colors[i], ls=lss[i], alpha=0.6)
        # ax.fill_between(r+i*0.03, -noise_err, noise_err, color=colors[i], alpha=0.2)
        ax.set_ylabel(r'$\left<\tau\right>$')
        ax.set_xlabel('Radius [arcmin]')
    ax.axhline(0, c='k', alpha=0.8, lw=1., zorder=30, ls=':')
    ax.plot(r, tau, c='k', label=r'$\tau$ signal', ls='-', zorder=2)
    # plot legend in pleasing order
    handles, labels = ax.get_legend_handles_labels()
    order = [3,4,2,0,1]
    #order = [1,3,4,0,2]
    ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
              frameon=False, framealpha=1., ncol=2, loc='upper right')
    # ax.legend(frameon=False, framealpha=1., ncol=2, loc='upper right')

    ax.grid(False)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-4, -4))
    ax.set_xlim((plot_r[0], plot_r[-1]))
    ax.set_ylim((-0.5e-4, 1e-4))
    plt.savefig(fn+'.pdf', metadata=dict(Subject=str(data_txt_list)), dpi=300, bbox_inches='tight')
    plt.savefig(fn+'.png', dpi=100, bbox_inches='tight')

def plot_noise_curves(data_txt_list, fn):
    pass

def plot_signal_ratios(data_txt_list, fn, title, labels, lss):
    fig, ax = plt.subplots(1,1, figsize=(10,5))

    # labels = [r'$\tau$ TI', 'tSZ TI', 'kSZ TI', 'CIB TI', 'AdvACT spectrum']
    # labels = [r'$\tau$ Sgn', 'tSZ Sgn', 'kSZ Sgn', 'CIB Sgn', 'AdvACT spectrum']
    colors = ['k', 'r', 'b', 'orange', 'c']
    # lss = ['-', '-', '-', '-', '--']

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
        plot_signal_curve(args.data_files, args.outfile, args.title, args.labels, args.ls)

    if args.plot_type == 'signal-band':
        plot_signal_with_noise_band(args.data_files, args.noise_files, args.outfile, args.title, args.labels, args.ls)

    elif args.plot_type == 'noise':
        plot_noise_curves(args.data_files, args.outfile)

    elif args.plot_type == 'tau_ratio':
        plot_signal_ratios(args.data_files, args.outfile, args.title, args.labels, args.ls)

    else:
        raise ValueError('Invalid plot type specified:', args.plot_type)

if __name__ == '__main__':
    main()
