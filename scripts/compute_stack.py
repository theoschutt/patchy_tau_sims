"""compute_stack.py: Compute stacked profile using arbitrary large- and small scale aperture photometry files."""

import os, sys
import numpy as np

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Compute stacked profile using arbitrary large- and small-scale aperture photometry files.')
    parser.add_argument(
        '--t_small',
        help='Full path to the small scale filtMap photometry file')
    parser.add_argument(
        '--t_large',
        help='Full path to the large scale filtNoiseStdDev photometry file (for weights, i.e. numerator of stack)')
    parser.add_argument(
        '--t_large_norm',
        help='Full path to large scale filtNoiseStdDev photometry file (for norm, i.e. denominator of stack)')
    parser.add_argument(
        '--balance_signs',
        action='store_const',
        const=True,
        default=False,
        help='Whether to apply sign balancing to the stacked profile')
    parser.add_argument(
        '--t_large_min',
        default=None,
        type=float,
        help='Minimum |T_large| to cut from the stack. Uses t_large_norm values')
    parser.add_argument(
        '--tag',
        help='unique string for directory/filename of output files')
    parser.add_argument(
        '--outpath',
        default='/home/theo/Documents/research/CMB/patchy_tau_sims/output/thumbstack',
        help='Non-default path to directory to save output files. subdir will be written named --tag.')

    args = parser.parse_args()

    return args

def compute_stack(t_small, t_large, t_large_norm, est, t_large_min, balance_signs):
    if t_large_min is not None:
        mask = np.abs(t_large_norm[:,0]) >= t_large_min
        print(len(t_small), sum(mask[:]), mask.shape)
        t_small = t_small[mask,:]
        t_large = t_large[mask,:]
        t_large_norm = t_large_norm[mask,:]
        print(t_small.shape)
    if est == 'ti':
        weights = -t_large
        norm = 1./np.sum(t_large_norm**2, axis=0)
    if est == 'sgn':
        weights = -np.sign(t_large)
        norm = 1./np.sum(np.abs(t_large_norm), axis=0)
    print(weights)

    if balance_signs:
        diff = np.sum(t_large[:,0] > 0) - np.sum(t_large[:,0] < 0)
        print('EqualSignedWeights diff:', diff)
        if diff == 0:
            pass # do nothing. The number of + and - weights are already equal
        else:
            # get the indices of the excess signed weights
            if diff > 0:
               idxExcess = np.where(t_large_norm[:,0] > 0)[0]
            else:
               idxExcess = np.where(t_large_norm[:,0] < 0)[0]
            # randomly select N=diff samples of the excess signed weights
            # TODO: We probably want this selection to be reproducible. How to pick seed?
            #       Do we want the same samples removed for all bootstrap iterations?
            rng = np.random.default_rng()
            idxToRemove = rng.choice(idxExcess, np.abs(diff), replace=False)
            # now set the weights of those randoms to zero so
            # they are effectively removed from the stack
            weights[idxToRemove,:] = 0
            # now update the norm
            if est=='ti':
               norm = 1./np.sum(np.delete(t_large_norm, idxToRemove, axis=0)**2, axis=0)
            elif est=='sgn':
               norm = 1./np.sum(np.abs(np.delete(t_large_norm, idxToRemove,
                                                 axis=0)), axis=0)

    print(est, weights, np.sum(weights), norm, norm.shape)
    s2Full = np.var(t_small, axis=0)
    stack = norm * np.sum(t_small * weights, axis=0)
    sStack = norm * np.sqrt(np.sum(s2Full * weights**2, axis=0))

    print('s2Full, stack, sStack:', s2Full, stack, sStack)
    print('s2Full, stack, sStack shape:', s2Full.shape, stack.shape, sStack.shape)
    return stack, sStack

def load_photometry(t_small_fn, t_large_fn, t_large_norm_fn):
    t_small = np.genfromtxt(t_small_fn)
    t_large = np.genfromtxt(t_large_fn)
    t_large_norm = np.genfromtxt(t_large_norm_fn)

    return t_small, t_large, t_large_norm

def write_log(args, pathOut, tag):
    # Write text file logging what command line args were used
    log_fn = os.path.join(pathOut, 'args_%s.log'%tag)
    print('Writing argument log file:', log_fn)
    if not isinstance(args, dict):
        arg_dict = vars(args)
    else:
        arg_dict = args
    with open(log_fn, 'w') as f:
        for arg in arg_dict:
            f.write('%s: %s\n'%(str(arg), str(arg_dict[arg])))

def save_stack(stack, sStack, est, pathOut):
    data = np.zeros((9, 3))
    rApArcmin = np.linspace(1., 6., 9)
    data[:,0] = rApArcmin
    data[:,1] = stack
    data[:,2] = sStack
    fn = os.path.join(pathOut, "tauring_tau_%s_uniformweight_measured.txt"%est)
    np.savetxt(fn, data)

def main():
    args = parse_args()

    pathOut = os.path.join(args.outpath, args.tag)
    if not os.path.exists(pathOut):
        os.makedirs(pathOut)

    write_log(args, pathOut, args.tag)

    t_small, t_large, t_large_norm = load_photometry(
        args.t_small, args.t_large, args.t_large_norm)

    for est in ['ti', 'sgn']:
        stack, sStack = compute_stack(
            t_small, t_large, t_large_norm, est, args.t_large_min, args.balance_signs)
        save_stack(stack, sStack, est, pathOut)

if __name__ == '__main__':
    main()
