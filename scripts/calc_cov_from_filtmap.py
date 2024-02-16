"""calc_cov_from_filtmap.py: afterburner calculation of the covariance matrix
from saved aperture photometry filter data. Just reset `outdir` in main()."""
import os
import numpy as np

def get_w_and_norm(t_large, est):
    """calculate weights and norm for given estimator."""
    if est=='ti':
        weights = -t_large
        norm = 1./np.sum(t_large**2, axis=0)
    # tau: weighted by sign of large scale temperature
    elif est=='sgn':
        weights = -np.sign(t_large)
        norm = 1./np.sum(np.abs(t_large), axis=0)

    return weights, norm

def equalize_signs(t_large, weights, est, i_bootstrap):
    """equalize positive and negative signed weights in stack."""
    # get the difference N_pos - N_neg.
    # Some weights might be exactly zero, so count both explicitly.
    # For tau estimators, weights for a stamp's apertures are all equal,
    # so only use the first aperture's weight for counting.
    diff = np.sum(weights[:,0] > 0) - np.sum(weights[:,0] < 0)
    # print('EqualSignedWeights diff:', diff)
    if diff == 0:
        pass # do nothing. The number of + and - weights are already equal
    else:
        # get the indices of the excess signed weights
        if diff > 0:
            idxExcess = np.where(weights[:,0] > 0)[0]
        else:
            idxExcess = np.where(weights[:,0] < 0)[0]
        # randomly select N=diff samples of the excess signed weights
        rng = np.random.default_rng(i_bootstrap)
        idxToRemove = rng.choice(idxExcess, np.abs(diff), replace=False)
        # now set the weights of those randoms to zero so
        # they are effectively removed from the stack
        weights[idxToRemove,:] = 0
        # now update the norm
        if est=='ti':
            norm = 1./np.sum(np.delete(t_large, idxToRemove, axis=0)**2, axis=0)
        elif est=='sgn':
            norm = 1./np.sum(np.abs(np.delete(t_large, idxToRemove,
                                              axis=0)), axis=0)

        # print('equalSignedWeights \ndiff:', diff)
        # print('N_removed:', len(idxToRemove))
        # print('New weights:', weights)
        # print('New norm:', norm)

    return weights, norm

def do_bootstrap(i_bootstrap, t_small, t_large):
    """do the bootstrap resampling and return the samples for ti and sgn."""

    print(f'Performing bootstrap sampling with {n_bootstrap} iterations.')
    # will have shape (n_bootstrap, n_apertures)
#     ti_boot_resamples = np.empty((n_bootstrap, t_small.shape[1]))
#     sgn_boot_resamples = np.empty((n_bootstrap, t_small.shape[1]))
    n_obj = len(t_small)
    all_idx = np.arange(n_obj)

#    for i in range(n_bootstrap):
   if i_bootstrap%1000 == 0:
       print('Bootstrap index:', i_bootstrap)
   rng = np.random.default_rng(i_bootstrap)
   idx = rng.choice(all_idx, size=n_obj, replace=True)

   ts = t_small[idx,:]
   tl = t_large[idx,:]
   boot_resamples = []
   for est, boot_resamples in zip(['ti', 'sgn'], [ti_boot_resamples, sgn_boot_resamples]):
       weights, norm = get_w_and_norm(tl, est)
       weights, norm = equalize_signs(tl, weights, est, i_bootstrap)

       stack = norm * np.sum(ts * weights, axis=0)
       # sStack = norm * np.sqrt(np.sum(np.var(ts, axis=0) * weights**2, axis=0))
       # print(stack)
       # boot_resamples[i_bootstrap,:] = stack
       boot_resamples.append(stack)

    # return ti_boot_resamples, sgn_boot_resamples
    return boot_resamples

# def do_bootstrap(n_bootstrap, t_small, t_large):
#     """do the bootstrap resampling and return the samples for ti and sgn."""
# 
#     print(f'Performing bootstrap sampling with {n_bootstrap} iterations.')
#     # will have shape (n_bootstrap, n_apertures)
#     ti_boot_resamples = np.empty((n_bootstrap, t_small.shape[1]))
#     sgn_boot_resamples = np.empty((n_bootstrap, t_small.shape[1]))
#     n_obj = len(t_small)
#     all_idx = np.arange(n_obj)
# 
#     for i in range(n_bootstrap):
#         if i%1000 == 0:
#             print('Bootstrap index:', i)
#         rng = np.random.default_rng(i)
#         idx = rng.choice(all_idx, size=n_obj, replace=True)
# 
#         ts = t_small[idx,:]
#         tl = t_large[idx,:]
# 
#         for est, boot_resamples in zip(['ti', 'sgn'], [ti_boot_resamples, sgn_boot_resamples]):
#             weights, norm = get_w_and_norm(tl, est)
#             weights, norm = equalize_signs(tl, weights, est, i)
# 
#             stack = norm * np.sum(ts * weights, axis=0)
#             # sStack = norm * np.sqrt(np.sum(np.var(ts, axis=0) * weights**2, axis=0))
#             # print(stack)
#             boot_resamples[i,:] = stack
# 
#     return ti_boot_resamples, sgn_boot_resamples

def calc_and_save_cov(boot_resamples, outpath):
    """what it says on the tin."""
    print('Calculating covariance matrix...')
    stack_cov = np.cov(boot_resamples, rowvar=False)
    print(f'Saving as {outpath}')
    np.savetxt(outpath, stack_cov)

    return stack_cov

def main():
    outdir = ('/home/theo/Documents/research/CMB/patchy_tau_sims/output/thumbstack/'
        'cmass_m_10x10_spt_ilc_beam1.6_final')
    print(f'Loading data from {outdir}')
    t_small = np.genfromtxt(os.path.join(outdir, 'tauring_filtmap.txt'))
    t_large = np.genfromtxt(os.path.join(outdir, 'tauring_filtnoisestddev.txt'))
    # TODO: add sharedmem stuff
    n_bootstrap = 10000
    ti_boot_resamples, sgn_boot_resamples = do_bootstrap(n_bootstrap, t_small, t_large)

    for est, boot_resamples in zip(['ti', 'sgn'], [ti_boot_resamples, sgn_boot_resamples]):
        outpath = os.path.join(outdir, f'cov_tauring_tau_{est}_uniformweight_bootstrap_manual.txt')
        calc_and_save_cov(boot_resamples, outpath)

if __name__ == '__main__':
    main()
