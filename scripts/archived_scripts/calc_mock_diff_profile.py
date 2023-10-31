#!/bin/bash/python
import os
import pickle
import numpy as np

mocks_dir = ''
dict_fn = os.path.join(mocks_dir, '')
with open(dict_fn, 'rb') as f:
    collated_dict = pickle.load(f)

# collated_dict has 8 sets of thumbstack outputs
# *-lpf_*-hpf x (ti or sgn)
corr_type = ['cmb-lpf_cmb-hpf',
             'cmb-lpf_cmb+lens-hpf',
             'cmb-lpf_tsz-hpf',
             'cmb+tsz-lpf_tsz-hpf']

# will want per-seed differences: CMBxCMB - CMBxCMB+lens, CMB+tszxtsz - CMBxtsz
diff_pairs = [(corr_type[0], corr_type[1]), # lensing bias
         (corr_type[3], corr_type[2])] # tSZ mean bias
diff_names = ['diff_lens', 'diff_tsz']
est_type = ['ti', 'sgn']

# to make output file same as thumbstack's include radial bin value as first col
r = np.linspace(1., 6., 9)

stack_list = []
sStack_list = []
for pair, name in zip(diff_pairs, diff_names):
    for est in est_type:
        for seed in range(len(collated_dict)):
            diff_stack = collated_dict[seed][pair[0]][est]['stack'] - collated_dict[seed][pair[1]][est]['stack']
            stack_list.append(diff_stack)

        stack_arr = np.vstack(stack_list)

        profile = np.average(stack_arr, axis=0)
        sProfile = np.std(stack_arr, axis=0) / np.sqrt(len(collated_dict))

        profile_info = np.vstack([r, profile, sProfile]).T

        np.savetxt('hi-tsz+lens_100mocks_%s_%s_measured.txt'%(name, est), profile_info, fmt='%.5e')
