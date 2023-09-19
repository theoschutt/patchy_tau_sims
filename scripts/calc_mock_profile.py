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
est_type = ['ti', 'sgn']

# to make output file same as thumbstack's include radial bin value as first col
r = np.linspace(1., 6., 9)

stack_list = []
sStack_list = []
for corr in corr_type:
    for est in est_type:
        for seed in range(len(collated_dict)):
            stack_list.append(collated_dict[seed][corr][est]['stack'])
            sStack_list.append(collated_dict[seed][corr][est]['sStack'])

        stack_arr = np.vstack(stack_list)
        sStack_arr = np.vstack(sStack_list)

        stack = np.average(stack_arr, axis=0)
        sStack = np.sqrt(len(collated_dict)) * np.average(sStack_arr, axis=0)

        stack_info = np.vstack([r, stack, sStack]).T

        np.savetxt('hi-tsz+lens_100mocks_%s_%s_measured.txt'%(corr, est), stack_info, fmt='%.5e')
