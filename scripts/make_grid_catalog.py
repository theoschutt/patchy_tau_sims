#!/bin/bash
import os
import numpy as np

n_profiles = 10
odir = '/home/theo/Documents/research/CMB/patchy_tau_sims/output/catalog/'
catdir = 'grid_10x10_%ix%isrc_periodic'%(n_profiles, n_profiles)
fn = os.path.join(odir, catdir, 'catalog.txt')
if not os.path.exists(os.path.join(odir, catdir)):
    os.makedirs(os.path.join(odir, catdir))

side_length = 10.
ramin = 200.
ramax = ramin + side_length
decmin = 10.
decmax = decmin + side_length
# want edges to be half the separation of other sources to maintain period BCs
edge_buffer = side_length/n_profiles/2.

ra_list = np.linspace(ramin+edge_buffer, ramax-edge_buffer, n_profiles)
dec_list = np.linspace(decmin+edge_buffer, decmax-edge_buffer, n_profiles)

print(ra_list)
print(dec_list)

ra_arr, dec_arr = np.meshgrid(ra_list, dec_list)

all_ras = ra_arr.flatten()
all_decs = dec_arr.flatten()

with open(fn, 'w') as f:
    for r,d in zip(all_ras, all_decs):
        f.write("%f %f\n"%(r,d))
