import sys
sys.path.append('../../ThumbStack')

import universe
from universe import *
import mass_conversion
from mass_conversion import *
import catalog
from catalog import *
import thumbstack
from thumbstack import *
import cmb
from cmb import *
# load catalogs
# 7642 entries over [200<RA<210, 10<DEC<20]
catpath = '/home/theo/Documents/research/CMB/patchy_tau_sims/output/catalog/grid_10x10_10x10src_periodic/catalog.txt'
# catpath = '/home/theo/Documents/research/CMB/patchy_tau_sims/output/catalog/cmass_m_10x10_v2/catalog.txt'
# randcatpath = '/home/theo/Documents/research/CMB/patchy_tau_sims/output/catalog/cmass_m_10x10_randradec_v2/catalog.txt'
catname = 'grid_10x10_10x10src_periodic'
# catname = 'cmass_m_10x10_v2_fwhm5'
# randname = '%s_randradec'%catname
cattype = 'radec'
#cattype = 'cmass'
# 10k entries, CMASS positions over [200<RA<211.5, 10<DEC<21.4]
# catpath = '/home/theo/Documents/research/CMB/tau_sims/output/catalog/cmass_m_10kgal/cmass_m_10kgal_catalog.txt'
# 10k entries, RADEC randomized over [200<RA<210, 10<DEC<20]
# randcatpath = '/home/theo/Documents/research/CMB/tau_sims/output/catalog/cmass_m_10kgal_randradec/catalog.txt'
# catname='cmass_m_10kgal_sig5'
# randname = '%s_randradec'%catname


# set up cosmology and mass conversion for catalog
u = Universe() # MarianaUniv() had a weird bug not inheriting Universe class attributes
massConversion = MassConversionKravtsov14()
# 10k entries, CMASS positions over [200<RA<210, 10<DEC<20]
# cmass10x10 = Catalog(
#     u,
#     massConversion,
#     name='cmass_m_10x10_randradec_v3',
#     pathInCatalog=catpath,
#     save=True
# )

cmass10x10 = Catalog(
    u,
    massConversion,
    name=catname,
    pathInCatalog=catpath,
    save=True,
    catType=cattype,
    workDir='..'
)
"""
cmass10x10rand = Catalog(
    u,
    massConversion,
    name=randname,
    pathInCatalog=randcatpath,
    save=True
)
"""
# calc mean tau over all galaxies
# mean_tau = np.mean(cmass10x10.integratedTau)
# print('Mean tau: ', mean_tau)

# Generate empty square map, then make a mock map
# Generate an empty square map with RA in [200., 210.] and DEC in [10., 20.]
# convention for defining box corners is [[dec_min, ra_max],[dec_max, ra_min]]
side_length = 10.
# side_length = 11.4
ra_min = 200.
ra_max = ra_min + side_length
dec_min = 10.
dec_max = dec_min + side_length

box = np.array([[dec_min, ra_max],[dec_max, ra_min]]) * utils.degree
resArcmin = 0.5 # map pixel size [arcmin]
shape,wcs = enmap.geometry(pos=box, res=resArcmin * utils.arcmin, proj='car')

# create a mask that keeps the whole area
boxMask = enmap.ones(shape, wcs=wcs)

# automatically writes maps
fwhm = 5.
sigma = fwhm / np.sqrt(8.*np.log(2))
test = False
cmass10x10.generateMockMaps(boxMask, sigma=sigma, test=test)
#cmass10x10rand.generateMockMaps(boxMask, sigma=sigma, test=test)
