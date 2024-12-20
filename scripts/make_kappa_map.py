"""make_kappa_map.py: makes a kappa map from a galaxy map of Dirac delta
functions and a given average redshift and halo mass for the galaxy sample.
Heavily borrowed from Abhi Maniyar's code.
"""
import os
import sys
import numpy as np
from pixell import enmap
# from astropy.cosmology import Planck18 as cosmo
from scipy import integrate, special
from make_noise_maps import make_map, save_map
sys.path.append("../../ThumbStack")
from universe import Universe
from flat_map import FlatMap

####################
# CMASS
# map_dir = ""
# z_gal = 0.55
# r_vir = 1.6 # [arcmin]
# M_gal = 2.e13

# unWISE blue
# map_dir = ""
# z_gal = 0.6
# M_gal = 1.4e13

# unWISE green
# map_dir = ""
# z_gal = 1.1
# M_gal = 1.3e13

####################
# CONSTANTS [in SI units]

Msun_to_Kg = 1.989e30
Mpc_to_m = 3.086e22
G_gravity = 6.67e-11
c_light = 3e8
z_cmb = 1100.
h = 0.7

# next four functions cribbed from Manu Schaan's HaloGen.profile.ProfNFW
def j0(x):
   """relevant for isotropic Fourier transform in 2d
   """
   return special.jn(0, x)

def rS_rhoS_c(U, m, z):
    """comoving scale radius for NFW profile
    in Mpc/h
    """
    Rvir = U.frvir(m, z)
    # m_nonlin = U.nonLinMass(0.)
    m_nonlin = 1e12
    #c = 10./(1.+z) * (m / m_nonlin)**(-0.2)   # from Takada & Jain 2002
    c = 9./(1.+z) * (m / m_nonlin)**(-0.13) # Takada & Jain 2003
    # scale radius
    rS = Rvir / c  # in Mpc/h
    # normalize the mass within rVir to be mVir
    rhoS = m / (4.*np.pi*rS**3)
    rhoS /= np.log(1.+c) - c/(1.+c)  # (Msun/h) / (Mpc/h)^3
    return rS, rhoS, c

def rho2d(R, U, m, z, trunc=1.):
    """comoving 2d density
    where the 3d NFW is truncated at trunc*rVir
    in (Msun/h) / (Mpc/h)^2
    R in Mpc/h
    x = R/RS = theta/thetaS, dimensionless
    """
    rS, rhoS, c = rS_rhoS_c(U, m, z)
    X = R / rS  # dimensionless
    # projected density in (Msun/h) / (Mpc/h)^2
    result = 2.*rhoS*rS
    # truncate the 3d profile at rvir
    if X>trunc*c:
       result *= 0.
    else:
       f = lambda x: x * 1./(x*(1.+x)**2) / np.sqrt(x**2 - X**2)
       # result *= integrate.quad(f, X, trunc*c, epsabs=0., epsrel=1.e-3)[0]
       result *= integrate.quad(f, X, trunc*c, epsrel=1.e-3)[0]
    return result

def rho2d_notrunc(R, U, m, z):
    """Compute the 2D comoving density with the truncation radius set to
    infinity. This has an analytic form.
    """
    # TODO
    pass

def rho2d_rescaled(R, U, m, z, trunc=1.):
    """Renormalize the nontruncated rho2d to contain the total mass
    within the given truncation radius.
    """
    # TODO: CHECK
    rVir = U.frvir(m, z)
    rS, rhoS, _ = rS_rhoS_c(U, m, z)
    # truncation radius over scale radius
    xMax = trunc * rVir/rS
    totMass = 4./3. * np.pi * rS**3 * rhoS
    totMass = xMax - np.log(1 + xMax)
    rho2d_rescaled = rho2d_notrunc(R, U, m, z) * totMass

    return rho2d_rescaled

# this needs to be a function of |ell| that can be passed to
# filterFourierIsotropic
def nfw2d_fourier(ell, U, m, z):
    """Compute the Fourier transform of the 2D projected NFW profile
    given the background cosmology, halo mass and object redshift.
    """
    chi = U.bg.comoving_distance(z)
    k = ell / chi
    integrand = lambda R: 4.*np.pi*R**2 * rho2d(R, U, m, z) * j0(k * R)
    # rho2d_fourier = integrate.quad(integrand, 0., np.inf, epsabs=0.,
    rho2d_fourier = integrate.quad(integrand, 0., np.inf,
                                   epsrel=1.e-2)[0]
    rho2d_fourier /= chi**2

    return rho2d_fourier

def get_rvir(U, m, z):
    """Returns virial radius in radians given mass and redshift.
    """
    chi_gal = U.bg.comoving_distance(z) # [Mpc/h/rad]
    Rvir = U.frvir(m, z) # [Mpc/h]
    Rvir_rad = Rvir / chi_gal # [rad]
    Rvir_amin = Rvir_rad * 180. * 60. / np.pi # [arcmin]
    print('Comoving Rvir [Mpc], Rvir [arcmin]:', Rvir, Rvir_amin)

    return Rvir_amin

def gauss2d_fourier(ell, sigma):
    """2D Gaussian in Fourier space given sigma in radians.
    """
    gauss2d_fourier =  np.exp(-0.5*ell**2* sigma**2)

    return gauss2d_fourier

def make_kappa_map(sample, profile='nfw', out_version='test', test=False):
    """Make kappa map from unWISE Dirac delta g map.
    sample = 'blue' or 'green'
    """
    U = Universe()
    unwise_dict = {
        'blue': {
            'z': 0.6,
            'm': 1.4e13, # [Msun]
            'ngal': 335070,
        },
        'green': {
            'z': 1.1,
            'm': 1.3e13, # [Msun]
            'ngal': 180025,
        }
    }
    m_gal = unwise_dict[sample]['m']
    m_gal /= 0.7
    z_gal = unwise_dict[sample]['z']

    # read map
    map_dir = f"../output/catalog/unwise_{sample}_10x10_nomask_correct_norm"
    map_fn = os.path.join(map_dir, "mock_count_dirac_car.fits")
    print('Reading counts map:', map_fn)
    dirac_enmap = enmap.read_map(map_fn) # [1/arcmin^2]
    print('Done.')

    # copy map to make base FlatMap
    baseMap = make_map()
    print('baseMap.sizeX:', baseMap.sizeX)
    print('baseMap.sizeY:', baseMap.sizeY)
    print('baseMap.nX:', baseMap.nX)
    print('baseMap.nY:', baseMap.nY)

    # multiply by map pixel areas [in arcmin^2] to get back Dirac delta map
    pixsizemap = dirac_enmap.pixsizemap() # [steradians]
    dirac_enmap *= pixsizemap * (180.*60./np.pi)**2 # [dimless]

    # approx. area of single pixel [str]
    pixarea = baseMap.sizeX/baseMap.nX * baseMap.sizeY/baseMap.nY

    print("np.sum(pixsizemap) [str]:", np.sum(pixsizemap))
    print("nX * nY * pixarea [str]:", pixarea * baseMap.nX * baseMap.nY)
    print("Should be equal to", 100. * (np.pi/180)**2)
    print("Sample ngal:", unwise_dict[sample]['ngal'])
    print("np.sum(dirac_enmap):", np.sum(dirac_enmap))

    # make kernel for 2D projected NFW profile
    if profile == 'nfw':
        filter_kernel = lambda l: nfw2d_fourier(l, U, m_gal, z_gal)

    # make kernel for 2D Gaussian
    elif profile == 'gauss':
        rvir_rad = get_rvir(U, m_gal, z_gal)
        sigma = rvir_rad/3 # [radians]
        filter_kernel = lambda l: gauss2d_fourier(l, sigma)
    else:
        print('Invalid filter type.')
        sys.exit()

    # Fourier transform delta map and multiply with profile kernel
    diracFourier = baseMap.fourier(dirac_enmap)
    print(f'Applying {profile} mass profile filter...')
    diracFourierFilt = baseMap.filterFourierIsotropic(filter_kernel,
                                                      dataFourier=diracFourier)
    print('Done.')

    # inverse Fourier transform to get the convolved real space map
    dirac_map_filt = baseMap.inverseFourier(diracFourierFilt)
    print("np.sum(dirac_map_filt):", np.sum(dirac_map_filt))
    
    print('Calculating kappa map...')
    # get comoving distance to galaxy mean z
    chi_gal = U.bg.comoving_distance(z_gal) # [Mpc/h/rad]
    chi_gal_m = chi_gal*Mpc_to_m

    # calc normalization such that
    # $\int d^2r \Sigma(r) = Ngal \times {\rm Mass}$
    # where $d^2r = {pix area (in steradians)} \times \chi_{gal}^2$
    M_gal_kg = m_gal * Msun_to_Kg
    d2r = pixarea * chi_gal_m**2
    Sigma_map = dirac_map_filt / d2r * M_gal_kg

    # check Sigma_map norm
    print("\int d^2r \Sigma(r) / M_gal_kg = Ngal? :",
        np.sum(Sigma_map) * pixarea * chi_gal_m**2 / M_gal_kg)

    # now make kappa map from Sigma map
    d_lens = U.bg.comoving_distance(z_gal)*Mpc_to_m
    d_source = U.bg.comoving_distance(z_cmb)*Mpc_to_m
    d_lenssource = d_source - d_lens
    a_gal = 1. / (1. + z_gal)
    constants = 4 * np.pi * G_gravity / c_light**2 / a_gal
    kappa_map = Sigma_map * constants * d_lens * d_lenssource / d_source
    print('Done.')

    # write map in enmap and flatmap format
    stem = os.path.basename(map_dir)
    out_fn = os.path.join(map_dir, f"{stem}_{profile}_kappa_map_{out_version}.fits")
    print('Writing kappa map:', out_fn)
    enmap.write_map(out_fn, kappa_map)


def main():
    profile = 'gauss'
    # profile = 'nfw'
    out_version = 'v1'

    for sample in ['blue', 'green']:
        print(f'Making kappa map for unWISE {sample} sample.')
        make_kappa_map(sample, profile, out_version)

if __name__ == '__main__':
    main()

