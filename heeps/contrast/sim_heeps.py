"""
ScopeSim processing of a HCI (RAVC) data cube produced by HEEPS

The cube contains the temporal evolution of the coronagraphic PSF using
the METIS RAVC in the L band. ScopeSim uses the cube as the input source
and produced detector images with background and detector noise components.
"""
import os
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

import scopesim as sim
#sim.download_package(['instruments/METIS', 'telescopes/ELT', 'locations/Armazones'])

from synphot import SourceSpectrum

# Change if you want to use another cube
#FNAME = "onaxis_PSF_L_RAVC.fits"

WORKDIR = '/home/tom/METIS/HEEPS/heeps_metis/'
OUTFILE = WORKDIR+"heeps_scopesimed.fits"

# Change FLUX_SCALE to vary the star's brightness relative to Vega
FLUX_SCALE = 1.

# This sets the path to the METIS configuration. Change the path to
# match your setup.
#PKGS_DIR = os.path.abspath(
#    os.path.join(
#        os.path.dirname("/home/tom/METIS")))
#sim.rc.__config__['!SIM.file.local_packages_path'] = PKGS_DIR


def sim_heeps(psf_ON, psf_OFF, header, **conf):
    """Do the simulation"""

    # Prepare the source
    cube = psf_ON #fits.open(psf_ON)
    imwcs = WCS(header).sub([1, 2])
    naxis1, naxis2, naxis3 = (header['NAXIS1'],
                              header['NAXIS2'],
                              header['NAXIS3'])

    # Set up the instrument
    cmd = sim.UserCommands(use_instrument='METIS', set_modes=['img_lm'])

    # Set up detector size to match the input cube - using the full
    # METIS detector of 2048 x 2048 pixels would require too much memory
    cmd["!DET.width"] = naxis1
    cmd["!DET.height"] = naxis2

    metis = sim.OpticalTrain(cmd)

    # Set included and excluded effects
    #metis['armazones_atmo_skycalc_ter_curve'].include = True
    #metis['armazones_atmo_default_ter_curve'].include = False

    #metis['scope_surface_list'].include = False
    #metis['eso_combined_reflection'].include = True

    #metis['detector_array_list'].include = False
    #metis['detector_window'].include = True

    #metis['scope_vibration'].include = False
    #metis['detector_linearity'].include = False
    #metis['metis_psf_img'].include = False

    # We attach the spectrum of Vega to the image
    spec = SourceSpectrum.from_file(WORKDIR+"alpha_lyr_stis_008.fits")


    # Exposure parameters
    dit = header['CDELT3']
    ndit = 1

    # Loop over the input cube
    nplanes = 10
    #nplanes = naxis3  # for full cube

    nplanes = min(nplanes, naxis3)

    for i in range(nplanes):
        print("Plane", i+1, "/", nplanes)
        imhdu = fits.ImageHDU(header=imwcs.to_header(),
                              data=cube[i, :, :] * FLUX_SCALE)
        src = sim.Source(spectra=[spec], image_hdu=imhdu)

        metis.observe(src)
        outhdu = metis.readout(dit=dit, ndit=ndit)[0]

        cube[i, :, :] = outhdu[1].data.astype(np.float32)[:naxis1,:naxis2] ### REMOVE THE SLICING AT THE END!

    # We have filled up the original cube with the simulated data
    # Write out only the simulated planes
    fits.writeto(OUTFILE, data=cube[:nplanes, :, :],
                 header=header,
                 overwrite=True)

    # TODO: Return something meaningful!
    return psf_ON, psf_OFF

if __name__ == "__main__":
    sim_heeps()
