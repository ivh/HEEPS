from heeps.util.img_processing import resize_cube
from heeps.util.coord import mas2rms
import os.path
from astropy.io import fits
import numpy as np

def load_errors(nframes=20, nstep=1, npupil=285, add_phase=True, add_amp=False, 
        f_phase='', f_amp='', add_point_err=False, f_point_err='', 
        add_apo_drift=False, apo_drift=0.02, verbose=False, **conf):
    
    ''' Load wavefront errors.
    
    nframes (int):
        number of frames to crop the input data
    nstep (int):
        take 1 frame every nstep (cubesize = nframes/nstep)
    npupil (int):
        size of pupil
    lam (float):
        wavelength in m
    add_phase (bool):
        true if adding phase screens
    add_amp (bool):
        true if adding amplitude screens
    f_phase (str):
        path to phase screens fits file
    f_amp (str):
        path to amplitude screens fits file
    add_point_err (bool):
        true if adding pointing errors
    f_point_err (str):
        path to pointing errors fits file
    add_apo_drift (bool):
        true if adding apodizer drift

    '''
    # size of cubes/arrays
    nscreens = int((nframes/nstep) + 0.5)

    # load phase screens
    if add_phase is True:
        assert(os.path.isfile(f_phase) and os.path.splitext(f_phase)[1] == '.fits'), \
            "'f_phase' must be a valid fits file."
        phase_screens = fits.getdata(f_phase)[:nframes][::nstep] # in meters
        phase_screens = np.nan_to_num(phase_screens)
        phase_screens = resize_cube(phase_screens, npupil)
        if verbose is True:
            print("Load phase screens from '%s'"%os.path.basename(f_phase))
            print('   nscreens=%s (nframes=%s, nstep=%s)'%(len(phase_screens), nframes, nstep))
    else:
        phase_screens = [None]*nscreens

    # load amp screens
    if add_amp is True:
        assert(os.path.isfile(f_amp) and os.path.splitext(f_amp)[1] == '.fits'), \
            "'f_amp' must be a valid fits file."
        amp_screens = np.array(fits.getdata(f_amp), ndmin=3)
        if len(amp_screens) > 1: # cube
            amp_screens = amp_screens[:nframes][::nstep]
        amp_screens = np.nan_to_num(amp_screens)
        amp_screens = resize_cube(amp_screens, npupil)
        if verbose is True:
            print("Load amp screens from '%s'"%os.path.basename(f_amp))
            print('   nscreens=%s'%(len(amp_screens)))
    else:
        amp_screens = np.array([None]*nscreens)
    
    # load pointing errors (in mas)
    if add_point_err is True:
        tiptilts = np.array(fits.getdata(f_point_err), ndmin= 2)
        if len(tiptilts) > 1: # cube
            tiptilts = tiptilts[:nframes][::nstep]
        # convert mas to rms phase error
        tiptilts = mas2rms(tiptilts, conf['diam_ext'])
        if verbose is True:
            print("Load pointing errors from '%s'"%os.path.basename(f_point_err))
            print('   nscreens=%s'%(len(tiptilts)))
    else:
        tiptilts = np.array([None]*nscreens)

    # load apodizer drift
    if add_apo_drift is True and 'RAVC' in conf['mode']:
        misaligns = np.array([[x,0,0,0,0,0] \
            for x in np.linspace(-apo_drift/2, apo_drift/2, 12000)])[:nframes][::nstep]
        if verbose is True:
            print('Load apodizer drit=%s %% ptv'%apo_drift)
    else:
        misaligns = np.array([None]*nscreens)

    
    return phase_screens, amp_screens, tiptilts, misaligns