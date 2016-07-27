import os
import re
import subprocess
import glob
import astropy.io.fits as fits

from astropy import wcs
import numpy as np
import scipy.ndimage as ndimage
import scipy.stats

import sys
from copy import copy
import configparser as ConfigParser

#from pyklip.instruments.P1640_support import P1640spots
#from pyklip.instruments.P1640_support import P1640utils
#from pyklip.instruments.P1640_support import P1640_cube_checker

from scipy.interpolate import interp1d

class MAGAOData(Data):
    
    """
    A sequence of P1640 Data. Each P1640Data object has the following fields and functions 
    Args:
        filepaths: list of filepaths to occulted files
        skipslices: a list of datacube slices to skip (supply index numbers e.g. [0,1,2,3])
        corefilepaths: a list of filepaths to core (i.e. unocculted) files, for contrast calc
        spot_directory: (None) path to the directory where the spot positions are stored. Defaults to P1640.ini val
    Attributes:
        input: Array of shape (N,y,x) for N images of shape (y,x)
        centers: Array of shape (N,2) for N centers in the format [x_cent, y_cent]
        filenums: Array of size N for the numerical index to map data to file that was passed in
        filenames: Array of size N for the actual filepath of the file that corresponds to the data
        PAs: Array of N for the parallactic angle rotation of the target (used for ADI) [in degrees]
        wvs: Array of N wavelengths of the images (used for SDI) [in microns]. For polarization data, defaults to "None"
        wcs: Array of N wcs astormetry headers for each image.
        IWA: a floating point scalar (not array). Specifies to inner working angle in pixels
        output: Array of shape (b, len(files), len(uniq_wvs), y, x) where b is the number of different KL basis cutoffs
        spot_flux: Array of N of average satellite spot flux for each frame
        contrast_scaling: Flux calibration factors (multiply by image to "calibrate" flux)
        flux_units: units of output data [DN, contrast]
        prihdrs: not used for P1640, set to None
        exthdrs: Array of N P1640 headers (these are written by the P1640 cube extraction pipeline)
    Methods:
        readdata(): reread in the data
        savedata(): save a specified data in the P1640 datacube format (in the 1st extension header)
        calibrate_output(): calibrates flux of self.output
    """

    ##########################
   ### Class Initialization ###
    ##########################
    #Some static variables to define the MAGAO instrument
    centralwave = {} #in microns
    fpm_diam = {} #in pixels
    flux_zeropt = {}
    spot_ratio = {} #w.r.t. central star
    lenslet_scale = 1.0 #arcseconds per pixel (pixel scale)
    ifs_rotation = 0.0 #degrees CCW from +x axis to zenith
    
    observatory_latitude = 0.0

    #read in MAGAO configuration file and set these static variables
    package_directory = os.path.dirname(os.path.abspath(__file__))
    configfile = package_directory + "/" + "MAGAO.ini"
    config = ConfigParser.ConfigParser()
    try:
        config.read(configfile)
        #get pixel scale
        lenselet_scale = float(config.get("instrument", "ifs_lenslet_scale")) #arcsecond/pixel
        #get IFS rotation
        ifs_rotation = float(config.get("instrument", "ifs_rotation")) #degrees
        
