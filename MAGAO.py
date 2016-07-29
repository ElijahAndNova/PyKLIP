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

    #I'm marking things that I'm not sure if we need with a "#!"

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
        lenselet_scale = float(config.get("instrument", "ifs_lenslet_scale")) #!
        #get IFS rotation
        ifs_rotation = float(config.get("instrument", "ifs_rotation"))
        bands = ['HA'. 'CONT']
        for band in bands:
            centralwave[band] = float(config.get("instrument", "cen_wave_{0}".format(band)))
            fpm_diam[band] = float(config.get("instrument", "fpm_dian_{0}".format(band))) #!
            flux_zeropt[band] = float(config.get("instrument", "zero_pt_flux_{0}".format(band))) #!
        observatory_latitude = float(vonfig.get("observatory", "observatory_lat"))
    except ConfigParser.Error as e:
        print("Error reading MAGAO configuration file: {0}".format(e.message))
        raise e
    
    #########################
   ###    Constructors     ###
    #########################
    def __init__(self, filepaths=None):
        """
        Initialization code for MAGAOData
        """
        super(MAGAOData, self).__init__()
        self._output = None
        if filepaths is None:
            self._input = None
            self._centers = None
            self._filenums = None
            self._filenames = None
            self._PAs = None
            self._wvs = None
            self._wcs = None
            self._IWA = None
            self.spot_flux = None #!
            self.star_flux = None
            self.contrast_scaling = None
            self.prihdrs = None
            self.exthdrs = None
        else:
            self.readdata(filepaths)
    
    ##############################
   ### Instance Required Fields ###
    ##############################
    @properrt
    def input(self):
        return self._input
    @input.setter
    def input(self, newval):
        self._input = newval
    
    @property
    def centers(self):
        return self._centers
    @centers.setter
    def centers(self, newval):
        self._centers = newval

    @property
    def PAs(self):
        return self._PAs
    @PAs.setter
    def PAs(self, newval):
        self._PAs = newval
    
    @property
    def wvs(self):
        return self._wvs
    @wvs.setter
    def wvs(self, newval):
        self._wvs = newval
    
    @property
    def wcs(self):
        return self._wcs
    @wcs.setter
    def wcs(self, newval):
        self._wcs = newval

    @property
    def IWA(self):
        return self._IWA
    @IWA.setter
    def IWA(self, newval):
        self._IWA = newval
    
    @property
    def output(self):
        return self._output
    @output.setter
    def output(self, newval):
        self._output = newval

    ###################
   ###    Methods    ###
    ###################
        
    def readdata(self, filepaths):
        """
        Method to open and read a list of MAGAO data
        """
        if isinstance(filepaths, str):
            filepaths = [filepaths]

        data = []
        filenums = []
        filenaes = []
        rot_angles = []
        wvs = []
        centers = []
        wcs_hdrs = []
        star_fluxes = []
        spot_fluxes = [] #!
        prihdrs = []
        
        for index, filepath in enumerate(filepaths):
            cube, center, pa, wv, astr_hdrs, filt_band, fpm_band, ppm_band, star_flux, spot_flux, prihdr, exthdr = _magao_process_file(filepath, index)
        
            data.append(cube)
            centers.append(center)
            star_fluxes.append(star_flux)
            spot_fluxes.append(spot_flux) #!
            rot_angles.append(pa)
            wvs.append(wv)
            filenums.append(np.ones(pa.shape[0]) * index)
            wcs_hdrs.append(astr_hdrs) #!
            prihdrs.append(prihdr)
            filenames.append([filepath for i in range(pa.shape[0])])
            
        data = np.array(data)
        dims = data.shape
        data = data.reshape([dims[0] * dims[1], dims[2], dims[3]])
        filenums = np.array(filenums).reshape([dims[0] * dims[1]])
        filenames = np.array(filenames).reshape([dims[0] * dims[1]])
        rot_angles = np.array(rot_angles).reshape([dims[0] * dims[1]])
        wvs = np.array(wvs).reshape([dims[0] * dims[1]])
        wcs_hdrs = np.array(wcs_hdrs).reshape([dims[0] * dims[1]])
        centers = np.array(centers).reshape([dims[0] * dims[1], 2])
        star_fluxes = np.array(star_fluxes).reshape([dims[0] * dims[1]])
        spot_fluxes = np.array(spot_fluxes).reshape([dims[0] * dims[1]]) #!

        self._input = data
        self._centers = centers
        self._filenums = filenums
        self._filenames = filenames
        self._PAs = rot_angles
        self._wvs = wvs
        self._wcs = None #wvs_hdrs
        self.spot_flux = spot_fluxes
        self.IWA = MAGAOData.fpm_diam[fpm_band] / 2.0 #!
        self.star_flux = star_fluxes
        self.contrast_scaling = 1./star_fluxes
        self.prihdrs = prihdrs

        
    def _magao_process_file(filepath, filetype):
        #filetype == 0 --> HA
        #filetype == 1 --> CONT
        print("Reading File: {0}".format(filepath))
        hdulist = fits.open(os.getcwd()+"/rotoff_preproc.fits")
        rotangles = hdulist[0].data
        rotangles = np.array(rotangles)
        hdulist.close()
        hdulist = fits.open(filepath)
        try:
            cube = hdulist[0].data
            exthdr = None
            prihdr = hdulist[0].header
            
            if filetype == 0:
                filt_band = "H-Alpha"
            elif filetype == 1:
                filt_band = "Continuum"
            
            fpm_band = filt_band
            ppm_band = None
            
            wvs = [1.0]
            center = [[225,225]]
            
            dims = cube.shape
            x, y = np.meshgrid(np.arange(dims[1], dtype=np.float32), np.arange(dims[0], dtype=np.float32))
            nx = center[0][0] - (x - center[0][0])
            minval = np.min([np.nanmin(cube), 0.0])
            flipped.cube = ndimage.map_coordinates(np.copy(cube), [y, nx], cval=minval * 5.0)
            
            star_flux = calc_starflux(flipped_cube, center) #WRITE THIS FUNCTION
            cube = flipped_cube.reshape([1, flipped_cube.shape[0], flipped_cube.shape[1]])
            parang = rotangles
            astr_hdrs = np.repeat(None, 1)
            spot_fluxes = [[1]] #!
        finally:
            hdulist.close()
        
        return cube, center, parang, wvs, astr_hdrs, filt_band, fpm_band, ppm_band, star_flux, spot_fluxes, prihdr, exthdr



