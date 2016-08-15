#slice up h-a data cube, feed that in and get it to run (incorporate rotoff into corresponding headers)

import glob
import pyklip.instruments.GPI as GPI
import pyklip.instruments.MAGAO as MAGAO
import pyklip.parallelized as parallelized
import numpy as np
import pyklip.klip as klip
from astropy.io import fits

#filelist = glob.glob("20141218_H_Spec/*.fits")
filelist = glob.glob("../HD142527/HD142527/8Apr14/MERGED_long_sets/sliced/*.fits")
#filelist = glob.glob("spiral/sliced/*.fits")
dataset = MAGAO.MAGAOData(filelist)
#dataset = GPI.GPIData(filelist)

outputFileName = "tutorialObject"

parallelized.klip_dataset(dataset, outputdir="", fileprefix=outputFileName, annuli=1, subsections=1, movement=1, numbasis=[1,2,3,4,5,10,20,50,100], calibrate_flux=False, mode="ADI")

print("Shape of dataset.output is " + str(dataset.output.shape))
print("Shape of dataset.output[1] is " + str(dataset.output[1].shape))
avgframe = np.nanmean(dataset.output[1], axis=(0,1))
print("Shape of avgframe is " + str(avgframe.shape))
calib_frame = dataset.calibrate_output(avgframe)

print("Shape of calib_frame: " + str(calib_frame.shape))
#seps, contrast = klip.meas_contrast(calib_frame, dataset.IWA, 1.1/GPI.GPIData.lenslet_scale, 3.5)

print("Completed klipping. Rotating images")
hdulist = fits.open(outputFileName+"-KLmodes-all.fits")
cube = hdulist[1].data
hdulist.close()
cube = cube[:,:,::-1]
"""
newCube = []
for i in range(len(cube)):
    newCube.append(cube[len(cube)-i-1])
"""
hdulist = fits.PrimaryHDU(cube)
hdulist.writeto("_"+outputFileName+"-KLmodes-all.fits")


print("Complete")
