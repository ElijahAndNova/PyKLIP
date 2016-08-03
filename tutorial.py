#slice up h-a data cube, feed that in and get it to run (incorporate rotoff into corresponding headers)

import glob
import pyklip.instruments.GPI as GPI
import pyklip.instruments.MAGAO as MAGAO
import pyklip.parallelized as parallelized
import numpy as np
import pyklip.klip as klip

#filelist = glob.glob("../../../Python_code/20141218_H_Spec/*.fits")
filelist = glob.glob("../HD142527/HD142527/8Apr14/MERGED_long_sets/sliced/*.fits")
dataset = MAGAO.MAGAOData(filelist)

parallelized.klip_dataset(dataset, outputdir="", fileprefix="tutorialObject", annuli=9, subsections=4, movement=3, numbasis=[1,20,100], calibrate_flux=True, mode="ADI")

avgframe = np.nanmean(dataset.output[1], axis=(0,1))
calib_frame = dataset.calibrate_output(avgframe)
seps, contrast = klip.meas_contrast(calib_frame, dataset.IWA, 1.1/GPI.GPIData.lenslet_scale, 3.5)

print("Complete")
