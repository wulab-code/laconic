README
actin_difference_MOC.m is an example of how to calculate PCC or MOC, assuming that you 
have an image mask of cells. This mask is applied to Lifeact images and LPR images

actin_LPR_moc.m is used to calculate the PCC or MOC between a lifeact difference image and
an LPR image

actin_LPR_randomizer.m takes as input the Lifeact difference image and the LPR image. This
file scrambles the LPR image and recalculates the PCC or MOC

analyze_10x_cells_LPR_example.m is an example of taking as input a directory of segmented images and produces
the LPR of all the segmented cells in the directory. 

deep_learning_training_validation_test.m was used to calculate the error in the calculation
of fluorescence and LPR

hires_LPR_example.m is an example of how the high magnification (per pixel LPR) was performed

segmentation_dl.m uses the deep learning network to segment images. 
extract_cells.m, filtercellareas.m, getcellinterior.m, getminicell.m, resizetophat.m are 
helper files. 

linkcells.m links together trajectories of segmented cells. 

The deep learning network is too large for github to host but is available on request