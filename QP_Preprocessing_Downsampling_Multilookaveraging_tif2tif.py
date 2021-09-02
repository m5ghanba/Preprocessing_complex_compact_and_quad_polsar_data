# Purpose: RCM complex CP preprocessing (Main program)             Mohsen Ghanbari May 2020

# First, we ingest one file C3.tif. So the p_directory should be modified in the program. Also, set the parameters:
# IS_THERE_BOXCAR_AVERAGING and IS_THERE_D_SAMPLING, as well as D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C, and BOXCAR_COEFF.
# Then, the program generates the elements of the CP coherence matrix: c11, c12_real, c22, c12_imag
# Then, if required, the box-car averaging and downsampling are performed.
# The output is the elements of the CP coherence matrix that need to be converted to .bin files using the corresponding MATLAB script.

import tifffile as tiff
import numpy as np
import scipy.io
from scipy import ndimage
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from skimage.measure import block_reduce
from tkinter import Tk
from tkinter.filedialog import askdirectory

IS_THERE_BOXCAR_AVERAGING = True
IS_THERE_D_SAMPLING = False


D_SAMPLING_COEFF_R = 2 # downsampling coefficient: the number of rows in the window in which the pixels are replaced with the average value
D_SAMPLING_COEFF_C = 2 # downsampling coefficient: the number of columns in the window in which the pixels are replaced with the average value

BOXCAR_COEFF = 5 # boxcar averaging coefficient: the size of the averaging window


# p_directory = 'C:/Users/vip-mohsen/Desktop/Results_realCp/Winnipeg/Calibrated_by_PCI/without boxcar/'
p_directory = askdirectory(title='Select Folder') # shows dialog box and return the path
p_directory = p_directory + "/"

C3 = tiff.imread(p_directory + 'C3.tif')

c11 = C3[:,:,0]
c12_real = C3[:,:,1]
c22 = C3[:,:,2]
c13_real = C3[:,:,3]
c23_real = C3[:,:,4]
c33 = C3[:,:,5]
c12_imag = C3[:,:,6]
c13_imag = C3[:,:,7]
c23_imag = C3[:,:,8]



#taking a subscene
# import cv2
# img = cv2.imread(p_directory + "image_for_labeling.png")
# row_start = 0
# row_end = 700
# col_start = 150
# col_end = 650
# img_sub = img[row_start:row_end, col_start:col_end]
# tiff.imwrite(p_directory + 'subscene.tif', img_sub)
# c11 = c11[row_start:row_end, col_start:col_end]
# c12_real = c12_real[row_start:row_end, col_start:col_end]
# c22 = c22[row_start:row_end, col_start:col_end]
# c12_imag = c12_imag[row_start:row_end, col_start:col_end]
# c13_real = c13_real[row_start:row_end, col_start:col_end]
# c33 = c33[row_start:row_end, col_start:col_end]
# c13_imag = c13_imag[row_start:row_end, col_start:col_end]
# c23_real = c23_real[row_start:row_end, col_start:col_end]
# c23_imag = c23_imag[row_start:row_end, col_start:col_end]

if IS_THERE_BOXCAR_AVERAGING:


    fig1 = plt.figure()
    fig2 = plt.figure()
    fig3 = plt.figure()

    ax_c11 = fig1.add_subplot(121)
    ax_c11.set_title('c11')
    ax_c22 = fig2.add_subplot(121)
    ax_c22.set_title('c22')
    ax_c33 = fig3.add_subplot(121)
    ax_c33.set_title('c33')

    ax_c11_filt = fig1.add_subplot(122)
    ax_c11_filt.set_title('c11-filtered')
    ax_c22_filt = fig2.add_subplot(122)
    ax_c22_filt.set_title('c22-filtered')
    ax_c33_filt = fig3.add_subplot(122)
    ax_c33_filt.set_title('c33-filtered')

    ax_c11.imshow(c11, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))
    ax_c22.imshow(c22, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))
    ax_c33.imshow(c33, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))

    # Save just the portion inside the axes
    # extent = ax_c11.get_window_extent().transformed(fig1.dpi_scale_trans.inverted())
    # fig1.savefig(p_directory + 'HH.png', bbox_inches=extent)
    # extent = ax_c22.get_window_extent().transformed(fig2.dpi_scale_trans.inverted())
    # fig2.savefig(p_directory + 'HV.png', bbox_inches=extent)

    my_dpi = 96
    fig_hh = plt.figure(frameon=False, figsize=(c11.shape[0]/my_dpi, c11.shape[0]/my_dpi), dpi=my_dpi)
    ax_hh = plt.Axes(fig_hh, [0., 0., 1., 1.])
    ax_hh.set_axis_off()
    fig_hh.add_axes(ax_hh)
    fig_hv = plt.figure(frameon=False, figsize=(c11.shape[0]/my_dpi, c11.shape[0]/my_dpi), dpi=my_dpi)
    ax_hv = plt.Axes(fig_hv, [0., 0., 1., 1.])
    ax_hv.set_axis_off()
    fig_hv.add_axes(ax_hv)
    ax_hh.imshow(c11, cmap='gray', norm=colors.Normalize(vmin=0, vmax=1))
    ax_hv.imshow(c22, cmap='gray', norm=colors.Normalize(vmin=0, vmax=1))
    fig_hh.savefig(p_directory + 'HH.png', bbox_inches='tight', pad_inches=0, dpi=my_dpi)
    fig_hv.savefig(p_directory + 'HV.png', bbox_inches='tight', pad_inches=0, dpi=my_dpi)

    c11 = ndimage.uniform_filter(c11, size = BOXCAR_COEFF)
    c22 = ndimage.uniform_filter(c22, size=BOXCAR_COEFF)
    c33 = ndimage.uniform_filter(c33, size=BOXCAR_COEFF)
    c12_real = ndimage.uniform_filter(c12_real, size=BOXCAR_COEFF)
    c12_imag = ndimage.uniform_filter(c12_imag, size=BOXCAR_COEFF)
    c13_real = ndimage.uniform_filter(c13_real, size=BOXCAR_COEFF)
    c13_imag = ndimage.uniform_filter(c13_imag, size=BOXCAR_COEFF)
    c23_real = ndimage.uniform_filter(c23_real, size=BOXCAR_COEFF)
    c23_imag = ndimage.uniform_filter(c23_imag, size=BOXCAR_COEFF)


    ax_c11_filt.imshow(c11, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))
    ax_c22_filt.imshow(c22, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))
    ax_c33_filt.imshow(c33, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))

    # Save just the portion inside the axes
    # extent = ax_c11_filt.get_window_extent().transformed(fig1.dpi_scale_trans.inverted())
    # fig1.savefig(p_directory + 'HH_avg.png', bbox_inches=extent)
    # extent = ax_c22_filt.get_window_extent().transformed(fig2.dpi_scale_trans.inverted())
    # fig2.savefig(p_directory + 'HV_avg.png', bbox_inches=extent)

    fig_hh_f = plt.figure(frameon=False, figsize=(c11.shape[0]/my_dpi, c11.shape[0]/my_dpi), dpi=my_dpi)
    ax_hh_f = plt.Axes(fig_hh_f, [0., 0., 1., 1.])
    ax_hh_f.set_axis_off()
    fig_hh_f.add_axes(ax_hh_f)
    fig_hv_f = plt.figure(frameon=False, figsize=(c11.shape[0]/my_dpi, c11.shape[0]/my_dpi), dpi=my_dpi)
    ax_hv_f = plt.Axes(fig_hv_f, [0., 0., 1., 1.])
    ax_hv_f.set_axis_off()
    fig_hv_f.add_axes(ax_hv_f)
    ax_hh_f.imshow(c11, cmap='gray', norm=colors.Normalize(vmin=0, vmax=1))
    ax_hv_f.imshow(c22, cmap='gray', norm=colors.Normalize(vmin=0, vmax=1))
    fig_hh_f.savefig(p_directory + 'HH_avg.png', bbox_inches='tight', pad_inches=0)
    fig_hv_f.savefig(p_directory + 'HV_avg.png', bbox_inches='tight', pad_inches=0)

if IS_THERE_D_SAMPLING:

    fig1 = plt.figure()
    fig2 = plt.figure()
    fig3 = plt.figure()
    ax_c11 = fig1.add_subplot(121)
    ax_c11.set_title('c11')
    ax_c22 = fig2.add_subplot(121)
    ax_c22.set_title('c22')
    ax_c33 = fig3.add_subplot(121)
    ax_c33.set_title('c33')
    ax_c11_filt = fig1.add_subplot(122)
    ax_c11_filt.set_title('c11-downsampled')
    ax_c22_filt = fig2.add_subplot(122)
    ax_c22_filt.set_title('c22-downsampled')
    ax_c33_filt = fig3.add_subplot(122)
    ax_c33_filt.set_title('c33-downsampled')
    ax_c11.imshow(c11, cmap='gray', norm=colors.Normalize(vmin=0, vmax=1))
    ax_c22.imshow(c22, cmap='gray', norm=colors.Normalize(vmin=0, vmax=1))
    ax_c33.imshow(c33, cmap='gray', norm=colors.Normalize(vmin=0, vmax=1))

    c11 = block_reduce(c11, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c22 =  block_reduce(c22, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c33 =  block_reduce(c33, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c12_real =  block_reduce(c12_real, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c12_imag = block_reduce(c12_imag, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c13_real =  block_reduce(c13_real, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c13_imag = block_reduce(c13_imag, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c23_real =  block_reduce(c23_real, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c23_imag = block_reduce(c23_imag, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)

    ax_c11_filt.imshow(c11, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))
    ax_c22_filt.imshow(c22, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))
    ax_c33_filt.imshow(c33, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))




C3 = np.zeros((c11.shape[0], c11.shape[1], 9))
C3[:,:,0] = c11
C3[:,:,1] = c12_real
C3[:,:,2] = c22
C3[:,:,3] = c13_real
C3[:,:,4] = c23_real
C3[:,:,5] = c33
C3[:,:,6] = c12_imag
C3[:,:,7] = c13_imag
C3[:,:,8] = c23_imag
tiff.imwrite(p_directory + 'C3.tif', C3)

# 0:C11, 1:C12_real, 2:C22, 3:C13_real, 4:C23_real, 5:C33, 6:C12_imag, 7:C13_imag, 8:C23_imag
covariance = {'c11':c11,'c12_real':c12_real, 'c22':c22, 'c13_real':c13_real, 'c23_real':c23_real, 'c33':c33,
             'c12_imag':c12_imag, 'c13_imag':c13_imag, 'c23_imag':c23_imag}

scipy.io.savemat(p_directory + 'CovarianceMatrixElements.mat', covariance)

