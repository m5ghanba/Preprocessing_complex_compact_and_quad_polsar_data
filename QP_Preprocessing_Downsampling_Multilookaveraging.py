# Purpose: RCM complex CP preprocessing (Main program)             Mohsen Ghanbari May 2020

# First, we ingest two files: CH and CV. So the p_directory should be modified in the program. Also, set the parameters:
# IS_THERE_BOXCAR_AVERAGING and IS_THERE_D_SAMPLING, as well as D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C, and BOXCAR_COEFF.
# Then, the program generates the elements of the CP coherence matrix: c11, c12_real, c22, c12_imag
# Then, if required, the box-car averaging and downsampling are performed.
# The output is the elements of the CP coherence matrix that need to be converted to .bin files using the corresponding MATLAB script.
# Subscene functionality added!
import tifffile as tiff
import numpy as np
import scipy.io
from scipy import ndimage
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

hh_tiff = tiff.imread(p_directory + '*HH.tif')
vv_tiff = tiff.imread(p_directory + '*VV.tif')
hv_tiff = tiff.imread(p_directory + '*HV.tif')

hh = hh_tiff[:,:,0] + 1j * hh_tiff[:,:,1]
vv = vv_tiff[:,:,0] + 1j * vv_tiff[:,:,1]
hv = hv_tiff[:,:,0] + 1j * hv_tiff[:,:,1]

c11 = np.real(np.multiply(hh, np.conj(hh)))
c22 = np.real(np.multiply(hv, np.conj(hv)))
c33 = np.real(np.multiply(vv, np.conj(vv)))
c12_real = np.real(np.multiply(hh, np.conj(hv)))
c12_imag = np.imag(np.multiply(hh, np.conj(hv)))
c13_real = np.real(np.multiply(hh, np.conj(vv)))
c13_imag = np.imag(np.multiply(hh, np.conj(vv)))
c23_real = np.real(np.multiply(hv, np.conj(vv)))
c23_imag = np.imag(np.multiply(hv, np.conj(vv)))



#taking a subscene
# import cv2
# img = cv2.imread(p_directory + "image_for_labeling.png")
# row_start = 5091
# row_end = 6091
# col_start = 1327
# col_end = 2127
#
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

