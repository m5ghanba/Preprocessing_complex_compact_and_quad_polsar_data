# Purpose: RCM complex CP preprocessing (Main program)             Mohsen Ghanbari May 2020

# First, we ingest two files: CH and CV. So the p_directory should be modified in the program. Also, set the parameters:
# IS_THERE_BOXCAR_AVERAGING and IS_THERE_D_SAMPLING, as well as D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C, and BOXCAR_COEFF.
# Then, the program generates the elements of the CP coherence matrix: c11, c12_real, c22, c12_imag
# Then, if required, the box-car averaging and downsampling are performed.
# The output is the elements of the CP coherence matrix that need to be converted to .bin files using the corresponding MATLAB script.

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

h = tiff.imread(p_directory + '*CH.tif')
v = tiff.imread(p_directory + '*CV.tif')

ch = h[:,:,0] + 1j * h[:,:,1]
cv = v[:,:,0] + 1j * v[:,:,1]

#ch = h
#cv = v

c11 = np.real(np.multiply(ch, np.conj(ch)))
c22 = np.real(np.multiply(cv, np.conj(cv)))
c12_real = np.real(np.multiply(ch, np.conj(cv)))
c12_imag = np.imag(np.multiply(ch, np.conj(cv)))


#taking a subscene
# import cv2
# img = cv2.imread(p_directory + "image_for_labeling.png")
# img_sub = img[10000:11000, 6000:7200]# 4826*2:5637*2, 288*2:924*2      10000:11000, 6000:7200
# tiff.imwrite(p_directory + 'subscene.tif', img_sub)
# c11 = c11[10000:11000, 6000:7200]
# c12_real = c12_real[10000:11000, 6000:7200]
# c22 = c22[10000:11000, 6000:7200]
# c12_imag = c12_imag[10000:11000, 6000:7200]

if IS_THERE_BOXCAR_AVERAGING:


    fig1 = plt.figure()
    fig2 = plt.figure()
    ax_c11 = fig1.add_subplot(121)
    ax_c11.set_title('c11')
    ax_c22 = fig2.add_subplot(121)
    ax_c22.set_title('c22')
    ax_c11_filt = fig1.add_subplot(122)
    ax_c11_filt.set_title('c11-filtered')
    ax_c22_filt = fig2.add_subplot(122)
    ax_c22_filt.set_title('c22-filtered')
    ax_c11.imshow(c11, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))
    ax_c22.imshow(c22, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))

    c11 = ndimage.uniform_filter(c11, size = BOXCAR_COEFF)
    c22 = ndimage.uniform_filter(c22, size=BOXCAR_COEFF)
    c12_real = ndimage.uniform_filter(c12_real, size=BOXCAR_COEFF)
    c12_imag = ndimage.uniform_filter(c12_imag, size=BOXCAR_COEFF)

    ax_c11_filt.imshow(c11, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))
    ax_c22_filt.imshow(c22, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))


if IS_THERE_D_SAMPLING:


    fig1 = plt.figure()
    fig2 = plt.figure()
    ax_c11 = fig1.add_subplot(121)
    ax_c11.set_title('c11')
    ax_c22 = fig2.add_subplot(121)
    ax_c22.set_title('c22')
    ax_c11_filt = fig1.add_subplot(122)
    ax_c11_filt.set_title('c11-downsampled')
    ax_c22_filt = fig2.add_subplot(122)
    ax_c22_filt.set_title('c22-downsampled')
    ax_c11.imshow(c11, cmap='gray', norm=colors.Normalize(vmin=0, vmax=1))
    ax_c22.imshow(c22, cmap='gray', norm=colors.Normalize(vmin=0, vmax=1))

    c11 = block_reduce(c11, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c22 =  block_reduce(c22, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c12_real =  block_reduce(c12_real, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)
    c12_imag = block_reduce(c12_imag, block_size=(D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C), func=np.mean)



    ax_c11_filt.imshow(c11, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))
    ax_c22_filt.imshow(c22, cmap = 'gray', norm = colors.Normalize(vmin = 0, vmax =1))


# tiff.imwrite(p_directory + 'SV.tif', SV)

coherence = {'c11':c11,'c12_real':c12_real,'c22':c22,'c12_imag':c12_imag}
scipy.io.savemat(p_directory + 'CoherenceMatrixElements.mat', coherence)

