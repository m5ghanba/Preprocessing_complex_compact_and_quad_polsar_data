# Preprocessing_complex_compact_and_quad_polsar_data
Deriving multilook complex data from SLC compact and quad PolSAR data with applying multilook averaging and down-sampling 



CP_Preprocessing_Downsampling_Multilookaveraging.py
Purpose: RCM complex CP preprocessing           

First, we ingest two files: CH and CV. So the p_directory should be modified in the program. Also, set the parameters:
IS_THERE_BOXCAR_AVERAGING and IS_THERE_D_SAMPLING, as well as D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C, and BOXCAR_COEFF.
Then, the program generates the elements of the CP coherence matrix: 
c11
c12_real
c22
c12_imag
Then, if required, the box-car averaging and downsampling are performed.
The output is the elements of the CP coherence matrix that need to be converted to .bin files using the corresponding MATLAB script.


QP_Preprocessing_Downsampling_Multilookaveraging.py
Purpose: RCM complex QP preprocessing          

First, we ingest three files: *HH.tif, *HV.tif and *VV.tif. So the p_directory should be modified in the program. Also, set the parameters:
IS_THERE_BOXCAR_AVERAGING and IS_THERE_D_SAMPLING, as well as D_SAMPLING_COEFF_R, D_SAMPLING_COEFF_C, and BOXCAR_COEFF.
Then, the program generates the elements of the QP covariance matrix: 
c11
c22 
c33
c12_real
c12_imag 
c13_real 
c13_imag 
c23_real 
c23_imag 
Then, if required, the box-car averaging and downsampling are performed.
The output is the elements of the QP covariance coherence matrix that need to be converted to .bin files using the corresponding MATLAB script.
Subscene functionality and tiff output option added!
