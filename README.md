# Subband Image Compression Coder


## Project Background
Final project for the class CPE 462 Intro to Image Processing and Coding which required a small group to develop an application involving image processing topics including; 
compression, enhancement, segmenatation, restoration, or 3-D data imaging.

The subband (or wavelet) image coding technique has become very popular during the last few years. A major reason is that it significantly outperforms the current JPEG image coding standard in most occasions. In fact, a subband based coding algorithm will become the baseline in the next generation JPEG2000 image coding.

In this project, we implemented a subband image coder in a single Matlab script. Our script performs the subband decomposition, scalar quantization and entropy encoding to a typical input image, and produce a coded bit stream stored as a data file. The decoder then reads in this coded file and performs entropy decoding and subband reconstruction, and finally produce a reconstructed image in the same format of the input image. It also calculates the Peak-Signal-to-Noise ratio of your reconstructed image to evaluate the performance of the image coder. the script takes one input that controls the quantization step size, this parameter will ultimately be used to control the size of the coded data file (or the compression ratio).

##File Breakdown
subband_encoding_decoding.m - Matlab script that takes in a single .png and produces a bitstream data file named, 'binary.txt' as well as an output image when it is reconstructed from the bitstream data file
Project Report.docx - Report detailing how the script works and breakdowns the steps involved with subband image coding techniques.
/Test_Images: Contains .png files used to test the scripts performance.

### Testing
The Test_Images folder contains samples to test the script, simply change the filename variable to the desired image in this folder and run the script. in Matlab, this will give a PSNR of the output image compared to the original.