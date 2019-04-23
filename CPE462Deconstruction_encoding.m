clc;
close all;
clear all; %#ok<*CLALL>

% Input the image file
X=uint8(imread('abstract.jpg'));

% Create new image to compare to output
imwrite(X,'in.jpg');
[r,c,p] = size(X);

% Seperate image by RGB colors
X1= X(:,:,1);
X2= X(:,:,2);
X3= X(:,:,3);

% Displaying the original image
figure; 
imshow(X);
title('Original image');

% define the wavelet type and how many level of decomposition performed
wavelet = 'haar';
N = 1;
[LoD,HiD] = wfilters (wavelet,'d'); 

% Perform a level 1 decomposition of the image along the RGB channels
% [cA1,cH1,cV1,cD1] = dwt2(X1,wavelet);
% [cA2,cH2,cV2,cD2] = dwt2(X2,wavelet);
% [cA3,cH3,cV3,cD3] = dwt2(X3,wavelet);
% [height, width] = size(cA1);
%  rX1 = uint8(idwt2(qA1,qH1,qV1,qD1,wavelet));
%  rX2 = uint8(idwt2(qA2,qH2,qV2,qD2,wavelet));
%  rX3 = uint8(idwt2(qA3,qH3,qV3,qD3,wavelet));
[C1,S1] = wavedec2(X1,N,LoD,HiD);
[C2,S2] = wavedec2(X2,N,LoD,HiD);
[C3,S3] = wavedec2(X3,N,LoD,HiD);

% Compress(quantize) each channel individually
[thr1,sorh1,keepapp1] = ddencmp('cmp','wv',C1);
[thr2,sorh2,keepapp2] = ddencmp('cmp','wv',C2);
[thr3,sorh3,keepapp3] = ddencmp('cmp','wv',C3);
[XC1,CXC1,LXC1,PERF1,PERFL1] = wdencmp('gbl',C1,S1,wavelet,N,thr1,sorh1,keepapp1);
[XC2,CXC2,LXC2,PERF2,PERFL2] = wdencmp('gbl',C2,S2,wavelet,N,thr2,sorh2,keepapp2);
[XC3,CXC3,LXC3,PERF3,PERFL3] = wdencmp('gbl',C3,S3,wavelet,N,thr3,sorh3,keepapp3);

% Break down the compressed wavlets by appox. horz. verit. & diag.
cA1 = appcoef2(CXC1,LXC1,wavlet,N);
[cH1,cV1,cD1] = detcoef2('all',CXC1,LXC1,N);
cA2 = appcoef2(CXC2,LXC2,wavlet,N);
[cH2,cV2,cD2] = detcoef2('all',CXC2,LXC2,N);
cA3 = appcoef2(CXC3,LXC3,wavlet,N);
[cH3,cV3,cD3] = detcoef2('all',CXC3,LXC3,N);

% Add encoding here also make the output something you can decode into
% other file           
