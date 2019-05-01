clc;
close all;
clear all; %#ok<*CLALL>

% Input the image file
X=uint8(imread('colorwheel.png'));

% Create new image to compare to output
imwrite(X,'in.png');
[r,c,p] = size(X);

% Seperate image by RGB colors
X1= X(:,:,1);
X2= X(:,:,2);
X3= X(:,:,3);

% Displaying the original image
figure; 
imshow(X);
title('Original image');

% define the decomposition wavelet type and how many levels thershhold
% levels to use for quantization
wavelet = 'haar'; %type of decompisition
N=10; %levels to quantize image by (can be 1 to 20)

% Perform a level 1 decomposition of the image along the RGB channels
[cA1,cH1,cV1,cD1] = dwt2(X1,wavelet);
[cA2,cH2,cV2,cD2] = dwt2(X2,wavelet);
[cA3,cH3,cV3,cD3] = dwt2(X3,wavelet);
[height, width] = size(cA1);

% Find N threshhold levels (only need to do one set since all sets share
% range of values)
threshA=multithresh(cA1,N);
threshH=multithresh(cH1,N);
threshV=multithresh(cV1,N);
threshD=multithresh(cD1,N);

% Contrusts a vector such that the min value in each quantization interval
% is assigned to the N levels of the output image
valuesMinA= [min(cA1(:)) threshA];
valuesMinH= [min(cH1(:)) threshH];
valuesMinV= [min(cV1(:)) threshD];
valuesMinD= [min(cD1(:)) threshV];

% Quantize the components of the decomposed image using N threshhold levels
%These quantitized components will be encoded later (qcA1, qcH1, qcV1, qcD1
%etc.)
[qcA1, indexA1] = imquantize(cA1,threshA,valuesMinA);
[qcH1, indexH1] = imquantize(cH1,threshH,valuesMinH);
[qcV1, indexV1] = imquantize(cV1,threshV,valuesMinV);
[qcD1, indexD1] = imquantize(cD1,threshD,valuesMinD);

[qcA2, indexA2] = imquantize(cA2,threshA,valuesMinA);
[qcH2, indexH2] = imquantize(cH2,threshH,valuesMinH);
[qcV2, indexV2] = imquantize(cV2,threshV,valuesMinV);
[qcD2, indexD2] = imquantize(cD2,threshD,valuesMinD);

[qcA3, indexA3] = imquantize(cA3,threshA,valuesMinA);
[qcH3, indexH3] = imquantize(cH3,threshH,valuesMinH);
[qcV3, indexV3] = imquantize(cV3,threshV,valuesMinV);
[qcD3, indexD3] = imquantize(cD3,threshD,valuesMinD);

% to be used in the reconstruction REMOVE IN FINAL SUBMISSION
%Recompose image into color channels using inverse wavelet transform and
%approximation, horizantal, vertical, and diagonal components
rX1 = uint8(idwt2(qcA1,qcH1,qcV1,qcD1,wavelet));
rX2 = uint8(idwt2(qcA2,qcH2,qcV2,qcD2,wavelet));
rX3 = uint8(idwt2(qcA3,qcH3,qcV3,qcD3,wavelet));

%Reconstruct image using the 3 color channels
Q(:,:,1)=rX1;
Q(:,:,2)=rX2;
Q(:,:,3)=rX3;

%Show reconstructed compressed image
figure;
imwrite(Q,'out.png');
imshow(Q);

% Add encoding here also make the output something you can decode into
% other file           

%Encode approximation of red channel component (qcA1)
%Find the probability of each of the symbols (pixel values 0-255)
A1Vec = qcA1(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(A1Vec,double(unique(A1Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[RedA1HuffDict,avglen]=huffmandict(symbols,prob);

%Encode the red approximation component using the dictionary for it
encodedA1Img = huffmanenco(A1Vec, RedA1HuffDict);

