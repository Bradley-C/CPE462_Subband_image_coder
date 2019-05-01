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

% to be used in the reconstruction REMOVE IN FINAL SUBMISSION ------%
%Recompose image into color channels using inverse wavelet transform and
%approximation, horizantal, vertical, and diagonal components
rX1 = uint8(idwt2(qcA1,qcH1,qcV1,qcD1,wavelet));
rX2 = uint8(idwt2(qcA2,qcH2,qcV2,qcD2,wavelet));
rX3 = uint8(idwt2(qcA3,qcH3,qcV3,qcD3,wavelet));
%--------------------------------------------------------------------%

%Reconstruct image using the 3 color channels
Q(:,:,1)=rX1;
Q(:,:,2)=rX2;
Q(:,:,3)=rX3;

%Show reconstructed compressed image
figure;
imwrite(Q,'out.png');
imshow(Q);

%---------------- HUFFMAN ENCODING ----------------------------%   

%-------------------- RED CHANNEL -----------------------------%
%Encode approximation red channel component (qcA1)
%Find the probability of each of the symbols
A1Vec = qcA1(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(A1Vec,double(unique(A1Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[RedA1HuffDict,avglenA1]=huffmandict(symbols,prob);

%Encode the red approximation component using the dictionary for it
encodedA1 = huffmanenco(A1Vec, RedA1HuffDict);

%THESE WILL BE MOVED TO THE DECODER AND ARE JUST FOR TESTING ------%
%Decode the component
decodedA1 = huffmandeco(encodedA1, RedA1HuffDict);

%Turn the component from a vector back into an array
remadeqcA1 = reshape(decodedA1, [256, 256]);
%-------------------------REMOVE LATER ----------------------------------%

%Encode horizantal red channel component (qcH1)
%Turn the component into a vector
H1Vec = qcH1(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(H1Vec,double(unique(H1Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[RedH1HuffDict,avglenH1]=huffmandict(symbols,prob);

%Encode the red horizantal component using the dictionary for it
encodedH1 = huffmanenco(H1Vec, RedH1HuffDict);


%Encode vertical red channel component (qcV1)
%Turn the component into a vector
V1Vec = qcV1(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(V1Vec,double(unique(V1Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[RedV1HuffDict,avglenV1]=huffmandict(symbols,prob);

%Encode the red vertical component using the dictionary for it
encodedV1 = huffmanenco(V1Vec, RedV1HuffDict);



%Encode diaganol red channel component (qcD1)
%Turn the component into a vector
D1Vec = qcD1(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(D1Vec,double(unique(D1Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[RedD1HuffDict,avglenD1]=huffmandict(symbols,prob);

%Encode the red diagonal component using the dictionary for it
encodedD1 = huffmanenco(D1Vec, RedD1HuffDict);
%-------------------- RED CHANNEL -----------------------------%



%-------------------- GREEN CHANNEL -----------------------------%
%Encode approximation green channel component (qcA2)
%Find the probability of each of the symbols
A2Vec = qcA2(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(A2Vec,double(unique(A2Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[GreenA2HuffDict,avglenA2]=huffmandict(symbols,prob);

%Encode the green approximation component using the dictionary for it
encodedA2 = huffmanenco(A2Vec, GreenA2HuffDict);



%Encode horizantal green channel component (qcH2)
%Turn the component into a vector
H2Vec = qcH2(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(H2Vec,double(unique(H2Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[GreenH2HuffDict,avglenH2]=huffmandict(symbols,prob);

%Encode the green horizantal component using the dictionary for it
encodedH2 = huffmanenco(H2Vec, GreenH2HuffDict);



%Encode vertical green channel component (qcV2)
%Turn the component into a vector
V2Vec = qcV2(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(V2Vec,double(unique(V2Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[GreenV2HuffDict,avglenV2]=huffmandict(symbols,prob);

%Encode the green vertical component using the dictionary for it
encodedV2 = huffmanenco(V2Vec, GreenV2HuffDict);



%Encode diaganol green channel component (qcD2)
%Turn the component into a vector
D2Vec = qcD2(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(D2Vec,double(unique(D2Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[GreenD2HuffDict,avglenD2]=huffmandict(symbols,prob);

%Encode the green diagonal component using the dictionary for it
encodedD2 = huffmanenco(D2Vec, GreenD2HuffDict);
%-------------------- GREEN CHANNEL -----------------------------%



%-------------------- BLUE CHANNEL -----------------------------%
%Encode approximation blue channel component (qcA3)
%Find the probability of each of the symbols
A3Vec = qcA3(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(A3Vec,double(unique(A3Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[BlueA3HuffDict,avglenA3]=huffmandict(symbols,prob);

%Encode the blue approximation component using the dictionary for it
encodedA3 = huffmanenco(A3Vec, BlueA3HuffDict);


%Encode horizantal blue channel component (qcH3)
%Turn the component into a vector
H3Vec = qcH3(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(H3Vec,double(unique(H3Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[BlueH3HuffDict,avglenH3]=huffmandict(symbols,prob);

%Encode the blue horizantal component using the dictionary for it
encodedH3 = huffmanenco(H3Vec, BlueH3HuffDict);


%Encode vertical blue channel component (qcV3)
%Turn the component into a vector
V3Vec = qcV3(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(V3Vec,double(unique(V3Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[BlueV3HuffDict,avglenV3]=huffmandict(symbols,prob);

%Encode the blue vertical component using the dictionary for it
encodedV3 = huffmanenco(V3Vec, BlueV3HuffDict);



%Encode diaganol blue channel component (qcD3)
%Turn the component into a vector
D3Vec = qcD3(:);
%Find unique symbols and how many of them occur in the image
[prob, symbols] = hist(D3Vec,double(unique(D3Vec)));
%Determine the probability of each symbol occuring
prob=prob/sum(prob);

%Create a Huffman code dictionary using the symbols and probabilities we
%just found
[BlueD3HuffDict,avglenD3]=huffmandict(symbols,prob);

%Encode the blue diagonal component using the dictionary for it
encodedD3 = huffmanenco(D3Vec, BlueD3HuffDict);
%-------------------- BLUE CHANNEL -----------------------------%
